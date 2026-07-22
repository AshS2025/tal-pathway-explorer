"""
pathway_scoring.py
==================

Wraps DORAnet's `pathway_ranking()` and parses its output into
structured Python data so the UI (or any post-processing step) can
display per-criterion scores, sort by composite score, etc.

DORAnet's `pathway_ranking()` writes a file `{job}_ranked_pathways.txt`
with this per-pathway block format:

    ranking N
    final score X.X
    atomic economy X.XX
    pathway by-product NNN
    intermediate by-product {smiles: count, ...}
    reaction SMILES, name, and enthalpy:
    <K reaction SMILES lines>
    <K reaction names>
    <K enthalpies>
    <blank line>

We parse each block into a `RankedPathway` dataclass.

The default weights match DORAnet's own defaults (steps=4, thermo=2,
by-products=2, atom_economy=1). Chemists can override them by passing
their own `weights` dict to `generate_base_rankings`.
"""

from __future__ import annotations

import ast
import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Descriptors

from doranet.modules.post_processing.post_processing import pathway_ranking


DEFAULT_WEIGHTS = {
    "reaction_thermo": 2,
    "number_of_steps": 4,
    "by_product_number": 2,
    "atom_economy": 1,
    "salt_score": 0,
    "in_reaxys": 0,
    "coolness": 0,
}


@dataclass
class RankedPathway:
    """One pathway's structured data after DORAnet's pathway_ranking."""

    rank: int
    final_score: float
    atomic_economy: float
    pathway_byproduct_count: int
    intermediate_byproducts: dict           # {smiles: count}
    reaction_smiles: list                   # ["reactants>>products", ...]
    reaction_names: list                    # ["Dehydration", ...]
    reaction_enthalpies: list               # [float | "No_Thermo", ...]
    # Equilibrator scoring — filled in by decorate_with_equilibrator()
    # after DORAnet's ranking finishes. `equilibrator_dgs` is a per-step
    # list of ΔG'° values (kJ/mol) or None where a step's compounds
    # aren't in equilibrator's database. max/avg summarise across the
    # scoreable steps only.
    equilibrator_dgs: list = field(default_factory=list)
    equilibrator_max_dg: Optional[float] = None
    equilibrator_avg_dg: Optional[float] = None
    equilibrator_coverage: float = 0.0        # fraction of steps we could score
    # Lemnisca blend — filled in by apply_lemnisca_blend().
    #   lemnisca_score  = the Lemnisca sub-score (geomean of the custom
    #                     criteria: stability, diversity, ...), 0–1.
    #   blended_score   = the FINAL ranking key = geomean(DORAnet, Lemnisca)
    #                     with the tier-2 layer weights, 0–1.
    #   lemnisca_components = raw 0–1 grade per criterion (incl. "doranet").
    lemnisca_score: Optional[float] = None
    blended_score: Optional[float] = None
    lemnisca_components: dict = field(default_factory=dict)

    @property
    def num_steps(self) -> int:
        return len(self.reaction_smiles)

    @property
    def max_dh(self):
        """Max reaction enthalpy across steps (kJ/mol). Returns None if
        the pathway has no thermo data at all (every step was
        'No_Thermo', meaning RMG wasn't enabled during expansion).

        DORAnet's `pathway_ranking` uses max dH — not average — as the
        thermo criterion. The chemist rationale: a pathway's
        feasibility is bottlenecked by its worst (most endothermic)
        step. One reaction at dH = +80 kJ/mol dominates the ranking
        even if the other 5 steps are all exothermic.
        """
        real = [dh for dh in self.reaction_enthalpies if dh is not None]
        return max(real) if real else None

    @property
    def avg_dh(self):
        """Average reaction enthalpy across steps with real thermo data
        (kJ/mol). Not used for ranking (DORAnet uses max_dh) but often
        useful for eyeballing overall pathway energetics."""
        real = [dh for dh in self.reaction_enthalpies if dh is not None]
        return sum(real) / len(real) if real else None


def generate_base_rankings(
    starter: str,
    target: str,
    helpers: list,
    job_name: str = "TAL",
    weights: Optional[dict] = None,
    molecule_thermo_calculator=None,
    max_rxn_thermo_change: float = 15,
) -> list:
    """
    Run DORAnet's pathway ranking on the pathways in
    `{job_name}_pathways.txt` and parse the resulting
    `{job_name}_ranked_pathways.txt` into a list of RankedPathway.

    Returns the list of ranked pathways (sorted highest score first).
    """
    if weights is None:
        weights = DEFAULT_WEIGHTS
    # num_process=1 is critical when running under Streamlit on Windows.
    # DORAnet's default multiprocessing spawns worker processes that
    # re-import the calling module — under Streamlit, that's
    # streamlit_app.py, which references variables (like `pathways`) that
    # only exist inside a Run handler, not at module scope. The workers
    # crash on import and Streamlit retries forever. Single-process is
    # only marginally slower for the ~50-pathway workloads we handle.
    pathway_ranking(
        starters=[starter],
        helpers=helpers,
        target=[target],
        weights=weights,
        num_process=1,
        job_name=job_name,
        molecule_thermo_calculator=molecule_thermo_calculator,
        max_rxn_thermo_change=max_rxn_thermo_change,
    )
    return parse_ranked_pathways(job_name)


import re as _re
_BIO_OP_PATTERN = _re.compile(r"^rule\d+")   # JN1224MIN: rule0087, rule1118, ...


def is_bio_op(op_name: str) -> bool:
    """
    Return True if this reaction was fired by a bio (enzymatic) operator.
    We rely on JN1224MIN's naming convention: bio operator names look
    like "rule0087", "rule1118", etc. Chem operators have human-readable
    names like "Dehydration of Alcohol". This heuristic matches every
    bio rule DORAnet ships without needing a separate lookup.
    """
    if not op_name:
        return False
    return bool(_BIO_OP_PATTERN.match(op_name.strip()))


def decorate_with_equilibrator(
    ranked_pathways: list,
    equilibrator_client,
    max_abs_dg_threshold: Optional[float] = None,
) -> list:
    """
    Compute equilibrator ΔG'° per reaction for each ranked pathway,
    populate the RankedPathway.equilibrator_* fields, and (optionally)
    prune pathways whose worst BIO step exceeds `max_abs_dg_threshold`.

    Routing: equilibrator is designed for enzymatic reactions in the
    KEGG universe. We only run it on reactions fired by a bio operator
    (identified via `is_bio_op(name)`). Chem reactions get None so RMG
    is the only thermo signal for them, matching each tool's domain of
    validity.

    Pruning rule: pathways whose worst BIO step |ΔG'°| > threshold are
    dropped. Chem-only pathways cannot be pruned by equilibrator (they
    have no bio steps to check). Reactions we couldn't score
    (bio step whose compounds aren't in equilibrator's DB) don't count
    against the threshold — conservative: we only drop what we can
    prove is bad.
    """
    if equilibrator_client is None:
        return ranked_pathways

    out = []
    for p in ranked_pathways:
        # Only score bio steps. Chem steps get None.
        dgs = []
        for smi, name in zip(p.reaction_smiles, p.reaction_names):
            if is_bio_op(name):
                dgs.append(equilibrator_client.dG_prime(smi))
            else:
                dgs.append(None)
        p.equilibrator_dgs = dgs

        # Aggregates only across bio reactions we actually scored.
        n_bio = sum(1 for name in p.reaction_names if is_bio_op(name))
        real = [d for d in dgs if d is not None]
        if real:
            p.equilibrator_max_dg = max(real, key=abs)   # worst |ΔG|
            p.equilibrator_avg_dg = sum(real) / len(real)
            # coverage = fraction of BIO steps we successfully scored
            # (not fraction of all steps, since chem steps aren't in
            # scope for equilibrator). 100% if this pathway has no
            # bio steps at all.
            p.equilibrator_coverage = (
                len(real) / n_bio if n_bio > 0 else 1.0
            )
        else:
            p.equilibrator_max_dg = None
            p.equilibrator_avg_dg = None
            p.equilibrator_coverage = 1.0 if n_bio == 0 else 0.0

        if max_abs_dg_threshold is not None:
            if (p.equilibrator_max_dg is not None
                    and abs(p.equilibrator_max_dg) > max_abs_dg_threshold):
                continue

        out.append(p)
    return out


def parse_ranked_pathways(job_name: str) -> list:
    """
    Read `{job_name}_ranked_pathways.txt` and return a list of
    `RankedPathway` objects.

    File format (one block per pathway, blank-line separated):
        ranking N
        final score X.X
        atomic economy X.XX
        pathway by-product NNN
        intermediate by-product {smiles: count, ...}
        reaction SMILES, name, and enthalpy:
        <K reaction SMILES>
        <K reaction names>
        <K enthalpies>
    """
    path = f"{job_name}_ranked_pathways.txt"
    with open(path, encoding="utf-8") as f:
        text = f.read()

    # Split on blank lines — each block is one pathway. Trailing/leading
    # blank lines are safely handled by filtering empties.
    blocks = [
        b.strip() for b in text.split("\n\n") if b.strip()
    ]
    pathways = []
    for block in blocks:
        pathway = _parse_one_block(block)
        if pathway is not None:
            pathways.append(pathway)
    return pathways


def _parse_one_block(block: str) -> Optional[RankedPathway]:
    """Parse a single pathway block into a RankedPathway. Returns None
    if the block looks malformed (e.g. missing the reactions section)."""
    lines = block.split("\n")

    # First few lines are labelled key: value pairs
    rank = None
    final_score = None
    atomic_economy = None
    pathway_byproduct_count = None
    intermediate_byproducts = {}
    header_idx = None

    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith("ranking "):
            rank = int(line.split(None, 1)[1])
        elif line.startswith("final score "):
            final_score = float(line[len("final score "):].strip())
        elif line.startswith("atomic economy "):
            atomic_economy = float(line[len("atomic economy "):].strip())
        elif line.startswith("pathway by-product "):
            pathway_byproduct_count = int(
                line[len("pathway by-product "):].strip()
            )
        elif line.startswith("intermediate by-product "):
            raw = line[len("intermediate by-product "):].strip()
            try:
                intermediate_byproducts = ast.literal_eval(raw)
            except (ValueError, SyntaxError):
                intermediate_byproducts = {}
        elif line.startswith("reaction SMILES"):
            header_idx = i
            break

    if header_idx is None or rank is None:
        return None

    # Everything after the header line is K reaction lines + K names +
    # K enthalpies. Total remaining lines = 3 * K.
    remaining = [ln for ln in lines[header_idx + 1:] if ln.strip()]
    if len(remaining) % 3 != 0:
        return None
    k = len(remaining) // 3
    smiles_block    = remaining[:k]
    names_block     = remaining[k:2 * k]
    enthalpies_raw  = remaining[2 * k:]
    enthalpies = []
    for e in enthalpies_raw:
        e = e.strip()
        if e == "No_Thermo":
            enthalpies.append(None)
        else:
            try:
                enthalpies.append(float(e))
            except ValueError:
                enthalpies.append(None)

    return RankedPathway(
        rank=rank,
        final_score=final_score if final_score is not None else 0.0,
        atomic_economy=atomic_economy if atomic_economy is not None else 0.0,
        pathway_byproduct_count=pathway_byproduct_count or 0,
        intermediate_byproducts=intermediate_byproducts,
        reaction_smiles=smiles_block,
        reaction_names=names_block,
        reaction_enthalpies=enthalpies,
    )


# ======================================================================
# "Lemnisca" blend — custom criteria layered ON TOP of DORAnet's score
# ======================================================================
#
# We do NOT reimplement DORAnet's scoring. DORAnet's whole composite
# (steps, thermo, atom economy, by-products — with the user's weights)
# enters the blend as ONE ingredient via DoranetScoreCriterion. On top
# of it we add only the criteria DORAnet lacks: intermediate stability
# and procedure diversity (cost / toxicity / DORA-XGB feasibility come
# later). apply_lemnisca_blend() combines them into `lemnisca_score` and
# re-orders the pathways.
#
# Convention (same as DORAnet): each criterion returns a list[float],
# one score per pathway in [0, 1], HIGHER = BETTER.


class PathwayCriterion(ABC):
    """A scoring criterion over a batch of RankedPathway objects.

    `floor` is where this criterion's grade bottoms out when blended into
    the geometric mean: 0.5 means "discounts but never gates"; 0.0 means
    "a zero here gates the whole pathway to 0" (use only for genuine
    dealbreakers, e.g. a catastrophic intermediate)."""

    name: str = "criterion"
    floor: float = 0.5

    @abstractmethod
    def score(self, pathways: list) -> list:
        ...


class DoranetScoreCriterion(PathwayCriterion):
    """DORAnet's WHOLE composite score as a single criterion, min-max
    normalized across the batch to [0, 1] (higher = better).

    This is the seam that keeps DORAnet's well-tuned, weight-adjustable
    score intact: its internal component weights are set upstream (in
    generate_base_rankings); here we just take the resulting final_score
    and rescale it so it can be blended with the custom [0,1] criteria.
    """

    name = "doranet"

    def score(self, pathways: list) -> list:
        vals = [
            (p.final_score if p.final_score is not None else 0.0)
            for p in pathways
        ]
        if not vals:
            return []
        lo, hi = min(vals), max(vals)
        if hi == lo:
            # All pathways scored identically — no discrimination here.
            return [1.0 for _ in vals]
        return [(v - lo) / (hi - lo) for v in vals]


class IntermediateStabilityCriterion(PathwayCriterion):
    """Penalize pathways whose INTERMEDIATES carry structural features
    that make them unstable / unsafe / hard to isolate at the bench.

    Hazard tiers (derivable from SMILES alone, no extra data needed):
      TIER A — catastrophic, one match zeroes the intermediate:
        peroxide, organic azide, diazo, diazonium, gem-dinitro.
      TIER B — major (0.3 each): acyl halide, radical center, MW > 800.
      TIER C — minor (0.1 each): net formal charge, 3-membered ring,
        gem-diol, 500 < MW <= 800, nitro group.

    Per-intermediate score:  max(0, 1 - 1.0*A - 0.3*B - 0.1*C)
    Per-pathway score:       MIN over intermediates (weakest-link).

    Absolute, NOT batch-normalized: 0.7 means "70% viable" on its own.
    `excluded_smiles` (starter, target, helpers) are skipped.
    """

    name = "stability"
    floor = 0.0   # a catastrophic intermediate (score 0) GATES the pathway

    _TIER_A_SMARTS = {
        "peroxide":     "[OX2]-[OX2]",
        "azide_a":      "[N-]=[N+]=N",
        "azide_b":      "N=[N+]=[N-]",
        "diazo":        "[CX3]=[N+]=[N-]",
        "diazonium":    "[#6][NX1]#[NX2+]",
        "gem_dinitro":  "[CX4]([NX3](=O)=O)[NX3](=O)=O",
    }
    _TIER_B_SMARTS = {
        "acyl_halide":  "[CX3](=O)[F,Cl,Br,I]",
    }
    _TIER_C_SMARTS = {
        "gem_diol":     "[CX4]([OX2H])[OX2H]",
        "nitro":        "[NX3](=O)=O",
    }
    _MW_TIER_B = 800.0
    _MW_TIER_C = 500.0

    def __init__(self, excluded_smiles=None):
        self._excluded = {
            c for s in (excluded_smiles or [])
            if (c := self._canonical(s)) is not None
        }
        self._compiled_a = {n: Chem.MolFromSmarts(s) for n, s in self._TIER_A_SMARTS.items()}
        self._compiled_b = {n: Chem.MolFromSmarts(s) for n, s in self._TIER_B_SMARTS.items()}
        self._compiled_c = {n: Chem.MolFromSmarts(s) for n, s in self._TIER_C_SMARTS.items()}
        self._score_cache: dict = {}

    @staticmethod
    def _canonical(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol is not None else None

    def _intermediates_in(self, p) -> set:
        """Unique molecules across the pathway's reactions, minus excluded.
        RankedPathway stores each step as 'reactants>>products'."""
        seen: set = set()
        for rxn in p.reaction_smiles:
            lhs, _, rhs = rxn.partition(">>")
            for s in lhs.split(".") + rhs.split("."):
                s = s.strip()
                if s and s not in self._excluded:
                    seen.add(s)
        return seen

    def _score_intermediate(self, smiles: str) -> float:
        cached = self._score_cache.get(smiles)
        if cached is not None:
            return cached
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            self._score_cache[smiles] = 0.0
            return 0.0
        tier_a = sum(len(mol.GetSubstructMatches(pat)) for pat in self._compiled_a.values())
        tier_b = sum(len(mol.GetSubstructMatches(pat)) for pat in self._compiled_b.values())
        tier_c = sum(len(mol.GetSubstructMatches(pat)) for pat in self._compiled_c.values())
        tier_b += sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        if Chem.GetFormalCharge(mol) != 0:
            tier_c += 1
        mw = Descriptors.MolWt(mol)
        if mw > self._MW_TIER_B:
            tier_b += 1
        elif mw > self._MW_TIER_C:
            tier_c += 1
        tier_c += sum(1 for r in mol.GetRingInfo().AtomRings() if len(r) == 3)
        s = max(0.0, 1.0 - 1.0 * tier_a - 0.3 * tier_b - 0.1 * tier_c)
        self._score_cache[smiles] = s
        return s

    def score(self, pathways: list) -> list:
        out = []
        for p in pathways:
            inters = self._intermediates_in(p)
            out.append(min((self._score_intermediate(s) for s in inters), default=1.0))
        return out


class ProcedureDiversityCriterion(PathwayCriterion):
    """Penalize routes that need many DISTINCT procedures (operator
    names). Each distinct procedure needs its own condition development,
    safety review, and analytical method, so a route repeating one
    reaction is cheaper to scale than one using N different reactions at
    the same step count.

        s = 1 - (distinct_names - 1) / max(1, total_steps - 1)

    All-same procedure → 1.0; all-different → 0.0; single step → 1.0.
    Absolute, not batch-normalized.
    """

    name = "diversity"

    def score(self, pathways: list) -> list:
        out = []
        for p in pathways:
            n = p.num_steps
            if n <= 1:
                out.append(1.0)
                continue
            distinct = len(set(p.reaction_names))
            out.append(1.0 - (distinct - 1) / (n - 1))
        return out


class FeasibilityCriterion(PathwayCriterion):
    """DORA-XGB enzymatic reaction-feasibility.

    DORA-XGB is trained on ENZYMATIC reactions, so it's only meaningful
    for bio steps — chem steps are treated as feasible (1.0) so they
    don't get penalised by an out-of-domain prediction. Each bio step is
    scored 0-1 (higher = more feasible) via the DoraXGBClient; the
    pathway takes its WEAKEST bio step (like stability's weakest-link).

    Soft (floor 0.5): it's a model prediction, so it discounts rather
    than hard-gates. Only active when a dora_client is supplied.
    """

    name = "feasibility"
    floor = 0.5

    def __init__(self, dora_client):
        self._client = dora_client

    def score(self, pathways: list) -> list:
        out = []
        for p in pathways:
            vals = []
            for smi, op_name in zip(p.reaction_smiles, p.reaction_names):
                if is_bio_op(op_name):
                    s = self._client.feasibility(smi)
                    if s is not None:
                        vals.append(s)
            out.append(min(vals) if vals else 1.0)  # chem-only → no penalty
        return out


class EnzymeLoadCriterion(PathwayCriterion):
    """Metabolic-engineering burden: the MINIMUM number of distinct enzymes
    you'd need to express to build this route in a host.

    The multiple UniProt hits per rule are ALTERNATIVE catalysts for that
    one step (same reaction, different organisms), so a step needs just one
    of them — not all. Better still: a multifunctional enzyme that appears
    in two different steps' candidate lists can catalyze BOTH, so the true
    minimum is a set-cover over the steps, which can be fewer than the
    number of bio steps. Chem steps and no-known-enzyme steps don't count.
    (This is an optimistic lower bound — the operators are somewhat
    promiscuous, so a shared enzyme *could* cover both steps, not a
    guarantee that it does.)

    Absolute (not batch-normalized): raw = DECAY ** min_enzymes, so each
    additional required enzyme discounts multiplicatively but never zeroes
    the score. Soft (floor 0.5): a preference, not a dealbreaker.
    """

    name = "enzyme_load"
    floor = 0.5
    _DECAY = 0.85

    def score(self, pathways: list) -> list:
        # lazy import avoids a module cycle (enzyme_info imports is_bio_op)
        from enzyme_info import minimum_enzyme_count
        return [self._DECAY ** minimum_enzyme_count(p.reaction_names)
                for p in pathways]


# Tier-2 (layer) default weights: DORAnet's chemistry score vs. the
# Lemnisca process-viability score. DORAnet counts double by default.
LAYER_DEFAULT_WEIGHTS = {
    "doranet":  2.0,
    "lemnisca": 1.0,
}
# Tier-1 (component) default weights inside the Lemnisca sub-score.
# (DORA-XGB feasibility is applied as a generation-phase PRUNE, not a
# ranking component — see run_pipeline's feasibility_prune_threshold.)
LEMNISCA_DEFAULT_WEIGHTS = {
    "stability":   1.0,
    "diversity":   1.0,
    "enzyme_load": 1.0,
}


def _floored(x: float, floor: float) -> float:
    """Remap a raw [0,1] grade into [floor, 1] linearly. floor=0 leaves
    it unchanged (so a 0 survives and can gate); floor=0.5 lifts the worst
    case to 0.5 so the criterion discounts but never gates."""
    return floor + (1.0 - floor) * x


def _weighted_geomean(pairs: list) -> Optional[float]:
    """Weighted geometric mean of (value, weight) pairs:
        (∏ vᵢ^wᵢ)^(1/Σwᵢ)
    Any value == 0 with weight > 0 makes the whole result 0 (the gate).
    Returns None if no pair has positive weight."""
    active = [(v, w) for v, w in pairs if w > 0]
    if not active:
        return None
    if any(v <= 0.0 for v, _ in active):
        return 0.0
    wsum = sum(w for _, w in active)
    log_sum = sum(w * math.log(v) for v, w in active)
    return math.exp(log_sum / wsum)


def apply_lemnisca_blend(
    ranked_pathways: list,
    layer_weights: Optional[dict] = None,
    lemnisca_weights: Optional[dict] = None,
    excluded_smiles=None,
) -> list:
    """
    Two-tier weighted GEOMETRIC-mean blend, then re-sort + re-rank.

    Tier 1 — Lemnisca sub-score: geometric mean of the custom criteria,
        each remapped into [floor, 1]. Stability (floor 0) can drive the
        sub-score to 0 → the gate; diversity (floor 0.5) only discounts.

            lemnisca = ( ∏ floored(cᵢ)^wᵢ )^(1/Σwᵢ)      (lemnisca_weights)

    Tier 2 — final blended score: geometric mean of DORAnet (remapped
        into [0.5, 1] so it never gates) and the Lemnisca sub-score
        (NOT floored, so a stability gate of 0 propagates through):

            blended = ( doranet^a · lemnisca^b )^(1/(a+b))   (layer_weights)

    Stores raw 0–1 grades on `lemnisca_components` (incl. "doranet"), the
    sub-score on `lemnisca_score`, and the ranking key on `blended_score`.
    `excluded_smiles` (starter, target, helpers) skip stability scoring.
    """
    if not ranked_pathways:
        return ranked_pathways
    layer_weights = layer_weights or LAYER_DEFAULT_WEIGHTS
    lemnisca_weights = lemnisca_weights or LEMNISCA_DEFAULT_WEIGHTS

    doranet_c = DoranetScoreCriterion()
    lem_criteria = [
        IntermediateStabilityCriterion(excluded_smiles),
        ProcedureDiversityCriterion(),
        EnzymeLoadCriterion(),
    ]
    doranet_grades = doranet_c.score(ranked_pathways)
    lem_grades = {c.name: c.score(ranked_pathways) for c in lem_criteria}

    for i, p in enumerate(ranked_pathways):
        comps: dict = {}
        # --- Tier 1: Lemnisca sub-score ---
        lem_pairs = []
        for c in lem_criteria:
            raw = lem_grades[c.name][i]
            comps[c.name] = raw
            w = float(lemnisca_weights.get(c.name, 0.0))
            lem_pairs.append((_floored(raw, c.floor), w))
        lemnisca_sub = _weighted_geomean(lem_pairs)
        if lemnisca_sub is None:        # no lemnisca weights → neutral
            lemnisca_sub = 1.0

        # --- Tier 2: final blend ---
        dora_raw = doranet_grades[i]
        comps["doranet"] = dora_raw
        final_pairs = [
            (_floored(dora_raw, doranet_c.floor),
             float(layer_weights.get("doranet", 0.0))),
            (lemnisca_sub,              # NOT floored: lets the gate through
             float(layer_weights.get("lemnisca", 0.0))),
        ]
        blended = _weighted_geomean(final_pairs)
        if blended is None:            # both layer weights 0 → fall back
            blended = _floored(dora_raw, doranet_c.floor)

        p.lemnisca_score = lemnisca_sub
        p.blended_score = blended
        p.lemnisca_components = comps

    reordered = sorted(
        ranked_pathways,
        key=lambda p: (p.blended_score if p.blended_score is not None else 0.0),
        reverse=True,
    )
    for new_rank, p in enumerate(reordered, 1):
        p.rank = new_rank
    return reordered
