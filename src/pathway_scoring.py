"""
pathway_scoring.py
Pathway-level scoring framework.

WHAT THIS IS (vs. recipe_rankers.py)
------------------------------------
`recipe_rankers.py` scores INDIVIDUAL RECIPES during network expansion
to steer the search — it picks which (operator, reactants) to expand
next. This file scores WHOLE PATHWAYS after `pathway_finder` has
enumerated candidate routes to a target. The two layers don't compete;
they work on different objects at different stages of the pipeline:

    network expansion        →  recipe ranker (per-recipe, soft order)
    candidate pathways list  →  pathway scorer (per-path, soft rank)

DESIGN PHILOSOPHY: COMPOSABILITY
--------------------------------
Different employees will weight pathway criteria differently — a
process engineer caring about manufacturing simplicity, a green-
chemistry reviewer caring about atom economy, etc. So each criterion
lives in its own small class, and `WeightedPathwayScorer` combines
any subset with user-chosen weights. Don't hard-code one objective.

WHY CRITERIA SCORE A BATCH, NOT ONE PATH AT A TIME
--------------------------------------------------
A pathway's raw signal (e.g. max ΔH across its reactions) is only
meaningful relative to the other candidates in the same job. DORAnet's
`pathway_ranking` normalizes each criterion to [0, 1] using the
min/max within the batch (post_processing.py:1569-1584, :1798-1806).
We follow the same convention: a criterion takes `list[Pathway]` and
returns `list[float]` aligned 1:1, each in [0, 1] with higher = better.

This keeps every criterion's contribution comparable when combined.

CURRENT SLICE
-------------
This first slice ships the two simplest criteria — `StepsCriterion`
and `ThermoCriterion` — so the pipeline can run end-to-end before we
add atom economy (which needs DORAnet's stoichiometric path-unrolling
math) or cost / sustainability / manufacturing (which need catalogs).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Descriptors

from pathway_tools import Pathway, parse_reaction_string


# =====================================================================
# Base class
# =====================================================================

class PathwayCriterion(ABC):
    """
    A scoring criterion over a batch of pathways.

    Convention (verified against DORAnet's pathway_ranking):
        - input:  list[Pathway]
        - output: list[float], one score per pathway, same order
        - each score lies in [0, 1]; HIGHER = BETTER
    """

    name: str = "criterion"

    @abstractmethod
    def score(self, pathways: list[Pathway]) -> list[float]:
        ...


# =====================================================================
# Criterion: number of reaction steps
# =====================================================================

class StepsCriterion(PathwayCriterion):
    """
    Fewer reaction steps = higher score. Normalized linearly so the
    shortest pathway in the batch gets 1.0 and the longest gets 0.0.

    Mirrors DORAnet's number-of-steps block (post_processing.py:
    1793-1806). Step count is the single strongest empirical predictor
    of process feasibility in DORAnet's default weights (weight=4 by
    default; everything else is 0-2).
    """

    name = "steps"

    def score(self, pathways: list[Pathway]) -> list[float]:
        if not pathways:
            return []
        steps = [p.num_steps for p in pathways]
        min_step, max_step = min(steps), max(steps)
        # DORAnet uses diff = -(max - min) so that (i - min)/diff + 1
        # gives shortest→1, longest→0. If all paths share a length the
        # divisor would be 0 — fall back to a tiny number so every
        # pathway scores ~1.0 (they're equally good on this axis).
        diff = -(max_step - min_step)
        if diff == 0:
            return [1.0] * len(pathways)
        return [(s - min_step) / diff + 1.0 for s in steps]


# =====================================================================
# Criterion: thermodynamic favorability
# =====================================================================

class ThermoCriterion(PathwayCriterion):
    """
    Lower max-ΔH-along-the-pathway = higher score. Normalized linearly
    so the most exothermic pathway in the batch gets 1.0 and the most
    endothermic gets 0.0.

    "Max" is the right aggregator because a pathway is gated by its
    WORST step: a single highly endothermic reaction makes the route
    impractical even if every other step is favorable. Same choice
    DORAnet makes (post_processing.py:1556-1567).

    Missing data handling: dH values arrive as strings in the reaction
    record, with the literal "No_Thermo" standing in for "calculator
    couldn't price this molecule." Following DORAnet's convention
    (post_processing.py:1554-1567): individual No_Thermo reactions are
    skipped when aggregating, but the pathway's max is still computed
    over its priced reactions. A pathway is assigned score 0 only when
    EVERY reaction is No_Thermo — in that case we genuinely have nothing
    to compare on.

    Chemistry note: max-over-available is a lower bound on the actual
    worst-step ΔH. If a real bottleneck happens to be the one unknown
    step, the bound is loose — but ignoring all priced data points
    because one is missing is strictly worse.
    """

    name = "thermo"

    def score(self, pathways: list[Pathway]) -> list[float]:
        if not pathways:
            return []

        max_H_list: list[float | None] = []
        for path in pathways:
            path_max = None
            for rxn in path.reactions:
                dH = parse_reaction_string(rxn)["dH"]
                if dH is None:
                    continue   # skip unknown steps, don't disqualify the path
                path_max = dH if path_max is None else max(path_max, dH)
            # path_max stays None only if EVERY reaction was No_Thermo
            max_H_list.append(path_max)

        numeric = [v for v in max_H_list if v is not None]
        if not numeric:
            # No pathway has complete thermo — nothing to compare.
            return [0.0] * len(pathways)

        min_H, max_H = min(numeric), max(numeric)
        diff = -(max_H - min_H)
        # If every priced pathway has the same max-dH, every score
        # collapses to 1.0; use a small divisor as a tiebreaker so the
        # math still runs (DORAnet does the same).
        if diff == 0:
            return [1.0 if v is not None else 0.0 for v in max_H_list]

        return [
            (v - min_H) / diff + 1.0 if v is not None else 0.0
            for v in max_H_list
        ]


# =====================================================================
# Criterion: intermediate stability / lab-handling safety
# =====================================================================

class IntermediateStabilityCriterion(PathwayCriterion):
    """
    Penalize pathways whose INTERMEDIATES carry structural features that
    make them unstable, unsafe, or hard to isolate at the bench.

    Why this matters (more than thermo, for slice-2 process work)
    ------------------------------------------------------------
    Process chemists kill more routes for safety/isolability than for
    thermodynamics. A 2-step pathway through a peroxide intermediate
    is *worse* than a 5-step pathway through stable solids — the short
    route can't be run at all. ΔH being favorable doesn't help if you
    can't put the intermediate in a flask overnight.

    Why intermediates, not reactions
    --------------------------------
    Reaction safety (conditions, exotherm risk, reagent toxicity) needs
    data DORAnet operators don't carry — temperature, solvent, catalyst.
    Intermediate stability is fully derivable from SMILES alone
    (substructure matches + RDKit atom properties + MW), so it ships now
    with no new infrastructure. Reaction-level hazard scoring belongs
    in a future `ReactionHazardCriterion` keyed on operator-name tags.

    Hazard tiers
    ------------
    TIER A — catastrophic (one match alone drops the intermediate to 0):
        • peroxide / hydroperoxide / acyl peroxide  (R-O-O-R)
        • organic azide   (R-N=N+=N-)
        • diazo (R2C=N+=N-)
        • diazonium (Ar-N#N+)
        • gem-dinitro (one carbon bearing two NO2 groups)
        Rationale: each is on industry "do not handle" lists. A real
        process plant won't permit them at scale.

    TIER B — major (0.3 score loss per match):
        • acyl halide                     (very moisture-sensitive)
        • free radical center             (won't survive isolation)
        • molecular weight > 800 Da       (purification headache)
        Rationale: workable in research with care, but adds enough
        process risk to bias against the route.

    TIER C — minor (0.1 score loss per match):
        • non-zwitterionic formal charge  (likely an ion / unstable carbocation)
        • 3-membered ring                 (epoxides/aziridines fine, but flagged)
        • gem-diol                        (equilibrium with carbonyl, hard to isolate)
        • 500 < MW <= 800 Da              (heavier intermediates)
        • single nitro group              (mild — many nitro compounds are stable)
        Rationale: yellow flags worth surfacing in the score but not
        decisive on their own.

    Score per intermediate:
        s = max(0, 1 - 1.0·A_count - 0.3·B_count - 0.1·C_count)

    Score per pathway:
        s = MIN(s over intermediates)
        — weakest-link convention. A pathway is as safe as its worst
        intermediate, mirroring the way max-ΔH gates thermo.

    Absolute, not batch-normalized
    ------------------------------
    Unlike Steps/Thermo, this score is meaningful in absolute terms:
    0.7 = "70% viable" without reference to the rest of the batch.
    So we DON'T normalize by min/max across pathways — a batch of
    uniformly hazardous routes shouldn't have one re-rescaled to 1.0
    just because it's "least bad."

    Parameters
    ----------
    excluded_smiles : set[str] | None
        SMILES of molecules to skip when collecting intermediates —
        typically your starter, helpers, and target (you've already
        decided to make/use them; their stability isn't in question).
        SMILES are canonicalized internally before comparison.
        If None, every molecule in the pathway is scored.
    """

    name = "stability"

    # SMARTS patterns are compiled once on first use and shared across
    # instances. Keeping them at class level documents the hazard list
    # and lets users override by subclassing if they want a different
    # set for a different chemistry domain.
    _TIER_A_SMARTS = {
        "peroxide":        "[OX2]-[OX2]",
        "azide_a":         "[N-]=[N+]=N",
        "azide_b":         "N=[N+]=[N-]",
        "diazo":           "[CX3]=[N+]=[N-]",
        "diazonium":       "[#6][NX1]#[NX2+]",
        "gem_dinitro":     "[CX4]([NX3](=O)=O)[NX3](=O)=O",
    }
    _TIER_B_SMARTS = {
        "acyl_halide":     "[CX3](=O)[F,Cl,Br,I]",
    }
    _TIER_C_SMARTS = {
        "gem_diol":        "[CX4]([OX2H])[OX2H]",
        "nitro":           "[NX3](=O)=O",
    }

    _MW_TIER_B = 800.0   # > 800 Da → tier B
    _MW_TIER_C = 500.0   # 500–800 → tier C

    def __init__(self, excluded_smiles=None):
        self._excluded = {
            self._canonical(s) for s in (excluded_smiles or [])
            if self._canonical(s) is not None
        }
        self._compiled_a = {n: Chem.MolFromSmarts(s) for n, s in self._TIER_A_SMARTS.items()}
        self._compiled_b = {n: Chem.MolFromSmarts(s) for n, s in self._TIER_B_SMARTS.items()}
        self._compiled_c = {n: Chem.MolFromSmarts(s) for n, s in self._TIER_C_SMARTS.items()}
        self._score_cache: dict[str, float] = {}

    @staticmethod
    def _canonical(smiles):
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol is not None else None

    def _intermediates_in(self, path: Pathway) -> set[str]:
        """All unique molecules in the pathway, minus excluded ones."""
        seen: set[str] = set()
        for rxn in path.reactions:
            parsed = parse_reaction_string(rxn)
            for s in parsed["reactants"] + parsed["products"]:
                if s not in self._excluded:
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

        # Substructure hits per tier
        tier_a = sum(
            len(mol.GetSubstructMatches(p)) for p in self._compiled_a.values()
        )
        tier_b = sum(
            len(mol.GetSubstructMatches(p)) for p in self._compiled_b.values()
        )
        tier_c = sum(
            len(mol.GetSubstructMatches(p)) for p in self._compiled_c.values()
        )

        # Radical centers — RDKit property, not a SMARTS match
        radicals = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        tier_b += radicals

        # Formal charge — only flag if molecule is a net-charged ion (not
        # a drawn-zwitterion like sulfobetaine, whose net charge is 0).
        net_charge = Chem.GetFormalCharge(mol)
        if net_charge != 0:
            tier_c += 1

        # Molecular weight bands
        mw = Descriptors.MolWt(mol)
        if mw > self._MW_TIER_B:
            tier_b += 1
        elif mw > self._MW_TIER_C:
            tier_c += 1

        # 3-membered rings — tier C (epoxides/aziridines are common but
        # 3MRs are worth surfacing as a yellow flag in general)
        small_rings = sum(1 for r in mol.GetRingInfo().AtomRings() if len(r) == 3)
        tier_c += small_rings

        s = max(0.0, 1.0 - 1.0 * tier_a - 0.3 * tier_b - 0.1 * tier_c)
        self._score_cache[smiles] = s
        return s

    def score(self, pathways: list[Pathway]) -> list[float]:
        out: list[float] = []
        for path in pathways:
            inters = self._intermediates_in(path)
            if not inters:
                # No assessable intermediates → pathway is just starter
                # → helpers → target with nothing to flag. Score 1.
                out.append(1.0)
                continue
            out.append(min(self._score_intermediate(s) for s in inters))
        return out


# =====================================================================
# Criterion: procedure (operator-name) diversity
# =====================================================================

class ProcedureDiversityCriterion(PathwayCriterion):
    """
    Penalize pathways that require many DISTINCT procedures.

    What we count
    -------------
    Distinct operator NAMES (not "reaction types" — DORAnet's
    `reaction_type` field is uniformly "Catalytic" across all 388
    forward operators, which is too coarse to drive diversity scoring).
    Operator names like "Aldol Condensation" or "Hydrogenolysis of
    Ethers" each correspond to a distinct DORAnet SMARTS recipe with
    its own conditions in practice.

    Why this matters chemically / process-wise
    ------------------------------------------
    Every distinct procedure used in a route requires its own:
        - reaction-condition development (T, solvent, catalyst loading)
        - safety review
        - in-process analytical method (HPLC method, NMR assignment)
        - workup protocol
    A route running three different alkylations is genuinely cheaper to
    take to scale than a route running an alkylation, a hydrogenation
    and an oxidation — even at the same step count. This is standard
    process-chemistry teaching: Anderson, *Practical Process Research
    and Development*, 2nd ed., Academic Press 2012, Ch. 2 ("Reactions
    you've never run before are slower and more expensive than reactions
    you've run a hundred times").

    Score
    -----
        s = 1 - (distinct_names - 1) / max(1, total_reactions - 1)

        - All steps use the same procedure: distinct = 1 → s = 1.0
        - Every step uses a different procedure: distinct = N → s = 0.0
        - Single-step pathway:                                   s = 1.0

    Absolute, not batch-normalized — the meaning ("fraction of step
    budget spent on procedure repeats") is the same regardless of the
    other pathways in the batch.

    What this DOESN'T capture
    -------------------------
    Two mechanistically-similar operators with different names (say,
    two aldol variants) count as distinct procedures. That's arguably
    correct for development-cost purposes — even mechanistic siblings
    typically need their own optimization — but a future criterion
    could group by mechanism family if we ever get a research-grounded
    taxonomy. DORAnet's current `reaction_type` field doesn't provide
    one.
    """

    name = "diversity"

    def score(self, pathways: list[Pathway]) -> list[float]:
        out: list[float] = []
        for path in pathways:
            n = len(path.reactions)
            if n <= 1:
                out.append(1.0)
                continue
            distinct = len({parse_reaction_string(r)["op_name"] for r in path.reactions})
            out.append(1.0 - (distinct - 1) / (n - 1))
        return out


# =====================================================================
# Criterion: chemical / biological regime mixing
# =====================================================================

class ChemBioSwitchCriterion(PathwayCriterion):
    """
    Penalize pathways that mix chemical and biological reactions.

    Why this matters at manufacturing scale
    ---------------------------------------
    A chemical step (organic solvent, often elevated T/P, metal/acid/
    base catalyst) and a biological step (aqueous buffer, ambient T,
    enzyme or cell catalyst) live on fundamentally different process
    platforms. Switching between them mid-route forces:
        - a full solvent swap and workup between regimes
        - separate sets of process equipment (reactor materials, GMP
          tier, containment)
        - independent QA / QC for each regime
        - typically a hold/transfer point with associated yield loss
    All-chem or all-bio routes inherit one process platform; mixed
    routes inherit two and pay the integration cost. Biocatalysis
    review literature is consistent on this point — e.g. Sheldon &
    Woodley, *Chem. Rev.* 118, 801 (2018), on the practicalities of
    integrating biocatalysis into chemical processes.

    Identifying bio reactions
    -------------------------
    We use DORAnet's curated `bio_rxn_names` set (loaded from
    `JN3604IMT_rules.tsv` — a published enzymatic reaction-rule library
    used in metabolic-route prediction). A reaction is "bio" iff its
    operator name is in that set; everything else is "chem". This
    avoids inventing our own classification — the labels come straight
    from the rule file DORAnet itself uses everywhere else (atom
    economy, by-product counting, cofactor handling).

    Score (regime purity)
    ---------------------
        s = max(n_chem, n_bio) / total_reactions

        - All-chem or all-bio:   s = 1.0
        - 80 / 20 split:         s = 0.8
        - 50 / 50 split:         s = 0.5
        - One bio + many chem:   s ≈ 1.0 − 1/N  (only mildly penalized)

    Why this aggregator (and not "count of regime switches"):
    pathway files store reactions as a SET, not a topologically-sorted
    sequence — so "number of switches in execution order" is undefined
    without first topologically sorting the pathway. Regime purity
    captures the same manufacturing intuition (mixing is bad,
    regardless of order) without requiring a sort.

    Current TAL caveat
    ------------------
    Until enzymatic operators are added to the TAL whitelist, every
    pathway is all-chem and this criterion returns 1.0 for everyone.
    It activates the moment the first bio operator gets pulled into
    the network.
    """

    name = "chem_bio"

    def __init__(self):
        # Defer the DORAnet import so users without enzymatic data
        # installed still get a clean import; only fail loudly when
        # the criterion is actually exercised.
        from doranet.modules.post_processing.post_processing import bio_rxn_names
        self._bio_names = bio_rxn_names

    def score(self, pathways: list[Pathway]) -> list[float]:
        out: list[float] = []
        for path in pathways:
            n = len(path.reactions)
            if n == 0:
                out.append(1.0)
                continue
            n_bio = sum(
                1 for r in path.reactions
                if parse_reaction_string(r)["op_name"] in self._bio_names
            )
            n_chem = n - n_bio
            out.append(max(n_chem, n_bio) / n)
        return out


# =====================================================================
# Composite scorer
# =====================================================================

@dataclass
class ScoredPathway:
    pathway: Pathway
    final_score: float
    components: dict[str, float]  # criterion name → raw [0,1] score


class WeightedPathwayScorer:
    """
    Combine any number of `PathwayCriterion`s into a single ranking via
    a user-defined weighted sum:

        final_score = Σ(weight_i · score_i) / Σ(weight_i)

    This is the seam where employee preference enters pathway ranking.
    Build the scorer at the call site with whatever weights make sense
    for the current task — e.g. a process engineer might weight steps
    high and atom economy low; a green-chemistry reviewer the opposite.

    Parameters
    ----------
    components : list[tuple[PathwayCriterion, float]]
        (criterion, weight). Weights need not sum to 1 — the scorer
        normalizes by the weight total.
    """

    def __init__(self, components: list[tuple[PathwayCriterion, float]]):
        if not components:
            raise ValueError("Need at least one (criterion, weight) component")
        weight_sum = sum(w for _, w in components)
        if weight_sum <= 0:
            raise ValueError("Sum of weights must be positive")
        self.components = tuple(components)
        self._weight_sum = float(weight_sum)

    def score(self, pathways: list[Pathway]) -> list[ScoredPathway]:
        if not pathways:
            return []

        per_criterion_scores: list[list[float]] = []
        for criterion, _ in self.components:
            per_criterion_scores.append(criterion.score(pathways))

        results: list[ScoredPathway] = []
        for i, path in enumerate(pathways):
            total = 0.0
            comp_dict: dict[str, float] = {}
            for (criterion, weight), scores in zip(self.components, per_criterion_scores):
                s = scores[i]
                comp_dict[criterion.name] = s
                total += weight * s
            results.append(
                ScoredPathway(
                    pathway=path,
                    final_score=total / self._weight_sum,
                    components=comp_dict,
                )
            )
        return results


# =====================================================================
# Entry point
# =====================================================================

def score_pathways_from_file(
    job_name: str,
    scorer: WeightedPathwayScorer,
) -> list[ScoredPathway]:
    """
    Load pathways from `{job_name}_pathways.txt`, run them through
    `scorer`, and return them sorted best-first.
    """
    from pathway_tools import load_pathways_from_file
    pathways = load_pathways_from_file(job_name)
    scored = scorer.score(pathways)
    scored.sort(key=lambda sp: sp.final_score, reverse=True)
    return scored
