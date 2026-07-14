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
from dataclasses import dataclass, field
from typing import Optional

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
