"""
config.py — central weight table for pathway scoring.

Every criterion's weight lives here so a single change rebalances
both the markdown ranking AND the PDF page order. Callers can pull
DEFAULT_WEIGHTS directly or use `build_unified_profile()` which
attaches the DORAnet-derived criteria when their precomputed values
are available.

Profiles can capture different employee/team priorities — a process
engineer caring about manufacturing simplicity, a green-chemistry
reviewer caring about atom economy, etc. Add new ones by appending
another constant below.
"""

from pathway_scoring import (
    StepsCriterion,
    ThermoCriterion,
    IntermediateStabilityCriterion,
    ProcedureDiversityCriterion,
    ChemBioSwitchCriterion,
    AtomEconomyCriterion,
    ByProductCountCriterion,
)


# ----------------------------------------------------------------------
# Central weight table — change numbers here and rerun.
# ----------------------------------------------------------------------
DEFAULT_WEIGHTS = {
    # Custom criteria
    "steps":              4.0,
    "thermo":             2.0,
    "stability":          2.0,
    "diversity":          2.0,
    "chem_bio":           2.0,
    # DORAnet-derived criteria (only attached when their precomputed
    # values are available — see build_unified_profile)
    "atom_economy":       1.0,
    "by_product":         2.0,
}


def build_unified_profile(
    atom_economy_by_index: dict | None = None,
    by_product_by_index: dict | None = None,
    weights: dict | None = None,
) -> list[tuple]:
    """
    Build a (criterion, weight) list for WeightedPathwayScorer that
    blends the custom criteria with DORAnet's atom-economy and
    by-product criteria. The DORAnet criteria are only included when
    their precomputed value dicts are passed in (typically from
    parse_doranet_ranked_file).
    """
    w = weights or DEFAULT_WEIGHTS
    profile = [
        (StepsCriterion(),                    w["steps"]),
        (ThermoCriterion(),                   w["thermo"]),
        (IntermediateStabilityCriterion(),    w["stability"]),
        (ProcedureDiversityCriterion(),       w["diversity"]),
        (ChemBioSwitchCriterion(),            w["chem_bio"]),
    ]
    if atom_economy_by_index is not None:
        profile.append(
            (AtomEconomyCriterion(atom_economy_by_index), w["atom_economy"])
        )
    if by_product_by_index is not None:
        profile.append(
            (ByProductCountCriterion(by_product_by_index), w["by_product"])
        )
    return profile


# Backwards-compat: a simple two-criterion profile (matches DORAnet's
# original weight emphasis: steps 4, thermo 2). Use build_unified_profile
# in new code.
DEFAULT_PROFILE = [
    (StepsCriterion(),  DEFAULT_WEIGHTS["steps"]),
    (ThermoCriterion(), DEFAULT_WEIGHTS["thermo"]),
]
