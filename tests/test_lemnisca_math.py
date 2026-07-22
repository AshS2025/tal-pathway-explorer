"""
Tests for the Lemnisca ranking blend — the product's custom IP, and until
now completely untested. Every expected number here is hand-computed from
the formulas in pathway_scoring.py so a silent regression in a floor, a
weight, or the geometric-mean gate gets caught immediately.
"""
import math
import pytest

from pathway_scoring import (
    _floored, _weighted_geomean,
    DoranetScoreCriterion, ProcedureDiversityCriterion,
    IntermediateStabilityCriterion, FeasibilityCriterion, EnzymeLoadCriterion,
    apply_lemnisca_blend, RankedPathway,
)


def mk(final_score, smiles, names):
    """Build a RankedPathway. num_steps == len(smiles)."""
    return RankedPathway(
        rank=0, final_score=final_score, atomic_economy=0.0,
        pathway_byproduct_count=0, intermediate_byproducts={},
        reaction_smiles=list(smiles), reaction_names=list(names),
        reaction_enthalpies=["No_Thermo"] * len(smiles),
    )


# ---------------------------------------------------------------- _floored
def test_floored_floor_zero_is_identity():
    assert _floored(0.0, 0.0) == 0.0
    assert _floored(1.0, 0.0) == 1.0
    assert _floored(0.4, 0.0) == pytest.approx(0.4)


def test_floored_lifts_worst_case_to_floor():
    # floor=0.5 remaps [0,1] -> [0.5,1]
    assert _floored(0.0, 0.5) == 0.5
    assert _floored(1.0, 0.5) == 1.0
    assert _floored(0.5, 0.5) == 0.75


# -------------------------------------------------------- _weighted_geomean
def test_weighted_geomean_equal_weights():
    # sqrt(0.8 * 0.5) = sqrt(0.4)
    assert _weighted_geomean([(0.8, 1), (0.5, 1)]) == pytest.approx(math.sqrt(0.4))


def test_weighted_geomean_gate_on_zero():
    # any zero value with positive weight zeroes the whole result (the gate)
    assert _weighted_geomean([(0.0, 1), (0.9, 1)]) == 0.0


def test_weighted_geomean_ignores_zero_weight_factors():
    # the 0.9 factor has weight 0, so only 0.5 counts
    assert _weighted_geomean([(0.9, 0), (0.5, 1)]) == pytest.approx(0.5)


def test_weighted_geomean_none_when_no_positive_weight():
    assert _weighted_geomean([]) is None
    assert _weighted_geomean([(0.5, 0), (0.6, 0)]) is None


def test_weighted_geomean_weighted():
    # exp((2*ln0.9 + 1*ln0.4) / 3) = 0.686817...
    assert _weighted_geomean([(0.9, 2), (0.4, 1)]) == pytest.approx(0.68682, abs=1e-4)


# ------------------------------------------------ DoranetScoreCriterion
def test_doranet_minmax_normalized():
    scores = DoranetScoreCriterion().score([mk(30, [], []), mk(10, [], []), mk(20, [], [])])
    assert scores == pytest.approx([1.0, 0.0, 0.5])


def test_doranet_all_equal_gives_ones():
    assert DoranetScoreCriterion().score([mk(5, [], []), mk(5, [], [])]) == [1.0, 1.0]


def test_doranet_none_score_treated_as_zero():
    assert DoranetScoreCriterion().score([mk(None, [], []), mk(10, [], [])]) == \
        pytest.approx([0.0, 1.0])


# --------------------------------------------- ProcedureDiversityCriterion
@pytest.mark.parametrize("names,expected", [
    (["a", "a"], 1.0),          # 2 steps, 1 distinct -> all same
    (["a", "b"], 0.0),          # 2 steps, 2 distinct -> all different
    (["a", "b", "c"], 0.0),     # 3 steps, 3 distinct
    (["a", "a", "b"], 0.5),     # 3 steps, 2 distinct
    (["a"], 1.0),               # single step
])
def test_diversity(names, expected):
    smiles = ["X>>Y"] * len(names)
    assert ProcedureDiversityCriterion().score([mk(1, smiles, names)])[0] == \
        pytest.approx(expected)


# ------------------------------------------- IntermediateStabilityCriterion
def test_stability_benign_pathway_is_one():
    p = mk(1, ["CCO>>CCO.O"], ["opA"])          # ethanol + water, both benign
    assert IntermediateStabilityCriterion().score([p])[0] == 1.0


def test_stability_peroxide_intermediate_gates_to_zero():
    p = mk(1, ["CCO>>OO"], ["opA"])             # OO = peroxide, Tier-A catastrophic
    assert IntermediateStabilityCriterion().score([p])[0] == 0.0


def test_stability_excluded_smiles_are_skipped():
    # same peroxide, but excluded (as starter/target/helper would be) -> not penalized
    p = mk(1, ["CCO>>OO"], ["opA"])
    assert IntermediateStabilityCriterion(excluded_smiles=["OO"]).score([p])[0] == 1.0


# ------------------------------------------------ FeasibilityCriterion (bio)
class _FakeDora:
    def __init__(self, table):
        self.table = table

    def feasibility(self, rxn):
        return self.table.get(rxn)


def test_feasibility_weakest_bio_step_and_chem_ignored():
    p = mk(1, ["A>>B", "B>>C"], ["rule0087", "Dehydration of Alcohol"])
    # only the bio step "A>>B" is scored; the chem step is ignored
    crit = FeasibilityCriterion(_FakeDora({"A>>B": 0.4}))
    assert crit.score([p])[0] == pytest.approx(0.4)


def test_feasibility_chem_only_pathway_is_one():
    p = mk(1, ["A>>B"], ["Dehydration of Alcohol"])   # no bio step at all
    assert FeasibilityCriterion(_FakeDora({})).score([p])[0] == 1.0


def test_feasibility_unscoreable_bio_step_not_penalized():
    # client returns None (couldn't score) -> no penalty, stays 1.0
    p = mk(1, ["A>>B"], ["rule0087"])
    assert FeasibilityCriterion(_FakeDora({"A>>B": None})).score([p])[0] == 1.0


# -------------------------------------------- EnzymeLoadCriterion
def test_enzyme_load_uses_minimum_enzymes():
    crit = EnzymeLoadCriterion()
    chem = mk(1, ["A>>B"], ["Dehydration of Alcohol"])   # 0 enzymes -> 1.0
    one = mk(1, ["A>>B"], ["rule0087"])                  # 1 enzyme  -> 0.85
    s = crit.score([chem, one])
    assert s[0] == pytest.approx(1.0)
    assert s[1] == pytest.approx(0.85)


# ================================================= apply_lemnisca_blend =====
# Two benign CHEM pathways under DEFAULT weights (layer doranet=2/lemnisca=1;
# lemnisca stability=1/diversity=1/enzyme_load=1). Benign intermediates ->
# stability 1.0; chem-only -> enzyme_load 1.0.
def test_blend_default_weights_numbers():
    p1 = mk(30, ["CCO>>CCO.O", "CCO>>CCO"], ["opA", "opA"])   # diversity 1.0
    p2 = mk(10, ["CCO>>CCO.O", "CCO>>CCO"], ["opA", "opB"])   # diversity 0.0
    out = apply_lemnisca_blend([p1, p2])

    # doranet min-max: p1 -> 1.0, p2 -> 0.0
    assert p1.lemnisca_components["doranet"] == pytest.approx(1.0)
    assert p2.lemnisca_components["doranet"] == pytest.approx(0.0)
    assert p1.lemnisca_components["stability"] == 1.0
    assert p2.lemnisca_components["diversity"] == pytest.approx(0.0)
    # chem-only routes carry no enzyme burden
    assert p1.lemnisca_components["enzyme_load"] == pytest.approx(1.0)
    assert p2.lemnisca_components["enzyme_load"] == pytest.approx(1.0)

    # p1: lemnisca = geomean(floored(1,0)=1, floored(1,0.5)=1, 1) = 1.0
    #     blended  = geomean(floored(1,0.5)=1 ^2, 1 ^1) = 1.0
    assert p1.lemnisca_score == pytest.approx(1.0)
    assert p1.blended_score == pytest.approx(1.0)

    # p2: lemnisca = geomean(1, floored(0,0.5)=0.5, 1) = 0.5 ** (1/3)
    #     blended  = geomean(floored(0,0.5)=0.5 ^2, lemnisca ^1) = 0.5833
    assert p2.lemnisca_score == pytest.approx(0.5 ** (1 / 3))
    assert p2.blended_score == pytest.approx(0.583, abs=1e-3)

    # sorted best-first, ranks reassigned
    assert out[0] is p1 and out[1] is p2
    assert [p.rank for p in out] == [1, 2]


def test_blend_catastrophic_intermediate_overrides_best_doranet():
    """The crown-jewel property: a pathway with the BEST DORAnet score but a
    catastrophic intermediate must be gated to 0 and sort LAST."""
    p_gate = mk(100, ["CCO>>OO"], ["opA"])   # highest DORAnet, but peroxide -> stability 0
    p_good = mk(1, ["CCO>>CCO"], ["opA"])    # lowest DORAnet, benign
    out = apply_lemnisca_blend([p_gate, p_good])

    assert p_gate.lemnisca_score == 0.0      # stability gate propagated
    assert p_gate.blended_score == 0.0       # ... all the way to the final score
    assert p_good.blended_score > 0.0
    assert out[0] is p_good and out[1] is p_gate   # gated one sorts last


def test_blend_zero_lemnisca_weights_is_neutral():
    # all lemnisca weights 0 -> sub-score falls back to neutral 1.0
    p = mk(10, ["CCO>>CCO"], ["opA"])
    apply_lemnisca_blend([p], lemnisca_weights={"stability": 0, "diversity": 0})
    assert p.lemnisca_score == pytest.approx(1.0)


def test_blend_zero_lemnisca_layer_weight_uses_only_doranet():
    # layer lemnisca weight 0 -> blended = floored(doranet, 0.5)
    p1 = mk(30, ["CCO>>CCO"], ["opA"])   # doranet -> 1.0 -> floored 1.0
    p2 = mk(10, ["CCO>>CCO"], ["opA"])   # doranet -> 0.0 -> floored 0.5
    apply_lemnisca_blend([p1, p2], layer_weights={"doranet": 1, "lemnisca": 0})
    assert p1.blended_score == pytest.approx(1.0)
    assert p2.blended_score == pytest.approx(0.5)
