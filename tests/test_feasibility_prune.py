"""DORA-XGB generation-phase prune. Uses a fake dora_client with hand-picked
scores so pathways can be forced above/below the threshold deterministically
(no real model needed). Exercises the same FeasibilityCriterion + filter that
run_pipeline uses."""
import pytest

from pathway_scoring import RankedPathway, FeasibilityCriterion


def mk(rank, steps):
    return RankedPathway(
        rank=rank, final_score=None, atomic_economy=0.0,
        pathway_byproduct_count=0, intermediate_byproducts={},
        reaction_smiles=[s for s, _ in steps],
        reaction_names=[n for _, n in steps],
        reaction_enthalpies=["No_Thermo"] * len(steps),
    )


class _FakeDora:
    def __init__(self, table):
        self.table = table

    def feasibility(self, rxn):
        return self.table.get(rxn)


# 5 pathways: bio steps mixed with one chem-only route.
P1 = mk(1, [("A>>B", "rule0087"), ("B>>C", "rule1118")])   # bio 0.90, 0.80
P2 = mk(2, [("A>>D", "rule0087"), ("D>>C", "rule0891")])   # bio 0.90, 0.20
P3 = mk(3, [("A>>E", "Dehydration"), ("E>>C", "Aldol")])   # chem-only
P4 = mk(4, [("A>>F", "rule1118"), ("F>>C", "rule0891")])   # bio 0.10, 0.15
P5 = mk(5, [("A>>G", "rule0087"), ("G>>C", "Dehydration")])  # bio 0.55 + chem
PATHS = [P1, P2, P3, P4, P5]

CLIENT = _FakeDora({
    "A>>B": 0.90, "B>>C": 0.80, "A>>D": 0.90, "D>>C": 0.20,
    "A>>F": 0.10, "F>>C": 0.15, "A>>G": 0.55,
})


def test_weakest_bio_step_and_chem_immune():
    scores = FeasibilityCriterion(CLIENT).score(PATHS)
    assert scores == pytest.approx([0.80, 0.20, 1.0, 0.10, 0.55])


def _prune(thr):
    scores = FeasibilityCriterion(CLIENT).score(PATHS)
    kept = [p for p, s in zip(PATHS, scores) if s >= thr]
    return {p.rank for p in kept}


@pytest.mark.parametrize("thr,kept", [
    (0.0,  {1, 2, 3, 4, 5}),   # nothing pruned
    (0.5,  {1, 3, 5}),         # drops P2 (0.20), P4 (0.10)
    (0.75, {1, 3}),            # also drops P5 (0.55)
    (0.85, {3}),               # only chem-only survives
])
def test_prune_thresholds(thr, kept):
    assert _prune(thr) == kept
