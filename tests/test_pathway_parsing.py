"""Parsing the raw DORAnet pathway file into RankedPathway objects, plus the
bio/chem operator classifier. Uses a real captured pathway file as a fixture
so the format assertions can't drift from reality."""
from conftest import FIXTURES

from pipeline import _parse_unranked_pathways
from pathway_scoring import is_bio_op

# fixture is tests/fixtures/smoke_test_pathways.txt (9 chem TAL->sorbic routes)
FIXTURE_JOB = str(FIXTURES / "smoke_test")


def test_parses_all_pathways():
    paths = _parse_unranked_pathways(FIXTURE_JOB)
    assert len(paths) == 9


def test_parsed_pathways_are_unranked_and_wellformed():
    paths = _parse_unranked_pathways(FIXTURE_JOB)
    for p in paths:
        assert p.final_score is None                       # unranked marker
        n = len(p.reaction_smiles)
        assert n == len(p.reaction_names) == len(p.reaction_enthalpies)
        assert all(">>" in s for s in p.reaction_smiles)   # reactants>>products


def test_parsed_pathways_sorted_by_step_count():
    steps = [p.num_steps for p in _parse_unranked_pathways(FIXTURE_JOB)]
    assert steps == sorted(steps)      # non-decreasing
    assert min(steps) == 5             # the two 5-step routes come first


def test_parsed_chem_pathways_have_named_operators():
    paths = _parse_unranked_pathways(FIXTURE_JOB)
    all_names = {n for p in paths for n in p.reaction_names}
    assert "Dehydration of Alcohol" in all_names
    # these are chem operators, so none should read as bio rules
    assert not any(is_bio_op(n) for n in all_names)


def test_missing_file_returns_empty_not_crash():
    assert _parse_unranked_pathways(str(FIXTURES / "does_not_exist")) == []


def test_is_bio_op_classification():
    assert is_bio_op("rule0087")
    assert is_bio_op("rule1118")
    assert not is_bio_op("Dehydration of Alcohol")
    assert not is_bio_op("Hydrogenation of Alkene")
