"""validate_config edge cases — the first line of defense against a bad run.
A UI shows these messages live, so each branch must fire correctly."""
import pytest

from pipeline import PipelineConfig, validate_config

TAL = "Cc1cc(O)cc(=O)o1"
SORBIC = "CC=CC=CC(=O)O"


def cfg(**over):
    base = dict(starter_smiles=TAL, target_smiles=SORBIC, domain="chem",
                direction="bidirectional", generations=3)
    base.update(over)
    return PipelineConfig.from_dict(base)


def test_valid_config_passes():
    assert validate_config(cfg()) is None


@pytest.mark.parametrize("over,needle", [
    ({"starter_smiles": ""},            "Starter"),
    ({"target_smiles": ""},             "Target"),
    ({"starter_smiles": "C(C"},         "Invalid starter"),
    ({"target_smiles": "C(C"},          "Invalid target"),
    ({"domain": "food"},                "Domain"),
    ({"direction": "sideways"},         "Direction"),
    ({"strategy": "guess"},             "Strategy"),
    ({"generations": 0},                "between"),
    ({"generations": 99},               "between"),
    ({"helpers": ["O", "C(C"]},         "helper"),
    # runaway-search guardrails
    ({"beam_size": 0},                          "Beam size"),
    ({"beam_size": 10_000_000},                 "too large"),
    ({"strategy": "cartesian", "generations": 5}, "cartesian"),
])
def test_invalid_configs_are_rejected(over, needle):
    msg = validate_config(cfg(**over))
    assert msg is not None, f"expected rejection for {over}"
    assert needle.lower() in msg.lower(), f"{needle!r} not in {msg!r}"


@pytest.mark.parametrize("over", [
    {},                                              # defaults (priority_queue, gen 3)
    {"generations": 8},                              # deep, but guided beam -> allowed
    {"strategy": "cartesian", "generations": 3},     # cartesian at the cap -> allowed
    {"beam_size": 100_000},                          # exactly the max -> allowed
])
def test_reasonable_configs_pass(over):
    assert validate_config(cfg(**over)) is None


def test_starter_equals_target_rejected():
    # same molecule written two ways -> canonicalization catches it
    msg = validate_config(cfg(starter_smiles=TAL, target_smiles="c1(C)cc(O)cc(=O)o1"))
    assert msg is not None and "same molecule" in msg.lower()
