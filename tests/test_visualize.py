"""Graph rendering — specifically that op_labels (enzyme names) reach the
edge tooltips. Uses the committed chem pathway fixture; no network."""
from conftest import FIXTURES

from visualize_pathways import visualize_pathways, _chain_endpoints

FIXTURE_JOB = str(FIXTURES / "smoke_test")
STARTER = "Cc1cc(O)cc(=O)o1"   # TAL
TARGET = "CC=CC=CC(=O)O"       # sorbic acid


def test_op_labels_appear_in_edge_tooltip(tmp_path):
    out = str(tmp_path / "g.html")
    path = visualize_pathways(
        job_name=FIXTURE_JOB, starter_smiles=STARTER, target_smiles=TARGET,
        pathway_filter="all", output_html=out,
        op_labels={"Dehydration of Alcohol": "SENTINEL-ENZYME"},
        top_n_threshold=10 ** 9,
    )
    html = open(path, encoding="utf-8").read()
    assert "SENTINEL-ENZYME" in html          # label injected into a tooltip


# The Node-B scenario: a step `Node B + acetyl-CoA -> triketide` where the
# feedstock (acetyl-CoA) is listed first. The arrow must start from the
# INTERMEDIATE (Node B), so Node B connects forward instead of dangling.
PRODUCED = {"nodeB", "triketide", "TAL"}
CONSUMED = {"acetyl", "malonyl", "nodeB", "triketide"}


def test_arrow_follows_intermediate_not_feedstock():
    src, dst = _chain_endpoints(
        ["acetyl", "nodeB"], ["triketide"], PRODUCED, CONSUMED,
        starter="acetyl", target="TAL",
    )
    assert src == "nodeB"        # the growing chain, not the re-used feedstock
    assert dst == "triketide"


def test_first_step_starts_from_the_starter():
    # acetyl + malonyl -> nodeB : neither reactant is a prior product, so the
    # arrow starts at the starter (acetyl), producing nodeB.
    src, dst = _chain_endpoints(
        ["acetyl", "malonyl"], ["nodeB"], PRODUCED, CONSUMED,
        starter="acetyl", target="TAL",
    )
    assert src == "acetyl"
    assert dst == "nodeB"


def test_last_step_ends_at_the_target():
    src, dst = _chain_endpoints(
        ["triketide"], ["TAL"], PRODUCED, CONSUMED, starter="acetyl", target="TAL",
    )
    assert src == "triketide" and dst == "TAL"


def test_no_intermediate_falls_back_to_first():
    # nothing produced/consumed matches -> deterministic fallback to firsts
    src, dst = _chain_endpoints(
        ["X"], ["Y"], set(), set(), starter="Z", target="T",
    )
    assert src == "X" and dst == "Y"


def test_renders_without_op_labels(tmp_path):
    # backwards-compatible: op_labels is optional
    out = str(tmp_path / "g2.html")
    path = visualize_pathways(
        job_name=FIXTURE_JOB, starter_smiles=STARTER, target_smiles=TARGET,
        pathway_filter="all", output_html=out, top_n_threshold=10 ** 9,
    )
    assert open(path, encoding="utf-8").read()  # non-empty HTML, no crash
