"""Graph rendering — specifically that op_labels (enzyme names) reach the
edge tooltips. Uses the committed chem pathway fixture; no network."""
from conftest import FIXTURES

from visualize_pathways import visualize_pathways

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


def test_renders_without_op_labels(tmp_path):
    # backwards-compatible: op_labels is optional
    out = str(tmp_path / "g2.html")
    path = visualize_pathways(
        job_name=FIXTURE_JOB, starter_smiles=STARTER, target_smiles=TARGET,
        pathway_filter="all", output_html=out, top_n_threshold=10 ** 9,
    )
    assert open(path, encoding="utf-8").read()  # non-empty HTML, no crash
