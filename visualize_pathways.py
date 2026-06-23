"""
visualize_pathways.py — example caller

This is a thin wrapper around the reusable function
`visualize_pathways()` in src/visualize_pathways.py. It renders the
TAL -> sorbic acid pathways from bidir_combined_pathways.txt.

Edit the call below for other targets, or call the function from any
script after running pathway_finder.

Usage:
    python visualize_pathways.py

Output:
    sorbic_pathways_graph.html
"""

import os
import sys

# Make src/ importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from visualize_pathways import visualize_pathways


if __name__ == "__main__":
    out_path = visualize_pathways(
        job_name="bidir_combined",
        starter_smiles="Cc1cc(O)cc(=O)o1",
        target_smiles="CC=CC=CC(=O)O",
        starter_label="TAL",
        target_label="sorbic acid",
        pathway_filter="all",          # or "shortest" or [1, 4]
        output_html="sorbic_pathways_graph.html",
    )
    print(f"Wrote {out_path}")
