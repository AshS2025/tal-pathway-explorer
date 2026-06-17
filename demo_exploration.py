"""
TAL Downstream Exploration — Demo

What this demonstrates
----------------------
The brief asks: "What products can feasibly be produced with TAL as
the starting point?"

This script answers that by running open-exploration mode:
  1. Expand chemistry forward from TAL with no preset target.
  2. Score every reachable molecule by an interestingness heuristic
     (carbon-count gate, Bertz complexity, functional-group
     diversity, aromatic bonus).
  3. Return the top N endpoints with ranked pathways to each.
  4. Cross-check against a curated list of literature TAL
     derivatives (sorbic acid, phloroglucinol, etc.) as a
     diagnostic — were the expected products reached?

How to run
----------
From the repo root:

    python demo_exploration.py

Runs in ~30-60 seconds on a laptop. Writes a markdown report you
can open in any viewer to walk through with stakeholders.
"""

import os
import sys
import time
import warnings

# Project src on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

# UTF-8 stdout so Greek letters in chemistry names don't crash
sys.stdout.reconfigure(encoding="utf-8")

# Silence the RDKit "Explicit valence" warning stream — it floods
# the console during expansion because many trial reactions produce
# valence-violating intermediates that are rejected by filters.
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

from pathway_tools import explore_downstream
from tal_downstream_derivatives import TAL_DOWNSTREAM_DERIVATIVES


# ---------------------------------------------------------------------
# Configuration — kept conservative so the demo finishes in under a
# minute on a laptop. For a deeper exploration, bump `gen` to 3 (the
# cartesian depth warning will fire — that's expected).
# ---------------------------------------------------------------------
JOB_NAME = "tal_demo"
STARTER_FILE = "demo_starter.smi"
HELPER_FILE = "demo_helpers.smi"

TAL_SMILES = "Cc1cc(O)cc(=O)o1"   # 4-hydroxy-6-methyl-2H-pyran-2-one
HELPERS = ["O", "[H][H]"]         # water, hydrogen

GENERATIONS = 2
TOP_N = 10
MAX_NUM_RXNS = 3


def banner(text):
    print()
    print("=" * 64)
    print(f" {text}")
    print("=" * 64)


def write_smiles_file(path, smiles_list):
    with open(path, "w") as f:
        for s in smiles_list:
            f.write(s + "\n")


def main():
    banner("TAL DOWNSTREAM EXPLORATION  —  LIVE DEMO")
    print()
    print("Question (from the project brief):")
    print('  "What products can feasibly be produced with TAL as the')
    print('   starting point?"')
    print()
    print("Approach:")
    print("  1. Expand chemistry forward from TAL with no preset target")
    print("  2. Score every reachable molecule for 'interestingness'")
    print("  3. Find pathways to the top endpoints")
    print("  4. Cross-check against literature TAL derivatives")
    print()
    print(f"Settings: starter=TAL, gen={GENERATIONS}, "
          f"strategy=cartesian, top_n={TOP_N}")
    print(f"Chem-only (bio disabled for demo speed)")

    write_smiles_file(STARTER_FILE, [TAL_SMILES])
    write_smiles_file(HELPER_FILE, HELPERS)

    t0 = time.time()

    result = explore_downstream(
        starter=STARTER_FILE,
        helpers=HELPER_FILE,
        job_name=JOB_NAME,
        top_n=TOP_N,
        max_num_rxns=MAX_NUM_RXNS,
        network_kwargs=dict(
            gen=GENERATIONS,
            direction="forward",
            molecule_thermo_calculator=None,
            max_rxn_thermo_change=15.0,
            max_atoms={"C": 50, "O": 8, "N": 0, "S": 0},
            max_molecular_weight=500,
            allow_multiple_reactants="default",
            strategy="cartesian",
            min_carbons=0,
            include_chem=True,
            include_bio=False,
        ),
        scoring_kwargs=dict(
            carbon_window=(3, 12),
            require_oxygen=True,
        ),
        derivatives_list=TAL_DOWNSTREAM_DERIVATIVES,
    )

    elapsed = time.time() - t0

    # ----------------------------------------------------------------
    # Final summary — the headline numbers a stakeholder wants to see.
    # ----------------------------------------------------------------
    banner("RESULTS")
    print()
    print(f"  Time elapsed:                {elapsed:.1f} seconds")
    print(f"  Top endpoints surfaced:      {len(result['top_endpoints'])}")
    print(f"  Endpoints with pathways:     "
          f"{sum(1 for v in result['pathways'].values() if v)} "
          f"/ {len(result['top_endpoints'])}")

    if result["derivative_matches"] is not None:
        n_match = len(result["derivative_matches"])
        n_total = n_match + len(result["derivative_missing"])
        print(f"  Literature derivatives hit:  {n_match} / {n_total}")
        if n_match:
            print()
            print("  Found in network:")
            for d in result["derivative_matches"]:
                print(f"    - {d['name']}")
        if result["derivative_missing"]:
            print()
            print(f"  Not reached this run (diagnostic):")
            for d in result["derivative_missing"][:5]:
                print(f"    - {d['name']}")
            if len(result["derivative_missing"]) > 5:
                print(f"    ... + {len(result['derivative_missing']) - 5} more")

    print()
    print(f"  Top 5 endpoints by interestingness score:")
    for i, es in enumerate(result["top_endpoints"][:5], 1):
        n_paths = len(result["pathways"].get(es.smiles, []))
        print(f"    {i}. {es.smiles}")
        print(f"       score={es.score:.3f}  C={es.carbons}  "
              f"Bertz={es.bertz:.0f}  FGs={es.n_functional_groups}  "
              f"pathways={n_paths}")

    banner("REPORT")
    print()
    print(f"  Full ranked report (every endpoint + every pathway")
    print(f"  + cross-check):")
    print()
    print(f"      {JOB_NAME}_exploration_report.md")
    print()
    print(f"  Open in any markdown viewer (VS Code: Ctrl+Shift+V).")
    print()
    print("=" * 64)


if __name__ == "__main__":
    main()
