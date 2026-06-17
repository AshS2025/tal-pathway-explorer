"""
Smoke test for open-exploration mode.

Runs a tiny gen=1 chem-only expansion from TAL and prints the top-N
endpoints, pathways found per endpoint, and the cross-check against
the literature derivatives list.

Deliberately small: gen=1, chem only, no bio (to dodge the cofactor
combinatorics issue), top_n=5. Should finish in tens of seconds.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

# Make stdout UTF-8 so chemistry names with greek letters don't crash
sys.stdout.reconfigure(encoding="utf-8")

from pathway_tools import explore_downstream
from tal_downstream_derivatives import TAL_DOWNSTREAM_DERIVATIVES


# Minimal starter + helper SMILES files
def _write_smiles_file(path, smiles_list):
    with open(path, "w") as f:
        for s in smiles_list:
            f.write(s + "\n")


_write_smiles_file("test_starter.smi", ["Cc1cc(O)cc(=O)o1"])  # TAL
_write_smiles_file("test_helpers.smi", ["O", "[H][H]"])       # H2O, H2


result = explore_downstream(
    starter="test_starter.smi",
    helpers="test_helpers.smi",
    job_name="tal_explore_smoke",
    top_n=5,
    max_num_rxns=2,
    network_kwargs=dict(
        gen=2,
        direction="forward",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        # Match the working settings from test_tal_pipeline.py
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

print("\n" + "=" * 60)
print("EXPLORATION SMOKE TEST RESULT")
print("=" * 60)
print(f"Top endpoints scored:         {len(result['top_endpoints'])}")
print(f"Endpoints with >=1 pathway:   "
      f"{sum(1 for v in result['pathways'].values() if v)}")
if result["derivative_matches"] is not None:
    print(f"Known derivatives matched:    "
          f"{len(result['derivative_matches'])}/"
          f"{len(TAL_DOWNSTREAM_DERIVATIVES)}")
print("\nReport written to: tal_explore_smoke_exploration_report.md")
