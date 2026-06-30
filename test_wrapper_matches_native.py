"""
test_wrapper_matches_native.py
==============================

Apples-to-apples diagnostic: call our refactored wrapper with EXACTLY
the same setup that the native test (test_one_bio_rule.py) used —
3 rules, all 41 cofactors, gen=3 — and time it.

If this finishes in ~30s, the wrapper is fine; the previous slowness
was the extra rules / cofactor restriction in test_forward_polyketide.py.
If this still hangs, our wrapper has a defect even with the same inputs
as native.
"""

import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from network_generation import generate_network_tal


ACETYL_COA = (
    "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)
MALONYL_COA = (
    "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)

# Same 3 rules as test_one_bio_rule.py
RULES_3 = frozenset({"rule1118", "rule0087", "rule0891"})

STARTER = "wrapper_matches_native_starters.smi"


def main():
    with open(STARTER, "w") as f:
        f.write(ACETYL_COA + "\n")
        f.write(MALONYL_COA + "\n")

    print("=" * 60)
    print(" WRAPPER == NATIVE DIAGNOSTIC")
    print("=" * 60)
    print(" Rules: 3, cofactors: all 41 (no restriction), gen: 3")
    print()

    t0 = time.time()
    network = generate_network_tal(
        job_name="wrapper_matches_native",
        starters=STARTER,
        gen=3,
        direction="forward",
        max_rxn_thermo_change=15.0,
        max_atoms=None,                       # match native default (no atom cap)
        allow_multiple_reactants="default",   # native default
        include_chem=False,
        include_bio=True,
        bio_allow_multiple_reactants=True,
        bio_whitelist=RULES_3,
        included_cofactors=None,              # no cofactor restriction
    )
    elapsed = time.time() - t0
    print(f"\n[DONE] wrapper finished in {elapsed:.1f}s")
    print(f"  mols={len(network.mols)}, rxns={len(network.rxns)}")


if __name__ == "__main__":
    main()
