"""
test_wrapper_bio_gen1.py
========================

Diagnostic: minimum-viable run of our wrapper for bio only, gen=1.

Goal: find out whether the slowness in `generate_network_tal` for bio
expansion is setup overhead (operator loading, strategy init) or
per-generation expansion cost. If gen=1 takes minutes the bottleneck
is upstream of expansion. If gen=1 is fast, the cost compounds with
generation count.

Same chemistry, rules, and cofactors as test_forward_polyketide.py,
but bumped down to gen=1 with the gating counter exposed.
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

POLYKETIDE_WHITELIST = frozenset({
    "rule1118", "rule0087", "rule0126", "rule0350", "rule0891",
})
COFACTORS = ("CoA", "CO2", "WATER", "H+")

STARTER = "wrapper_bio_gen1_starter.smi"
HELPERS = "wrapper_bio_gen1_helpers.smi"


def main():
    with open(STARTER, "w") as f:
        f.write(ACETYL_COA + "\n")
        f.write(MALONYL_COA + "\n")
    with open(HELPERS, "w") as f:
        f.write("O\n[H][H]\n")

    print("=" * 60)
    print(" WRAPPER BIO DIAGNOSTIC — gen=1, bio-only")
    print("=" * 60)

    t0 = time.time()
    network = generate_network_tal(
        job_name="wrapper_bio_gen1",
        starters=STARTER,
        helpers=HELPERS,
        gen=1,
        direction="forward",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        max_atoms={"C": 35, "O": 20, "N": 8, "S": 2, "P": 5},
        max_molecular_weight=1100,
        allow_multiple_reactants="default",
        strategy="cartesian",
        min_carbons=0,
        include_chem=False,
        include_bio=True,
        bio_allow_multiple_reactants=True,
        bio_whitelist=POLYKETIDE_WHITELIST,
        included_cofactors=COFACTORS,
    )
    elapsed = time.time() - t0
    print(f"\n[DONE] gen=1 finished in {elapsed:.1f}s")
    print(f"  mols={len(network.mols)}, rxns={len(network.rxns)}")


if __name__ == "__main__":
    main()
