"""
test_one_bio_rule.py
====================

ABSOLUTE MINIMUM bio test. Goal: prove that even ONE bio rule can
fire in our setup, going from acetyl-CoA + malonyl-CoA to
acetoacetyl-CoA in a single generation.

If this doesn't finish in seconds, the bio expansion is fundamentally
broken in DORAnet on this hardware — not just slow with many rules.

Setup:
  - Rule: rule1118 only  (Claisen condensation 1, makes acetoacetyl-CoA)
  - Starters: acetyl-CoA, malonyl-CoA
  - Cofactors: a minimal pool (CoA, CO2, water, H+)
  - gen: 1
  - allow_multiple_reactants: True (required for Claisens)
  - Calls DORAnet's NATIVE enzymatic.generate_network — no project wrapper

Success criterion:
  - finishes in < 30s
  - acetoacetyl-CoA appears in the network
"""

import csv
import os
import sys
import time
import importlib
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

from rdkit import RDLogger, Chem
RDLogger.DisableLog("rdApp.*")


ACETYL_COA = (
    "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)
MALONYL_COA = (
    "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)
EXPECTED_ACETOACETYL_COA = (
    "CC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)

RULES_TO_KEEP = {
    "rule1118",   # Claisen 1
    "rule0087",   # Claisen 2
    "rule0126",   # Claisen 2 variant  — does this explode native too?
    "rule0350",   # Claisen 2 variant
    "rule0891",   # cyclization
}
STARTER_FILE = "one_rule_starters.smi"
FILTERED_TSV = Path("one_rule_filtered.tsv").resolve()


def main():
    print("=" * 64)
    print(f" SCALED BIO TEST  —  {len(RULES_TO_KEEP)} rules, gen=1")
    print("=" * 64)

    # Step 1: filter JN1224MIN to just the chosen rules
    import doranet
    src = (Path(doranet.__file__).parent
           / "modules" / "enzymatic" / "JN1224MIN_rules.tsv")
    n_kept = 0
    with open(src, "r", encoding="utf-8") as fin, \
         open(FILTERED_TSV, "w", encoding="utf-8", newline="") as fout:
        r = csv.DictReader(fin, delimiter="\t")
        w = csv.DictWriter(fout, fieldnames=r.fieldnames, delimiter="\t")
        w.writeheader()
        for row in r:
            if row["Name"] in RULES_TO_KEEP:
                w.writerow(row)
                n_kept += 1
    print(f"\nFiltered ruleset to {sorted(RULES_TO_KEEP)} — "
          f"{n_kept} row(s) in TSV")

    # Step 2: patch DORAnet's registry
    gn_module = importlib.import_module(
        "doranet.modules.enzymatic.generate_network"
    )
    gn_module.AVAILABLE_RULESETS["SUBSET_ONLY"] = FILTERED_TSV

    # Step 3: write starter SMILES
    with open(STARTER_FILE, "w") as f:
        f.write(ACETYL_COA + "\n")
        f.write(MALONYL_COA + "\n")

    print(f"Starters: acetyl-CoA + malonyl-CoA")
    print(f"Rules:    {sorted(RULES_TO_KEEP)}")
    print(f"Gen:      3")
    print(f"Multi-substrate: True")
    print(f"\nLaunching DORAnet native expansion...")

    t0 = time.time()
    network = gn_module.generate_network(
        job_name="one_rule_bio",
        starters=STARTER_FILE,
        gen=3,
        direction="forward",
        allow_multiple_reactants=True,
        targets=EXPECTED_ACETOACETYL_COA,
        ruleset="SUBSET_ONLY",
    )
    elapsed = time.time() - t0
    print(f"\n[DONE] Finished in {elapsed:.1f}s")
    print(f"  Molecules: {len(network.mols)}")
    print(f"  Reactions: {len(network.rxns)}")

    canon_expected = Chem.MolToSmiles(Chem.MolFromSmiles(EXPECTED_ACETOACETYL_COA))
    network_smiles = {
        Chem.MolToSmiles(Chem.MolFromSmiles(m.smiles))
        for m in network.mols
        if Chem.MolFromSmiles(m.smiles) is not None
    }
    found = canon_expected in network_smiles
    print(f"\nacetoacetyl-CoA produced? {found}")

    TAL = "Cc1cc(O)cc(=O)o1"
    DIOXOHEXANOYL_COA = (
        "CC(=O)CC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)"
        "COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)"
        "[C@H](O)[C@@H]1OP(=O)(O)O"
    )
    tal_canon = Chem.MolToSmiles(Chem.MolFromSmiles(TAL))
    dx_canon  = Chem.MolToSmiles(Chem.MolFromSmiles(DIOXOHEXANOYL_COA))
    print(f"3,5-dioxohexanoyl-CoA in network? {dx_canon in network_smiles}")
    print(f"TAL in network? {tal_canon in network_smiles}")


if __name__ == "__main__":
    main()
