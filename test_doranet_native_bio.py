"""
test_doranet_native_bio.py
==========================

DIAGNOSTIC: bypasses our generate_network_tal wrapper entirely.
Calls DORAnet's own enzymatic.generate_network directly with our
5-rule polyketide whitelist (we filter JN1224MIN down to just those
5 rules, write a temp TSV, and patch DORAnet's ruleset registry to
load it).

If this finishes fast, the bottleneck is OUR wrapper (specifically
the chem-filter pipeline we apply to bio reactions).
If this also hangs, DORAnet's own bio expansion can't handle this
chemistry on this hardware regardless of wrapper.

Apples-to-apples comparison: same 5 rules as our wrapper, same
chemistry, just DORAnet's native expansion + filter pipeline.
"""

import csv
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from pathlib import Path

# Same 5 polyketide rules as our wrapper's MINIMAL_POLYKETIDE_WHITELIST
POLYKETIDE_RULES = {
    "rule1118",
    "rule0087",
    "rule0126",
    "rule0350",
    "rule0891",
}

ACETYL_COA = (
    "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)
MALONYL_COA = (
    "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)
TAL = "Cc1cc(O)cc(=O)o1"

STARTER_FILE = "doranet_native_bio_starters.smi"


def make_filtered_tsv(source_tsv, output_tsv, keep_rule_names):
    """
    Write a new TSV containing only the rows whose Name column is in
    keep_rule_names. Preserves header and all columns.
    """
    n_kept = 0
    with open(source_tsv, "r", encoding="utf-8") as fin, \
         open(output_tsv, "w", encoding="utf-8", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames,
                                delimiter="\t")
        writer.writeheader()
        for row in reader:
            if row["Name"] in keep_rule_names:
                writer.writerow(row)
                n_kept += 1
    return n_kept


def main():
    print("=" * 64)
    print(" DORANET NATIVE BIO  —  5-rule polyketide subset")
    print("=" * 64)

    # Step 1: write a filtered TSV with only our 5 polyketide rules
    import doranet
    doranet_enzymatic_dir = (
        Path(doranet.__file__).parent / "modules" / "enzymatic"
    )
    source_tsv = doranet_enzymatic_dir / "JN1224MIN_rules.tsv"
    filtered_tsv = Path("doranet_native_bio_filtered_rules.tsv").resolve()

    print(f"\nFiltering {source_tsv.name} to {len(POLYKETIDE_RULES)} rules...")
    n_kept = make_filtered_tsv(source_tsv, filtered_tsv, POLYKETIDE_RULES)
    print(f"  wrote {filtered_tsv.name} ({n_kept} rules)")

    if n_kept != len(POLYKETIDE_RULES):
        print(f"  WARNING: expected {len(POLYKETIDE_RULES)} rules, got {n_kept}")

    # Step 2: patch DORAnet's ruleset registry to recognise our
    # filtered tsv. Note: enzymatic/__init__.py shadows the submodule
    # name with the function name, so we have to import via sys.modules
    # to actually get the module object (where AVAILABLE_RULESETS lives).
    import importlib
    gn_module = importlib.import_module(
        "doranet.modules.enzymatic.generate_network"
    )
    gn_module.AVAILABLE_RULESETS["POLYKETIDE_ONLY"] = filtered_tsv

    # Step 3: write starter SMILES file
    with open(STARTER_FILE, "w") as f:
        f.write(ACETYL_COA + "\n")
        f.write(MALONYL_COA + "\n")

    print()
    print("Settings:")
    print(f"  ruleset:                   POLYKETIDE_ONLY ({n_kept} rules)")
    print(f"  gen:                       3")
    print(f"  allow_multiple_reactants:  True")
    print(f"  target:                    TAL")
    print()
    print("Launching DORAnet native expansion...")
    t0 = time.time()

    try:
        network = gn_module.generate_network(
            job_name="doranet_native_bio",
            starters=STARTER_FILE,
            gen=3,
            direction="forward",
            allow_multiple_reactants=True,
            targets=TAL,
            ruleset="POLYKETIDE_ONLY",
        )
    except Exception as exc:
        elapsed = time.time() - t0
        print(f"\n[ERROR] DORAnet native raised after {elapsed:.1f}s:")
        print(f"  {type(exc).__name__}: {exc}")
        return

    elapsed = time.time() - t0
    print(f"\n[DONE] DORAnet native finished in {elapsed:.1f}s")
    print(f"  Molecules: {len(network.mols)}")
    print(f"  Reactions: {len(network.rxns)}")

    from rdkit import Chem
    target_canon = Chem.MolToSmiles(Chem.MolFromSmiles(TAL))
    network_smiles = {
        Chem.MolToSmiles(Chem.MolFromSmiles(m.smiles))
        for m in network.mols
        if Chem.MolFromSmiles(m.smiles) is not None
    }
    print(f"\nTAL in network? {target_canon in network_smiles}")


if __name__ == "__main__":
    main()
