"""
Novel-chemistry test: TAL -> 2,4-pentanedione (acetylacetone) + CO2.

This transformation was reported in the literature (ring-opening +
decarboxylation in water, no catalyst) but isn't a canonical TAL
derivative — the question is whether the network finds it on its own.

  TAL  +  H2O  -->  triacetic acid (open chain)         [hydrolysis]
  triacetic acid -->  acetylacetone  +  CO2             [decarboxylation]

Both operators ('Hydrolysis of Ethers, Esters, Anhydrides' and
'Carboxylic Acids Decarboxylation') are in the TAL chem whitelist —
so if the route exists structurally in DORAnet's SMARTS, the
pathway finder should surface it at gen=3.
"""

import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from network_generation import generate_network_tal
from pathway_tools import find_pathways_to_target, load_pathways_from_file
from recipe_rankers import ForwardProductTanimotoRanker


STARTER_FILE = "test_novel_starter.smi"
HELPER_FILE = "test_novel_helpers.smi"

TAL = "CC1=CC(O)=CC(=O)O1"  # KEKULIZED form — non-aromatic lactone
TARGET = "CC(=O)CC(=O)C"   # 2,4-pentanedione (acetylacetone)


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" NOVEL CHEMISTRY TEST  —  TAL -> 2,4-pentanedione + CO2")
    print("=" * 64)

    write_smi(STARTER_FILE, [TAL])
    write_smi(HELPER_FILE, ["O", "[H][H]"])

    t0 = time.time()
    network = generate_network_tal(
        job_name="novel_test",
        starters=STARTER_FILE,
        helpers=HELPER_FILE,
        gen=4,
        direction="forward",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        # Tight atom limits keep cartesian tractable at gen=4 — the
        # full TAL+CO2+H2O+acetylacetone chain stays well within these
        # bounds (TAL=6C, triacetic=6C, acetylacetone=5C).
        max_atoms={"C": 8, "O": 5, "N": 0, "S": 0},
        max_molecular_weight=200,
        allow_multiple_reactants="default",
        strategy="cartesian",
        min_carbons=0,
        include_chem=True,
        include_bio=False,
    )

    # Check if target is in network at all
    from rdkit import Chem
    target_canon = Chem.MolToSmiles(Chem.MolFromSmiles(TARGET))
    network_smiles = {
        Chem.MolToSmiles(Chem.MolFromSmiles(m.smiles))
        for m in network.mols
        if Chem.MolFromSmiles(m.smiles) is not None
    }
    in_network = target_canon in network_smiles
    print(f"\nTarget canonical SMILES: {target_canon}")
    print(f"Target in network?       {in_network}")
    print(f"Network has {len(network_smiles)} molecules total.")

    if not in_network:
        print("\nTarget not reached. Probable causes:")
        print("  - decarboxylation SMARTS doesn't match the open-chain")
        print("    triacetic acid intermediate")
        print("  - the open-chain intermediate itself isn't being formed")
        print("    (i.e., the hydrolysis operator doesn't fire on TAL's")
        print("    aromatic 2-pyranone lactone).")
        # Show which TAL-derivatives ARE in the network for context
        candidates = [
            ("triacetic acid (open chain)", "CC(=O)CC(=O)CC(=O)O"),
            ("acetic acid",                  "CC(=O)O"),
            ("parasorbic acid",              "CC1CC=CC(=O)O1"),
            ("2H-pyran-2-one (no methyl)",   "O=C1OC=CC=C1"),
        ]
        print("\nDiagnostic — what IS in the network:")
        for name, smi in candidates:
            canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
            mark = "FOUND" if canon in network_smiles else "missing"
            print(f"  [{mark}] {name}  ({canon})")
        return

    # Find pathways
    print(f"\nSearching for TAL -> acetylacetone pathways...")
    find_pathways_to_target(
        network=network,
        starter=TAL,
        target=TARGET,
        helpers=["O", "[H][H]"],
        generations=4,
        max_num_rxns=10,
        job_name="novel_test",
    )

    pathways = load_pathways_from_file("novel_test")
    if not pathways:
        print("\nTarget was in the network, but no pathway was traced "
              "to it. May indicate disconnected route or filter pruning.")
        return

    print(f"\n[SUCCESS] Found {len(pathways)} pathway(s) to acetylacetone.")
    for i, p in enumerate(pathways, 1):
        print(f"\n  Pathway {i}  ({p.num_steps} steps):")
        for j, rxn in enumerate(p.reactions, 1):
            from pathway_tools import parse_reaction_string
            parsed = parse_reaction_string(rxn)
            arrow = (
                " + ".join(parsed["reactants"]) + "  ->  "
                + " + ".join(parsed["products"])
            )
            print(f"    Step {j}. [{parsed['op_name']}]")
            print(f"             {arrow}")

    print(f"\nTotal runtime: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
