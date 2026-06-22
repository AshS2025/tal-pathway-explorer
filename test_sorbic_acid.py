"""
Forward chem-only test: TAL -> sorbic acid via the literature
3-step route (which expands to ~6 elementary steps in DORAnet's
SMARTS-level operators).

Strategy: tighten the atom budget to EXACTLY what the literature
path requires. Sorbic acid is C6/O2, TAL is C6/O3, every intermediate
is at most C6/O3 (water briefly enters a ring as O before
dehydration removes it). Anything bigger — Diels-Alder dimers, ether
adducts with the helpers — is structurally blocked by the cap.

Expected literature path (mapped to our operator names):
  1.  TAL  --Keto-enol Tautomerization Reverse-->  non-aromatic TAL
  2.       --Hydrogenation of C=C-->  partially-reduced TAL
  3.       --Reduction of ketone-->  HMP
  4.  HMP  --Dehydration of Alcohol-->  PSA
  5.  PSA  --Hydrolysis of Esters, Intramolecular-->  open-chain hydroxy-acid
  6.       --Dehydration of Alcohol-->  sorbic acid
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


STARTER_FILE = "test_sorbic_starter.smi"
HELPER_FILE = "test_sorbic_helpers.smi"

TAL = "Cc1cc(O)cc(=O)o1"
# Stereo-less SMILES — DORAnet's Dehydration of Alcohol operator
# produces the alkene without E/Z. The stereo-specific form silently
# fails to match what the network produces. Hidden bug found via
# test_psa_to_sorbic.py.
SORBIC_ACID = "CC=CC=CC(=O)O"            # 2,4-hexadienoic acid (no stereo)


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" TAL -> SORBIC ACID  —  forward chem, tight atom budget")
    print("=" * 64)

    write_smi(STARTER_FILE, [TAL])
    write_smi(HELPER_FILE, ["O", "[H][H]"])

    t0 = time.time()
    network = generate_network_tal(
        job_name="sorbic_test",
        starters=STARTER_FILE,
        helpers=HELPER_FILE,
        gen=8,
        direction="forward",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        # EXACTLY what the literature path needs:
        #   sorbic acid     = C6/O2
        #   TAL             = C6/O3
        #   all intermediates fit within (6, 3)
        # Anything bigger (dimers, ether adducts) gets blocked by
        # the cap, forcing the search through the real route.
        max_atoms={"C": 6, "O": 3, "N": 0, "S": 0},
        max_molecular_weight=150,
        allow_multiple_reactants="default",
        # Priority-queue beam search ranked by Tanimoto similarity to
        # sorbic acid. Each intermediate in the literature path
        # (non-aromatic TAL, HMP, PSA, open-chain hexenoate) shares
        # the 6-carbon + carboxyl/lactone motif with sorbic acid, so
        # the ranker has a real signal to follow (unlike forward
        # Tanimoto for tiny early fragments).
        strategy="priority_queue",
        targets=SORBIC_ACID,
        recipe_ranker=ForwardProductTanimotoRanker(SORBIC_ACID),
        beam_size=1000,
        min_carbons=0,
        include_chem=True,
        include_bio=False,
    )

    elapsed = time.time() - t0
    print(f"\nExpansion finished in {elapsed:.1f}s.")

    from rdkit import Chem
    target_canon = Chem.MolToSmiles(Chem.MolFromSmiles(SORBIC_ACID))
    network_smiles = {
        Chem.MolToSmiles(Chem.MolFromSmiles(m.smiles))
        for m in network.mols
        if Chem.MolFromSmiles(m.smiles) is not None
    }
    in_network = target_canon in network_smiles
    print(f"\nSorbic acid canonical: {target_canon}")
    print(f"In network?            {in_network}")
    print(f"Network has {len(network_smiles)} molecules total.")

    # Diagnostic — check intermediates from the literature path
    print("\nLiterature-path intermediates check:")
    intermediates = [
        ("Non-aromatic TAL (diketo)", "CC1=CC(=O)CC(=O)O1"),
        ("HMP (saturated lactone)",    "CC1CC(O)CC(=O)O1"),
        ("PSA (parasorbic acid)",      "CC1CC=CC(=O)O1"),
        ("Open-chain hydroxy-hexenoate", "CC(O)CC=CC(=O)O"),
        ("Sorbic acid (no-stereo)",    "CC=CC=CC(=O)O"),
    ]
    for name, smi in intermediates:
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        mark = "FOUND" if canon in network_smiles else "miss "
        print(f"  [{mark}] {name:35s} {canon}")

    if not in_network:
        print("\nSorbic acid not reached. See intermediate check above")
        print("for the first step where the route breaks down.")
        return

    # Found it — trace the pathway
    print(f"\nSearching for TAL -> sorbic acid pathways...")
    find_pathways_to_target(
        network=network,
        starter=TAL,
        target=SORBIC_ACID,
        helpers=["O", "[H][H]"],
        generations=8,
        max_num_rxns=15,
        job_name="sorbic_test",
    )
    try:
        pathways = load_pathways_from_file("sorbic_test")
    except FileNotFoundError:
        print("\nNo pathway file. Pathway_finder rejected all routes.")
        return

    if not pathways:
        print("\nSorbic acid in network but no traceable pathway.")
        return

    print(f"\n[SUCCESS] Found {len(pathways)} pathway(s) TAL -> sorbic acid.")
    from pathway_tools import parse_reaction_string
    for i, p in enumerate(pathways[:3], 1):
        print(f"\n  Pathway {i}  ({p.num_steps} steps):")
        for j, rxn in enumerate(p.reactions, 1):
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
