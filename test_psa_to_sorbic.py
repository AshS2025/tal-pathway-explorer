"""
Focused chemistry test: PSA -> sorbic acid in 2 steps.

If this works, we know:
  - the network builder + filter pipeline can find the route
  - the issue isolating TAL -> sorbic acid is navigation from TAL
    to PSA (the hard part with the aromatic ring), not navigation
    from PSA to sorbic acid

If this DOESN'T work, the problem is deeper — something in the
pipeline rejects intermediates the hand-tested SMARTS produce.

Expected path:
  PSA --[Hydrolysis of Esters, Intramolecular + H2O]--> open-chain hexenoate
  open-chain hexenoate --[Dehydration of Alcohol]--> sorbic acid
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
from rmg_thermo import RMGThermoClient


STARTER_FILE = "test_psa_starter.smi"
HELPER_FILE = "test_psa_helpers.smi"

PSA = "CC1CC=CC(=O)O1"                    # parasorbic acid
# Sorbic acid SMILES — using NO-stereo form because DORAnet's
# Dehydration of Alcohol operator produces the alkene without
# specifying E/Z. The stereo-specific SMILES "C/C=C/C=C/C(=O)O"
# canonicalizes differently and silently fails to match what the
# network produces. This was hiding TAL -> sorbic acid all along.
SORBIC_ACID = "CC=CC=CC(=O)O"             # 2,4-hexadienoic acid (no stereo)


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" PSA -> SORBIC ACID  —  short forward chem test")
    print("=" * 64)

    write_smi(STARTER_FILE, [PSA])
    write_smi(HELPER_FILE, ["O", "[H][H]"])

    t0 = time.time()
    # Spawn the RMG thermo bridge. Inside this `with` block:
    #   - molecule_thermo_calculator=calc passes RMG into DORAnet's
    #     existing Chem_Rxn_dH_Calculator + Rxn_dH_Filter pipeline,
    #   - reactions whose computed dH exceeds max_rxn_thermo_change
    #     (in kJ/mol — RMG's unit) are HARD-REJECTED, never enter
    #     the network. That cuts ~30% of candidate reactions on
    #     this small expansion (see test_rmg_pruning.py for the
    #     before/after numbers).
    print("Spawning RMG thermo server (~3-10s for first load)...")
    with RMGThermoClient() as calc:
        network = generate_network_tal(
            job_name="psa_test",
            starters=STARTER_FILE,
            helpers=HELPER_FILE,
            gen=4,
            direction="forward",
            molecule_thermo_calculator=calc,    # ← was None
            max_rxn_thermo_change=15.0,         # kJ/mol — matches RMG units
            # Same tight cap as the TAL->sorbic test
            max_atoms={"C": 6, "O": 3, "N": 0, "S": 0},
            max_molecular_weight=150,
            allow_multiple_reactants="default",
            strategy="priority_queue",
            targets=SORBIC_ACID,
            recipe_ranker=ForwardProductTanimotoRanker(SORBIC_ACID),
            beam_size=200,
            min_carbons=0,
            include_chem=True,
            include_bio=False,
        )

        elapsed = time.time() - t0
        print(f"\nExpansion finished in {elapsed:.1f}s.")
        print(f"  RMG queries: {len(calc._cache)} unique molecules looked up")

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

    print("\nIntermediate check:")
    for name, smi in [
        ("PSA (starter)",                   PSA),
        ("Open-chain hydroxy-hexenoate",    "CC(O)CC=CC(=O)O"),
        ("Sorbic acid",                     SORBIC_ACID),
    ]:
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        mark = "FOUND" if canon in network_smiles else "miss "
        print(f"  [{mark}] {name:33s} {canon}")

    if not in_network:
        print("\nSorbic acid not reached. If the open-chain hexenoate")
        print("IS in the network but sorbic acid isn't, the dehydration")
        print("step isn't firing under DORAnet's full pipeline (even")
        print("though it fires when the SMARTS is run directly).")
        return

    print(f"\nSearching for PSA -> sorbic acid pathways...")
    find_pathways_to_target(
        network=network,
        starter=PSA,
        target=SORBIC_ACID,
        helpers=["O", "[H][H]"],
        generations=4,
        max_num_rxns=6,
        job_name="psa_test",
    )
    try:
        pathways = load_pathways_from_file("psa_test")
    except FileNotFoundError:
        print("\nNo pathway file. pathway_finder rejected all routes.")
        return

    if not pathways:
        print("\nSorbic acid in network but no traceable pathway.")
        return

    print(f"\n[SUCCESS] Found {len(pathways)} pathway(s) PSA -> sorbic acid.")
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
