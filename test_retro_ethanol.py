"""
Clean retro bio test: ethanol <- pyruvate (fermentation, retro).

The textbook fermentation chain is:
    pyruvate --[pyruvate decarboxylase]--> acetaldehyde --[ADH]--> ethanol

In retro, the search walks backward from ethanol toward pyruvate.
Both steps are single-substrate single-enzyme reactions — the
non-pathological case bio expansion is designed for. With our new
default (bio_allow_multiple_reactants=False, i.e. single-substrate
filter ON), this should run cleanly with the full bio whitelist.

Purpose: validate that the retro-bio pipeline works end-to-end on
NORMAL biology, separate from the polyketide complications around TAL.
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
from recipe_rankers import FeedstockProximityRanker


STARTER_FILE = "test_ethanol_starter.smi"
HELPER_FILE  = "test_ethanol_helpers.smi"

ETHANOL  = "CCO"
PYRUVATE = "CC(=O)C(=O)O"
ACETALDEHYDE = "CC=O"


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" RETRO BIO TEST  —  ethanol  <--  pyruvate")
    print("=" * 64)
    print()
    print(" Expected route (forward biology):")
    print("   pyruvate --[pyruvate decarboxylase]--> acetaldehyde")
    print("   acetaldehyde --[alcohol dehydrogenase]--> ethanol")
    print()
    print(" Retro search walks: ethanol -> acetaldehyde -> pyruvate")
    print()

    write_smi(STARTER_FILE, [ETHANOL])
    write_smi(HELPER_FILE,  ["O", "[H][H]"])

    t0 = time.time()
    network = generate_network_tal(
        job_name="retro_ethanol",
        starters=STARTER_FILE,
        helpers=HELPER_FILE,
        gen=3,                          # only 2 steps needed; buffer of 1
        direction="retro",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        max_atoms={"C": 6, "O": 4, "N": 0, "S": 0},
        max_molecular_weight=200,
        # Mentor's scope: retro = bio only
        include_chem=False,
        include_bio=True,
        # NEW default: single-substrate filter ON. Fermentation steps
        # are single-substrate, so this is the right regime.
        # (Don't set bio_allow_multiple_reactants — let it default to False)
        strategy="priority_queue",
        targets=PYRUVATE,
        recipe_ranker=FeedstockProximityRanker([PYRUVATE]),
        beam_size=50,
        min_carbons=0,
    )

    elapsed = time.time() - t0
    print(f"\nNetwork expansion finished in {elapsed:.1f}s.")

    from rdkit import Chem
    network_smiles = {
        Chem.MolToSmiles(Chem.MolFromSmiles(m.smiles))
        for m in network.mols
        if Chem.MolFromSmiles(m.smiles) is not None
    }
    print(f"Network has {len(network_smiles)} molecules.")

    # Check intermediates
    print("\nIntermediate check:")
    for name, smi in [
        ("Ethanol (starter)", ETHANOL),
        ("Acetaldehyde",      ACETALDEHYDE),
        ("Pyruvate",          PYRUVATE),
    ]:
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        mark = "FOUND" if canon in network_smiles else "miss "
        print(f"  [{mark}] {name:25s} {canon}")

    pyruvate_canon = Chem.MolToSmiles(Chem.MolFromSmiles(PYRUVATE))
    if pyruvate_canon not in network_smiles:
        print("\nPyruvate not reached. Diagnostic — 2-3 carbon molecules"
              " in the network:")
        for smi in sorted(network_smiles):
            m = Chem.MolFromSmiles(smi)
            if m is None: continue
            n_c = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")
            if 2 <= n_c <= 3:
                print(f"  {smi}")
        return

    # Trace the pathway
    print(f"\nSearching ethanol <- pyruvate pathways...")
    find_pathways_to_target(
        network=network,
        starter=ETHANOL,
        target=PYRUVATE,
        helpers=["O", "[H][H]"],
        generations=3,
        max_num_rxns=8,
        job_name="retro_ethanol",
    )
    try:
        pathways = load_pathways_from_file("retro_ethanol")
    except FileNotFoundError:
        print("\nNo pathway file written.")
        return
    if not pathways:
        print("\nPyruvate reached but no pathway traced.")
        return

    print(f"\n[SUCCESS] Found {len(pathways)} pathway(s) ethanol <- pyruvate.")
    from pathway_tools import parse_reaction_string
    for i, p in enumerate(pathways[:3], 1):
        print(f"\n  Pathway {i}  ({p.num_steps} enzyme steps):")
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
