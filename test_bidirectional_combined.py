"""
test_bidirectional_combined.py
==============================

BIDIRECTIONAL TAL -> sorbic acid with full pathway extraction.

The key insight (buried in DORAnet's source): pretreat_networks
accepts a LIST of networks and automatically flips any reaction
flagged with Reaction_direction="retro" back to forward direction.
That means:

  1. Build a forward network from TAL.
  2. Build a retro network from sorbic acid (reactions stored in
     retro direction with Reaction_direction="retro" meta).
  3. Pass BOTH to pretreat_networks(networks=[fwd, retro], ...).
     The retro reactions get auto-flipped to forward direction and
     merged with the forward network's reactions. Duplicates are
     dropped.
  4. Run pathway_finder on the unified pretreated file. It traces
     continuous TAL -> ... -> sorbic_acid routes through whatever
     intermediates connect the two halves.

This is the proper bidirectional pathway extraction — what
test_meet_in_middle.py was trying to do manually, but using
DORAnet's built-in machinery.
"""

import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from network_generation import generate_network_tal
from pathway_tools import load_pathways_from_file
from recipe_rankers import FeedstockProximityRanker
from doranet.modules.post_processing.post_processing import (
    pretreat_networks,
    pathway_finder,
)


TAL         = "Cc1cc(O)cc(=O)o1"
SORBIC_ACID = "CC=CC=CC(=O)O"

PSA       = "CC1CC=CC(=O)O1"
HEXENOATE = "CC(O)CC=CC(=O)O"

JOB_NAME = "bidir_combined"
HALF_GEN = 4


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" BIDIRECTIONAL TAL -> sorbic acid via combined network")
    print("=" * 64)

    write_smi("bidir_fwd_starter.smi", [TAL])
    write_smi("bidir_fwd_helpers.smi", ["O", "[H][H]"])
    write_smi("bidir_retro_starter.smi", [SORBIC_ACID])
    write_smi("bidir_retro_helpers.smi", ["O", "[H][H]"])

    # ----- STEP 1: build the forward network from TAL --------------
    print(f"\n[1/4] Forward expansion from TAL (gen={HALF_GEN}, cartesian)")
    t0 = time.time()
    fwd_net = generate_network_tal(
        job_name=f"{JOB_NAME}_fwd",
        starters="bidir_fwd_starter.smi",
        helpers="bidir_fwd_helpers.smi",
        gen=HALF_GEN,
        direction="forward",
        molecule_thermo_calculator=None,
        max_atoms={"C": 6, "O": 3, "N": 0, "S": 0},
        max_molecular_weight=150,
        strategy="cartesian",
        include_chem=True,
        include_bio=False,
    )
    print(f"      finished in {time.time()-t0:.1f}s  |  "
          f"{len(fwd_net.mols)} mols, {len(fwd_net.rxns)} reactions")

    # ----- STEP 2: build the retro network from sorbic acid --------
    print(f"\n[2/4] Retro expansion from sorbic acid (gen={HALF_GEN})")
    t0 = time.time()
    retro_net = generate_network_tal(
        job_name=f"{JOB_NAME}_retro",
        starters="bidir_retro_starter.smi",
        helpers="bidir_retro_helpers.smi",
        gen=HALF_GEN,
        direction="retro",
        molecule_thermo_calculator=None,
        max_atoms={"C": 6, "O": 3, "N": 0, "S": 0},
        max_molecular_weight=150,
        strategy="priority_queue",
        targets=TAL,
        recipe_ranker=FeedstockProximityRanker([TAL]),
        beam_size=200,
        include_chem=True,
        include_bio=False,
    )
    print(f"      finished in {time.time()-t0:.1f}s  |  "
          f"{len(retro_net.mols)} mols, {len(retro_net.rxns)} reactions")

    # ----- STEP 3: combine via pretreat_networks -------------------
    # DORAnet's pretreat_networks auto-flips retro reactions to
    # forward direction before merging. The unified output is a
    # single JSON of forward-direction reactions covering both
    # halves. Starter must be TAL since the combined network is
    # presented as "TAL -> ... -> sorbic_acid".
    print(f"\n[3/4] Combining forward + retro into a unified network...")
    t0 = time.time()
    pretreat_networks(
        networks=[fwd_net, retro_net],
        starters=[TAL],
        helpers=["O", "[H][H]"],
        total_generations=HALF_GEN * 2,    # generous bound
        job_name=JOB_NAME,
    )
    print(f"      combined in {time.time()-t0:.1f}s")

    # ----- STEP 4: trace pathways TAL -> sorbic_acid ---------------
    print(f"\n[4/4] Searching for TAL -> sorbic acid pathways through the"
          f" combined network...")
    t0 = time.time()
    pathway_finder(
        starters=[TAL],
        helpers=["O", "[H][H]"],
        target=[SORBIC_ACID],
        search_depth=HALF_GEN * 2,        # allow full bidirectional depth
        max_num_rxns=12,
        job_name=JOB_NAME,
    )
    print(f"      pathway_finder finished in {time.time()-t0:.1f}s")

    # Read the result
    try:
        pathways = load_pathways_from_file(JOB_NAME)
    except FileNotFoundError:
        print("\n  No pathway file written.")
        return
    if not pathways:
        print("\n  No pathways found in the combined network.")
        return

    print(f"\n[SUCCESS] {len(pathways)} pathway(s) TAL -> sorbic acid")
    print("           via bidirectional combined network.\n")

    from pathway_tools import parse_reaction_string
    for i, p in enumerate(pathways[:5], 1):
        print(f"  Pathway {i}  ({p.num_steps} steps):")
        for j, rxn in enumerate(p.reactions, 1):
            parsed = parse_reaction_string(rxn)
            arrow = (
                " + ".join(parsed["reactants"]) + "  ->  "
                + " + ".join(parsed["products"])
            )
            print(f"    Step {j}. [{parsed['op_name']}]")
            print(f"             {arrow}")
        print()


if __name__ == "__main__":
    main()
