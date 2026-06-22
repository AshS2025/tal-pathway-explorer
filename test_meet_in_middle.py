"""
test_meet_in_middle.py
======================

WHAT THIS IS
------------
Bidirectional (meet-in-the-middle) search for TAL -> sorbic acid.

Why bidirectional?
------------------
The full TAL -> sorbic acid route is ~6 elementary chemistry steps:
    TAL -> non-aromatic TAL -> hydrogenated -> HMP -> PSA -> hexenoate -> sorbic
A single forward expansion at gen=6 explodes on a laptop. But we
don't need to walk all 6 steps from one side. If we walk 3 steps
forward from TAL AND 3 steps backward from sorbic acid, both halves
should reach roughly the same middle molecule (likely HMP or PSA).
Where they overlap, we have a valid full pathway.

Mathematically: two gen=3 searches cost ~2 × 3^k work. One gen=6
search costs ~6^k. For k≈10 the bidirectional is ~30x cheaper.

How this script implements it
-----------------------------
1. Build a FORWARD network from TAL using DORAnet's forward chem.
   Targets sorbic acid with the Tanimoto ranker (steers chemistry
   toward sorbic-like fragments).
2. Build a RETRO network from sorbic acid using DORAnet's retro
   chem. Targets TAL with FeedstockProximityRanker (steers
   backward walking toward TAL-like fragments).
3. Compute the intersection of the two molecule sets. Any molecule
   that appears in BOTH is a "meeting point" — a valid intermediate
   on the path TAL ↔ sorbic acid.
4. For each meeting point, trace the forward pathway TAL ->
   meeting_point and report it. The user can read off the half-routes
   on each side.

CAVEAT — retro chem operators are NOT whitelist-filtered
--------------------------------------------------------
generate_network_tal applies TAL_REACTION_WHITELIST only to forward
chem ops. The retro side loads all 388 retro SMARTS, which means the
retro search space is broader than the forward search space. That's
usually fine for "did we meet?" diagnostics but worth knowing if you
want a strict apples-to-apples comparison.
"""

import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

from rdkit import RDLogger, Chem
RDLogger.DisableLog("rdApp.*")

from network_generation import generate_network_tal
from pathway_tools import find_pathways_to_target, load_pathways_from_file
from recipe_rankers import FeedstockProximityRanker


# ---- starter / target ------------------------------------------------
TAL = "Cc1cc(O)cc(=O)o1"
SORBIC_ACID = "CC=CC=CC(=O)O"           # no-stereo form

# Waypoint pool for the forward ranker. Including PSA and the
# open-chain hexenoate ensures the forward beam doesn't prune
# HMP -> PSA -> hexenoate recipes just because their products
# still look only modestly similar to sorbic acid. Any molecule
# in this pool counts as a "good direction" for the search.
PSA       = "CC1CC=CC(=O)O1"
HEXENOATE = "CC(O)CC=CC(=O)O"
FORWARD_WAYPOINT_POOL = [SORBIC_ACID, HEXENOATE, PSA]

FWD_STARTER = "test_mim_fwd_starter.smi"
FWD_HELPERS = "test_mim_fwd_helpers.smi"
RETRO_STARTER = "test_mim_retro_starter.smi"
RETRO_HELPERS = "test_mim_retro_helpers.smi"

# Depth of each half. The full path is 6 steps, so 3+3 should be
# enough to meet — and either side staying at gen=3 stays tractable.
HALF_GEN = 4


# ---- helper to extract canonical SMILES from a network ---------------
def canonical_smiles_set(network):
    """Return the canonical SMILES of every molecule in the network."""
    out = set()
    for mol in network.mols:
        try:
            rd = Chem.MolFromSmiles(mol.smiles)
            if rd is not None:
                out.add(Chem.MolToSmiles(rd))
        except Exception:
            continue
    return out


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" MEET-IN-THE-MIDDLE  —  TAL  ↔  sorbic acid")
    print("=" * 64)
    print(f" Each half runs at gen={HALF_GEN}. Looking for overlap molecules.")

    # ----- write SMILES files -----
    write_smi(FWD_STARTER, [TAL])
    write_smi(FWD_HELPERS, ["O", "[H][H]"])
    write_smi(RETRO_STARTER, [SORBIC_ACID])
    write_smi(RETRO_HELPERS, ["O", "[H][H]"])

    # ----- FORWARD half: TAL -> middle -----
    print()
    print(f"[1/3] Forward expansion from TAL, gen={HALF_GEN}, target=sorbic acid")
    t0 = time.time()
    fwd_network = generate_network_tal(
        job_name="mim_fwd",
        starters=FWD_STARTER,
        helpers=FWD_HELPERS,
        gen=HALF_GEN,
        direction="forward",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        # Tight atom cap matched to sorbic acid chemistry.
        max_atoms={"C": 6, "O": 3, "N": 0, "S": 0},
        max_molecular_weight=150,
        allow_multiple_reactants="default",
        # Cartesian for the forward half — gives the search breadth
        # so the HMP -> PSA dehydration actually fires. With the
        # C=6/O=3 atom cap the search stays small (~100-300 mols
        # at gen=4) because there's no room for Diels-Alder dimers.
        # No ranker needed for cartesian.
        strategy="cartesian",
        min_carbons=0,
        include_chem=True,
        include_bio=False,
    )
    fwd_t = time.time() - t0
    fwd_mols = canonical_smiles_set(fwd_network)
    print(f"     finished in {fwd_t:.1f}s   |  {len(fwd_mols)} mols, "
          f"{len(fwd_network.rxns)} reactions")

    # ----- RETRO half: sorbic acid -> middle -----
    print()
    print(f"[2/3] Retro expansion from sorbic acid, gen={HALF_GEN}, "
          f"target=TAL (via feedstock proximity)")
    t0 = time.time()
    retro_network = generate_network_tal(
        job_name="mim_retro",
        starters=RETRO_STARTER,
        helpers=RETRO_HELPERS,
        gen=HALF_GEN,
        direction="retro",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        max_atoms={"C": 6, "O": 3, "N": 0, "S": 0},
        max_molecular_weight=150,
        allow_multiple_reactants="default",
        strategy="priority_queue",
        targets=TAL,
        # FeedstockProximityRanker takes a POOL of "destination" SMILES
        # and scores recipes by their substrates' max similarity to
        # any of them. For our purposes the pool is just TAL — we
        # want retro to walk toward TAL-like fragments.
        recipe_ranker=FeedstockProximityRanker([TAL]),
        beam_size=200,
        min_carbons=0,
        include_chem=True,
        include_bio=False,
    )
    retro_t = time.time() - t0
    retro_mols = canonical_smiles_set(retro_network)
    print(f"     finished in {retro_t:.1f}s   |  {len(retro_mols)} mols, "
          f"{len(retro_network.rxns)} reactions")

    # ----- INTERSECTION: meeting points -----
    print()
    print("[3/3] Computing meeting points (forward mols ∩ retro mols)")

    # Exclude trivial overlaps: starter, target, helpers — these are
    # in both networks by definition and don't count as "meeting in
    # the middle."
    excluded = {
        Chem.MolToSmiles(Chem.MolFromSmiles(TAL)),
        Chem.MolToSmiles(Chem.MolFromSmiles(SORBIC_ACID)),
        Chem.MolToSmiles(Chem.MolFromSmiles("O")),
        Chem.MolToSmiles(Chem.MolFromSmiles("[H][H]")),
    }
    meeting_points = (fwd_mols & retro_mols) - excluded

    if not meeting_points:
        print("     No overlap. Either:")
        print("       - HALF_GEN too small (try 4)")
        print("       - retro chem operators producing fragments the")
        print("         forward search doesn't see (whitelist mismatch)")
        return

    print(f"     {len(meeting_points)} meeting molecules found:")
    for smi in sorted(meeting_points):
        # Recognize known intermediates by their canonical SMILES
        labels = []
        named = {
            "CC1=CC(=O)CC(=O)O1":  "non-aromatic TAL",
            "CC1CC(O)CC(=O)O1":    "HMP",
            "CC1CC=CC(=O)O1":      "PSA",
            "CC(O)CC=CC(=O)O":     "open-chain hexenoate",
        }
        canon_named = {Chem.MolToSmiles(Chem.MolFromSmiles(s)): n
                       for s, n in named.items()}
        label = canon_named.get(smi, "")
        print(f"     • {smi}{('  (' + label + ')') if label else ''}")

    # ----- TRACE pathways through each meeting point -----
    print()
    print("Tracing TAL -> meeting -> sorbic_acid for each meeting point...")
    for i, meet_smi in enumerate(sorted(meeting_points), 1):
        print(f"\n--- Meeting point {i}: {meet_smi} ---")

        # Forward half: TAL -> meet_smi in fwd_network
        try:
            find_pathways_to_target(
                network=fwd_network,
                starter=TAL,
                target=meet_smi,
                helpers=["O", "[H][H]"],
                generations=HALF_GEN,
                max_num_rxns=6,
                job_name=f"mim_fwd_{i}",
            )
            fwd_paths = load_pathways_from_file(f"mim_fwd_{i}")
        except (FileNotFoundError, Exception) as exc:
            fwd_paths = []
            print(f"  forward trace failed: {exc}")

        # Retro half: sorbic_acid (starter) -> meet_smi in retro_network
        try:
            find_pathways_to_target(
                network=retro_network,
                starter=SORBIC_ACID,
                target=meet_smi,
                helpers=["O", "[H][H]"],
                generations=HALF_GEN,
                max_num_rxns=6,
                job_name=f"mim_retro_{i}",
            )
            retro_paths = load_pathways_from_file(f"mim_retro_{i}")
        except (FileNotFoundError, Exception) as exc:
            retro_paths = []
            print(f"  retro trace failed: {exc}")

        if fwd_paths and retro_paths:
            fp = fwd_paths[0]
            rp = retro_paths[0]
            print(f"  TAL  →  {meet_smi}  ({fp.num_steps} forward steps)")
            print(f"  {meet_smi}  →  sorbic acid  ({rp.num_steps} retro steps)")
            print(f"  TOTAL: {fp.num_steps + rp.num_steps} steps via this "
                  f"meeting point.")
        else:
            print(f"  forward paths: {len(fwd_paths)}, "
                  f"retro paths: {len(retro_paths)} "
                  f"(couldn't reconstruct full route through this point)")


if __name__ == "__main__":
    main()
