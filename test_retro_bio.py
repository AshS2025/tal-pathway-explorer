"""
End-to-end retro bio test: glucose -> TAL via enzymes.

We expand BACKWARD from TAL using flipped bio SMARTS, then check
whether glucose appears in the resulting network. If it does, we
trace the pathway from TAL back to glucose — which is the
biosynthetic route, read in the natural direction.

Mentor's scope: forward = chem only, retro = bio only. This is the
first end-to-end probe of that retro-bio path.
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


STARTER_FILE = "test_retro_starter.smi"
HELPER_FILE = "test_retro_helpers.smi"

TAL = "Cc1cc(O)cc(=O)o1"

# Pool of candidate feedstocks. Biology's actual TAL biosynthesis runs
# glucose -> pyruvate -> acetyl-CoA -> malonyl-CoA -> 3x condensation
# -> TAL. Reaching glucose from TAL at gen=2 is unrealistic (~6 steps);
# the closer intermediates (acetyl-CoA, pyruvate) should be much
# fewer steps away and provide a cleaner end-to-end validation.
#
# Note: acetyl-CoA + many other cofactors are auto-preloaded as
# coreactant helpers by generate_network_tal, so they exist in the
# network from t=0. The real question is whether pathway_finder can
# trace TAL -> feedstock via actual retro reactions.
FEEDSTOCKS = {
    "glucose":     "OCC(O)C(O)C(O)C(O)C=O",
    "pyruvate":    "CC(=O)C(=O)O",
    "acetyl-CoA":  "CC(=O)S",                  # simplified thiol-only form
    "acetate":     "CC(=O)O",
    "malonate":    "OC(=O)CC(=O)O",
}


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" RETRO BIO TEST  —  TAL  <--  feedstock pool")
    print("=" * 64)
    print()
    print(" Feedstock pool:")
    for name, smi in FEEDSTOCKS.items():
        print(f"   {name:14s} {smi}")
    print()

    # In retro mode the "starter" is the molecule we walk backward FROM.
    # We want to walk back from TAL, so TAL goes in starters.
    write_smi(STARTER_FILE, [TAL])
    # Helpers: water + H2 mirror what's used forward. Bio cofactors are
    # auto-preloaded by generate_network_tal when include_bio=True.
    write_smi(HELPER_FILE, ["O", "[H][H]"])

    t0 = time.time()
    network = generate_network_tal(
        job_name="retro_bio_test",
        starters=STARTER_FILE,
        helpers=HELPER_FILE,
        gen=2,
        direction="retro",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        # Glucose is C6/O6, TAL is C6/O3 — both fit comfortably.
        # Tight limits help with the bio combinatorial pressure.
        max_atoms={"C": 10, "O": 8, "N": 4, "S": 0},
        max_molecular_weight=300,
        allow_multiple_reactants="default",
        # Priority-queue beam search with a feedstock-proximity ranker.
        # Without this, the 342 bio rules x 41 cofactors combinatorial
        # explosion makes even gen=2 intractable on a laptop. Ranker
        # is pointed at the whole feedstock pool — any match counts.
        strategy="priority_queue",
        targets=list(FEEDSTOCKS.values())[0],   # post-expansion check target
        recipe_ranker=FeedstockProximityRanker(list(FEEDSTOCKS.values())),
        beam_size=100,
        min_carbons=0,
        # Mentor's scope: retro = bio only
        include_chem=False,
        include_bio=True,
        bio_allow_multiple_reactants=True,
    )

    elapsed_expand = time.time() - t0
    print(f"\nExpansion finished in {elapsed_expand:.1f}s.")

    # Which feedstocks landed in the network?
    from rdkit import Chem
    network_smiles = {
        Chem.MolToSmiles(Chem.MolFromSmiles(m.smiles))
        for m in network.mols
        if Chem.MolFromSmiles(m.smiles) is not None
    }
    print(f"\nNetwork has {len(network_smiles)} molecules total.")
    print()
    print("Feedstock check:")
    hits = []
    misses = []
    for name, smi in FEEDSTOCKS.items():
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        present = canon in network_smiles
        mark = "FOUND" if present else "miss "
        print(f"  [{mark}] {name:14s} {canon}")
        if present:
            hits.append((name, smi))
        else:
            misses.append(name)

    if not hits:
        print("\nNo feedstocks reached. Diagnostic — molecules with"
              " 6 carbons + multiple OH groups in the network:")
        sugars = []
        for smi in network_smiles:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            n_c = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")
            n_oh = len(m.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
            if n_c == 6 and n_oh >= 2:
                sugars.append(smi)
        for s in sugars[:10]:
            print(f"  {s}")
        if not sugars:
            print("  (none — network doesn't reach sugar-like molecules)")
        return

    # Trace pathways to each feedstock that was reached.
    print(f"\nTracing pathways TAL <- {len(hits)} feedstock(s)...")
    from pathway_tools import parse_reaction_string

    for name, smi in hits:
        sub_job = f"retro_bio_{name}"
        print(f"\n--- TAL <- {name} ---")
        try:
            find_pathways_to_target(
                network=network,
                starter=TAL,
                target=smi,
                helpers=["O", "[H][H]"],
                generations=2,
                max_num_rxns=10,
                job_name=sub_job,
            )
            pathways = load_pathways_from_file(sub_job)
        except FileNotFoundError:
            print("  No pathway file written.")
            continue
        if not pathways:
            print("  Feedstock in network but no pathway traced.")
            continue
        print(f"  Found {len(pathways)} pathway(s).")
        for i, p in enumerate(pathways[:3], 1):
            print(f"    Pathway {i} ({p.num_steps} enzyme steps):")
            for j, rxn in enumerate(p.reactions, 1):
                parsed = parse_reaction_string(rxn)
                arrow = (
                    " + ".join(parsed["reactants"]) + "  ->  "
                    + " + ".join(parsed["products"])
                )
                print(f"      Step {j}. [{parsed['op_name']}]")
                print(f"               {arrow}")

    print(f"\nTotal runtime: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
