"""
test_tal_centered_combined.py
=============================

TAL-centered combined graph. TAL sits at the centre of the network
with two branches:

  LEFT  : TAL → acetyl-CoA     (retro bio, 3 enzyme steps)
  RIGHT : TAL → sorbic acid    (bidirectional chem, ~6 chemistry steps)

Three networks expand outward:

  1. BIO RETRO    : TAL → ... → acetyl-CoA
  2. CHEM FORWARD : TAL → ... → (halfway toward sorbic)
  3. CHEM RETRO   : sorbic acid → ... → (halfway back toward TAL)

All three are merged with starter=TAL. pathway_finder runs TWICE:
  - target = acetyl-CoA      → produces the left branch pathways
  - target = sorbic acid     → produces the right branch pathways

The two pathway files are concatenated. visualize_pathways then renders
TAL as the starter on the left, with branches diverging to acetyl-CoA
and sorbic acid as terminal leaves.
"""

import os
import re
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from network_generation import generate_network_tal
from pathway_tools import load_pathways_from_file, parse_reaction_string
from recipe_rankers import FeedstockProximityRanker
from visualize_pathways import visualize_pathways
from doranet.modules.post_processing.post_processing import (
    pretreat_networks,
    pathway_finder,
)


ACETYL_COA = (
    "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)
MALONYL_COA = (
    "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)
TAL         = "Cc1cc(O)cc(=O)o1"
SORBIC_ACID = "CC=CC=CC(=O)O"

POLYKETIDE_WHITELIST = frozenset({
    "rule1118",   # Claisen 1
    "rule0087",   # Claisen 2
    "rule0891",   # cyclization
})

JOB_NAME = "tal_centered_combined"
BIO_GEN        = 3
CHEM_FWD_GEN   = 4
CHEM_RETRO_GEN = 4


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def _merge_pathway_files(left_file, right_file, out_file):
    """
    Concatenate two `{job_name}_pathways.txt` files into a unified
    file, renumbering the "pathway number X" markers so the combined
    file has a single contiguous sequence.

    The file format is the same one DORAnet's pathway_finder produces
    and load_pathways_from_file consumes — load_pathways_from_file just
    splits on lines containing "pathway number", so renumbering keeps
    it parseable.
    """
    def read_pathways_section(path):
        with open(path, encoding="utf-8") as f:
            text = f.read()
        # split on the marker line (keep markers)
        pieces = re.split(r"(?=^pathway number \d+)", text, flags=re.MULTILINE)
        # piece[0] is whatever leads up to the first marker (header / summary)
        header = pieces[0]
        bodies = pieces[1:]
        return header, bodies

    left_header, left_bodies = read_pathways_section(left_file)
    _,           right_bodies = read_pathways_section(right_file)

    with open(out_file, "w", encoding="utf-8") as out:
        out.write(left_header)
        idx = 1
        for body in left_bodies + right_bodies:
            renumbered = re.sub(
                r"^pathway number \d+",
                f"pathway number {idx}",
                body,
                count=1,
                flags=re.MULTILINE,
            )
            out.write(renumbered)
            idx += 1
    return idx - 1


def main():
    print("=" * 70)
    print(" TAL-CENTERED COMBINED  —  3 networks, 2 branches from TAL")
    print("=" * 70)

    write_smi("talcentered_bio_retro_starter.smi",   [TAL])
    write_smi("talcentered_chem_fwd_starter.smi",    [TAL])
    write_smi("talcentered_chem_fwd_helpers.smi",    ["O", "[H][H]"])
    write_smi("talcentered_chem_retro_starter.smi",  [SORBIC_ACID])
    write_smi("talcentered_chem_retro_helpers.smi",  ["O", "[H][H]"])

    overall_t0 = time.time()

    # ----- 1. bio retro: TAL → acetyl-CoA ------------------------------
    print(f"\n[1/6] Bio retro:    TAL → acetyl-CoA  (polyketide chain reversed)")
    t0 = time.time()
    bio_retro_net = generate_network_tal(
        job_name=f"{JOB_NAME}_bio_retro",
        starters="talcentered_bio_retro_starter.smi",
        gen=BIO_GEN,
        direction="retro",
        max_rxn_thermo_change=15.0,
        allow_multiple_reactants="default",
        include_chem=False,
        include_bio=True,
        bio_allow_multiple_reactants=True,
        bio_whitelist=POLYKETIDE_WHITELIST,
        targets=ACETYL_COA,
    )
    print(f"      finished in {time.time()-t0:.1f}s  |  "
          f"{len(bio_retro_net.mols)} mols, {len(bio_retro_net.rxns)} reactions")

    # ----- 2. chem forward: TAL → halfway ------------------------------
    print(f"\n[2/6] Chem forward: TAL → (halfway toward sorbic acid)")
    t0 = time.time()
    chem_fwd_net = generate_network_tal(
        job_name=f"{JOB_NAME}_chem_fwd",
        starters="talcentered_chem_fwd_starter.smi",
        helpers="talcentered_chem_fwd_helpers.smi",
        gen=CHEM_FWD_GEN,
        direction="forward",
        max_rxn_thermo_change=15.0,
        max_atoms={"C": 6, "O": 3, "N": 0, "S": 0},
        max_molecular_weight=150,
        strategy="cartesian",
        include_chem=True,
        include_bio=False,
    )
    print(f"      finished in {time.time()-t0:.1f}s  |  "
          f"{len(chem_fwd_net.mols)} mols, {len(chem_fwd_net.rxns)} reactions")

    # ----- 3. chem retro: sorbic acid ← halfway ------------------------
    print(f"\n[3/6] Chem retro:   sorbic acid ← (halfway back toward TAL)")
    t0 = time.time()
    chem_retro_net = generate_network_tal(
        job_name=f"{JOB_NAME}_chem_retro",
        starters="talcentered_chem_retro_starter.smi",
        helpers="talcentered_chem_retro_helpers.smi",
        gen=CHEM_RETRO_GEN,
        direction="retro",
        max_rxn_thermo_change=15.0,
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
          f"{len(chem_retro_net.mols)} mols, {len(chem_retro_net.rxns)} reactions")

    # ----- 4. merge all 3 with TAL as starter --------------------------
    print(f"\n[4/6] Combining all 3 networks (pretreat_networks, starter=TAL)")
    t0 = time.time()
    pretreat_networks(
        networks=[bio_retro_net, chem_fwd_net, chem_retro_net],
        starters=[TAL],
        helpers=[MALONYL_COA, "O", "[H][H]"],
        total_generations=BIO_GEN + CHEM_FWD_GEN + CHEM_RETRO_GEN,
        job_name=JOB_NAME,
    )
    print(f"      combined in {time.time()-t0:.1f}s")

    # ----- 5. two pathway searches, both starting from TAL -------------
    print(f"\n[5/6] Pathway search A: TAL → acetyl-CoA")
    job_left = f"{JOB_NAME}_left"
    # pathway_finder needs its own pretreated file under job_left;
    # easiest is to re-pretreat into that job name (cheap — just copies).
    pretreat_networks(
        networks=[bio_retro_net, chem_fwd_net, chem_retro_net],
        starters=[TAL],
        helpers=[MALONYL_COA, "O", "[H][H]"],
        total_generations=BIO_GEN + CHEM_FWD_GEN + CHEM_RETRO_GEN,
        job_name=job_left,
    )
    pathway_finder(
        starters=[TAL],
        helpers=[MALONYL_COA, "O", "[H][H]"],
        target=[ACETYL_COA],
        search_depth=BIO_GEN,             # tight: only bio depth matters for this branch
        max_num_rxns=6,
        job_name=job_left,
    )

    print(f"\n      Pathway search B: TAL → sorbic acid")
    job_right = f"{JOB_NAME}_right"
    pretreat_networks(
        networks=[bio_retro_net, chem_fwd_net, chem_retro_net],
        starters=[TAL],
        helpers=[MALONYL_COA, "O", "[H][H]"],
        total_generations=BIO_GEN + CHEM_FWD_GEN + CHEM_RETRO_GEN,
        job_name=job_right,
    )
    pathway_finder(
        starters=[TAL],
        helpers=[MALONYL_COA, "O", "[H][H]"],
        target=[SORBIC_ACID],
        # match the proven-clean bidirectional test (~6 pathways).
        # search_depth = chem half-gens × 2  (the bio depth does not
        # contribute to chem-side reorderings on this branch).
        search_depth=CHEM_FWD_GEN + CHEM_RETRO_GEN,
        max_num_rxns=12,
        job_name=job_right,
    )

    # Merge the two pathway files
    left_file  = f"{job_left}_pathways.txt"
    right_file = f"{job_right}_pathways.txt"
    out_file   = f"{JOB_NAME}_pathways.txt"
    left_exists  = os.path.exists(left_file)
    right_exists = os.path.exists(right_file)
    if not left_exists and not right_exists:
        print("\nNo pathway files written for either branch.")
        return
    if left_exists and right_exists:
        n_merged = _merge_pathway_files(left_file, right_file, out_file)
        print(f"\n      merged: {n_merged} total pathways  →  {out_file}")
    else:
        single = left_file if left_exists else right_file
        with open(single, encoding="utf-8") as src, \
             open(out_file, "w", encoding="utf-8") as dst:
            dst.write(src.read())
        print(f"\n      only one branch produced pathways → using {single}")

    pathways = load_pathways_from_file(JOB_NAME)
    print(f"\n[SUCCESS] {len(pathways)} pathway(s) from TAL")
    for i, p in enumerate(pathways[:4], 1):
        print(f"\n  Pathway {i}  ({p.num_steps} steps):")
        for j, rxn in enumerate(p.reactions, 1):
            parsed = parse_reaction_string(rxn)
            short_reas = [r if len(r) < 40 else r[:37] + "..." for r in parsed["reactants"]]
            short_pros = [pp if len(pp) < 40 else pp[:37] + "..." for pp in parsed["products"]]
            arrow = " + ".join(short_reas) + "  ->  " + " + ".join(short_pros)
            print(f"    Step {j}. [{parsed['op_name']}]")
            print(f"             {arrow}")

    # ----- 6. visualize ------------------------------------------------
    # starter = TAL (centre), target = sorbic acid (one of the two
    # endpoints — visualize_pathways uses it for label positioning).
    # acetyl-CoA appears as the leaf of the other branch.
    print(f"\n[6/6] Rendering TAL-centered graph")
    html_path = visualize_pathways(
        job_name=JOB_NAME,
        starter_smiles=TAL,
        target_smiles=SORBIC_ACID,
        starter_label="TAL",
        target_label="sorbic acid",
        helpers=[MALONYL_COA, "O", "[H][H]"],
        pathway_filter="all",
    )
    print(f"      {html_path}")
    print(f"\nTotal runtime: {time.time() - overall_t0:.1f}s")


if __name__ == "__main__":
    main()
