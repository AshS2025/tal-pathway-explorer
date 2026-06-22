"""
ONE-SHOT TAL -> sorbic acid expansion. Heavy compute — intended for
Google Cloud (or any machine with > 16 GB RAM and a few cores).

Why this exists
---------------
On a laptop, the full TAL -> sorbic acid route (6 elementary steps
through DORAnet's SMARTS) consistently runs the process to 700-800 MB
before the user has to kill it. We've verified the chemistry works
end-to-end by splitting into TAL -> HMP (3 steps) and PSA -> sorbic
(2 steps) on the laptop. This script is the single-network version
for a beefier machine.

What it does
------------
Forward chem-only expansion starting from TAL, with the atom budget
sized exactly for the literature route. The expected literature chain
(translated to DORAnet operators):

  1. TAL  --[Keto-enol Tautomerization Reverse]--> non-aromatic TAL
  2.      --[Hydrogenation of C=C]--> partially-reduced TAL
  3.      --[Reduction of ketone -> alcohol]--> HMP
  4. HMP  --[Dehydration of Alcohol]--> PSA (parasorbic acid)
  5. PSA  --[Hydrolysis of Esters, Intramolecular + H2O]-->
                                          open-chain hydroxy-hexenoate
  6.      --[Dehydration of Alcohol]--> sorbic acid

Two strategies — pick one
-------------------------
This file has both knobs ready. Default is CARTESIAN (exhaustive,
needs RAM) which is the right thing on cloud. To swap, change
STRATEGY below to "priority_queue".

Cloud resource notes
--------------------
- Expected RAM: 4-16 GB (cartesian at gen=8 with our tight atom cap)
- Expected runtime: 5-30 minutes (varies with machine)
- Output: stdout + a markdown summary at tal_to_sorbic_full_report.md
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


STARTER_FILE = "test_tal_sorbic_starter.smi"
HELPER_FILE  = "test_tal_sorbic_helpers.smi"

TAL         = "Cc1cc(O)cc(=O)o1"
SORBIC_ACID = "CC=CC=CC(=O)O"     # no-stereo form (matches what
                                  # DORAnet's Dehydration operator
                                  # actually produces — the E,E
                                  # canonical SMILES silently mismatches)

# Strategy knob — change to "priority_queue" for laptop runs.
STRATEGY = "cartesian"
GENERATIONS = 8


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" TAL -> SORBIC ACID  —  full single-network forward expansion")
    print(" (heavy compute — intended for cloud)")
    print("=" * 64)
    print(f" strategy = {STRATEGY}, gen = {GENERATIONS}")
    print(f" target   = {SORBIC_ACID}")

    write_smi(STARTER_FILE, [TAL])
    write_smi(HELPER_FILE,  ["O", "[H][H]"])

    t0 = time.time()

    # Common settings shared by both strategies
    common = dict(
        job_name="tal_to_sorbic_full",
        starters=STARTER_FILE,
        helpers=HELPER_FILE,
        gen=GENERATIONS,
        direction="forward",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        # Atom budget sized EXACTLY for the literature route.
        # TAL=C6/O3, intermediates ≤ C6/O3, sorbic=C6/O2.
        # The C=6 cap blocks Diels-Alder dimers (would need C=12).
        max_atoms={"C": 6, "O": 3, "N": 0, "S": 0},
        max_molecular_weight=150,
        min_carbons=0,
        include_chem=True,
        include_bio=False,
    )

    if STRATEGY == "cartesian":
        network = generate_network_tal(strategy="cartesian", **common)
    elif STRATEGY == "priority_queue":
        network = generate_network_tal(
            strategy="priority_queue",
            targets=SORBIC_ACID,
            recipe_ranker=ForwardProductTanimotoRanker(SORBIC_ACID),
            beam_size=500,
            **common,
        )
    else:
        raise ValueError(f"Unknown strategy: {STRATEGY!r}")

    elapsed_expand = time.time() - t0
    print(f"\nNetwork expansion finished in {elapsed_expand:.1f}s.")

    # Check intermediates
    from rdkit import Chem
    network_smiles = {
        Chem.MolToSmiles(Chem.MolFromSmiles(m.smiles))
        for m in network.mols
        if Chem.MolFromSmiles(m.smiles) is not None
    }
    target_canon = Chem.MolToSmiles(Chem.MolFromSmiles(SORBIC_ACID))
    in_network = target_canon in network_smiles

    print(f"\nSorbic acid canonical: {target_canon}")
    print(f"In network?            {in_network}")
    print(f"Network has {len(network_smiles)} molecules total.")

    print("\nLiterature-path intermediates check:")
    intermediates = [
        ("TAL (starter)",                "Cc1cc(O)cc(=O)o1"),
        ("Non-aromatic TAL (diketo)",   "CC1=CC(=O)CC(=O)O1"),
        ("HMP (saturated lactone)",      "CC1CC(O)CC(=O)O1"),
        ("PSA (parasorbic acid)",        "CC1CC=CC(=O)O1"),
        ("Open-chain hydroxy-hexenoate", "CC(O)CC=CC(=O)O"),
        ("Sorbic acid (no-stereo)",      SORBIC_ACID),
    ]
    for name, smi in intermediates:
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        mark = "FOUND" if canon in network_smiles else "miss "
        print(f"  [{mark}] {name:35s} {canon}")

    if not in_network:
        print("\nSorbic acid not reached at gen=" + str(GENERATIONS) + ".")
        print("First missing intermediate (above) shows where the route")
        print("breaks down. If only sorbic acid itself is missing but the")
        print("open-chain hexenoate is present, the final dehydration step")
        print("was pruned by the search; bump GENERATIONS or beam_size.")
        return

    # Trace the pathway
    print(f"\nTracing TAL -> sorbic acid pathway(s)...")
    find_pathways_to_target(
        network=network,
        starter=TAL,
        target=SORBIC_ACID,
        helpers=["O", "[H][H]"],
        generations=GENERATIONS,
        max_num_rxns=15,
        job_name="tal_to_sorbic_full",
    )
    try:
        pathways = load_pathways_from_file("tal_to_sorbic_full")
    except FileNotFoundError:
        print("\nNo pathway file written.")
        return
    if not pathways:
        print("\nSorbic acid in network but no traceable pathway.")
        return

    print(f"\n[SUCCESS] Found {len(pathways)} pathway(s) TAL -> sorbic acid.")
    from pathway_tools import parse_reaction_string

    # Also write a markdown report for sharing
    report_lines = [
        "# TAL -> Sorbic Acid Pathway Report",
        "",
        f"Found **{len(pathways)} pathway(s)** TAL -> sorbic acid via "
        f"{STRATEGY} expansion, gen={GENERATIONS}.",
        "",
    ]
    for i, p in enumerate(pathways[:5], 1):
        print(f"\n  Pathway {i}  ({p.num_steps} steps):")
        report_lines.append(f"## Pathway {i} ({p.num_steps} steps)")
        for j, rxn in enumerate(p.reactions, 1):
            parsed = parse_reaction_string(rxn)
            arrow = (
                " + ".join(parsed["reactants"]) + "  ->  "
                + " + ".join(parsed["products"])
            )
            print(f"    Step {j}. [{parsed['op_name']}]")
            print(f"             {arrow}")
            report_lines.append(f"{j}. `{parsed['op_name']}` — {arrow}")
        report_lines.append("")

    with open("tal_to_sorbic_full_report.md", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))
    print("\nReport written to tal_to_sorbic_full_report.md")
    print(f"\nTotal runtime: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
