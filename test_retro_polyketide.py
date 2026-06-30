"""
test_retro_polyketide.py
========================

The CEO-scoped retro bio search: TAL <- acetyl-CoA via the polyketide
synthase pathway. Three enzyme steps in the literature mechanism
(two Claisen condensations + a cyclization).

Four levers, all together, are what make this tractable:

  1. TAL_BIO_POLYKETIDE_WHITELIST — only 7 rules instead of 348.
     Cuts recipe generation cost ~50x relative to the standard bio
     whitelist.

  2. included_cofactors — only the 8 cofactors that actually appear
     in polyketide chemistry, instead of DORAnet's full 41. Each
     "Any" cofactor slot now has 8 candidates instead of 41, so a
     3-slot rule generates 8**3 = 512 candidate recipes instead of
     41**3 = ~68,000. Roughly 130x cheaper per rule per generation.

  3. bio_allow_multiple_reactants=True — required because the
     Claisen condensations consume two real substrates per step
     (acetyl-CoA + malonyl-CoA, then acetoacetyl-CoA + malonyl-CoA).
     Without this, the polyketide rules silently filter out and the
     route is structurally unreachable.

  4. priority_queue + FeedstockProximityRanker(acetyl-CoA) — beam
     search bounded by acetyl-CoA similarity. Keeps the network from
     compounding across generations.

OUTPUT
------
  test_retro_polyketide_pathways.txt  — raw pathway file from DORAnet
  test_retro_polyketide_graph.html    — interactive DAG (one HTML file
                                          we can open in a browser)
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
from visualize_pathways import visualize_pathways


# Minimal whitelist for the TAL <- acetyl-CoA retro chain ONLY. Drops
# the acetyl-CoA carboxylase rules (rule0023, rule0730) because they
# describe acetyl-CoA -> malonyl-CoA, which isn't on the path we
# want to walk back: we just need the two Claisens + the cyclization.
MINIMAL_POLYKETIDE_WHITELIST = frozenset({
    "rule1118",   # malonyl-CoA -> acetoacetyl-CoA           (Claisen 1)
    "rule0087",   # acetoacetyl-CoA -> 3,5-dioxohexanoyl-CoA (Claisen 2)
    "rule0891",   # 3,5-dioxohexanoyl-CoA -> TAL             (cyclization)
})
# Dropped rule0126 and rule0350 (additional Claisen 2 variants) — they
# overlap with rule0087 and explode the search combinatorially when
# multi-substrate is enabled. The 3-rule chain is sufficient for the
# acetyl-CoA <- TAL retro walk.


# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
TAL = "Cc1cc(O)cc(=O)o1"

# DORAnet's canonical acetyl-CoA SMILES (matches the entry in
# all_cofactors.tsv). Stereochemistry-explicit form.
ACETYL_COA = (
    "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)

# Minimal cofactor pool — only the 4 cofactors the two Claisens
# and the cyclization actually use. The ATP/ADP/Pi cofactors went
# away with the carboxylase rules.
POLYKETIDE_COFACTORS = (
    "CoA",         # released by Claisens and cyclization
    "CO2",         # released by Claisens
    "WATER",       # generic
    "H+",          # generic proton donor/acceptor
)

STARTER_FILE = "test_retro_polyketide_starter.smi"
HELPER_FILE  = "test_retro_polyketide_helpers.smi"
JOB_NAME = "test_retro_polyketide"


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" RETRO POLYKETIDE  —  TAL  <--  acetyl-CoA")
    print("=" * 64)
    print(f" Rules:      {len(MINIMAL_POLYKETIDE_WHITELIST)} (Claisen 1 + Claisen 2 + cyclization)")
    print(f" Cofactors:  {len(POLYKETIDE_COFACTORS)} (CoA, CO2, water, H+)")
    print(f" Strategy:   cartesian (no ranker overhead)")
    print(f" Multi-sub:  True (required for Claisen condensations)")

    # In retro mode, the starter is what we walk backward FROM (TAL).
    # The "target" we want to discover is acetyl-CoA, supplied as the
    # priority-queue's ranker target.
    write_smi(STARTER_FILE, [TAL])
    write_smi(HELPER_FILE, ["O", "[H][H]"])

    t0 = time.time()
    network = generate_network_tal(
        job_name=JOB_NAME,
        starters=STARTER_FILE,
        helpers=HELPER_FILE,
        gen=3,
        direction="retro",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        # Atoms are wide because acetyl-CoA / malonyl-CoA are big
        # (~50 atoms). The polyketide intermediates fit well within
        # these bounds.
        max_atoms={"C": 35, "O": 20, "N": 8, "S": 2, "P": 5},
        max_molecular_weight=1100,
        allow_multiple_reactants="default",
        # Cartesian instead of priority_queue — with only 5 rules the
        # combinatorial fan-out is small enough that "try everything"
        # is faster than running the FeedstockProximityRanker on every
        # candidate (ranker is expensive on big CoA-containing mols).
        strategy="cartesian",
        min_carbons=0,
        include_chem=False,
        include_bio=True,
        bio_allow_multiple_reactants=True,
        bio_whitelist=MINIMAL_POLYKETIDE_WHITELIST,
        # included_cofactors NOT passed — it silently drops whitelisted
        # rules whose Reactants/Products columns reference cofactor
        # codes outside the included set. Full DORAnet cofactor pool
        # is fine and finishes in seconds.
    )
    elapsed_exp = time.time() - t0
    print(f"\nNetwork built in {elapsed_exp:.1f}s")
    print(f"  {len(network.mols)} mols, {len(network.rxns)} reactions")

    # Was acetyl-CoA reached?
    from rdkit import Chem
    target_canon = Chem.MolToSmiles(Chem.MolFromSmiles(ACETYL_COA))
    network_smiles = {
        Chem.MolToSmiles(Chem.MolFromSmiles(m.smiles))
        for m in network.mols
        if Chem.MolFromSmiles(m.smiles) is not None
    }
    in_network = target_canon in network_smiles
    print(f"\nacetyl-CoA in network?  {in_network}")

    if not in_network:
        print("\n(acetyl-CoA wasn't reached — try bumping gen or beam_size)")
        return

    # Trace pathways
    print(f"\nTracing TAL <- acetyl-CoA pathways...")
    find_pathways_to_target(
        network=network,
        starter=TAL,
        target=ACETYL_COA,
        helpers=["O", "[H][H]"],
        generations=3,
        max_num_rxns=6,
        job_name=JOB_NAME,
    )
    try:
        pathways = load_pathways_from_file(JOB_NAME)
    except FileNotFoundError:
        print("No pathway file written.")
        return

    if not pathways:
        print("acetyl-CoA in network but no pathway traced.")
        return

    print(f"\n[SUCCESS] Found {len(pathways)} pathway(s) TAL <- acetyl-CoA")
    from pathway_tools import parse_reaction_string
    for i, p in enumerate(pathways[:5], 1):
        print(f"\n  Pathway {i}  ({p.num_steps} enzyme steps):")
        for j, rxn in enumerate(p.reactions, 1):
            parsed = parse_reaction_string(rxn)
            short_reas = [r if len(r) < 30 else r[:27] + "..." for r in parsed["reactants"]]
            short_pros = [p_ if len(p_) < 30 else p_[:27] + "..." for p_ in parsed["products"]]
            arrow = " + ".join(short_reas) + "  ->  " + " + ".join(short_pros)
            print(f"    Step {j}. [{parsed['op_name']}]")
            print(f"             {arrow}")

    # Visualize
    print(f"\nRendering interactive graph...")
    html_path = visualize_pathways(
        job_name=JOB_NAME,
        starter_smiles=TAL,
        target_smiles=ACETYL_COA,
        starter_label="TAL",
        target_label="acetyl-CoA",
        pathway_filter="all",
    )
    print(f"  {html_path}")
    print(f"\nTotal: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
