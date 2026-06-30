"""
test_forward_polyketide.py
==========================

The CEO-scoped synthesis: acetyl-CoA + malonyl-CoA --> TAL via the
2-pyrone synthase polyketide pathway. Three enzyme steps.

WHY FORWARD AND NOT RETRO
-------------------------
Retro from TAL (test_retro_polyketide.py) is computationally
intractable on a laptop because:

  - bio rules with multi-substrate enabled try every combination of
    network molecules in their slots
  - polyketide SMARTS encode the ~50-atom CoA tail; every RunReactants
    call against these big patterns is ~50-100 ms in RDKit
  - retro starts with TAL (small) and grows outward into BIG
    intermediates — each new big molecule then has all 5 rules
    trying to fire on it, exploding the per-trial cost

Forward sidesteps all this by starting from the big molecules
(acetyl-CoA + malonyl-CoA) and applying highly selective polyketide
rules that only fire on specific patterns. Each generation adds 1-2
molecules. The whole search finishes in seconds.

Same chemistry, same operators, same answer — just a search direction
that matches how the chemistry actually proceeds in nature.

PATHWAY (forward, gen=3)
------------------------
  gen 1:  acetyl-CoA + malonyl-CoA  →  acetoacetyl-CoA + CO2 + CoA
                                       (Claisen condensation 1)
  gen 2:  acetoacetyl-CoA + malonyl-CoA  →  3,5-dioxohexanoyl-CoA + ...
                                       (Claisen condensation 2)
  gen 3:  3,5-dioxohexanoyl-CoA  →  TAL + CoA
                                       (intramolecular cyclization)
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


# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
TAL = "Cc1cc(O)cc(=O)o1"

# DORAnet's canonical acetyl-CoA SMILES (matches all_cofactors.tsv).
ACETYL_COA = (
    "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)

# Malonyl-CoA — same CoA tail with a carboxylated acetyl head.
# Not in DORAnet's cofactor table, so we supply it explicitly.
MALONYL_COA = (
    "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)

# Same 5-rule minimal whitelist as the retro test
MINIMAL_POLYKETIDE_WHITELIST = frozenset({
    "rule1118",   # Claisen 1
    "rule0087",   # Claisen 2
    "rule0126",   # Claisen 2 variant
    "rule0350",   # Claisen 2 variant
    "rule0891",   # cyclization
})

# Same 4-cofactor pool
POLYKETIDE_COFACTORS = ("CoA", "CO2", "WATER", "H+")

STARTER_FILE = "test_forward_polyketide_starter.smi"
HELPER_FILE  = "test_forward_polyketide_helpers.smi"
JOB_NAME = "test_forward_polyketide"


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def main():
    print("=" * 64)
    print(" FORWARD POLYKETIDE  —  acetyl-CoA + malonyl-CoA  →  TAL")
    print("=" * 64)
    print(f" Rules:      5 (Claisen 1 + Claisen 2 + cyclization)")
    print(f" Cofactors:  4 (CoA, CO2, water, H+)")
    print(f" Strategy:   cartesian")
    print(f" Multi-sub:  True (required for Claisen condensations)")
    print(f" Direction:  forward")

    # Both acetyl-CoA AND malonyl-CoA are starters — the search builds
    # from them.
    write_smi(STARTER_FILE, [ACETYL_COA, MALONYL_COA])
    write_smi(HELPER_FILE, ["O", "[H][H]"])

    t0 = time.time()
    network = generate_network_tal(
        job_name=JOB_NAME,
        starters=STARTER_FILE,
        helpers=HELPER_FILE,
        gen=3,
        direction="forward",
        molecule_thermo_calculator=None,
        max_rxn_thermo_change=15.0,
        max_atoms={"C": 35, "O": 20, "N": 8, "S": 2, "P": 5},
        max_molecular_weight=1100,
        allow_multiple_reactants="default",
        strategy="cartesian",
        min_carbons=0,
        include_chem=False,
        include_bio=True,
        bio_allow_multiple_reactants=True,
        bio_whitelist=MINIMAL_POLYKETIDE_WHITELIST,
        included_cofactors=POLYKETIDE_COFACTORS,
    )
    elapsed = time.time() - t0
    print(f"\nNetwork built in {elapsed:.1f}s")
    print(f"  {len(network.mols)} mols, {len(network.rxns)} reactions")

    # Was TAL reached?
    from rdkit import Chem
    target_canon = Chem.MolToSmiles(Chem.MolFromSmiles(TAL))
    network_smiles = {
        Chem.MolToSmiles(Chem.MolFromSmiles(m.smiles))
        for m in network.mols
        if Chem.MolFromSmiles(m.smiles) is not None
    }
    in_network = target_canon in network_smiles
    print(f"\nTAL in network? {in_network}")

    if not in_network:
        print("\n(TAL wasn't reached at gen=3 — likely a rule SMARTS"
              " mismatch. Check intermediate molecules in the network.)")
        return

    # Trace pathways: starter = acetyl-CoA, target = TAL
    print(f"\nTracing acetyl-CoA → TAL pathways...")
    find_pathways_to_target(
        network=network,
        starter=ACETYL_COA,
        target=TAL,
        helpers=["O", "[H][H]", MALONYL_COA],
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
        print("TAL in network but no pathway traced.")
        return

    print(f"\n[SUCCESS] Found {len(pathways)} pathway(s) acetyl-CoA → TAL")
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

    print(f"\nRendering interactive graph...")
    html_path = visualize_pathways(
        job_name=JOB_NAME,
        starter_smiles=ACETYL_COA,
        target_smiles=TAL,
        starter_label="acetyl-CoA",
        target_label="TAL",
        pathway_filter="all",
    )
    print(f"  {html_path}")
    print(f"\nTotal: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
