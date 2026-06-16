"""
test_tal_pipeline.py
End-to-end smoke test for the TAL pathway-explorer pipeline.

What this exercises
-------------------
1. Network generation: forward expansion from triacetic acid lactone
   (TAL, SMILES `CC1=CC(O)=CC(=O)O1`) using `generate_network_tal`.
2. Pathway finding: pretreat the network and search for pathways to a
   target picked from the generated network.
3. Pathway scoring: load pathways from the pathway file and score them
   with `WeightedPathwayScorer` + the default profile.

Run from project root:
    python test_tal_pipeline.py

The job_name "tal_test" controls all output file names so it's easy to
clean up between runs.
"""

import os
import sys
import shutil

# --- make src/ importable -------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(HERE, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

# All output files land next to this script.
os.chdir(HERE)

from network_generation import generate_network_tal
from pathway_tools import find_pathways_to_target, load_pathways_from_file
from pathway_scoring import (
    WeightedPathwayScorer,
    score_pathways_from_file,
    StepsCriterion,
    ThermoCriterion,
    IntermediateStabilityCriterion,
    ProcedureDiversityCriterion,
    ChemBioSwitchCriterion,
)


JOB_NAME = "tal_test"
STARTER = "CC1=CC(O)=CC(=O)O1"            # Triacetic acid lactone
HELPERS = ["[H][H]", "O", "CCO"]          # Hydrogen, water, ethanol
TARGET = "CCc1ccc(C)oc1=O"                # A known reachable molecule
                                          # (same target the notebook used)
GENERATIONS = 3
MAX_ATOMS = {
    "C": 25, "O": 4, "N": 2,
    "S": 0, "P": 0, "F": 0, "Cl": 0, "Br": 0, "I": 0,
}


# =====================================================================
# >>> TUNE PATHWAY-SCORING WEIGHTS HERE <<<
# =====================================================================
# Each criterion contributes  (weight_i * score_i) / sum(weights)  to the
# final ranking, so absolute magnitudes don't matter — only the RATIOS.
# Set a weight to 0 to disable a criterion.
#
# Examples:
#   steps=4, thermo=2  → DORAnet's default; favor short routes
#   steps=1, thermo=5  → "I care most about thermodynamic feasibility"
#   steps=0, thermo=1  → "Ignore step count; rank purely on thermo"
#   steps=3, thermo=3  → equal weight
#
# Add more (criterion, weight) tuples to SCORING_WEIGHTS as new criteria
# come online (AtomEconomyCriterion, CostCriterion, etc.).
# =====================================================================
SCORING_WEIGHTS = [
    (StepsCriterion(),  4.0),
    (ThermoCriterion(), 2.0),
    # `excluded_smiles` skips starter/helpers/target — we already trust
    # those, and assessing them would just drag the pathway score down
    # for things we've chosen to use anyway.
    (IntermediateStabilityCriterion(
        excluded_smiles=[STARTER, TARGET] + HELPERS,
    ), 3.0),
    # ProcedureDiversityCriterion: penalizes pathways that mix many distinct
    # operator families (each new procedure ≈ new optimization burden).
    # Source: Anderson, *Practical Process R&D*, 2nd ed., Ch. 2.
    (ProcedureDiversityCriterion(), 2.0),
    # ChemBioSwitchCriterion: penalizes chem⇄bio regime switching (different
    # solvents, equipment, vendors). Uses DORAnet's `bio_rxn_names`
    # (3604 enzymatic rules from JN3604IMT_rules.tsv) as ground truth.
    # Source: Sheldon & Woodley, *Chem. Rev.* 118, 801 (2018).
    (ChemBioSwitchCriterion(), 2.0),
]


def banner(text):
    print("\n" + "=" * 70)
    print(text)
    print("=" * 70)


def step_1_generate_network():
    banner("STEP 1: NETWORK GENERATION")
    network = generate_network_tal(
        job_name=JOB_NAME,
        starters=[STARTER],
        helpers=HELPERS,
        max_atoms=MAX_ATOMS,
        gen=GENERATIONS,
        direction="forward",
        min_carbons=0,
    )
    print(f"\nGenerated {len(network.mols)} molecules, {len(network.rxns)} reactions.")
    return network


def step_2_find_pathways(network):
    banner(f"STEP 2: PATHWAY SEARCH (target = {TARGET})")
    find_pathways_to_target(
        network=network,
        starter=STARTER,
        target=TARGET,
        helpers=HELPERS,
        generations=GENERATIONS,
        max_num_rxns=20,
        job_name=JOB_NAME,
    )
    # Sanity-check the artifact exists and is non-empty.
    pw_file = f"{JOB_NAME}_pathways.txt"
    assert os.path.exists(pw_file), f"Expected pathway file {pw_file} not found"
    size = os.path.getsize(pw_file)
    print(f"\n{pw_file} written ({size} bytes).")


def step_3_score_pathways():
    banner("STEP 3: PATHWAY SCORING")
    scorer = WeightedPathwayScorer(SCORING_WEIGHTS)
    print("Weights in use:")
    for criterion, weight in SCORING_WEIGHTS:
        print(f"  {criterion.name:8s} = {weight}")
    scored = score_pathways_from_file(JOB_NAME, scorer)
    print(f"\nLoaded and scored {len(scored)} pathways.")
    for i, sp in enumerate(scored, 1):
        print(f"\n--- ranking {i} ---")
        print(f"  final score   : {sp.final_score:.4f}")
        print(f"  components    : {sp.components}")
        print(f"  num_steps     : {sp.pathway.num_steps}")
        print(f"  reactions:")
        for rxn in sp.pathway.reactions:
            print(f"    {rxn}")
    return scored


def main():
    network = step_1_generate_network()
    step_2_find_pathways(network)
    scored = step_3_score_pathways()
    banner("ALL STEPS COMPLETED")
    print(f"Pathways scored: {len(scored)}")
    if scored:
        top = scored[0]
        print(f"Best pathway: {top.pathway.num_steps} steps, score {top.final_score:.4f}")


if __name__ == "__main__":
    main()
