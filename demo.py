"""
TAL Pathway Explorer — Comprehensive Demo

What this demonstrates
----------------------
Everything currently working end-to-end on this branch:

  1. ARCHITECTURE
       - Chem operators curated for TAL chemistry
       - Bio operators carrying enzyme catalysis data
       - 5 pluggable pathway-ranking criteria

  2. OPEN EXPLORATION (the brief's "what products from TAL?")
       - Forward expansion with no preset target
       - Endpoint scoring by interestingness
       - Cross-check against literature TAL derivatives

  3. DIRECTED SEARCH + RANKING
       - Use the top exploration endpoint as a target
       - Find pathways and rank with all 5 criteria
       - Show score breakdown per criterion

  4. WHAT'S NEXT
       - Roadmap framing

How to run
----------
From the repo root:

    python demo.py

Runs in well under a minute on a laptop.
"""

import glob
import os
import sys
import time
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

# Quiet RDKit's valence warnings and any deprecation chatter
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

import pandas as pd

from pathway_tools import explore_downstream
from pathway_scoring import (
    StepsCriterion,
    ThermoCriterion,
    IntermediateStabilityCriterion,
    ProcedureDiversityCriterion,
    ChemBioSwitchCriterion,
    WeightedPathwayScorer,
    score_pathways_from_file,
)
from tal_downstream_derivatives import TAL_DOWNSTREAM_DERIVATIVES
from tal_reaction_whitelist import TAL_REACTION_WHITELIST
from tal_bio_reaction_whitelist import TAL_BIO_REACTION_WHITELIST


JOB_NAME = "tal_demo"
STARTER_FILE = "demo_starter.smi"
HELPER_FILE = "demo_helpers.smi"

TAL_SMILES = "Cc1cc(O)cc(=O)o1"
HELPERS = ["O", "[H][H]"]

GENERATIONS = 2
TOP_N = 10
MAX_NUM_RXNS = 3


def banner(text, char="="):
    print()
    print(char * 70)
    print(f" {text}")
    print(char * 70)


def write_smiles_file(path, smiles_list):
    with open(path, "w") as f:
        for s in smiles_list:
            f.write(s + "\n")


# =====================================================================
# Scene 1 — architecture: operators + scoring framework
# =====================================================================
def scene_1_architecture():
    banner("[1/4]  ARCHITECTURE  —  what the tool knows about")

    print("\nReaction operators (the universe of moves the tool can make):")
    print(f"  Chem operators (curated TAL whitelist):  {len(TAL_REACTION_WHITELIST):>4}")
    print(f"  Bio operators  (TAL-filtered JN1224MIN): {len(TAL_BIO_REACTION_WHITELIST):>4}")
    print(f"  ── Total: {len(TAL_REACTION_WHITELIST) + len(TAL_BIO_REACTION_WHITELIST)} operators")

    # Show that every bio operator carries enzyme catalysis info.
    print("\nEvery bio operator carries enzyme catalysis data.")
    print("Sampling one rule from the JN1224MIN ruleset...")
    try:
        from network_generation import AVAILABLE_RULESETS
        rules_df = pd.read_csv(AVAILABLE_RULESETS["JN1224MIN"], sep="\t")
        # Pick the first rule in the bio whitelist that we can find in the DF
        sample_name = next(iter(TAL_BIO_REACTION_WHITELIST))
        row = rules_df[rules_df["Name"] == sample_name]
        if not row.empty:
            comments = row.iloc[0]["Comments"]
            enzymes = [e.strip() for e in str(comments).split(";") if e.strip()]
            reactants = row.iloc[0]["Reactants"]
            print(f"  Rule:       {sample_name}")
            print(f"  Reactants:  {reactants}")
            print(f"  Enzymes:    {len(enzymes)} UniProt IDs — "
                  f"{', '.join(enzymes[:5])}"
                  f"{'...' if len(enzymes) > 5 else ''}")
    except Exception as exc:
        print(f"  (could not load ruleset sample: {exc})")

    print("\nPathway-ranking criteria (all pluggable, all weight-tunable):")
    for cls in (
        StepsCriterion,
        ThermoCriterion,
        IntermediateStabilityCriterion,
        ProcedureDiversityCriterion,
        ChemBioSwitchCriterion,
    ):
        print(f"  - {cls.__name__}")


# =====================================================================
# Scene 2 — open exploration (the brief's downstream-applications Q)
# =====================================================================
def scene_2_exploration():
    banner("[2/4]  OPEN EXPLORATION  —  'what can TAL become?'")
    print()
    print("Brief: \"What products can feasibly be produced with TAL as")
    print(" the starting point?\"")
    print()
    print(f"Settings: starter=TAL, gen={GENERATIONS}, "
          f"chem-only (bio off for demo speed), top_n={TOP_N}")

    write_smiles_file(STARTER_FILE, [TAL_SMILES])
    write_smiles_file(HELPER_FILE, HELPERS)

    t0 = time.time()
    result = explore_downstream(
        starter=STARTER_FILE,
        helpers=HELPER_FILE,
        job_name=JOB_NAME,
        top_n=TOP_N,
        max_num_rxns=MAX_NUM_RXNS,
        network_kwargs=dict(
            gen=GENERATIONS,
            direction="forward",
            molecule_thermo_calculator=None,
            max_rxn_thermo_change=15.0,
            max_atoms={"C": 50, "O": 8, "N": 0, "S": 0},
            max_molecular_weight=500,
            allow_multiple_reactants="default",
            strategy="cartesian",
            min_carbons=0,
            include_chem=True,
            include_bio=False,
        ),
        scoring_kwargs=dict(carbon_window=(3, 12), require_oxygen=True),
        derivatives_list=TAL_DOWNSTREAM_DERIVATIVES,
    )
    elapsed = time.time() - t0

    print()
    print(f"Exploration finished in {elapsed:.1f}s.")
    print(f"  Top endpoints surfaced:  {len(result['top_endpoints'])}")
    print(f"  Endpoints with paths:    "
          f"{sum(1 for v in result['pathways'].values() if v)}/"
          f"{len(result['top_endpoints'])}")
    if result["derivative_matches"] is not None:
        n_match = len(result["derivative_matches"])
        n_total = n_match + len(result["derivative_missing"])
        print(f"  Literature derivatives:  {n_match}/{n_total} found")
        print(f"  (0/13 at gen=2 is expected — literature TAL→sorbic")
        print(f"   acid is 3 steps; the cross-check surfaces that gap.)")

    print("\nTop 5 endpoints by interestingness score:")
    for i, es in enumerate(result["top_endpoints"][:5], 1):
        n_paths = len(result["pathways"].get(es.smiles, []))
        print(f"  {i:2d}. {es.smiles}")
        print(f"      score={es.score:.3f}  C={es.carbons}  "
              f"Bertz={es.bertz:.0f}  FGs={es.n_functional_groups}  "
              f"pathways={n_paths}")

    return result


# =====================================================================
# Scene 3 — directed search + full 5-criteria ranking
# =====================================================================
def scene_3_ranking(exploration_result):
    banner("[3/4]  RANKING  —  full 5-criteria scoring on a chosen target")

    top = exploration_result["top_endpoints"]
    if not top:
        print("\n  (no endpoints from exploration to rank — skipping)")
        return
    target = top[0]
    target_job = f"{JOB_NAME}_ep01"

    print()
    print(f"Target (top exploration hit): {target.smiles}")
    print()
    print("Pathway pool already built in Scene 2 — now we apply the")
    print("full 5-criteria scorer with the default weights.")

    scorer = WeightedPathwayScorer([
        (StepsCriterion(),                    4.0),
        (ThermoCriterion(),                   2.0),
        (IntermediateStabilityCriterion(),    2.0),
        (ProcedureDiversityCriterion(),       2.0),
        (ChemBioSwitchCriterion(),            2.0),
    ])

    try:
        scored = score_pathways_from_file(target_job, scorer)
    except FileNotFoundError:
        print(f"\n  (pathway file {target_job}_pathways.txt not found — "
              f"Scene 2 must have failed)")
        return

    print(f"\nScored {len(scored)} pathway(s). Top 3:")
    for i, sp in enumerate(scored[:3], 1):
        print(f"\n  Pathway {i}  —  final score {sp.final_score:.3f}  "
              f"({sp.pathway.num_steps} steps)")
        for name, val in sp.components.items():
            print(f"    {name:32s} = {val:.3f}")


# =====================================================================
# Scene 4 — what's next
# =====================================================================
def scene_4_next():
    banner("[4/4]  WHAT'S NEXT")
    print("""
Near-term:
  1. Bio combinatorial fix
       - Default to priority-queue strategy when bio is enabled
       - Use TAL-relevant cofactor subset (~10 instead of 41)
       - Will let bio enzymes run end-to-end in seconds

  2. Retrosynthesis mode (feedstocks → TAL)
       - New ranker: FeedstockProximityRanker
       - Same pathway-finding primitive, flipped direction

  3. Streamlit UI
       - SMILES input + weight sliders + ranked output
       - Lets non-developers use the tool

Medium-term:
  4. Thermodynamics
       - Wire eQuilibrator (works for cataloged molecules today)
       - Decide on ChemAxon access vs RMG vs Joback fallback

  5. Domain-specific criterion packs
       - Agrochemical / cosmetic / polymer-grade filters
       - "Amazon-style search filters" with weight sliders
""")


# =====================================================================
def cleanup_artifacts():
    """
    Remove DORAnet's per-job intermediate files so the project root
    stays uncluttered. Keeps the consolidated exploration report.
    """
    patterns = [
        f"{JOB_NAME}_ep*_pathways.txt",
        f"{JOB_NAME}_ep*_network_pretreated.json",
        f"{JOB_NAME}_ep*_reaxys_batch_*.txt",
        f"{JOB_NAME}_ep*_reaxys_batch_*.csv",
        f"{JOB_NAME}_*_saved_network*",
        f"{JOB_NAME}_network_pretreated.json",
        f"{JOB_NAME}_pathways.txt",
        f"{JOB_NAME}_reaxys_batch_*.txt",
        f"{JOB_NAME}_reaxys_batch_*.csv",
    ]
    n_removed = 0
    for pat in patterns:
        for path in glob.glob(pat):
            try:
                os.remove(path)
                n_removed += 1
            except OSError:
                pass
    return n_removed


def main():
    banner("TAL PATHWAY EXPLORER  —  COMPREHENSIVE DEMO", char="#")
    t_total = time.time()
    scene_1_architecture()
    exploration_result = scene_2_exploration()
    scene_3_ranking(exploration_result)
    scene_4_next()

    n_cleaned = cleanup_artifacts()

    banner(f"DEMO COMPLETE  ({time.time() - t_total:.1f}s total)", char="#")
    print(f"\nFull markdown report: {JOB_NAME}_exploration_report.md")
    print("Open in any markdown viewer (VS Code: Ctrl+Shift+V).")
    if n_cleaned:
        print(f"\n[housekeeping] Cleaned up {n_cleaned} intermediate artifact files.\n")
    else:
        print()


if __name__ == "__main__":
    main()
