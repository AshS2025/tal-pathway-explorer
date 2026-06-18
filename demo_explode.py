"""
TAL Pathway Explorer — EXPLODE Demo

What this is for
----------------
Same pipeline as demo.py — architecture, open exploration, directed
search with full 5+2 criteria ranking, PDF output — but the search
space is cranked up to a scale that local hardware will NOT handle
gracefully. This file exists to *demonstrate the scaling limit*, not
to actually finish. Talking-point material when your boss asks "what
would it take to run this at production scale?"

What's been turned up vs demo.py
--------------------------------
- gen:                  3   -> 5    (cartesian: ~exponential growth)
- helpers:              3   -> 6    (water, H2, ethanol, methanol, propanol, acetic acid)
- max_atoms (C):        25  -> 80
- max_atoms (O):        4   -> 15
- max_molecular_weight: 500 -> 1000
- top_n (exploration):  10  -> 50
- include_bio:          False -> True (loads 348 bio ops + 41 cofactors)
- bio_allow_multiple_reactants: True (no single-substrate gate on bio)
- multiple precoded targets for directed search

Why this should crush a laptop
------------------------------
Cartesian expansion attempts EVERY (operator x molecule combination)
each generation. At gen=3 with chem-only and 3 helpers, demo.py
generates ~150 molecules. Here:

  ~500 operators x ~50 starting molecules x 5 generations
  with each gen ~10x molecule count growth
  plus 41 cofactors creating cofactor x cofactor combinatorial soup

Expected behavior: memory pressure within minutes, full crash or
unbounded runtime within 10-15 minutes on a laptop. If you ever
want a 'runs in finite time' version with bio enabled, switch
strategy to 'priority_queue' (beam search) — which is the actual
scaling fix.

How to run
----------
From the repo root:

    python demo_explode.py

You almost certainly want to kill it with Ctrl+C after a few minutes
rather than wait. The point is the scaling story, not the result.
"""

import glob
import os
import sys
import time
import warnings

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")

import pandas as pd

from pathway_tools import (
    explore_downstream,
    find_pathways_to_target,
    load_pathways_from_file,
)
from pathway_scoring import (
    StepsCriterion,
    ThermoCriterion,
    IntermediateStabilityCriterion,
    ProcedureDiversityCriterion,
    ChemBioSwitchCriterion,
    AtomEconomyCriterion,
    ByProductCountCriterion,
    WeightedPathwayScorer,
    parse_doranet_ranked_file,
    rewrite_ranked_file_with_unified_order,
)
from config import DEFAULT_WEIGHTS, build_unified_profile
from tal_downstream_derivatives import TAL_DOWNSTREAM_DERIVATIVES
from tal_reaction_whitelist import TAL_REACTION_WHITELIST
from tal_bio_reaction_whitelist import TAL_BIO_REACTION_WHITELIST

from doranet.modules.post_processing.post_processing import (
    pathway_ranking,
    pathway_visualization,
)


JOB_NAME = "tal_explode"
STARTER_FILE = "explode_starter.smi"
HELPER_FILE = "explode_helpers.smi"

TAL_SMILES = "Cc1cc(O)cc(=O)o1"
# Bigger helper set — more starting chemistry, more combinatorial spread
HELPERS = ["O", "[H][H]", "CCO", "CO", "CCCO", "CC(=O)O"]

# Cranked parameters that local hardware will struggle with
GENERATIONS = 5
TOP_N = 50
MAX_NUM_RXNS = 20

# Multiple precoded targets — each runs its own directed search
DIRECTED_TARGETS = [
    ("Sorbic acid (literature TAL derivative)",          "C/C=C/C=C/C(=O)O"),
    ("Phloroglucinol (literature TAL derivative)",       "Oc1cc(O)cc(O)c1"),
    ("Acetylacetone (literature ring-open/decarbox)",    "CC(=O)CC(=O)C"),
    ("Ethyl-methyl-pyranone (smaller cross-check)",      "CCc1ccc(C)oc1=O"),
]


def banner(text, char="="):
    print()
    print(char * 70)
    print(f" {text}")
    print(char * 70)


def write_smi(path, smiles_list):
    with open(path, "w") as f:
        for s in smiles_list:
            f.write(s + "\n")


def cleanup_artifacts():
    patterns = [
        f"{JOB_NAME}_*pathways.txt",
        f"{JOB_NAME}_*pretreated.json",
        f"{JOB_NAME}_*reaxys_batch_*",
        f"{JOB_NAME}_*_saved_network*",
        f"{JOB_NAME}_*ranked_pathways*.txt",
        "explode_*.smi",
    ]
    n = 0
    for pat in patterns:
        for p in glob.glob(pat):
            try:
                os.remove(p)
                n += 1
            except OSError:
                pass
    return n


def scene_1_architecture():
    banner("[1/4]  ARCHITECTURE  —  what the tool knows about (chem + BIO)")
    print()
    print(f"  Chem operators (curated TAL whitelist):  {len(TAL_REACTION_WHITELIST):>4}")
    print(f"  Bio operators  (TAL-filtered JN1224MIN): {len(TAL_BIO_REACTION_WHITELIST):>4}")
    print(f"  Bio cofactors auto-loaded as coreactants:  ~41")
    print(f"  ── Total ops:  ~{len(TAL_REACTION_WHITELIST) + len(TAL_BIO_REACTION_WHITELIST)}")
    print()
    print("  Sampling one bio rule's enzyme catalysis data...")
    try:
        from network_generation import AVAILABLE_RULESETS
        rules_df = pd.read_csv(AVAILABLE_RULESETS["JN1224MIN"], sep="\t")
        sample_name = next(iter(TAL_BIO_REACTION_WHITELIST))
        row = rules_df[rules_df["Name"] == sample_name]
        if not row.empty:
            enzymes = [e.strip() for e in str(row.iloc[0]["Comments"]).split(";")]
            print(f"  Rule: {sample_name}  |  enzymes: {len(enzymes)} UniProt IDs")
            print(f"  First 5 enzyme IDs: {', '.join(enzymes[:5])}")
    except Exception:
        pass
    print()
    print("  Pathway-ranking criteria active in this run (from config.py):")
    for name in ("steps", "thermo", "stability", "diversity",
                 "chem_bio", "atom_economy", "by_product"):
        print(f"    - {name:14s} weight={DEFAULT_WEIGHTS.get(name)}")


def scene_2_exploration():
    banner("[2/4]  OPEN EXPLORATION  —  bio + chem at gen=5  (EXPECT SLOW)")
    print()
    print(f"  Settings: gen={GENERATIONS}, chem+bio enabled, top_n={TOP_N}")
    print(f"  Helpers : {HELPERS}")
    print(f"  max_atoms: C=80, O=15  |  MW cap: 1000 Da")
    print()
    print("  At this scale, cartesian expansion will attempt millions of")
    print("  (operator x molecule combinations) per generation. If your")
    print("  laptop survives this scene, the directed-search scene below")
    print("  will likely finish it off. Press Ctrl+C any time.")

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
            max_atoms={"C": 80, "O": 15, "N": 4, "S": 0},
            max_molecular_weight=1000,
            allow_multiple_reactants="default",
            strategy="cartesian",
            min_carbons=0,
            include_chem=True,
            include_bio=True,                  # <-- the big multiplier
            bio_allow_multiple_reactants=True,
        ),
        scoring_kwargs=dict(carbon_window=(3, 20), require_oxygen=True),
        derivatives_list=TAL_DOWNSTREAM_DERIVATIVES,
    )
    elapsed = time.time() - t0
    print(f"\n  Exploration finished in {elapsed:.1f}s.")
    print(f"  Endpoints surfaced: {len(result['top_endpoints'])}")
    return result


def scene_3_directed_search(exploration_result):
    banner("[3/4]  DIRECTED SEARCH  —  4 precoded targets, full ranking")

    network = exploration_result["network"]
    starter_smi = exploration_result["starter_smiles"]
    helper_smiles = exploration_result["helper_smiles"]

    for label, target_smi in DIRECTED_TARGETS:
        print()
        print(f"--- Target: {label} ---")
        print(f"    SMILES: {target_smi}")

        sub_job = f"{JOB_NAME}_directed_{abs(hash(target_smi)) % 10000:04d}"
        try:
            find_pathways_to_target(
                network=network,
                starter=starter_smi,
                target=target_smi,
                helpers=helper_smiles,
                generations=GENERATIONS,
                max_num_rxns=MAX_NUM_RXNS,
                job_name=sub_job,
            )
            pathways = load_pathways_from_file(sub_job)
        except Exception as exc:
            print(f"    [skip] {exc}")
            continue

        if not pathways:
            print(f"    No pathway found (target not reachable).")
            continue

        # Wire DORAnet criteria into our unified scorer
        target_file = f"explode_target_{abs(hash(target_smi)) % 10000:04d}.smi"
        write_smi(target_file, [target_smi])
        try:
            pathway_ranking(
                starters=STARTER_FILE, helpers=HELPER_FILE,
                target=target_file, job_name=sub_job,
            )
            atom_econ, by_prod = parse_doranet_ranked_file(sub_job, pathways)
            profile = build_unified_profile(atom_econ, by_prod)
            scorer = WeightedPathwayScorer(profile)
            scored = scorer.score(pathways)
            scored.sort(key=lambda sp: sp.final_score, reverse=True)
            rewrite_ranked_file_with_unified_order(sub_job, scored)

            print(f"    Found {len(scored)} pathway(s). Top 3 ranked:")
            for i, sp in enumerate(scored[:3], 1):
                print(f"      Route {i}: final_score={sp.final_score:.3f} "
                      f"({sp.pathway.num_steps} steps)")
            pathway_visualization(
                starters=STARTER_FILE, helpers=HELPER_FILE,
                job_name=sub_job,
            )
            pdf = f"{sub_job}_pathways_visualized.pdf"
            if os.path.exists(pdf):
                print(f"    PDF: {pdf}")
        except Exception as exc:
            print(f"    [scoring/PDF error] {exc}")


def scene_4_what_youre_looking_at():
    banner("[4/4]  WHY THIS DEMO IS DELIBERATELY HARD")
    print("""
  This run reproduces demo.py's pipeline at a scale that exposes the
  scaling cliff of cartesian expansion:

    cartesian:  every operator x every molecule combination, every gen
    bio adds:   41 cofactor molecules whose pairwise combos explode
                inside multi-substrate bio operators
    gen=5:      molecule count compounds ~10x per generation

  Fixes (in priority order):

    1. strategy='priority_queue'
         beam search with the Tanimoto ranker. Caps expansion to
         beam_size recipes per iteration. Order-of-magnitude speedup.

    2. TAL-relevant cofactor subset
         instead of all 41 cofactors, inject only ~10 relevant ones
         (acetyl-CoA, CoA, ATP, NADH, etc.). 4x reduction at every
         bio operator firing.

    3. cloud / cluster compute
         if cartesian must be used (e.g., to guarantee enumeration
         completeness for a publication), run on a workstation or
         cloud VM. Local laptops are not the right hardware.

  This demo is the 'before' picture. The fixes above are the 'after'.
""")


def main():
    banner("TAL PATHWAY EXPLORER  —  EXPLODE DEMO  (expect to Ctrl+C)", char="#")
    t_total = time.time()
    print()
    print("This demo intentionally overwhelms a laptop. It exists to")
    print("show the scaling limit, not to actually finish. Press Ctrl+C")
    print("whenever you've seen enough.")
    print()

    write_smi(STARTER_FILE, [TAL_SMILES])
    write_smi(HELPER_FILE, HELPERS)

    try:
        scene_1_architecture()
        exploration_result = scene_2_exploration()
        scene_3_directed_search(exploration_result)
        scene_4_what_youre_looking_at()
    except KeyboardInterrupt:
        print("\n\n[INTERRUPTED]  Demo killed by user — expected behavior.")
        print("Look at the demo file's WHY section for fixes.")

    n_cleaned = cleanup_artifacts()
    banner(f"END  ({time.time() - t_total:.1f}s)", char="#")
    if n_cleaned:
        print(f"\n[housekeeping] Cleaned up {n_cleaned} intermediate files.\n")


if __name__ == "__main__":
    main()
