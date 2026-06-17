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

from pathway_tools import (
    explore_downstream,
    find_pathways_to_target,
    load_pathways_from_file,
)
from doranet.modules.post_processing.post_processing import (
    pathway_ranking,
    pathway_visualization,
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


JOB_NAME = "tal_demo"
STARTER_FILE = "demo_starter.smi"
HELPER_FILE = "demo_helpers.smi"

TAL_SMILES = "Cc1cc(O)cc(=O)o1"
# H2, water, ethanol. Ethanol unlocks alkylation chemistry,
# enabling multiple routes to ethyl-substituted derivatives —
# without it, the directed-search PDF would only have one page.
HELPERS = ["[H][H]", "O", "CCO"]

GENERATIONS = 3
TOP_N = 10
MAX_NUM_RXNS = 20

# Precoded target: 3-ethyl-6-methyl-2H-pyran-2-one. A TAL
# hydroalkylation derivative — multiple distinct routes exist
# (direct alkylation, hydrogenation-first, ether-intermediate) so
# we get a multi-page ranked PDF.
DIRECTED_TARGET = "CCc1ccc(C)oc1=O"
DIRECTED_TARGET_LABEL = "Ethyl-methyl-pyranone (TAL hydroalkylation derivative)"


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

    print("\nPathway-ranking criteria (all pluggable, all weight-tunable")
    print("from src/config.py:DEFAULT_WEIGHTS):")
    for name, cls in (
        ("steps",        StepsCriterion),
        ("thermo",       ThermoCriterion),
        ("stability",    IntermediateStabilityCriterion),
        ("diversity",    ProcedureDiversityCriterion),
        ("chem_bio",     ChemBioSwitchCriterion),
        ("atom_economy", AtomEconomyCriterion),
        ("by_product",   ByProductCountCriterion),
    ):
        weight = DEFAULT_WEIGHTS.get(name, "—")
        print(f"  - {cls.__name__:32s}  weight={weight}")


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
            # Tighter limits keep gen=3 tractable (matches the notebook
            # config that produced 5 ranked routes in ~10s).
            max_atoms={"C": 25, "O": 4, "N": 2, "S": 0},
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
# Scene 3 — directed search: precoded target, find + rank routes
# =====================================================================
def scene_3_directed_search(exploration_result):
    banner("[3/4]  DIRECTED SEARCH  —  precoded target, ranked routes")
    print()
    print("Question (from the project brief):")
    print('  "I want to make a specific molecule from TAL. What routes')
    print('   exist, and which is the best one?"')
    print()
    print(f"Precoded target:  {DIRECTED_TARGET_LABEL}")
    print(f"                  {DIRECTED_TARGET}")
    print()
    print("Reusing the network built in Scene 2 (no re-expansion).")

    network = exploration_result["network"]
    starter_smi = exploration_result["starter_smiles"]
    helper_smiles = exploration_result["helper_smiles"]

    directed_job = f"{JOB_NAME}_directed"
    try:
        find_pathways_to_target(
            network=network,
            starter=starter_smi,
            target=DIRECTED_TARGET,
            helpers=helper_smiles,
            generations=GENERATIONS,
            max_num_rxns=4,
            job_name=directed_job,
        )
        pathways = load_pathways_from_file(directed_job)
    except Exception as exc:
        print(f"\n  Pathway search failed: {exc}")
        return

    if not pathways:
        print(f"\n  No pathways to the target at gen={GENERATIONS}.")
        print(f"  Either the target isn't reachable in this network, or")
        print(f"  the route exceeds max_num_rxns. Try gen=3 for")
        print(f"  literature targets like sorbic acid.")
        return

    # ----------------------------------------------------------------
    # Run DORAnet's pathway_ranking FIRST so we can pull its per-pathway
    # atom_economy and by_product values into our unified scorer. The
    # ranked file it writes will be rewritten in our order below so the
    # PDF reflects the unified scoring.
    # ----------------------------------------------------------------
    target_file = "demo_directed_target.smi"
    write_smiles_file(target_file, [DIRECTED_TARGET])
    pathway_ranking(
        starters=STARTER_FILE,
        helpers=HELPER_FILE,
        target=target_file,
        job_name=directed_job,
    )
    atom_econ_by_idx, by_prod_by_idx = parse_doranet_ranked_file(
        directed_job, pathways
    )

    profile = build_unified_profile(
        atom_economy_by_index=atom_econ_by_idx,
        by_product_by_index=by_prod_by_idx,
    )
    scorer = WeightedPathwayScorer(profile)

    print(f"\n{len(pathways)} pathway(s) found. Scoring with the unified "
          f"{len(profile)}-criteria scorer...")
    scored = scorer.score(pathways)
    scored.sort(key=lambda sp: sp.final_score, reverse=True)

    # Rewrite the ranked file in OUR order so the PDF visualization
    # reflects the unified scoring, not DORAnet's defaults.
    rewrite_ranked_file_with_unified_order(directed_job, scored)

    print(f"\nAll {len(scored)} ranked routes:")
    for i, sp in enumerate(scored, 1):
        print()
        print(f"  Route {i}  —  final score {sp.final_score:.3f}  "
              f"({sp.pathway.num_steps} steps)")
        print(f"  Criterion breakdown:")
        for name, val in sp.components.items():
            print(f"    {name:32s} = {val:.3f}")
        print(f"  Reaction steps:")
        for j, rxn in enumerate(sp.pathway.reactions, 1):
            try:
                reactants_part, op_name, _, products_part = rxn.split(">")
                reactants = reactants_part.split(".")
                products = products_part.split(".")
                print(f"    {j}. [{op_name}]")
                print(f"       {' + '.join(reactants)}  →  "
                      f"{' + '.join(products)}")
            except Exception:
                print(f"    {j}. {rxn}")

    # Append the directed-search results to the exploration markdown
    # report so both modes live in one artifact.
    _append_directed_search_to_report(scored)

    # Generate the PDF directly from the rewritten ranked file —
    # pages will appear in OUR scorer's order.
    print()
    print("Generating pathway visualization PDF...")
    try:
        pathway_visualization(
            starters=STARTER_FILE,
            helpers=HELPER_FILE,
            job_name=directed_job,
        )
        pdf_path = f"{directed_job}_pathways_visualized.pdf"
        if os.path.exists(pdf_path):
            print(f"  PDF written: {pdf_path}")
        else:
            print("  PDF generation completed but file not found at "
                  f"{pdf_path}")
    except Exception as exc:
        print(f"  (PDF generation skipped: {exc})")


def _append_directed_search_to_report(scored):
    """Append a directed-search section to the exploration report."""
    report_path = f"{JOB_NAME}_exploration_report.md"
    lines = ["", "---", "", "# Directed Search — precoded target", ""]
    lines.append(f"**Target:** {DIRECTED_TARGET_LABEL}")
    lines.append("")
    lines.append(f"`{DIRECTED_TARGET}`")
    lines.append("")
    lines.append(f"## All {len(scored)} ranked routes")
    for i, sp in enumerate(scored, 1):
        lines.append("")
        lines.append(
            f"### Route {i} — final score {sp.final_score:.3f} "
            f"({sp.pathway.num_steps} steps)"
        )
        lines.append("")
        lines.append("**Criterion breakdown:**")
        for name, val in sp.components.items():
            lines.append(f"- `{name}`: {val:.3f}")
        lines.append("")
        lines.append("**Reaction steps:**")
        for j, rxn in enumerate(sp.pathway.reactions, 1):
            try:
                reactants_part, op_name, _, products_part = rxn.split(">")
                arrow = f"{reactants_part} → {products_part}"
                lines.append(f"{j}. `{op_name}` — {arrow}")
            except Exception:
                lines.append(f"{j}. {rxn}")
    lines.append("")
    with open(report_path, "a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


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
        f"{JOB_NAME}_directed_pathways.txt",
        f"{JOB_NAME}_directed_network_pretreated.json",
        f"{JOB_NAME}_directed_reaxys_batch_*.txt",
        f"{JOB_NAME}_directed_reaxys_batch_*.csv",
        f"{JOB_NAME}_directed_ranked_pathways*.txt",
        f"{JOB_NAME}_*_saved_network*",
        f"{JOB_NAME}_network_pretreated.json",
        f"{JOB_NAME}_pathways.txt",
        f"{JOB_NAME}_reaxys_batch_*.txt",
        f"{JOB_NAME}_reaxys_batch_*.csv",
        # Demo input files (regenerated each run)
        "demo_directed_target.smi",
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
    scene_3_directed_search(exploration_result)
    scene_4_next()

    n_cleaned = cleanup_artifacts()

    banner(f"DEMO COMPLETE  ({time.time() - t_total:.1f}s total)", char="#")
    print()
    print("Output files:")
    print(f"  Markdown report (exploration + directed search):")
    print(f"    {JOB_NAME}_exploration_report.md")
    print(f"      (in VS Code, open and press Ctrl+Shift+V for rendered view)")
    pdf_path = f"{JOB_NAME}_directed_pathways_visualized.pdf"
    if os.path.exists(pdf_path):
        print(f"  PDF visualization of the ranked directed-search route:")
        print(f"    {pdf_path}")
        print(f"      (double-click to open in your default PDF viewer)")
    if n_cleaned:
        print(f"\n[housekeeping] Cleaned up {n_cleaned} intermediate files.\n")
    else:
        print()


if __name__ == "__main__":
    main()
