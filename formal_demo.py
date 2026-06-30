"""
formal_demo.py
==============

End-to-end demo that exercises the three main workflows of the tool:

  PHASE 1 — Exploration (gen=2 cartesian)
    Forward expand from TAL with no target. Score every reachable
    molecule by interestingness. Write a markdown report listing the
    top endpoints.

  PHASE 2 — Bidirectional search for SORBIC ACID
    Forward expansion from TAL + retro expansion from sorbic acid.
    Combine via DORAnet's pretreat_networks (auto-flips retro
    reactions). Extract pathways. Write a pathways report + render
    an interactive HTML DAG.

  PHASE 3 — Bidirectional search for ACETYLACETONE
    Same workflow as Phase 2, different target.

OUTPUT FILES (every name contains "formal" + "demo")
----------------------------------------------------
  Phase 1
    formal_demo_explore_exploration_report.md
  Phase 2
    formal_demo_sorbic_pathways_report.md
    formal_demo_sorbic_graph.html
  Phase 3
    formal_demo_acetylacetone_pathways_report.md
    formal_demo_acetylacetone_graph.html

RUN
---
From the project root:

    python formal_demo.py
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

from rdkit import Chem

from network_generation import generate_network_tal
from pathway_tools import (
    explore_downstream,
    find_pathways_to_target,
    load_pathways_from_file,
    parse_reaction_string,
)
from recipe_rankers import FeedstockProximityRanker
from tal_downstream_derivatives import TAL_DOWNSTREAM_DERIVATIVES
from visualize_pathways import visualize_pathways

from doranet.modules.post_processing.post_processing import (
    pretreat_networks,
    pathway_finder,
    pathway_ranking,
)

from pathway_scoring import (
    WeightedPathwayScorer,
    parse_doranet_ranked_file,
)
from config import DEFAULT_WEIGHTS, build_unified_profile


# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------
TAL_SMILES = "Cc1cc(O)cc(=O)o1"
SORBIC_SMILES = "CC=CC=CC(=O)O"                # no-stereo form
ACETYLACETONE_SMILES = "CC(=O)CC(C)=O"         # canonical form

HELPERS = ["O", "[H][H]"]


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def banner(text, char="="):
    print()
    print(char * 70)
    print(f" {text}")
    print(char * 70)


def write_smi(path, smiles_list):
    with open(path, "w") as f:
        for s in smiles_list:
            f.write(s + "\n")


def write_pathways_report(pathways, starter_label, target_label, output_md):
    """Write a markdown report of the pathways."""
    lines = [
        f"# {starter_label} → {target_label} Pathways",
        "",
        f"Found **{len(pathways)} pathway(s)**.",
        "",
    ]
    for i, p in enumerate(pathways, 1):
        lines.append(f"## Pathway {i} — {p.num_steps} steps")
        lines.append("")
        for j, rxn in enumerate(p.reactions, 1):
            try:
                parsed = parse_reaction_string(rxn)
                arrow = (
                    " + ".join(parsed["reactants"])
                    + " → "
                    + " + ".join(parsed["products"])
                )
                dH = parsed["dH"]
                dH_s = f"{dH:.2f}" if dH is not None else "No_Thermo"
                lines.append(
                    f"{j}. **{parsed['op_name']}** "
                    f"(ΔH={dH_s})"
                )
                lines.append(f"   `{arrow}`")
            except Exception:
                lines.append(f"{j}. {rxn}")
        lines.append("")
    with open(output_md, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def run_bidirectional_search(
    job_name, starter_smiles, target_smiles, gen, max_atoms, max_mw,
):
    """
    Run a bidirectional pathway search:
      1. Forward expansion from starter (cartesian).
      2. Retro expansion from target (priority queue with
         FeedstockProximityRanker pointed back at starter).
      3. Combine networks with DORAnet's pretreat_networks (auto
         flips retro reactions).
      4. Run pathway_finder on the combined network.

    Returns the list of found pathways.
    """
    starter_file = f"{job_name}_starter.smi"
    helpers_file = f"{job_name}_helpers.smi"
    retro_starter_file = f"{job_name}_retro_starter.smi"

    write_smi(starter_file, [starter_smiles])
    write_smi(helpers_file, HELPERS)
    write_smi(retro_starter_file, [target_smiles])

    # 1. Forward expansion
    print(f"  [forward] expansion from starter (gen={gen}, cartesian)...")
    fwd_net = generate_network_tal(
        job_name=f"{job_name}_fwd",
        starters=starter_file,
        helpers=helpers_file,
        gen=gen,
        direction="forward",
        molecule_thermo_calculator=None,
        max_atoms=max_atoms,
        max_molecular_weight=max_mw,
        strategy="cartesian",
        include_chem=True,
        include_bio=False,
    )
    print(f"           {len(fwd_net.mols)} mols, {len(fwd_net.rxns)} rxns")

    # 2. Retro expansion (priority queue pulled toward the starter)
    print(f"  [retro]   expansion from target (gen={gen}, priority_queue)...")
    retro_net = generate_network_tal(
        job_name=f"{job_name}_retro",
        starters=retro_starter_file,
        helpers=helpers_file,
        gen=gen,
        direction="retro",
        molecule_thermo_calculator=None,
        max_atoms=max_atoms,
        max_molecular_weight=max_mw,
        strategy="priority_queue",
        targets=starter_smiles,
        recipe_ranker=FeedstockProximityRanker([starter_smiles]),
        beam_size=200,
        include_chem=True,
        include_bio=False,
    )
    print(f"           {len(retro_net.mols)} mols, {len(retro_net.rxns)} rxns")

    # 3. Combine into a unified network
    print(f"  [combine] pretreat_networks(forward + retro)...")
    pretreat_networks(
        networks=[fwd_net, retro_net],
        starters=[starter_smiles],
        helpers=HELPERS,
        total_generations=gen * 2,
        job_name=job_name,
    )

    # 4. Trace continuous pathways through the unified network
    print(f"  [trace]   pathway_finder on combined network...")
    pathway_finder(
        starters=[starter_smiles],
        helpers=HELPERS,
        target=[target_smiles],
        search_depth=gen * 2,
        max_num_rxns=15,
        job_name=job_name,
    )

    try:
        return load_pathways_from_file(job_name)
    except FileNotFoundError:
        return []


def rank_pathways_and_report(
    pathways, job_name, starter_smiles, target_smiles,
    starter_label, target_label, output_md,
):
    """
    Rank pathways with WeightedPathwayScorer (default weights from
    config.py) and write a "top N" markdown report.

    Calls DORAnet's pathway_ranking first to get atom-economy and
    by-product values; combines those with our custom criteria via
    build_unified_profile.
    """
    if not pathways:
        return None

    # 1. Run DORAnet pathway_ranking to get atom_economy + by_product
    starter_file = f"{job_name}_starter.smi"
    helpers_file = f"{job_name}_helpers.smi"
    target_file  = f"{job_name}_target.smi"
    write_smi(target_file, [target_smiles])

    try:
        pathway_ranking(
            starters=starter_file,
            helpers=helpers_file,
            target=target_file,
            job_name=job_name,
        )
        atom_econ, by_prod = parse_doranet_ranked_file(job_name, pathways)
    except Exception as exc:
        print(f"  [rank]    pathway_ranking failed ({exc}); "
              f"falling back to custom-only scoring")
        atom_econ, by_prod = {}, {}

    # 2. Build unified scorer
    profile = build_unified_profile(
        atom_economy_by_index=atom_econ if atom_econ else None,
        by_product_by_index=by_prod if by_prod else None,
    )
    scorer = WeightedPathwayScorer(profile)

    # 3. Score and sort
    scored = scorer.score(pathways)
    scored.sort(key=lambda sp: sp.final_score, reverse=True)

    # 4. Write markdown report
    n_show = min(5, len(scored))
    lines = [
        f"# Ranked Pathways — {starter_label} → {target_label}",
        "",
        f"**{len(scored)} pathways** scored with WeightedPathwayScorer "
        f"(default weights from `src/config.py`).",
        "",
        "## Default weights",
        "",
    ]
    for name, weight in DEFAULT_WEIGHTS.items():
        lines.append(f"- `{name}`: {weight}")
    lines.extend([
        "",
        f"## Top {n_show} ranked pathways",
        "",
    ])
    for rank, sp in enumerate(scored[:n_show], 1):
        lines.append(
            f"### Rank {rank} — final score {sp.final_score:.3f} "
            f"({sp.pathway.num_steps} steps)"
        )
        lines.append("")
        lines.append("**Criterion breakdown:**")
        lines.append("")
        for cname, cval in sp.components.items():
            lines.append(f"- `{cname}`: {cval:.3f}")
        lines.append("")
        lines.append("**Reaction steps:**")
        lines.append("")
        for j, rxn in enumerate(sp.pathway.reactions, 1):
            try:
                parsed = parse_reaction_string(rxn)
                arrow = (
                    " + ".join(parsed["reactants"])
                    + " → "
                    + " + ".join(parsed["products"])
                )
                lines.append(f"{j}. **{parsed['op_name']}** — `{arrow}`")
            except Exception:
                lines.append(f"{j}. {rxn}")
        lines.append("")

    with open(output_md, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    return scored


def run_forward_search(
    job_name, starter_smiles, target_smiles, gen, max_atoms, max_mw,
    max_num_rxns=10,
):
    """
    Forward-only directed search:
      1. Forward cartesian expansion from starter
      2. find_pathways_to_target traces routes starter -> target

    Returns the list of found pathways.
    """
    starter_file = f"{job_name}_starter.smi"
    helpers_file = f"{job_name}_helpers.smi"
    write_smi(starter_file, [starter_smiles])
    write_smi(helpers_file, HELPERS)

    print(f"  [forward] cartesian expansion (gen={gen})...")
    network = generate_network_tal(
        job_name=job_name,
        starters=starter_file,
        helpers=helpers_file,
        gen=gen,
        direction="forward",
        molecule_thermo_calculator=None,
        max_atoms=max_atoms,
        max_molecular_weight=max_mw,
        strategy="cartesian",
        include_chem=True,
        include_bio=False,
    )
    print(f"           {len(network.mols)} mols, {len(network.rxns)} rxns")

    print(f"  [trace]   pathway_finder (max_num_rxns={max_num_rxns})...")
    find_pathways_to_target(
        network=network,
        starter=starter_smiles,
        target=target_smiles,
        helpers=HELPERS,
        generations=gen,
        max_num_rxns=max_num_rxns,
        job_name=job_name,
    )
    try:
        return load_pathways_from_file(job_name)
    except FileNotFoundError:
        return []


def cleanup_intermediates():
    """Remove intermediate DORAnet files but keep the user-facing
    reports and HTML graphs."""
    patterns = [
        "formal_demo_*pathways.txt",
        "formal_demo_*_pretreated.json",
        "formal_demo_*_reaxys_batch_*",
        "formal_demo_*_saved_network*",
        "formal_demo_*_ranked_pathways*",
        "formal_demo_*_starter.smi",
        "formal_demo_*_helpers.smi",
        "formal_demo_*_retro_starter.smi",
        "formal_demo_*_target.smi",
        "formal_demo_*_ranked_pathways.txt",
        # exploration mode per-endpoint intermediates
        "formal_demo_explore_ep*_pathways.txt",
        "formal_demo_explore_ep*_pretreated.json",
        "formal_demo_explore_ep*_reaxys_batch_*",
        "formal_demo_explore*_saved_network*",
    ]
    n = 0
    for pat in patterns:
        for path in glob.glob(pat):
            try:
                os.remove(path)
                n += 1
            except OSError:
                pass
    return n


# ----------------------------------------------------------------------
# Phases
# ----------------------------------------------------------------------
def phase_1_exploration():
    banner("PHASE 1  —  Exploration (gen=2 cartesian from TAL)", char="#")
    print("Forward expansion from TAL with no target. Top reachable")
    print("molecules ranked by interestingness, with literature cross-check.")
    print()

    starter_file = "formal_demo_explore_starter.smi"
    helpers_file = "formal_demo_explore_helpers.smi"
    write_smi(starter_file, [TAL_SMILES])
    write_smi(helpers_file, HELPERS)

    t0 = time.time()
    result = explore_downstream(
        starter=starter_file,
        helpers=helpers_file,
        job_name="formal_demo_explore",
        top_n=10,
        max_num_rxns=3,
        network_kwargs=dict(
            gen=2,
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

    print(f"\nPhase 1 done in {elapsed:.1f}s.")
    print(f"  Top endpoints surfaced:  {len(result['top_endpoints'])}")
    if result['derivative_matches'] is not None:
        n = len(result['derivative_matches'])
        tot = n + len(result['derivative_missing'])
        print(f"  Literature derivatives:  {n}/{tot} reached")
    print(f"  Report: formal_demo_explore_exploration_report.md")


def phase_search(
    label, target_smi, target_label, mode, gen, max_atoms, max_mw,
    max_num_rxns=10,
):
    """
    Run a targeted search (forward-only OR bidirectional) and write
    the report + visualization.

    mode = "forward" or "bidirectional"
    """
    if mode == "bidirectional":
        banner(f"BIDIRECTIONAL SEARCH  —  TAL  ↔  {label}", char="#")
        print(f"Forward TAL expansion + retro {label} expansion,")
        print(f"combined via DORAnet's pretreat_networks.")
    else:
        banner(f"FORWARD SEARCH  —  TAL  →  {label}", char="#")
        print(f"Forward cartesian expansion + pathway_finder to {label}.")
    print()

    job_name = f"formal_demo_{label.lower().replace(' ', '_')}"

    t0 = time.time()
    if mode == "bidirectional":
        pathways = run_bidirectional_search(
            job_name=job_name,
            starter_smiles=TAL_SMILES,
            target_smiles=target_smi,
            gen=gen,
            max_atoms=max_atoms,
            max_mw=max_mw,
        )
    elif mode == "forward":
        pathways = run_forward_search(
            job_name=job_name,
            starter_smiles=TAL_SMILES,
            target_smiles=target_smi,
            gen=gen,
            max_atoms=max_atoms,
            max_mw=max_mw,
            max_num_rxns=max_num_rxns,
        )
    else:
        raise ValueError(f"Unknown mode {mode!r}")
    elapsed = time.time() - t0

    print(f"\nSearch done in {elapsed:.1f}s.")
    print(f"  Pathways found: {len(pathways)}")

    if not pathways:
        print(f"  (No pathways — skipping report + visualization)")
        return

    # Write all-pathways markdown report
    report_path = f"{job_name}_pathways_report.md"
    write_pathways_report(
        pathways, starter_label="TAL", target_label=target_label,
        output_md=report_path,
    )
    print(f"  Report:        {report_path}")

    # Rank pathways and write top-5 markdown report
    ranked_path = f"{job_name}_ranked_report.md"
    scored = rank_pathways_and_report(
        pathways=pathways,
        job_name=job_name,
        starter_smiles=TAL_SMILES,
        target_smiles=target_smi,
        starter_label="TAL",
        target_label=target_label,
        output_md=ranked_path,
    )
    if scored:
        print(f"  Ranked top 5:  {ranked_path}")
        top = scored[0]
        print(f"    Best route: {top.pathway.num_steps} steps, "
              f"final score = {top.final_score:.3f}")

    # Render interactive HTML graph
    html_path = visualize_pathways(
        job_name=job_name,
        starter_smiles=TAL_SMILES,
        target_smiles=target_smi,
        starter_label="TAL",
        target_label=target_label,
        pathway_filter="all",
        output_html=f"{job_name}_graph.html",
    )
    print(f"  Visualization: {html_path}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    banner("FORMAL DEMO — TAL pathway explorer", char="#")
    t_total = time.time()

    # Phase 1
    phase_1_exploration()

    # Phase 2 — sorbic acid via BIDIRECTIONAL search
    # (a long route — bidirectional cuts the effective depth in half)
    phase_search(
        label="sorbic acid",
        target_smi=SORBIC_SMILES,
        target_label="sorbic acid",
        mode="bidirectional",
        gen=4,
        max_atoms={"C": 6, "O": 3, "N": 0, "S": 0},
        max_mw=150,
    )

    # Phase 3 — acetylacetone via FORWARD-ONLY search
    # (a short route — bidirectional finds dozens of mechanistic
    # permutations; forward-only stays at the canonical 3 routes.)
    phase_search(
        label="acetylacetone",
        target_smi=ACETYLACETONE_SMILES,
        target_label="acetylacetone",
        mode="forward",
        gen=4,
        max_atoms={"C": 8, "O": 5, "N": 0, "S": 0},
        max_mw=200,
        max_num_rxns=10,
    )

    n_cleaned = cleanup_intermediates()
    total = time.time() - t_total
    banner(f"DEMO COMPLETE  ({total:.1f}s total)", char="#")
    print()
    print("User-facing outputs (the rest is cleaned up):")
    print("  formal_demo_explore_exploration_report.md")
    print("  formal_demo_sorbic_acid_pathways_report.md")
    print("  formal_demo_sorbic_acid_ranked_report.md")
    print("  formal_demo_sorbic_acid_graph.html")
    print("  formal_demo_acetylacetone_pathways_report.md")
    print("  formal_demo_acetylacetone_ranked_report.md")
    print("  formal_demo_acetylacetone_graph.html")
    if n_cleaned:
        print(f"\n[housekeeping] cleaned {n_cleaned} intermediate files.")
    print()


if __name__ == "__main__":
    main()
