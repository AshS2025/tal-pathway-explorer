"""
streamlit_app.py — TAL Pathway Explorer (v0.1, text-only)

Scope for this version:
  - Take a starter SMILES + target SMILES
  - Run a bidirectional chem search (forward from starter + retro from
    target) using our existing generate_network_tal wrapper
  - Trace pathways with DORAnet's pathway_finder
  - Show each pathway's reactions as text

Deliberately NOT in v0.1:
  - Interactive graph visualization
  - Bio expansion
  - Custom whitelists
  - Advanced tuning knobs (MW, atoms, thermo)

Once this works reliably, we add each piece back one at a time.
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

# Force UTF-8 on stdout/stderr so backend print statements with unicode
# arrows don't crash under Windows' default cp1252 encoding.
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except (AttributeError, ValueError):
    pass

import streamlit as st

# Make src/ importable.
_ROOT = Path(__file__).parent
sys.path.insert(0, str(_ROOT / "src"))

from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

from network_generation import generate_network_tal
from pathway_tools import (
    load_pathways_from_file,
    parse_reaction_string,
)
from visualize_pathways import visualize_pathways
from doranet.modules.post_processing.post_processing import (
    pretreat_networks,
    pathway_finder,
)


# --------------------------------------------------------------
# Page config
# --------------------------------------------------------------
st.set_page_config(
    page_title="TAL Pathway Explorer",
    page_icon="⚗️",
    layout="wide",
)

st.title("⚗️ TAL Pathway Explorer")
st.caption("v0.1 — text pathway output only. Graph view coming after this works.")


# --------------------------------------------------------------
# Sidebar — the essential inputs only
# --------------------------------------------------------------
with st.sidebar:
    st.header("Inputs")

    starter_smiles = st.text_input(
        "Starter SMILES",
        value="",
        placeholder="Cc1cc(O)cc(=O)o1",
        help="The molecule you want to start from.",
    )
    target_smiles = st.text_input(
        "Target SMILES",
        value="",
        placeholder="CC=CC=CC(=O)O",
        help="The molecule you want to reach.",
    )
    gen = st.slider(
        "Generations per side",
        min_value=1, max_value=6, value=3,
        help=(
            "Bidirectional search: this many steps forward from starter "
            "AND this many steps retro from target. Slider=3 finds "
            "pathways up to 6 steps long."
        ),
    )

    with st.expander("⚙️  Limits (prevents runaway expansion)", expanded=True):
        max_mw = st.number_input(
            "Max molecular weight (Da)",
            min_value=50, max_value=2000, value=200, step=10,
            help="Products above this weight are rejected.",
        )
        c1, c2, c3 = st.columns(3)
        with c1:
            max_c = st.number_input("Max C atoms", 1, 100, 10)
        with c2:
            max_o = st.number_input("Max O atoms", 0, 30, 5)
        with c3:
            max_n = st.number_input("Max N atoms", 0, 20, 2)
        max_dh = st.number_input(
            "Max |ΔH| per reaction (kJ/mol)",
            min_value=1.0, max_value=200.0, value=15.0, step=1.0,
            help=(
                "Reactions with an absolute enthalpy change larger than "
                "this are rejected as thermodynamically infeasible."
            ),
        )

    st.markdown("---")
    run_button = st.button("Run", type="primary", use_container_width=True)


# --------------------------------------------------------------
# Backend pipeline
# --------------------------------------------------------------
JOB_NAME = "streamlit_job"


def _cleanup():
    """Remove leftover job files so a stale one can't be shown as fresh."""
    for pattern in [
        f"{JOB_NAME}_pathways.txt",
        f"{JOB_NAME}_network_pretreated.json",
        f"{JOB_NAME}_reaxys_batch_query.txt",
        f"{JOB_NAME}_reaxys_batch_result.csv",
    ]:
        try:
            os.remove(pattern)
        except FileNotFoundError:
            pass


def _write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def _validate(starter, target):
    starter = starter.strip()
    target  = target.strip()
    if not starter:
        return None, "Please enter a starter SMILES."
    if not target:
        return None, "Please enter a target SMILES."
    if Chem.MolFromSmiles(starter) is None:
        return None, f"Invalid starter SMILES: `{starter}`"
    if Chem.MolFromSmiles(target) is None:
        return None, f"Invalid target SMILES: `{target}`"
    if Chem.MolToSmiles(Chem.MolFromSmiles(starter)) == \
       Chem.MolToSmiles(Chem.MolFromSmiles(target)):
        return None, (
            "Starter and target are the same molecule — nothing to "
            "search for."
        )
    return {"starter": starter, "target": target}, None


def _run(params, gen_count, limits):
    """Bidirectional chem expansion + pathway trace. Returns pathways list."""
    helpers = ["O", "[H][H]"]

    # File setup
    starter_file = f"{JOB_NAME}_starter.smi"
    helper_file  = f"{JOB_NAME}_helpers.smi"
    target_file  = f"{JOB_NAME}_target.smi"
    _write_smi(starter_file, [params["starter"]])
    _write_smi(helper_file, helpers)
    _write_smi(target_file, [params["target"]])

    # Forward expansion from starter
    net_fwd = generate_network_tal(
        job_name=f"{JOB_NAME}_fwd",
        starters=starter_file,
        helpers=helper_file,
        gen=gen_count,
        direction="forward",
        include_chem=True,
        include_bio=False,
        strategy="cartesian",
        max_atoms=limits["max_atoms"],
        max_molecular_weight=limits["max_mw"],
        max_rxn_thermo_change=limits["max_dh"],
    )

    # Retro expansion from target
    net_retro = generate_network_tal(
        job_name=f"{JOB_NAME}_retro",
        starters=target_file,
        helpers=helper_file,
        gen=gen_count,
        direction="retro",
        include_chem=True,
        include_bio=False,
        strategy="cartesian",
        max_atoms=limits["max_atoms"],
        max_molecular_weight=limits["max_mw"],
        max_rxn_thermo_change=limits["max_dh"],
    )

    # Merge, trace
    reach = gen_count * 2
    pretreat_networks(
        networks=[net_fwd, net_retro],
        starters=[params["starter"]],
        helpers=helpers,
        total_generations=reach,
        job_name=JOB_NAME,
    )
    pathway_finder(
        starters=[params["starter"]],
        helpers=helpers,
        target=[params["target"]],
        search_depth=reach,
        max_num_rxns=reach + 3,
        job_name=JOB_NAME,
    )

    pathway_file = f"{JOB_NAME}_pathways.txt"
    if not os.path.exists(pathway_file):
        return []
    return load_pathways_from_file(JOB_NAME)


def _truncate(s, n=60):
    return s if len(s) <= n else s[:n-3] + "..."


# --------------------------------------------------------------
# Main area — results
# --------------------------------------------------------------
if not run_button:
    st.info(
        "👈 Enter a starter and target on the left, choose the number of "
        "generations, and press **Run**."
    )
    st.markdown(
        """
        **Example**
        - Starter: `Cc1cc(O)cc(=O)o1` (TAL)
        - Target: `CC=CC=CC(=O)O` (sorbic acid)
        - Generations: **3** (finds pathways up to 6 steps)
        """
    )
    st.stop()

# --- Validate ---
params, err = _validate(starter_smiles, target_smiles)
if err:
    st.error(err)
    st.stop()

# --- Run ---
_cleanup()
t0 = time.time()
limits = {
    "max_mw": float(max_mw),
    "max_atoms": {"C": int(max_c), "O": int(max_o), "N": int(max_n)},
    "max_dh": float(max_dh),
}

try:
    with st.spinner(
        f"Running bidirectional search "
        f"({gen} generations forward + {gen} retro)…"
    ):
        pathways = _run(params, gen, limits)
    elapsed = time.time() - t0
except Exception as e:
    st.error(f"Pipeline error: `{type(e).__name__}: {e}`")
    st.stop()

# --- Results ---
if not pathways:
    st.warning(
        f"Network built in **{elapsed:.1f}s** but **no pathways found** "
        f"connecting the starter to the target within {gen*2} steps.\n\n"
        "Try increasing generations, or check that the starter and "
        "target are connectable through our default chem whitelist."
    )
    st.stop()

# Sort shortest first
pathways.sort(key=lambda p: p.num_steps)

st.success(
    f"✓ Found **{len(pathways)}** pathway(s) in **{elapsed:.1f}s**  •  "
    f"shortest: {pathways[0].num_steps} steps  •  "
    f"longest: {pathways[-1].num_steps} steps."
)

# ---- Two tabs: text pathways (guaranteed) + interactive graph (best-effort)
tab_pathways, tab_graph = st.tabs(["📋 Pathways", "🕸️ Graph"])

with tab_pathways:
    st.caption("Shortest first. Each pathway is listed step by step.")
    for p in pathways:
        with st.expander(
            f"Pathway {p.index} — {p.num_steps} step(s)",
            expanded=(p is pathways[0]),
        ):
            for i, rxn_str in enumerate(p.reactions, 1):
                parsed = parse_reaction_string(rxn_str)
                reactants = " + ".join(_truncate(r) for r in parsed["reactants"])
                products  = " + ".join(_truncate(pp) for pp in parsed["products"])
                dh = parsed["dH"]
                dh_str = f"  •  ΔH = {dh:.1f}" if dh is not None else ""
                st.markdown(
                    f"**Step {i}** — `{parsed['op_name']}`{dh_str}  \n"
                    f"{reactants} **→** {products}"
                )

with tab_graph:
    # Wrap-not-reimplement: call visualize_pathways to write the HTML,
    # then read that file and embed it. If anything goes wrong (bad
    # chemistry, RDKit hiccup, whatever) we surface the error and offer
    # the raw HTML as a download so the user isn't stuck.
    st.caption(
        "Interactive DAG. Hover a reaction arrow to see its name. "
        "Click an edge or node to highlight the pathway(s) it belongs to."
    )
    graph_path = None
    try:
        with st.spinner("Rendering interactive graph…"):
            graph_path = visualize_pathways(
                job_name=JOB_NAME,
                starter_smiles=params["starter"],
                target_smiles=params["target"],
                starter_label="starter",
                target_label="target",
                helpers=["O", "[H][H]"],
                pathway_filter="all",
            )
        with open(graph_path, "r", encoding="utf-8") as f:
            html_content = f.read()
        st.components.v1.html(html_content, height=720, scrolling=True)
    except Exception as e:
        st.warning(
            f"Interactive graph render failed: `{type(e).__name__}: {e}`.\n\n"
            "You can still download the pathway file below."
        )

# ---- Shared downloads (visible under both tabs)
st.markdown("---")
st.subheader("Downloads")
col_a, col_b = st.columns(2)
with col_a:
    try:
        with open(f"{JOB_NAME}_pathways.txt", "rb") as f:
            st.download_button(
                "📄 Pathway file (.txt)",
                f.read(),
                file_name=f"{JOB_NAME}_pathways.txt",
                mime="text/plain",
                use_container_width=True,
            )
    except FileNotFoundError:
        st.caption("(No pathway file to download.)")
with col_b:
    graph_html_path = f"{JOB_NAME}_graph.html"
    if os.path.exists(graph_html_path):
        with open(graph_html_path, "rb") as f:
            st.download_button(
                "🌐 Interactive graph (.html)",
                f.read(),
                file_name="pathway_graph.html",
                mime="text/html",
                use_container_width=True,
            )
    else:
        st.caption("(Graph HTML not generated.)")
