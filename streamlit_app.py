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

# ==== BACKEND-AGNOSTIC PIPELINE (also used by future React/FastAPI) ====
from pipeline import (
    PipelineConfig,
    run_pipeline,
    validate_config,
    cleanup_job_files,
)

# UI still needs these directly for display (not the pipeline itself)
from pathway_tools import parse_reaction_string           # for the rare fallback branch
from pathway_scoring import RankedPathway                 # for type hints in display
from visualize_pathways import visualize_pathways         # for the Graph tab

# ==== Thermo client factories (Streamlit-scoped — cached per session) ====
try:
    from equilibrator_client import EquilibratorClient
    _EQUILIBRATOR_AVAILABLE = True
except Exception:
    EquilibratorClient = None
    _EQUILIBRATOR_AVAILABLE = False
try:
    from rmg_thermo import RMGThermoClient
    _RMG_AVAILABLE = True
except Exception:
    RMGThermoClient = None
    _RMG_AVAILABLE = False


# --------------------------------------------------------------
# Windows-multiprocessing safety
# --------------------------------------------------------------
# DORAnet's pathway_ranking uses multiprocessing.Pool internally. On
# Windows, the 'spawn' method re-imports this file in each worker
# subprocess — which then executes Streamlit UI code with no widget
# state, and crashes with NameError forever. Instead of refactoring
# the whole app into a main() function, we replace multiprocessing.Pool
# with an in-process shim: same API, runs work synchronously in this
# process. DORAnet's pool.map calls execute the same logic without ever
# spawning a subprocess. Ranking is slightly slower single-threaded but
# our pathway counts are small (~50) so this is inconsequential.
import multiprocessing as _mp


class _InProcessAsyncResult:
    """Mimics multiprocessing.pool.AsyncResult but runs synchronously."""
    def __init__(self, value): self._value = value
    def get(self, timeout=None): return self._value
    def wait(self, timeout=None): pass
    def ready(self): return True
    def successful(self): return True


class _InProcessPool:
    """Drop-in replacement for multiprocessing.Pool that runs work in
    the calling process. Windows spawn is the actual failure mode
    for us — avoiding subprocess creation entirely is safer than
    trying to make the Streamlit script safe to re-import."""
    def __init__(self, *args, **kwargs): pass
    def __enter__(self): return self
    def __exit__(self, *args): pass
    def map(self, func, iterable):
        return [func(x) for x in iterable]
    def imap(self, func, iterable):
        for x in iterable:
            yield func(x)
    def imap_unordered(self, func, iterable):
        return self.imap(func, iterable)
    def apply(self, func, args=(), kwds=None):
        return func(*args, **(kwds or {}))
    def apply_async(self, func, args=(), kwds=None, callback=None,
                    error_callback=None):
        try:
            result = func(*args, **(kwds or {}))
            if callback is not None:
                callback(result)
            return _InProcessAsyncResult(result)
        except Exception as e:
            if error_callback is not None:
                error_callback(e)
            raise
    def starmap(self, func, iterable):
        return [func(*args) for args in iterable]
    def close(self): pass
    def join(self): pass
    def terminate(self): pass


# Patch both the multiprocessing namespace (for any future imports)
# AND DORAnet's already-cached Pool reference (imported at DORAnet's
# module-load time before we could patch multiprocessing globally).
_mp.Pool = _InProcessPool
from doranet.modules.post_processing import post_processing as _dpp
_dpp.Pool = _InProcessPool


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
    # Domain and direction are INDEPENDENT axes:
    #   Domain    = which operator set fires (chem / bio / both)
    #   Direction = which way we expand (forward / retro / bidirectional)
    # No starter/target swap is ever applied — the molecules are used
    # exactly as entered.
    domain_choice = st.radio(
        "Domain (operators)",
        options=["chem", "bio", "both"],
        index=0,
        format_func=lambda d: {
            "chem": "Chem — synthetic organic operators",
            "bio":  "Bio — enzymatic (JN1224MIN) operators",
            "both": "Both — chem + bio in one network",
        }[d],
        help=(
            "**Chem**: synthetic organic chemistry operators.\n\n"
            "**Bio**: enzymatic operators from the JN1224MIN rule set "
            "(whitelist below controls which rules fire).\n\n"
            "**Both**: merge chem + bio operators into a single network "
            "so pathways can mix enzymatic and synthetic steps."
        ),
    )
    direction_choice = st.radio(
        "Search direction",
        options=["bidirectional", "forward", "retro"],
        index=0,
        format_func=lambda d: {
            "bidirectional": "Bidirectional — expand from both ends (meet in middle)",
            "forward":       "Forward — expand from starter → target",
            "retro":         "Retro — expand back from target → starter",
        }[d],
        help=(
            "**Bidirectional**: expand forward from the starter AND "
            "backward from the target, then look for pathways where the "
            "two frontiers meet. Best coverage; each side runs to the "
            "generation depth below (so total path length can be 2×).\n\n"
            "**Forward**: only expand from the starter toward the target.\n\n"
            "**Retro**: only expand backward from the target toward the "
            "starter.\n\n"
            "Direction does NOT swap the molecules — starter is always "
            "what you start from, target is always what you want to reach."
        ),
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

    # Derived domain flags — driven by the Domain selector, not direction.
    include_chem = domain_choice in ("chem", "both")
    include_bio  = domain_choice in ("bio", "both")

    # Bio whitelist textarea, visible whenever bio operators are enabled.
    bio_whitelist_text = ""
    if include_bio:
        with st.expander("🧬  Bio rule whitelist", expanded=True):
            bio_whitelist_text = st.text_area(
                "One JN1224MIN rule per line (e.g. rule1118)",
                value="rule1118\nrule0087\nrule0891",
                height=120,
                help=(
                    "Which JN1224MIN operators to enable. The default 3 "
                    "rules are the polyketide chain: Claisen 1, Claisen 2, "
                    "cyclization → sufficient for acetyl-CoA → TAL. Add "
                    "more rules to explore broader biosynthesis. Blank = "
                    "fall back to the built-in TAL bio whitelist."
                ),
            )

    with st.expander("🎯  Search strategy", expanded=False):
        strategy = st.radio(
            "Strategy",
            options=["priority_queue", "cartesian"],
            index=0,
            format_func=lambda s: {
                "priority_queue": "Priority queue  (fast, target-guided)",
                "cartesian":       "Cartesian  (exhaustive, slower)",
            }[s],
            help=(
                "**Priority queue** uses a target-similarity ranker (Tanimoto "
                "for forward, feedstock-proximity for retro) to prune branches. "
                "Finds the most direct routes fast, but may miss some pathways "
                "the beam prunes.\n\n"
                "**Cartesian** expands every candidate reaction with no ranking. "
                "Complete — finds every pathway the limits allow — but slow "
                "and produces more permutation-noise results."
            ),
        )
        beam_size = st.number_input(
            "Beam size (priority queue only)",
            min_value=50, max_value=10000, value=1000, step=100,
            disabled=(strategy != "priority_queue"),
            help=(
                "How many top-ranked candidate reactions the priority queue "
                "expands per iteration. Bigger = more pathway coverage, "
                "slower runtime. 1000 is a good starting point; crank up if "
                "you want IP-diversification-style breadth."
            ),
        )

    with st.expander("🧪  Thermodynamics", expanded=False):
        thermo_enabled = st.checkbox(
            "Enable RMG thermodynamics",
            value=False,
            disabled=(not _RMG_AVAILABLE),
            help=(
                "When on, the ranker uses real per-reaction enthalpies "
                "(ΔH) computed by RMG. Startup cost: ~60 seconds on the "
                "first run of the session. Best suited for CHEM pathways."
            ),
        )
        if not _RMG_AVAILABLE:
            st.caption("⚠️  RMG env not detected. Install rmg_env to enable.")

        equilibrator_enabled = st.checkbox(
            "Enable equilibrator (biochemistry ΔG'° at pH 7)",
            value=False,
            help=(
                "When on, computes standard reaction free energies at "
                "physiological conditions using equilibrator_api. Best "
                "suited for BIO pathways with KEGG-listed cofactors. "
                "Startup cost: ~20 seconds on first run. Reactions whose "
                "compounds aren't in equilibrator's database are marked "
                "'—' rather than scored."
            ),
        )
        equilibrator_max_dg = st.number_input(
            "Prune pathways with any step |ΔG'°| > this value (kJ/mol)",
            min_value=0.0, max_value=500.0, value=100.0, step=10.0,
            disabled=(not equilibrator_enabled),
            help=(
                "Post-hoc filter: pathways whose worst step exceeds "
                "this |ΔG'°| are dropped from the results. Reactions "
                "we couldn't score (compound not in DB) don't count "
                "against the threshold. Set high (200+) to disable "
                "pruning and use equilibrator purely for ranking."
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


@st.cache_resource
def _get_equilibrator_client():
    """Spawn (or reuse a cached) equilibrator client. Init takes ~20s
    on first call — subsequent Runs in the same session reuse it."""
    if not _EQUILIBRATOR_AVAILABLE:
        return None
    try:
        return EquilibratorClient()
    except Exception as e:
        st.warning(
            f"Could not initialize equilibrator: "
            f"`{type(e).__name__}: {e}`. Ranking will run without ΔG'°."
        )
        return None


@st.cache_resource
def _get_rmg_client():
    """
    Spawn (or reuse a cached) RMG Python-2.7 subprocess for thermo
    calculations. The @st.cache_resource decorator ensures we pay the
    ~60-second RMG database load exactly ONCE per Streamlit session,
    not on every Run click. Subsequent queries against the running
    server are ~milliseconds thanks to its internal SMILES→Hf cache.

    Returns None if the RMG env isn't installed on this machine.
    """
    if not _RMG_AVAILABLE:
        return None
    try:
        return RMGThermoClient()
    except Exception as e:
        st.warning(
            f"Could not start RMG thermo server: `{type(e).__name__}: {e}`."
            " Ranking will run without thermodynamics."
        )
        return None


# Note: _cleanup, _validate, _run were deleted here. Their logic now
# lives in src/pipeline.py as `cleanup_job_files`, `validate_config`,
# and `run_pipeline` — used by both this Streamlit UI and any future
# React/FastAPI frontend.


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

# ==============================================================
# Build the PipelineConfig from widget state, hand it to the
# backend module. Validation + bio whitelist parsing
# + expansion + ranking + equilibrator decoration all live in
# `pipeline.run_pipeline` — this file only knows about widgets
# and display.
# ==============================================================
bio_whitelist_lines = None
if include_bio and bio_whitelist_text.strip():
    bio_whitelist_lines = [
        ln.strip() for ln in bio_whitelist_text.split("\n") if ln.strip()
    ]

config = PipelineConfig(
    starter_smiles=starter_smiles,
    target_smiles=target_smiles,
    domain=domain_choice,
    direction=direction_choice,
    generations=int(gen),
    strategy=strategy,
    beam_size=int(beam_size),
    max_molecular_weight=float(max_mw),
    max_atoms_c=int(max_c),
    max_atoms_o=int(max_o),
    max_atoms_n=int(max_n),
    max_rxn_dh=float(max_dh),
    bio_whitelist=bio_whitelist_lines,
    enable_rmg=bool(thermo_enabled),
    enable_equilibrator=bool(equilibrator_enabled),
    equilibrator_prune_max_abs_dg=float(equilibrator_max_dg),
    job_name=JOB_NAME,
)

err = validate_config(config)
if err:
    st.error(err)
    st.stop()

cleanup_job_files(JOB_NAME)

# Fetch (Streamlit-cached) clients — passed IN to run_pipeline so it
# stays UI-agnostic. A FastAPI backend would keep these as module
# singletons instead.
thermo_calc = None
if config.enable_rmg:
    with st.spinner(
        "Spawning RMG thermo server (first run of session takes ~60s)…"
    ):
        thermo_calc = _get_rmg_client()

equilibrator_client = None
if config.enable_equilibrator:
    with st.spinner(
        "Initializing equilibrator (first run of session takes ~20s)…"
    ):
        equilibrator_client = _get_equilibrator_client()

_domain_parts = []
if config.include_chem:
    _domain_parts.append("chem")
if config.include_bio:
    _domain_parts.append("bio")
domain_label = "+".join(_domain_parts) or "no-domain"
strategy_label = (
    f"priority queue (beam={config.beam_size})"
    if config.strategy == "priority_queue"
    else "cartesian (exhaustive)"
)
thermo_label = "with RMG thermo" if thermo_calc is not None else "no RMG"

_dir_label = {
    "bidirectional": f"bidirectional ({config.generations} fwd + {config.generations} retro)",
    "forward":       f"forward ({config.generations} generations)",
    "retro":         f"retro ({config.generations} generations)",
}[config.direction]
with st.spinner(
    f"Running {domain_label} search, {_dir_label}, "
    f"{strategy_label}, {thermo_label}…"
):
    result = run_pipeline(
        config,
        thermo_calc=thermo_calc,
        equilibrator_client=equilibrator_client,
    )

if not result.ok:
    st.error(f"Pipeline error: {result.error}")
    st.stop()

if result.n_pathways == 0:
    reason = result.diagnostics.get("reason", "unknown")
    st.warning(
        f"Network built in **{result.elapsed_seconds:.1f}s** but "
        f"**no pathways found** ({reason}).\n\n"
        "Try increasing generations, loosening the atom / MW / thermo "
        "limits, or expanding the bio whitelist if you're using bio operators."
    )
    st.stop()

ranked_pathways = result.ranked_pathways
n_eq_pruned = result.diagnostics.get("equilibrator_pruned", 0)
elapsed = result.elapsed_seconds

if n_eq_pruned:
    st.info(
        f"Equilibrator pruned **{n_eq_pruned}** pathway(s) with a step "
        f"|ΔG'°| > {config.equilibrator_prune_max_abs_dg} kJ/mol."
    )

st.success(
    f"✓ Found **{result.n_pathways}** pathway(s) in **{elapsed:.1f}s**  •  "
    f"ranked by composite score (steps, thermo, atom-economy, byproducts)."
)

# ---- Two tabs: ranked pathways table + interactive graph
tab_pathways, tab_graph = st.tabs(["📋 Pathways", "🕸️ Graph"])


def _fmt_dh(dh):
    return f"{dh:+.1f}" if dh is not None else "—"


def _fmt_cov(c):
    return f"{int(round(c*100))}%"


with tab_pathways:
    st.caption(
        "Ranked by DORAnet's 7-criterion composite score (steps, thermo, "
        "atom economy, byproducts). Higher = better trade-off."
    )
    summary_rows = [
        {
            "Rank": p.rank,
            "Score": round(p.final_score, 2),
            "Steps": p.num_steps,
            "Max ΔH (kJ/mol)": _fmt_dh(p.max_dh),
            "Max ΔG'° (kJ/mol)": _fmt_dh(p.equilibrator_max_dg),
            "ΔG'° coverage": _fmt_cov(p.equilibrator_coverage),
            "Atom econ.": round(p.atomic_economy, 2),
            "Byproducts": p.pathway_byproduct_count,
        }
        for p in ranked_pathways
    ]
    st.dataframe(
        summary_rows, use_container_width=True, hide_index=True,
    )

    st.markdown("### Details per pathway")
    for p in ranked_pathways:
        with st.expander(
            f"Rank {p.rank} — score {p.final_score:.2f} — "
            f"{p.num_steps} step(s)",
            expanded=(p is ranked_pathways[0]),
        ):
            max_dh_str  = _fmt_dh(p.max_dh)
            avg_dh_str  = _fmt_dh(p.avg_dh)
            max_dg_str  = _fmt_dh(p.equilibrator_max_dg)
            avg_dg_str  = _fmt_dh(p.equilibrator_avg_dg)
            st.markdown(
                f"**Atom economy:** {p.atomic_economy:.2f}  •  "
                f"**Byproducts:** {p.pathway_byproduct_count}  \n"
                f"**ΔH** (RMG): max {max_dh_str}, avg {avg_dh_str} kJ/mol  •  "
                f"**ΔG'°** (equilibrator): max {max_dg_str}, "
                f"avg {avg_dg_str} kJ/mol "
                f"(coverage {_fmt_cov(p.equilibrator_coverage)})"
            )
            for i, (smi, name, dh) in enumerate(zip(
                p.reaction_smiles, p.reaction_names,
                p.reaction_enthalpies,
            ), 1):
                lhs, rhs = smi.split(">>", 1)
                reactants = " + ".join(_truncate(r) for r in lhs.split("."))
                products  = " + ".join(_truncate(pp) for pp in rhs.split("."))
                dh_str = f"  •  ΔH = {dh:.1f}" if dh is not None else ""
                st.markdown(
                    f"**Step {i}** — `{name}`{dh_str}  \n"
                    f"{reactants} **→** {products}"
                )

with tab_graph:
    # No swap anymore: starter is always on the left, target on the
    # right, matching exactly what the user entered.
    graph_starter, graph_target = config.starter_smiles, config.target_smiles

    st.caption(
        "Interactive DAG. Hover a reaction arrow to see its name. "
        "Click an edge or node to highlight the pathway(s) it belongs to."
    )
    try:
        with st.spinner("Rendering interactive graph…"):
            graph_path = visualize_pathways(
                job_name=JOB_NAME,
                starter_smiles=graph_starter,
                target_smiles=graph_target,
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
