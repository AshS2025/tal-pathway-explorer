"""
src/visualize_pathways.py
=========================

A callable function that renders pathways from any
{job_name}_pathways.txt as an interactive layered DAG (HTML).

DESIGN
------
- Build the graph in NetworkX (pure data structure)
- Compute node positions in NetworkX (multipartite_layout by level)
- Hand to PyVis purely as renderer (positions locked, no physics)

Why split: PyVis's built-in hierarchical layout has quirks that
produced label overlap and misaligned columns. Computing positions
in NetworkX first sidesteps them entirely.

USAGE
-----
    from visualize_pathways import visualize_pathways

    visualize_pathways(
        job_name="bidir_combined",
        starter_smiles="Cc1cc(O)cc(=O)o1",
        target_smiles="CC=CC=CC(=O)O",
        starter_label="TAL",
        target_label="sorbic acid",
        pathway_filter="all",   # "shortest", "all", or [1, 4]
    )

Intermediate molecule labels default to (truncated) SMILES.
Full SMILES always appears on hover.
"""

from __future__ import annotations

import base64
import hashlib
import json
import os
from io import BytesIO
from typing import Iterable, Optional, Union

import networkx as nx
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Draw, rdMolDescriptors
from pyvis.network import Network

from pathway_tools import load_pathways_from_file, parse_reaction_string

RDLogger.DisableLog("rdApp.*")


# Bio cofactors auto-hide. CoA, CO2, NAD(H), NADP(H), ATP, etc. are
# "released" or "consumed" in enzymatic steps — they appear as bulky
# extra molecule nodes on every bio reaction, which buries the actual
# substrate→intermediate chain in cofactor clutter. We import DORAnet's
# canonical cofactor SMILES set and union it with the user-supplied
# helpers list so bio cofactors are hidden by default. Starter and
# target are never hidden, even if they happen to be in the cofactor
# set (e.g. acetyl-CoA can be the starter).
try:
    from doranet.modules.enzymatic.generate_network import (
        cofactors_clean as _BIO_COFACTOR_SMILES,
    )
except Exception:
    _BIO_COFACTOR_SMILES = frozenset()


# Categorical palette for per-pathway edge colors. 20 visually-distinct
# hues; cycles if more pathways exist (rare in practice — most demos
# show ≤10 pathways).
_PATHWAY_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#bcbd22", "#17becf", "#7f7f7f",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#c7c7c7",
]
_DEFAULT_EDGE_COLOR = "#555555"
_DIMMED_EDGE_COLOR  = "#e6e6e6"


def _pathway_color_map(pathway_indices):
    """Assign each pathway index a stable color from the palette."""
    return {
        idx: _PATHWAY_PALETTE[i % len(_PATHWAY_PALETTE)]
        for i, idx in enumerate(sorted(pathway_indices))
    }


def _molecule_metadata_html(smi, pathways_member, img_uri=None):
    """
    Build an HTML tooltip for a node: the molecule's chemical structure
    (rendered image) on top, followed by SMILES, formula, MW, heavy-atom
    count, ring count, and the pathways that include this molecule.

    img_uri, if given, is a self-contained base64 PNG data URI of the
    structure (see _render_molecule_data_uris) — embedded so the drawing
    shows on hover even inside Streamlit's sandboxed iframe, where local
    file paths would not load.

    pathways_member is a tuple/list of pathway indices for which this
    molecule appears in at least one reaction (computed in
    _build_networkx_graph).
    """
    img_html = (
        f'<img src="{img_uri}" width="200" '
        f'style="display:block;margin:0 auto 6px;background:#fff;'
        f'border:1px solid #ddd;border-radius:4px;">'
        if img_uri else ""
    )
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return f"{img_html}<b>SMILES:</b> {smi}"
    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    heavy = mol.GetNumHeavyAtoms()
    rings = rdMolDescriptors.CalcNumRings(mol)
    pw_str = (
        ", ".join(str(i) for i in pathways_member)
        if pathways_member else "—"
    )
    return (
        f"{img_html}"
        f"<b>SMILES:</b> {smi}<br>"
        f"<b>Formula:</b> {formula}<br>"
        f"<b>MW:</b> {mw:.2f} g/mol<br>"
        f"<b>Heavy atoms:</b> {heavy}<br>"
        f"<b>Rings:</b> {rings}<br>"
        f"<b>In pathway(s):</b> {pw_str}"
    )


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------
def visualize_pathways(
    job_name: str,
    starter_smiles: str,
    target_smiles: str,
    starter_label: str = "starter",
    target_label: str = "target",
    helpers: Iterable[str] = ("O", "[H][H]"),
    pathway_filter: Union[str, list] = "shortest",
    output_html: Optional[str] = None,
    img_dir: str = "pathway_images",
    extra_labels: Optional[dict] = None,
    top_n_threshold: int = 20,
    top_n: int = 5,
) -> str:
    """
    Render the pathways in `{job_name}_pathways.txt` as an interactive
    layered HTML DAG.

    Parameters
    ----------
    job_name : str
        Prefix for both the input pathway file and the default output.
    starter_smiles, target_smiles : str
        SMILES of the starter (left side) and target (right side).
    starter_label, target_label : str
        Human-readable labels shown on the starter and target nodes.
        Intermediate molecules are labeled with their SMILES (truncated
        for display; full SMILES on hover).
    helpers : iterable of str
        SMILES of helper molecules (water, H2, etc.) — these are
        omitted from the visualization to reduce clutter. Their role
        in each step still appears in the edge tooltip.
    pathway_filter : "shortest" | "all" | list[int]
        Which pathways from the file to include.
    output_html : str, optional
        Output path. Defaults to "{job_name}_graph.html".
    img_dir : str
        Directory where per-molecule PNGs are written.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    helpers_set = set(helpers)
    if output_html is None:
        output_html = f"{job_name}_graph.html"

    # Canonicalize the user-supplied SMILES so they match exactly what
    # the pathway file stores. Different SMILES strings can refer to
    # the same molecule (CC(=O)CC(=O)C vs CC(=O)CC(C)=O for acetyl-
    # acetone, for example). RDKit canonicalization makes them match.
    starter_smiles = _canonicalize(starter_smiles, "starter_smiles")
    target_smiles = _canonicalize(target_smiles, "target_smiles")
    helpers_set = {_canonicalize(h, f"helper {h!r}") for h in helpers_set}

    # extra_labels maps SMILES -> human label. Canonicalize the keys so
    # lookup against pathway-file SMILES is reliable. Anything in this
    # dict is treated as "important and visible" — never hidden, even
    # if it would otherwise match the cofactor set (e.g. malonyl-CoA).
    if extra_labels is None:
        extra_labels = {}
    extra_labels_canon = {
        _canonicalize(s, "extra_labels key"): lbl
        for s, lbl in extra_labels.items()
    }

    # Extend helpers with DORAnet's bio cofactor SMILES — these (free
    # CoA, CO2, NAD(H), ATP, etc.) get released or consumed in every
    # bio reaction and clutter the graph as fake "destination" nodes
    # if not hidden. We exclude the starter and target so they remain
    # visible even when they are themselves cofactors (e.g. when the
    # starter IS acetyl-CoA).
    helpers_set |= set(_BIO_COFACTOR_SMILES)
    helpers_set.discard(starter_smiles)
    helpers_set.discard(target_smiles)
    for s in extra_labels_canon:
        helpers_set.discard(s)

    # Build a stereo-free fingerprint set of the helpers so we match
    # molecules whose stereo annotations differ from the cofactor
    # table's flat SMILES (the bio network propagates stereo from
    # stereo-bearing starters; DORAnet's cofactor table is stereo-free).
    # Starter, target, and any explicitly-labelled molecules are removed
    # so they stay visible even when their flat skeleton matches a
    # cofactor.
    helper_fingerprints = {_stereo_free_canon(h) for h in helpers_set}
    helper_fingerprints.discard(_stereo_free_canon(starter_smiles))
    helper_fingerprints.discard(_stereo_free_canon(target_smiles))
    for s in extra_labels_canon:
        helper_fingerprints.discard(_stereo_free_canon(s))

    # Load and filter pathways
    pathways = load_pathways_from_file(job_name)
    selected = _select_pathways(pathways, pathway_filter)

    # Step 1: build the graph in NetworkX
    G = _build_networkx_graph(
        selected, starter_smiles, helpers_set,
        target_smiles=target_smiles,
        helper_fingerprints=helper_fingerprints,
    )

    # Step 2: compute positions in NetworkX
    positions = _compute_positions(G)

    # Render molecule structures as self-contained base64 PNG data URIs.
    # Used for BOTH the node thumbnails and the hover tooltips, so both
    # render reliably inside Streamlit's sandboxed iframe.
    img_data_uris = _render_molecule_data_uris(set(G.nodes))

    # Step 3: hand to PyVis purely as renderer
    _render_to_pyvis(
        G,
        positions,
        img_data_uris,
        starter_smiles=starter_smiles,
        target_smiles=target_smiles,
        starter_label=starter_label,
        target_label=target_label,
        output_html=output_html,
        extra_labels=extra_labels_canon,
    )

    # ── supplementary "top N" graph ────────────────────────────────────
    # When the full pathway file has more than `top_n_threshold` paths,
    # the all-pathways graph becomes visually overwhelming. We also
    # render a stripped-down version showing only the `top_n` shortest
    # pathways — easier to scan, useful for quick-look summaries and
    # CEO/mentor demos. The supplementary file lives alongside the main
    # output as `{output_html}` → `{output_html-base}_top{N}.html`.
    #
    # Ranking: shortest path first (by step count). This matches what a
    # chemist would call "the most direct route" without requiring the
    # user to configure weight values. Swap in a more sophisticated
    # ranker later if needed (pathway_scoring.WeightedPathwayScorer).
    if len(pathways) > top_n_threshold:
        ranked = sorted(
            enumerate(pathways, 1), key=lambda x: x[1].num_steps
        )
        top_indices = [idx for idx, _ in ranked[:top_n]]
        base, ext = os.path.splitext(output_html)
        top_html = f"{base}_top{top_n}{ext}"
        print(
            f"  [supplementary] {len(pathways)} pathways exceeds "
            f"threshold {top_n_threshold}; also rendering top {top_n} "
            f"(shortest) → {top_html}"
        )
        # Recursive call with disabled threshold to prevent infinite
        # recursion, and `pathway_filter` set to the top-N indices.
        visualize_pathways(
            job_name=job_name,
            starter_smiles=starter_smiles,
            target_smiles=target_smiles,
            starter_label=starter_label,
            target_label=target_label,
            helpers=helpers,
            pathway_filter=top_indices,
            output_html=top_html,
            img_dir=img_dir,
            extra_labels=extra_labels,
            top_n_threshold=10**9,   # don't recurse further
            top_n=top_n,
        )

    return os.path.abspath(output_html)


# ----------------------------------------------------------------------
# Internals
# ----------------------------------------------------------------------
def _canonicalize(smiles, label="smiles"):
    """
    Convert any SMILES into RDKit's canonical form.  Two SMILES that
    refer to the same molecule (CC(=O)CC(=O)C vs CC(=O)CC(C)=O) will
    produce the same canonical string, so node lookup against the
    pathway file's reactions becomes reliable.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"{label}: could not parse SMILES {smiles!r}")
    return Chem.MolToSmiles(mol)


def _stereo_free_canon(smiles):
    """
    Canonical SMILES with stereochemistry stripped.

    Use case: DORAnet's cofactor table stores stereo-free SMILES
    (`...OCC1OC(n2cnc3...`), while reaction networks built FROM real
    stereo-bearing starters (e.g. acetyl-CoA) emit stereo-annotated
    SMILES (`...OC[C@H]1O[C@@H](n2cnc3...`). Stripping stereo on both
    sides lets us match "this is the same molecule" without tracking
    explicit stereo variants of every cofactor.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)


def _select_pathways(pathways, filter_spec):
    if filter_spec == "all":
        return list(enumerate(pathways, 1))
    if filter_spec == "shortest":
        if not pathways:
            return []
        min_steps = min(p.num_steps for p in pathways)
        return [(i, p) for i, p in enumerate(pathways, 1)
                if p.num_steps == min_steps]
    if isinstance(filter_spec, list):
        return [(i, p) for i, p in enumerate(pathways, 1)
                if i in filter_spec]
    raise ValueError(f"Unknown pathway_filter: {filter_spec!r}")


def _build_networkx_graph(selected_pathways, starter, helpers,
                           *, target_smiles=None,
                           helper_fingerprints=None):
    """Build a DiGraph with unique molecules + unique reaction edges.

    Each edge carries:
        op          : reaction name
        pathways    : tuple of pathway indices using this edge
        helper_in   : tuple of helper SMILES consumed (visualised in tooltip)
        helper_out  : tuple of helper SMILES produced

    Each node carries:
        smiles      : the molecule SMILES (same as node id)
        level       : longest-path depth from starter (column index in layout)
        pathways    : tuple of pathway indices this molecule appears in,
                       computed as the union over its incident edges'
                       `pathways`. Used for click-highlight and tooltips.

    target_smiles, if provided, is bumped to a level greater than every
    other node so it always renders as the rightmost column.

    Note on multigraph: we use MultiDiGraph so two *different* reactions
    that happen to connect the same (r, p) pair stay as parallel edges
    rather than being collapsed. Collapsing them would force their
    pathway memberships to be unioned, which made click-highlight
    appear to light up unrelated pathways (an edge labelled with
    reaction X carried memberships from reactions X, Y, Z if X, Y, Z
    all coincidentally had the same r→p endpoints).
    """
    G = nx.MultiDiGraph()

    # A molecule SMILES is a helper if its raw form appears in the
    # explicit helpers set OR if its stereo-stripped canonical form
    # matches one of the helper fingerprints (covers cofactors emitted
    # by bio reactions with stereo annotations that differ from
    # DORAnet's flat table form).
    _fp_set = set(helper_fingerprints or ())
    def is_helper(smi):
        if smi in helpers:
            return True
        if _fp_set and _stereo_free_canon(smi) in _fp_set:
            return True
        return False

    # ONE EDGE PER REACTION, not one per (reactant × product) pair.
    # A Claisen condensation with 2 real reactants and 2 real products
    # otherwise fans out into 4 arrows for a single chemistry step,
    # and clicking any one of them lights up the same pathway across
    # all of them — making it look like the pathway "starts from
    # multiple molecules" or "has more arrows than steps." We pick the
    # FIRST non-helper reactant as the primary substrate and the FIRST
    # non-helper product as the primary product (DORAnet's SMARTS rules
    # are written so the main reacting atoms come first on each side).
    # Co-substrates and byproducts are recorded for the tooltip but do
    # not get their own arrows in the graph.
    edge_records = {}    # (u, v, op) -> {pathways:set, co_subs, byproducts}
    visible_mols = set()
    for idx, p in selected_pathways:
        for rxn_str in p.reactions:
            parsed = parse_reaction_string(rxn_str)
            reactants_real = [r for r in parsed["reactants"] if not is_helper(r)]
            products_real  = [pp for pp in parsed["products"] if not is_helper(pp)]
            reactants_help = [r for r in parsed["reactants"] if is_helper(r)]
            products_help  = [pp for pp in parsed["products"] if is_helper(pp)]
            if not reactants_real or not products_real:
                # Reaction where every reactant or every product is a
                # cofactor: skip to avoid an orphan edge.
                continue
            primary_in  = reactants_real[0]
            primary_out = products_real[0]
            visible_mols.add(primary_in)
            visible_mols.add(primary_out)
            key = (primary_in, primary_out, parsed["op_name"])
            rec = edge_records.get(key)
            if rec is None:
                edge_records[key] = {
                    "pathways":   {idx},
                    "co_subs":    tuple(reactants_real[1:] + reactants_help),
                    "byproducts": tuple(products_real[1:]  + products_help),
                }
            else:
                rec["pathways"].add(idx)

    for smi in visible_mols:
        G.add_node(smi, smiles=smi)

    for (u, v, op), data in edge_records.items():
        G.add_edge(
            u, v,
            op=op,
            pathways=tuple(sorted(data["pathways"])),
            helper_in=data["co_subs"],
            helper_out=data["byproducts"],
        )

    # 'level' attribute = longest path from starter (for layered layout)
    levels = _compute_node_levels(G, starter)
    for node, lvl in levels.items():
        G.nodes[node]["level"] = lvl
    for node in G.nodes:
        if "level" not in G.nodes[node]:
            G.nodes[node]["level"] = 0

    # Pin target to a level above every other node so it always lands
    # in the rightmost column. Without this, off-path chemistry branches
    # can end up at deeper levels than the target (e.g. when the network
    # includes byproduct chains that extend past the target's depth),
    # which pushes the target inward.
    if target_smiles is not None and target_smiles in G.nodes:
        other_max = max(
            (G.nodes[n]["level"] for n in G.nodes if n != target_smiles),
            default=0,
        )
        G.nodes[target_smiles]["level"] = other_max + 1

    # Per-node pathway membership: union of `pathways` across incident
    # edges. Used by the click-highlight handler in the rendered HTML.
    for node in G.nodes:
        member = set()
        for _, _, d in G.in_edges(node, data=True):
            member.update(d.get("pathways", ()))
        for _, _, d in G.out_edges(node, data=True):
            member.update(d.get("pathways", ()))
        G.nodes[node]["pathways"] = tuple(sorted(member))

    return G


def _compute_node_levels(G, starter):
    """Longest-path levels from the starter via topological walk."""
    levels = {starter: 0}
    try:
        for node in nx.topological_sort(G):
            if node == starter:
                continue
            preds = list(G.predecessors(node))
            if not preds:
                levels.setdefault(node, 0)
                continue
            pred_levels = [levels.get(p, 0) for p in preds]
            levels[node] = max(pred_levels) + 1
    except nx.NetworkXUnfeasible:
        # Cycle in the DAG (shouldn't happen for reaction pathways but
        # just in case) -- fall back to BFS shortest path
        for node in G.nodes:
            if nx.has_path(G, starter, node):
                levels[node] = nx.shortest_path_length(G, starter, node)
            else:
                levels[node] = 0
    return levels


def _compute_positions(G, scale_x=320, scale_y=180):
    pos = nx.multipartite_layout(G, subset_key="level", align="vertical")
    return {
        node: (x * scale_x * 4, y * scale_y * 4)
        for node, (x, y) in pos.items()
    }


def _safe_name(smiles):
    return hashlib.md5(smiles.encode("utf-8")).hexdigest()[:10] + ".png"


def _render_molecule_pngs(smiles_set, out_dir, size=(160, 160)):
    os.makedirs(out_dir, exist_ok=True)
    paths = {}
    for smi in smiles_set:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        full = os.path.join(out_dir, _safe_name(smi))
        try:
            Draw.MolToImage(mol, size=size).save(full)
            paths[smi] = full.replace("\\", "/")
        except Exception:
            continue
    return paths


def _render_molecule_data_uris(smiles_set, size=(200, 200)):
    """Render each molecule to a base64-encoded PNG *data URI*.

    Unlike _render_molecule_pngs (which writes files and returns paths),
    these URIs are self-contained, so the images render both as node
    thumbnails AND inside hover tooltips — and they survive being
    embedded in Streamlit's sandboxed iframe, where relative file paths
    silently fail to load. Also makes the downloaded HTML fully portable
    (no pathway_images/ folder needed).
    """
    uris = {}
    for smi in smiles_set:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        try:
            img = Draw.MolToImage(mol, size=size)
            buf = BytesIO()
            img.save(buf, format="PNG")
            b64 = base64.b64encode(buf.getvalue()).decode("ascii")
            uris[smi] = "data:image/png;base64," + b64
        except Exception:
            continue
    return uris


def _short_smiles(smi, max_len=22):
    """Truncate long SMILES so they fit cleanly under a node."""
    return smi if len(smi) <= max_len else smi[:max_len - 3] + "..."


def _render_to_pyvis(
    G,
    positions,
    img_data_uris,
    *,
    starter_smiles,
    target_smiles,
    starter_label,
    target_label,
    output_html,
    extra_labels=None,
):
    extra_labels = extra_labels or {}
    net = Network(
        height="800px",
        width="100%",
        directed=True,
        bgcolor="#ffffff",
        font_color="#222",
        notebook=False,
    )

    # Collect every pathway index present in the graph so we can assign
    # each one a distinct palette color. Both edges and nodes contribute
    # (every pathway must touch at least one edge, but a node-only walk
    # is cheap and future-proof).
    all_pathways = set()
    for _, _, d in G.edges(data=True):
        all_pathways.update(d.get("pathways", ()))
    for n in G.nodes:
        all_pathways.update(G.nodes[n].get("pathways", ()))
    pathway_colors = _pathway_color_map(all_pathways)

    for node in G.nodes:
        smi = node
        x, y = positions.get(smi, (0, 0))
        node_pathways = G.nodes[smi].get("pathways", ())

        # Decide label + style based on node role
        if smi == starter_smiles:
            label = f"{starter_label}\n(start)"
            border_color = "#28a745"
            border_width = 8
            size = 55
            font_color = "#0d5320"
        elif smi == target_smiles:
            label = f"{target_label}\n(END)"
            border_color = "#dc3545"
            border_width = 8
            size = 55
            font_color = "#7a1c25"
        elif smi in extra_labels:
            # Explicitly-named co-substrate (e.g. malonyl-CoA). Distinct
            # styling so the user knows it's a named/important input
            # rather than an unknown intermediate.
            label = extra_labels[smi]
            border_color = "#f0ad4e"          # amber border
            border_width = 5
            size = 45
            font_color = "#7a4d05"
        else:
            label = _short_smiles(smi)
            border_color = "#3a76c2"
            border_width = 2
            size = 35
            font_color = "#333"

        kwargs = dict(
            label=label,
            title=_molecule_metadata_html(
                smi, node_pathways, img_uri=img_data_uris.get(smi)
            ),
            x=float(x),
            y=float(y),
            physics=False,
            fixed={"x": True, "y": True},
            color={"border": border_color, "background": "#ffffff"},
            borderWidth=border_width,
            size=size,
            font={"size": 16, "color": font_color, "strokeWidth": 0,
                  "multi": True},
        )
        if smi in img_data_uris:
            kwargs["shape"] = "image"
            kwargs["image"] = img_data_uris[smi]
        net.add_node(smi, **kwargs)

    # Build edges with stable IDs so the JS click handler can correlate
    # an edge-click event to our (edge_id -> pathways) lookup table.
    # G is a MultiDiGraph so parallel edges between the same (u, v) pair
    # are walked separately. Edge IDs are simple counters — unique per
    # render, which is all the JS needs.
    #
    # PARALLEL EDGES: when more than one reaction connects the same
    # molecule pair, those edges share the same (u, v) endpoints in
    # vis.js. With the default cubicBezier smoothing they would draw
    # exactly on top of each other, so the user can't tell which is
    # which AND vis.js's hit-detection picks an arbitrary parallel on
    # click — which is what made the click-highlight look "wrong"
    # (clicked-on edge wasn't the one whose pathway lit up). We fix it
    # by assigning each parallel a different curve direction + roundness
    # so they fan out into visibly distinct arcs.
    from collections import defaultdict
    parallels_for_pair = defaultdict(int)
    for u, v, _k in G.edges(keys=True):
        parallels_for_pair[(u, v)] += 1
    parallel_index = defaultdict(int)

    edge_pathway_map = {}      # edge_id -> [pathway indices]
    for edge_idx, (u, v, _k, data) in enumerate(
        G.edges(keys=True, data=True)
    ):
        edge_id = f"e{edge_idx}"
        op = data.get("op", "?")
        pathways = data.get("pathways", ())
        helper_in = data.get("helper_in", ())
        helper_out = data.get("helper_out", ())
        helper_note = ""
        if helper_in:
            helper_note += " (+ " + ", ".join(helper_in) + ")"
        if helper_out:
            helper_note += " (- " + ", ".join(helper_out) + ")"
        pathway_label = ", ".join(str(i) for i in pathways)
        tooltip = (
            f"<b>{op}</b>{helper_note}<br>"
            f"In pathway(s): {pathway_label}"
        )
        edge_pathway_map[edge_id] = list(pathways)

        n_parallels = parallels_for_pair[(u, v)]
        i_in_group = parallel_index[(u, v)]
        parallel_index[(u, v)] += 1

        edge_kwargs = dict(
            id=edge_id,
            # NO label= — reaction name is hover-only via title=
            title=tooltip,
            arrows="to",
            color={"color": _DEFAULT_EDGE_COLOR,
                   "highlight": _DEFAULT_EDGE_COLOR},
            width=2,
        )
        # For singletons, use the global cubicBezier smoothing. For
        # parallel siblings, override with alternating CW/CCW curves at
        # increasing roundness so they fan apart visibly.
        if n_parallels > 1:
            direction = "curvedCW" if i_in_group % 2 == 0 else "curvedCCW"
            roundness = 0.20 + 0.18 * (i_in_group // 2)  # 0.20, 0.20, 0.38, 0.38, 0.56, ...
            edge_kwargs["smooth"] = {
                "enabled": True,
                "type": direction,
                "roundness": roundness,
            }
        net.add_edge(u, v, **edge_kwargs)

    # Layout disabled; positions are explicit. Physics off.
    # selectConnectedEdges = false: critical for our click-highlight to
    # look right. With it true, clicking an edge causes vis.js to also
    # mark every edge sharing an endpoint with the clicked edge as
    # "selected", and the selection styling visually re-colors those
    # edges. They then look "lit up" alongside the pathway we actually
    # want to highlight, which made the colors look inverted/random.
    net.set_options("""
    {
      "layout": { "improvedLayout": false },
      "physics": { "enabled": false },
      "edges": {
        "smooth": { "type": "cubicBezier", "forceDirection": "horizontal",
                     "roundness": 0.35 },
        "selectionWidth": 0,
        "chosen": false
      },
      "nodes": {
        "borderWidthSelected": 8,
        "shapeProperties": { "useImageSize": false, "interpolation": false }
      },
      "interaction": {
        "dragNodes": false,
        "dragView": true,
        "zoomView": true,
        "hover": true,
        "tooltipDelay": 100,
        "navigationButtons": true,
        "keyboard": true,
        "selectConnectedEdges": false
      }
    }
    """)

    net.write_html(output_html, notebook=False)

    # Inject the click-highlight JS. We append to the HTML after PyVis
    # has written it, so the script has access to the `network` and
    # `edges`/`nodes` DataSet objects defined by PyVis's template.
    node_pathway_map = {
        n: list(G.nodes[n].get("pathways", ())) for n in G.nodes
    }
    _inject_click_highlight(
        output_html,
        edge_pathway_map=edge_pathway_map,
        node_pathway_map=node_pathway_map,
        pathway_colors=pathway_colors,
    )


def _inject_click_highlight(
    output_html, *, edge_pathway_map, node_pathway_map, pathway_colors,
):
    """
    Append a <script> block to the pyvis HTML that wires up click-to-
    highlight. Behaviour:

      * Click an edge → highlight every pathway that edge belongs to.
        Edges in exactly one of the highlighted pathways take that
        pathway's color. Edges in multiple of the highlighted pathways
        take the average (RGB mean) of those colors — so a step shared
        by two pathways visibly carries both. Edges not in any of the
        clicked element's pathways dim to light grey.

      * Click a node → same logic, using the union of pathways the node
        appears in.

      * Click empty canvas → reset all edges to the default color.

    The data needed (which pathways each edge/node is in, the color
    assigned to each pathway) is baked into the HTML at render time, so
    no runtime computation is required beyond intersection + averaging.
    """
    payload = {
        "edges":   edge_pathway_map,
        "nodes":   node_pathway_map,
        "colors":  {str(k): v for k, v in pathway_colors.items()},
        "default": _DEFAULT_EDGE_COLOR,
        "dim":     _DIMMED_EDGE_COLOR,
    }
    payload_json = json.dumps(payload)
    script = f"""
<script type="text/javascript">
(function() {{
  // Wait for pyvis to finish constructing `network`. The pyvis template
  // assigns it as a global, but it may take a tick after DOMContentLoaded.
  function ready(cb) {{
    if (typeof network !== 'undefined' && network && network.body) {{
      cb();
    }} else {{
      setTimeout(function() {{ ready(cb); }}, 50);
    }}
  }}

  ready(function() {{
    var PAYLOAD = {payload_json};

    function mixHex(hexes) {{
      if (hexes.length === 1) return hexes[0];
      var r = 0, g = 0, b = 0;
      for (var i = 0; i < hexes.length; i++) {{
        var h = hexes[i];
        r += parseInt(h.substr(1, 2), 16);
        g += parseInt(h.substr(3, 2), 16);
        b += parseInt(h.substr(5, 2), 16);
      }}
      r = Math.round(r / hexes.length);
      g = Math.round(g / hexes.length);
      b = Math.round(b / hexes.length);
      function pad(x) {{ var s = x.toString(16); return s.length < 2 ? '0' + s : s; }}
      return '#' + pad(r) + pad(g) + pad(b);
    }}

    function colorEdgeForPathways(edgePathways, highlightSet) {{
      // Intersect the edge's pathways with the click's highlight set.
      // Pure intersection = pathways the user explicitly asked about.
      var hits = edgePathways.filter(function(p) {{
        return highlightSet.indexOf(p) >= 0;
      }});
      if (hits.length === 0) return PAYLOAD.dim;
      if (hits.length === 1) return PAYLOAD.colors[String(hits[0])];
      var hexes = hits.map(function(p) {{ return PAYLOAD.colors[String(p)]; }});
      return mixHex(hexes);
    }}

    function applyHighlight(highlightPathways) {{
      // highlightPathways: array of pathway indices to light up.
      // null/empty = reset everything to default.
      var updates = [];
      Object.keys(PAYLOAD.edges).forEach(function(eid) {{
        if (!highlightPathways || highlightPathways.length === 0) {{
          updates.push({{ id: eid, color: {{ color: PAYLOAD.default }} }});
        }} else {{
          var c = colorEdgeForPathways(PAYLOAD.edges[eid], highlightPathways);
          updates.push({{ id: eid, color: {{ color: c }} }});
        }}
      }});
      network.body.data.edges.update(updates);
    }}

    // Add a visible diagnostic panel so we don't need dev tools.
    var infoDiv = document.createElement('div');
    infoDiv.id = 'click-info';
    infoDiv.style.cssText = 'position:fixed;top:10px;right:10px;'
      + 'background:#fff;padding:10px 12px;border:1px solid #888;'
      + 'border-radius:4px;font-family:monospace;font-size:12px;'
      + 'max-width:380px;z-index:9999;box-shadow:0 2px 6px rgba(0,0,0,0.15);';
    infoDiv.innerHTML = '<b>Click info</b><br><span style="color:#888">'
      + '(click any edge / node / blank space to inspect)</span>';
    document.body.appendChild(infoDiv);

    function describeClick(kind, id, pathways) {{
      var pwStr = pathways.length ? pathways.join(', ') : '(none)';
      var html = '<b>Click info</b><br>'
        + '<b>kind:</b> ' + kind + '<br>'
        + '<b>id:</b> ' + (id.length > 50 ? id.substring(0, 47) + '...' : id) + '<br>'
        + '<b>pathways being highlighted:</b> [' + pwStr + ']';
      infoDiv.innerHTML = html;
    }}

    network.on('click', function(params) {{
      console.log('[click] params:', params);
      if (params.nodes && params.nodes.length > 0) {{
        var nid = params.nodes[0];
        var pw = PAYLOAD.nodes[nid] || [];
        console.log('[click] node:', nid, '-> pathways:', pw);
        describeClick('NODE', nid, pw);
        applyHighlight(pw);
      }} else if (params.edges && params.edges.length > 0) {{
        var eid = params.edges[0];
        var pw = PAYLOAD.edges[eid] || [];
        console.log('[click] edge:', eid, '-> pathways:', pw);
        describeClick('EDGE', eid, pw);
        applyHighlight(pw);
      }} else {{
        console.log('[click] blank canvas - reset');
        infoDiv.innerHTML = '<b>Click info</b><br><span style="color:#888">'
          + '(reset)</span>';
        applyHighlight(null);
      }}
    }});
  }});
}})();
</script>
"""
    with open(output_html, "r", encoding="utf-8") as f:
        html = f.read()
    if "</body>" in html:
        html = html.replace("</body>", script + "\n</body>")
    else:
        html = html + script
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html)
