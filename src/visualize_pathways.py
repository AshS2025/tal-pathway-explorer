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

import hashlib
import os
from typing import Iterable, Optional, Union

import networkx as nx
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw
from pyvis.network import Network

from pathway_tools import load_pathways_from_file, parse_reaction_string

RDLogger.DisableLog("rdApp.*")


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

    # Load and filter pathways
    pathways = load_pathways_from_file(job_name)
    selected = _select_pathways(pathways, pathway_filter)

    # Step 1: build the graph in NetworkX
    G = _build_networkx_graph(selected, starter_smiles, helpers_set)

    # Step 2: compute positions in NetworkX
    positions = _compute_positions(G)

    # Render molecule images (RDKit -> PNG)
    img_paths = _render_molecule_pngs(set(G.nodes), img_dir)

    # Step 3: hand to PyVis purely as renderer
    _render_to_pyvis(
        G,
        positions,
        img_paths,
        starter_smiles=starter_smiles,
        target_smiles=target_smiles,
        starter_label=starter_label,
        target_label=target_label,
        output_html=output_html,
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


def _build_networkx_graph(selected_pathways, starter, helpers):
    """Build a DiGraph with unique molecules + unique reaction edges."""
    G = nx.DiGraph()

    reactions = {}
    pathway_membership = {}
    for idx, p in selected_pathways:
        for rxn_str in p.reactions:
            parsed = parse_reaction_string(rxn_str)
            key = (
                frozenset(parsed["reactants"]),
                parsed["op_name"],
                frozenset(parsed["products"]),
            )
            reactions.setdefault(key, parsed)
            pathway_membership.setdefault(key, set()).add(idx)

    # Visible molecules = everything in any reaction, minus helpers
    visible_mols = set()
    for parsed in reactions.values():
        for m in parsed["reactants"]:
            if m not in helpers:
                visible_mols.add(m)
        for m in parsed["products"]:
            if m not in helpers:
                visible_mols.add(m)
    for smi in visible_mols:
        G.add_node(smi, smiles=smi)

    # Edges (helpers omitted from the graph but recorded in tooltips)
    for key, parsed in reactions.items():
        op = parsed["op_name"]
        path_indices = tuple(sorted(pathway_membership[key]))
        helper_in = [r for r in parsed["reactants"] if r in helpers]
        helper_out = [p for p in parsed["products"] if p in helpers]
        reactants = [r for r in parsed["reactants"] if r not in helpers]
        products = [p for p in parsed["products"] if p not in helpers]
        for r in reactants:
            for p in products:
                if G.has_edge(r, p):
                    G[r][p]["pathways"] = tuple(
                        sorted(set(G[r][p]["pathways"]) | set(path_indices))
                    )
                else:
                    G.add_edge(
                        r, p,
                        op=op,
                        pathways=path_indices,
                        helper_in=tuple(helper_in),
                        helper_out=tuple(helper_out),
                    )

    # 'level' attribute = longest path from starter (for layered layout)
    levels = _compute_node_levels(G, starter)
    for node, lvl in levels.items():
        G.nodes[node]["level"] = lvl
    for node in G.nodes:
        if "level" not in G.nodes[node]:
            G.nodes[node]["level"] = 0

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


def _short_smiles(smi, max_len=22):
    """Truncate long SMILES so they fit cleanly under a node."""
    return smi if len(smi) <= max_len else smi[:max_len - 3] + "..."


def _render_to_pyvis(
    G,
    positions,
    img_paths,
    *,
    starter_smiles,
    target_smiles,
    starter_label,
    target_label,
    output_html,
):
    net = Network(
        height="800px",
        width="100%",
        directed=True,
        bgcolor="#ffffff",
        font_color="#222",
        notebook=False,
    )

    for node in G.nodes:
        smi = node
        x, y = positions.get(smi, (0, 0))

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
        else:
            label = _short_smiles(smi)
            border_color = "#3a76c2"
            border_width = 2
            size = 35
            font_color = "#333"

        kwargs = dict(
            label=label,
            title=smi,                # full SMILES on hover
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
        if smi in img_paths:
            kwargs["shape"] = "image"
            kwargs["image"] = img_paths[smi]
        net.add_node(smi, **kwargs)

    for u, v, data in G.edges(data=True):
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
        net.add_edge(
            u, v,
            label=op if len(op) <= 32 else op[:30] + "...",
            title=tooltip,
            arrows="to",
            color={"color": "#555", "highlight": "#dc3545"},
            font={"size": 11, "color": "#444", "align": "middle",
                  "background": "white"},
            width=2,
        )

    # Layout disabled; positions are explicit. Physics off.
    net.set_options("""
    {
      "layout": { "improvedLayout": false },
      "physics": { "enabled": false },
      "edges": {
        "smooth": { "type": "cubicBezier", "forceDirection": "horizontal",
                     "roundness": 0.35 }
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
        "selectConnectedEdges": true
      }
    }
    """)

    net.write_html(output_html, notebook=False)
