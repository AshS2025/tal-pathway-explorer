# src/pathway_tools.py

from dataclasses import dataclass
from typing import Optional

from doranet.modules.post_processing.post_processing import (
    pretreat_networks,
    pathway_finder,
)
from doranet.modules.post_processing import post_processing as _doranet_pp
from doranet.modules.enzymatic.generate_network import AVAILABLE_RULESETS as _BIO_RULESETS

from rdkit import Chem
import json
import pandas as _pd


# ---------------------------------------------------------------------
# Bio rule name patch
# ---------------------------------------------------------------------
# DORAnet's pathway_finder detects bio reactions by checking if each
# reaction's operator name is in `post_processing.bio_rxn_names`, which
# at import time is populated ONLY from JN3604IMT_rules.tsv. Our wrapper
# expands using JN1224MIN, whose rule names (`rule0087`, `rule1118`,
# `rule0891`, …) are largely absent from JN3604IMT. Without this patch
# the pathway finder sees our bio reactions as chem reactions, never
# auto-adds cofactors to the helper set, and counts cofactor mass in
# atom economy — which kills the 0.3 atom-economy threshold for every
# CoA-tethered step and returns "No pathway found" even when the target
# is unambiguously in the network.
#
# Fix: union JN1224MIN's rule name set into `bio_rxn_names` once at
# import time. Idempotent: re-running just no-ops on already-present names.
def _patch_bio_rxn_names() -> int:
    added = 0
    for ruleset_path in _BIO_RULESETS.values():
        try:
            df = _pd.read_csv(ruleset_path, sep="\t")
        except FileNotFoundError:
            continue
        for name in df["Name"]:
            if name not in _doranet_pp.bio_rxn_names:
                _doranet_pp.bio_rxn_names.add(name)
                added += 1
    return added


_BIO_RXN_NAMES_ADDED = _patch_bio_rxn_names()


# =====================================================================
# Pathway data model
# =====================================================================

@dataclass
class Pathway:
    """
    A single candidate pathway parsed out of `{job_name}_pathways.txt`.

    `reactions` follows the same string convention DORAnet uses
    internally after `pathway_finder` runs:

        "reactants>op_name>dH$rea_stoi$pro_stoi>products"

    where reactants/products are dot-joined canonical SMILES, op_name is
    the operator that produced the reaction, dH is either a float (as a
    string) or the literal "No_Thermo", and the stoichiometry pieces are
    Python list literals (e.g. "[1, 1]") to be parsed with `eval`.

    Keeping the raw strings (rather than fully parsed structs) means a
    Pathway can be passed straight to anything that already speaks
    DORAnet's pathway format — including DORAnet's own helpers.
    """

    index: int                 # 1-based pathway number from the txt file
    reactions: list[str]

    @property
    def num_steps(self) -> int:
        return len(self.reactions)


def parse_reaction_string(rxn: str) -> dict:
    """
    Split a pathway reaction string into its pieces. Returns a dict with:
        reactants : list[str]   (canonical SMILES)
        op_name   : str
        dH        : float | None  (None when "No_Thermo")
        rea_stoi  : list[int]
        pro_stoi  : list[int]
        products  : list[str]
    """
    reactants_part, op_name, meta_part, products_part = rxn.split(">")
    dH_str, rea_stoi_str, pro_stoi_str = meta_part.split("$")
    dH = None if dH_str == "No_Thermo" else float(dH_str)
    return {
        "reactants": reactants_part.split("."),
        "op_name": op_name,
        "dH": dH,
        "rea_stoi": eval(rea_stoi_str),
        "pro_stoi": eval(pro_stoi_str),
        "products": products_part.split("."),
    }


def load_pathways_from_file(job_name: str) -> list[Pathway]:
    """
    Read `{job_name}_pathways.txt` (produced by `pathway_finder`) and
    rebuild a list of `Pathway` objects.

    Mirrors the parsing block inside DORAnet's `pathway_ranking`
    (post_processing.py:1465-1527) — same offsets, same reaction string
    reconstruction — so any pathway DORAnet recognizes, we recognize.
    """
    path = f"{job_name}_pathways.txt"
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()

    clean = [ln.strip() for ln in lines if ln != "\n"]

    markers = [idx for idx, ln in enumerate(clean) if "pathway number" in ln]

    pathways: list[Pathway] = []
    for path_idx, marker in enumerate(markers):
        # Stoichiometry line lives 4 lines below the marker; the first
        # 30 chars are a fixed label that we strip off before eval.
        stoi_list = eval(clean[marker + 4][30:])  # list[str]

        next_marker = markers[path_idx + 1] if path_idx + 1 < len(markers) else len(clean)
        block = next_marker - (marker + 6)
        step = block // 3  # SMILES, names, enthalpies each occupy `step` lines

        smiles_block = clean[marker + 6 : marker + 6 + step]
        name_block = clean[marker + 6 + step : marker + 6 + 2 * step]
        enthalpy_block = clean[marker + 6 + 2 * step : marker + 6 + 3 * step]

        reactions = []
        for i, smi_line in enumerate(smiles_block):
            reactants, products = smi_line.split(">>")
            rxn_string = (
                f"{reactants}>{name_block[i]}>{enthalpy_block[i]}"
                f"${stoi_list[i]}>{products}"
            )
            reactions.append(rxn_string)

        pathways.append(Pathway(index=path_idx + 1, reactions=reactions))

    return pathways


# =====================================================================
# Original pathway-generation entry point (unchanged)
# =====================================================================

def find_pathways_to_target(
    network,
    starter,
    target,
    helpers,
    generations,
    max_num_rxns,
    job_name,
):

    pretreat_networks(
        networks=[network],
        starters=[starter],
        helpers=helpers,
        total_generations=generations,
        job_name=job_name,
    )

    pretreated = json.load(
        open(f"{job_name}_network_pretreated.json")
    )

    print(f"Reactions in pretreated network: {len(pretreated)}")

    target_canonical = Chem.MolToSmiles(
        Chem.MolFromSmiles(target)
    )

    all_products = set()
    for rxn in pretreated:
        for p in rxn.split(">")[3].split("."):
            all_products.add(p)

    print(
        f"Target in pretreated network products: "
        f"{target_canonical in all_products}"
    )

    pathway_finder(
        starters=[starter],
        helpers=helpers,
        target=[target_canonical],
        search_depth=generations,
        max_num_rxns=max_num_rxns,
        job_name=job_name,
    )


# =====================================================================
# Open-exploration entry point
# =====================================================================

def explore_downstream(
    *,
    starter: str,
    helpers: str,
    job_name: str,
    network_kwargs: dict,
    top_n: int = 20,
    max_num_rxns: int = 4,
    scoring_kwargs: dict | None = None,
    derivatives_list: list[dict] | None = None,
):
    """
    Open-exploration entry point. Answers the brief's question:
    "What products can be feasibly produced with TAL as starting point?"

    Pipeline:
      1. Build a forward network from `starter` via generate_network_tal
         (no Product_Tanimoto_Filter — chemistry roams).
      2. Score every endpoint molecule with endpoint_scoring.
      3. Pick top_n endpoints.
      4. Run find_pathways_to_target for each, with a per-endpoint
         sub-job-name so the pretreated/pathway artifacts don't collide.
      5. (Optional) cross-check against `derivatives_list`: report which
         known literature derivatives appeared in the network.
      6. Write a consolidated markdown report.

    Parameters
    ----------
    starter, helpers : str
        Paths to SMILES files (same convention as generate_network_tal).
    job_name : str
        Prefix for all output artifacts. Per-endpoint runs use
        f"{job_name}_ep{idx:02d}".
    network_kwargs : dict
        Kwargs forwarded to generate_network_tal. Must include
        `gen`, `direction`, etc. `min_product_similarity` will be
        forced to None — exploration must not bias toward a target.
    top_n : int
        Number of top endpoints to pull pathways for.
    max_num_rxns : int
        Max steps in any single pathway returned by pathway_finder.
    scoring_kwargs : dict or None
        Forwarded to score_interestingness (carbon_window etc.).
    derivatives_list : list of dict or None
        Optional cross-check list; each dict needs "name" and "smiles".
        TAL_DOWNSTREAM_DERIVATIVES from
        src/tal_downstream_derivatives.py is the default canonical set.

    Returns
    -------
    dict with keys:
        top_endpoints : list[EndpointScore]
        pathways      : {smiles: list[Pathway]}
        derivative_matches  : list[dict] | None
        derivative_missing  : list[dict] | None
    """
    from network_generation import generate_network_tal, get_smiles_from_file
    from endpoint_scoring import rank_network_endpoints

    # Force-disable target similarity gating — exploration must not
    # bias toward any specific molecule.
    network_kwargs = dict(network_kwargs)
    network_kwargs["min_product_similarity"] = None
    network_kwargs.setdefault("targets", "")  # no target file needed

    print(f"\n{'='*60}")
    print(f"EXPLORATION MODE")
    print(f"{'='*60}")
    print(f"Job:     {job_name}")
    print(f"Starter: {starter}")
    print(f"Top N:   {top_n}")

    # --- 1. Build network -------------------------------------------------
    network = generate_network_tal(
        job_name=job_name,
        starters=starter,
        helpers=helpers,
        **network_kwargs,
    )

    # generate_network_tal takes file paths; find_pathways_to_target
    # takes raw SMILES. Resolve both formats here so per-endpoint
    # pathway runs get the right type.
    starter_smiles_list = get_smiles_from_file(starter)
    helper_smiles_list = get_smiles_from_file(helpers)
    starter_smi = starter_smiles_list[0] if starter_smiles_list else starter

    # --- 2. Score & rank --------------------------------------------------
    # Build exclusion set: starter SMILES + helper SMILES so they don't
    # rank as endpoints of themselves.
    excluded = set()
    for path in (starter, helpers):
        try:
            with open(path) as f:
                for ln in f:
                    s = ln.strip()
                    if s:
                        excluded.add(s)
        except (OSError, TypeError):
            pass

    print(f"\nScoring {len(network.mols)} molecules in network...")
    top_endpoints = rank_network_endpoints(
        network,
        exclude_smiles=excluded,
        top_n=top_n,
        **(scoring_kwargs or {}),
    )
    print(f"Top {len(top_endpoints)} endpoints by interestingness:")
    for i, es in enumerate(top_endpoints, 1):
        print(f"  {i:2d}. {es.smiles:40s} score={es.score:.3f}  "
              f"C={es.carbons} Bertz={es.bertz:.0f} FGs={es.n_functional_groups}")

    # --- 3. Find pathways to each top endpoint ----------------------------
    pathways_by_endpoint: dict[str, list[Pathway]] = {}
    generations = network_kwargs.get("gen", 2)
    for idx, es in enumerate(top_endpoints, 1):
        sub_job = f"{job_name}_ep{idx:02d}"
        print(f"\n--- pathways to endpoint {idx}/{len(top_endpoints)}: "
              f"{es.smiles} ---")
        try:
            find_pathways_to_target(
                network=network,
                starter=starter_smi,
                target=es.smiles,
                helpers=helper_smiles_list,
                generations=generations,
                max_num_rxns=max_num_rxns,
                job_name=sub_job,
            )
            pathways = load_pathways_from_file(sub_job)
        except Exception as exc:
            print(f"  [skip] {exc}")
            pathways = []
        pathways_by_endpoint[es.smiles] = pathways
        print(f"  → {len(pathways)} pathway(s) found")

    # --- 4. Cross-check against known-derivatives list --------------------
    matches: list[dict] | None = None
    missing: list[dict] | None = None
    if derivatives_list is not None:
        # Build canonical SMILES set of all molecules in the network.
        from rdkit import Chem
        network_canonical = set()
        for mol in network.mols:
            smi = getattr(mol, "smiles", None) or getattr(mol, "uid", None)
            if not smi:
                continue
            rd = Chem.MolFromSmiles(smi)
            if rd is not None:
                network_canonical.add(Chem.MolToSmiles(rd))

        matches, missing = [], []
        for entry in derivatives_list:
            rd = Chem.MolFromSmiles(entry["smiles"])
            if rd is None:
                missing.append(entry)
                continue
            canon = Chem.MolToSmiles(rd)
            (matches if canon in network_canonical else missing).append(entry)

        print(f"\nKnown derivatives in network: "
              f"{len(matches)}/{len(derivatives_list)}")

    result = {
        "top_endpoints": top_endpoints,
        "pathways": pathways_by_endpoint,
        "derivative_matches": matches,
        "derivative_missing": missing,
        "network": network,
        "starter_smiles": starter_smi,
        "helper_smiles": helper_smiles_list,
    }

    _write_exploration_report(result, job_name)
    return result


def _write_exploration_report(result: dict, job_name: str) -> None:
    """Write the consolidated exploration report as markdown."""
    out_path = f"{job_name}_exploration_report.md"
    lines: list[str] = []
    lines.append(f"# Exploration Report — {job_name}")
    lines.append("")
    lines.append(
        "Open exploration ranks every molecule reachable from the "
        "starter by an interestingness heuristic (carbon-count gate, "
        "Bertz complexity, functional-group diversity, aromatic bonus), "
        "then finds pathways to each top endpoint."
    )
    lines.append("")

    top = result["top_endpoints"]
    lines.append(f"## Top {len(top)} endpoints by interestingness")
    lines.append("")
    for i, es in enumerate(top, 1):
        lines.append(f"### {i}. `{es.smiles}` — score {es.score:.3f}")
        lines.append(
            f"- C={es.carbons}, Bertz={es.bertz:.1f}, "
            f"FGs={es.n_functional_groups}, "
            f"aromatic rings={es.n_aromatic_rings}"
        )
        lines.append(f"- breakdown: {es.breakdown}")
        pathways = result["pathways"].get(es.smiles, [])
        if not pathways:
            lines.append("- no pathways found (endpoint unreachable "
                         "within max_num_rxns)")
            lines.append("")
            continue
        lines.append(f"- pathways found: {len(pathways)} (top 3 shown)")
        for j, p in enumerate(pathways[:3], 1):
            lines.append(f"  - **pathway {j}** ({p.num_steps} steps)")
            for rxn in p.reactions:
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
                        f"    - `{parsed['op_name']}` "
                        f"(dH={dH_s}): {arrow}"
                    )
                except Exception:
                    lines.append(f"    - {rxn}")
        lines.append("")

    if result.get("derivative_matches") is not None:
        matches = result["derivative_matches"]
        missing = result["derivative_missing"]
        total = len(matches) + len(missing)
        lines.append("## Known-derivatives cross-check")
        lines.append("")
        lines.append(
            f"Of {total} literature TAL derivatives, "
            f"**{len(matches)}** appeared in the network, "
            f"{len(missing)} did not."
        )
        lines.append("")
        lines.append(f"### Found ({len(matches)})")
        for m in matches:
            lines.append(f"- **{m['name']}** — `{m['smiles']}`")
        lines.append("")
        lines.append(f"### Missing ({len(missing)})")
        for m in missing:
            lines.append(f"- **{m['name']}** — `{m['smiles']}`")
            if m.get("rationale"):
                lines.append(f"  - {m['rationale']}")
        lines.append("")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\nExploration report → {out_path}")