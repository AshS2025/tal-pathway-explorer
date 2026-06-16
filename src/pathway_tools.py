# src/pathway_tools.py

from dataclasses import dataclass
from typing import Optional

from doranet.modules.post_processing.post_processing import (
    pretreat_networks,
    pathway_finder,
)

from rdkit import Chem
import json


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