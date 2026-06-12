# src/pathway_tools.py

from doranet.modules.post_processing.post_processing import (
    pretreat_networks,
    pathway_finder,
)

from rdkit import Chem
import json


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