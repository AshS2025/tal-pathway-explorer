# src/tal_pathway_scoring.py

from doranet.modules.post_processing.post_processing import pathway_ranking


def generate_base_rankings(
    starter,
    target,
    helpers,
    job_name="TAL",
):
    """
    Run DORAnet's built-in pathway ranking.

    This creates:
        TAL_ranked_pathways.txt

    containing all pathways ranked by DORAnet's default metrics:
        - reaction thermodynamics
        - number of steps
        - by-product formation
        - atom economy
    """

    weights = {
        "reaction_thermo": 2,
        "number_of_steps": 4,
        "by_product_number": 2,
        "atom_economy": 1,
        "salt_score": 0,
        "in_reaxys": 0,
        "coolness": 0,
    }

    pathway_ranking(
        starters=[starter],
        helpers=helpers,
        target=[target],
        weights=weights,
        job_name=job_name,
    )

    print(f"Generated {job_name}_ranked_pathways.txt")


    # PARSE THROUGH THE OUTPUT AND EXTRACT DATA 

    def parse_ranked_pathways(job_name):
   
    # """
    # Read TAL_ranked_pathways.txt

    # Convert it into a Python structure like:

    # [
    #     {
    #         "rank": 1,
    #         "score": 7.4,
    #         "atom_economy": 0.81,
    #         "pathway_byproduct": 4,
    #         "reactions": [...]
    #     },
    #     ...
    # ]
    # """

    # TODO
        pass

# ON TOP OF THAT FRAMEWORK IS WHERE WE'D ADD ADDITIONAL ANALYSIS FILTERS AND THEN 
# RERANK THE PATHWAYS FOR OUR PURPOSES

# additional analysis metrics I'm going to add 
# cost - for retro it's feedstock price, for forward it's helper cost 
# bio-chem switches - operational complexity - (grounded in 2025 paper)
# intermediate analysis - toxicity (PubChem GHS flags or in silico Tox21) and 
# stability (InReaxys or DORA-XGB?)


# VERY IMPORTANT THING TO CONSIDER AS I MOVE AHEAD WITH ANALYSIS METRICS 
 # Weights should not be chosen to produce a single universal ranking. Instead, 
 # the framework should support multiple decision scenarios corresponding to different
 # process-development objectives - for example feasibility, sustainability, economics
# etc. -- this should be chosen by the user. -- actually this is going to get messy quick
# i'm going to first work on developing a 'lemnisca' score. and then if we want to add
# modifications to weights later based on company goals then we can do that later. 
