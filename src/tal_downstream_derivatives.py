"""
Known TAL (triacetic acid lactone) downstream derivatives from the
literature, used as a cross-check for open-exploration mode.

Open exploration ranks every endpoint in the network by interestingness
— this list answers a separate question: "Did we recover the products
the literature says TAL should be able to make?"

A miss is informative:
  - whitelist too narrow (chem or bio rule missing)
  - network didn't expand deep enough
  - synthesis route is non-trivial and the interestingness ranker
    pruned a critical intermediate

References:
  Chia, M., Schwartz, T. J., Shanks, B. H., & Dumesic, J. A. (2012).
    "Triacetic acid lactone as a potential biorenewable platform
    chemical." Green Chemistry 14, 1850-1853.
  Cardenas, J., & Da Silva, N. A. (2014). "Metabolic engineering of
    Saccharomyces cerevisiae for the overproduction of triacetic acid
    lactone." Metabolic Engineering 25, 194-203.
  Saunders, L. P., et al. (2015). "Triacetic acid lactone production in
    industrial Saccharomyces yeast strains." Journal of Industrial
    Microbiology & Biotechnology 42, 711-721.
  Markham, K. A., et al. (2018). "Rewiring Yarrowia lipolytica toward
    triacetic acid lactone for materials generation." PNAS 115, 2096-2101.
"""

TAL_DOWNSTREAM_DERIVATIVES = [
    {
        "name": "Sorbic acid",
        "smiles": "C/C=C/C=C/C(=O)O",
        "rationale": (
            "Major commercial food preservative (>30 kt/yr). TAL → "
            "sorbic acid via ring opening + dehydration is the "
            "headline commercial target for TAL platforms."
        ),
    },
    {
        "name": "Phloroglucinol (1,3,5-trihydroxybenzene)",
        "smiles": "Oc1cc(O)cc(O)c1",
        "rationale": (
            "Explosive precursor and dye intermediate. Accessed from "
            "TAL via aromatization of two condensed TAL units."
        ),
    },
    {
        "name": "Orcinol (5-methylresorcinol)",
        "smiles": "Cc1cc(O)cc(O)c1",
        "rationale": (
            "Dye and polymer monomer. Reported TAL derivative via "
            "Diels-Alder / aromatization pathway."
        ),
    },
    {
        "name": "Parasorbic acid (5,6-dihydro-6-methyl-2H-pyran-2-one)",
        "smiles": "CC1CC=CC(=O)O1",
        "rationale": (
            "Natural flavor compound. Direct partial hydrogenation "
            "of TAL."
        ),
    },
    {
        "name": "Triacetic acid (open-chain form)",
        "smiles": "CC(=O)CC(=O)CC(=O)O",
        "rationale": (
            "TAL is its lactone — opening the ring gives the linear "
            "triketoacid."
        ),
    },
    {
        "name": "2,4-dihydroxy-6-methylbenzoic acid (orsellinic acid)",
        "smiles": "Cc1cc(O)cc(O)c1C(=O)O",
        "rationale": (
            "Natural product intermediate (orsellinic acid family). "
            "Polyketide-derived aromatic analog of TAL — same C6+CO2H "
            "skeleton."
        ),
    },
    {
        "name": "Sorbyl alcohol (2,4-hexadien-1-ol)",
        "smiles": "C/C=C/C=C/CO",
        "rationale": (
            "Reduction product of sorbic acid. Fragrance/flavor "
            "intermediate; downstream of the headline TAL→sorbic route."
        ),
    },
    {
        "name": "δ-valerolactone",
        "smiles": "O=C1CCCCO1",
        "rationale": (
            "Bioplastic monomer precursor. Fully saturated "
            "6-membered lactone analog of TAL."
        ),
    },
    {
        "name": "4-hydroxy-2H-pyran-2-one",
        "smiles": "OC1=CC(=O)OC=C1",
        "rationale": (
            "De-methyl TAL — minimal pyranone scaffold; useful for "
            "studying whether the methyl-removal route is reachable."
        ),
    },
    {
        "name": "2,4-pentanedione (acetylacetone)",
        "smiles": "CC(=O)CC(=O)C",
        "rationale": (
            "Hydrolysis / retro-aldol product of TAL. Industrial "
            "solvent and ligand chemistry."
        ),
    },
    {
        "name": "Methyl sorbate",
        "smiles": "C/C=C/C=C/C(=O)OC",
        "rationale": (
            "Sorbic acid methyl ester. Common food-industry "
            "derivative."
        ),
    },
    {
        "name": "1,3,5-trihydroxymethylbenzene (5-methylphloroglucinol)",
        "smiles": "Cc1c(O)cc(O)cc1O",
        "rationale": (
            "Phloroglucinol with a methyl handle. Polymer/resin "
            "intermediate."
        ),
    },
    {
        "name": "4-hydroxy-6-methyltetrahydro-2H-pyran-2-one (DHMP)",
        "smiles": "CC1CC(O)CC(=O)O1",
        "rationale": (
            "Fully saturated TAL — hydrogenation product. Polymer "
            "precursor."
        ),
    },
]
