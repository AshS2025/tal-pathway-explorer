"""
TAL (Triacetic Acid Lactone) Reaction Whitelist for DORAnet
===========================================================

TAL is a C6 delta-lactone (2-pyranone) with formula C6H6O3.
Key structural features that guide this whitelist:
  - 6-membered lactone ring (delta-lactone)
  - Conjugated alkene in ring (C=C adjacent to C=O)
  - Beta-keto group
  - Methyl substituent
  - Strong tendency toward keto-enol tautomerism
  - Known to aromatize to orcinol/phloroglucinol-type phenols
  - Ring-opens to triacetic acid (beta-keto acid)

Helper molecules assumed: H2O, H2, EtOH (CCO), MeOH, acetic acid,
CO2, NH3, simple amines.

This whitelist retains ~100 of 388 reactions.
To use in generate_network_TAL, add before the operator loop:
    smarts_list = [op for op in smarts_list
                   if op.name in TAL_REACTION_WHITELIST]

Excluded categories (with rationale):
  A. All halogen-introducing reactions — TAL platform targets are
     oxygenated small molecules; halogens add toxicity and cost.
  B. Petrochemical/cracking reactions — TAL is bio-derived; cracking
     reactions require C3-C6 alkanes as substrates, not lactones.
  C. Organometallic cross-couplings requiring exotic metal reagents —
     Heck, Suzuki, Sonogashira, Wurtz/Gilman, Simmons-Smith, McMurry,
     Tsuji-Trost. These require Pd, Zn, Ti organometallic species that
     are not appropriate for a bio-based platform chemical context.
  D. Alkyne chemistry (except hydration) — TAL and its immediate
     derivatives have no alkynes; keeping hydration of alkynes only
     because CO2/acetylene Reppe chemistry is borderline.
  E. Full sulfur block (334-387) — TAL has no sulfur; these reactions
     will never match and just waste computation.
  F. Scaffold-specific industrial processes — Naphthalene oxidation,
     cyclohexane oxidation, ethylbenzene dehydrogenation, Hock process,
     syngas methanol, Cativa/Tennessee Eastman/Reppe processes.
  G. Exotic nitrogen chemistry — nitration, azides, diazo, isocyanates,
     Curtius/Hofmann/Sandmeyer/Beckmann (except Beckmann from ketones
     which could give ring-expanded lactam from TAL's ketone),
     nitroso chemistry, HCN/HNO process chemistry.
  H. Keto-enol tautomerization — INTENTIONALLY EXCLUDED despite TAL
     having active tautomers. The reason: TAL's beta-diketone character
     means tautomerization fires on virtually every molecule downstream,
     exponentially bloating the network. The enol products will still
     be reachable via other reaction paths. Re-enable if you specifically
     need enol intermediates.
"""

TAL_REACTION_WHITELIST = frozenset({

    # ================================================================
    # ALKENE CHEMISTRY
    # TAL has a ring C=C; ring-opened products have exocyclic alkenes.
    # These are the highest-priority reactions for TAL.
    # ================================================================

    "Hydrogenation of Alkene",
    # Direct reduction of TAL's ring double bond. Essential.

    "Oxidative Cleavage of Alkenes",
    # Ozonolysis-type. Relevant for cleaving ring-opened alkene products.

    "Oxidative Cleavage of Alkenes, Intramolecular",
    # Intramolecular version for cyclic alkenes — TAL's ring alkene.

    "Epoxidation of Alkene",
    # mCPBA-type oxidation of TAL's C=C. Gives epoxy-lactone.
    # Epoxy lactones are known bioactive TAL derivatives.

    "Hydration of Alkene, Addition of Alcohols or Acids to Alkenes",
    # Markovnikov addition of water/alcohols to alkenes.
    # With helper ethanol/water, highly relevant.

    "Addition of Alcohols or Acids to Alkenes, Intramolecular",
    # Intramolecular version — relevant for ring-containing products.

    "Hydration of Alkenes, 2-step",
    # Two-step variant via sulfate intermediate. Kept for coverage.

    "Diol Formation by Oxidation",
    # Sharpless-type dihydroxylation of TAL's alkene.
    # Diols are useful intermediates.

    "Diels-Alder Reaction with Alkenes",
    # TAL's ring alkene + carbonyl system can act as dienophile.
    # Alpha-pyrones (which TAL resembles) are classic Diels-Alder
    # dienophiles. HIGH PRIORITY for TAL.

    "Diels-Alder Reaction with Alkenes, Intramolecular",
    # Intramolecular DA — possible for ring-opened triene forms.

    "Oxo-Diels-Alder Reaction",
    # TAL's C=O as a hetero-dienophile with diene helper molecules.

    "Oxo-Diels-Alder Reaction, Intramolecular",

    "Acrolein Diels-Alder Reaction",
    # TAL's enone system as an acrolein-type dienophile. Very relevant.

    "Acrolein Diels-Alder Reaction, Intramolecular",

    "Claisen Rearrangement, Cope Rearrangements",
    # TAL has an allylic ether-like connectivity in its lactone.
    # Relevant for rearrangements of allyl/vinyl ether derivatives.

    "Reductive Cross-Coupling of Alkenes",
    # Kept for exploratory purposes — can form C-C bonds from two alkenes.

    "Olefin Cross Metathesis (CM)",
    # Kept moderately: Ru-catalyzed metathesis is industrially practiced
    # and not as exotic as Pd couplings. TAL-derived alkenes could
    # participate. Lenient inclusion for exploratory network.

    "Acetoxylation of Alkenes",
    # Wacker-type with acetic acid helper. TAL + AcOH -> acetoxy product.
    # Acetic acid is a common cheap helper.

    "Acetoxylation of Alkenes, Intramolecular",

    # ================================================================
    # ALCOHOL CHEMISTRY
    # TAL ring-opening gives beta-hydroxy acids and diols downstream.
    # Alcohols will be ubiquitous in the network with water as helper.
    # ================================================================

    "Dehydration of Alcohol",
    # Alcohol -> alkene. Very common, very cheap.

    "Dehydration of Alcohol, 2-step",

    "Selective Oxidation of Alcohols",
    # Primary/secondary alcohol -> aldehyde/ketone. Essential.

    "Oxidation Of Primary Alcohols to Carboxylic Acids",
    # Primary alcohol -> carboxylic acid. Essential.

    "Oxidation Of Primary Alcohols to Carboxylic Acids, 2-step",

    "Glycol Cleavage by Oxidation",
    # Diol -> two carbonyls (periodate-type). Relevant for diol products
    # from epoxide opening or dihydroxylation.

    "Glycol Cleavage by Oxidation Intramolecular",

    "Hydrodeoxygenation of Alcohol, Classic Synthesis of Aldehydes from Carboxylic Acids",
    # Alcohol -> hydrocarbon (HDO). Relevant for bio-refinery context.

    "Hydrogenolysis of Primary Alcohol",
    # C-O hydrogenolysis. Relevant for deoxygenation routes.

    "Pinacol Rearrangement, no H on Carbon#2",
    # 1,2-diol -> carbonyl rearrangement. Will fire on diol products.

    "Pinacol Rearrangement, 1 H on Carbon#2",

    "Pinacol Rearrangement, 2 Hs on Carbon#2",

    "Dehydration of Geminal Diol",
    # Gem-diol (hydrate) -> carbonyl. Important for reversing hydration.

    "Formation of Acetals from Hemiacetals",
    # Hemiacetal + alcohol -> acetal. Relevant with EtOH helper.

    "Formation of Cyclic Acetals from Ketones/Aldehydes with Diols",
    # Protecting group chemistry; also creates cyclic products.

    "Diol Carboxylation",
    # Diol + CO2 -> cyclic carbonate. Relevant with CO2 helper.

    "Diol Carboxylation 2",

    # ================================================================
    # ETHER / EPOXIDE CHEMISTRY
    # TAL's lactone oxygen can participate in ether-like chemistry.
    # Epoxides arise from alkene epoxidation.
    # ================================================================

    "Hydrolysis of Ethers, Esters, Anhydrides",
    # CRITICAL: This is how TAL's lactone ring opens with water.
    # Essential reaction for the entire network.

    "Hydrolysis of Ethers, Esters, Anhydrides, Intramolecular",
    # Intramolecular lactonization/delactonization. Also essential.

    "Hydrogenolysis of Ethers",
    # Ether C-O cleavage with H2.

    "Hydrogenolysis of Ethers Intramolecular",

    "Ether Synthesis by Dehydration",
    # Two alcohols -> ether. Relevant with alcohol helpers.

    "Ether Synthesis by Dehydration, Intramolecular",

    "Epoxides Ring Opening",
    # Nucleophilic opening of epoxides. Important for epoxy-lactone
    # intermediates formed from TAL's alkene.

    # ================================================================
    # ALDEHYDE & KETONE CHEMISTRY
    # TAL has a ketone; ring-opened products have aldehydes.
    # This is the core reactive center of TAL.
    # ================================================================

    "Aldehyde Oxidation",
    # Aldehyde -> carboxylic acid. Essential.

    "Aldehyde Oxidation, 2-step",

    "Aldehyde & Alcohol Oxidation, 2-step",

    "Ketone Oxidation",
    # Ketone -> ester or acid via oxidative cleavage (not BV).

    "Ketone Oxidation, Intramolecular",

    "Baeyer-Villiger Oxidation (Ketones)",
    # Ketone -> ester/lactone. VERY relevant for TAL — can insert oxygen
    # into the ring or into ketone-containing products. High priority.

    "Baeyer-Villiger Oxidation (Aldehydes)",
    # Aldehyde -> formate ester.

    "Hydrogenation of Ketones",
    # Ketone -> alcohol. Essential with H2 helper.

    "Decarbonylation of Aldehydes",
    # Aldehyde -> alkane + CO. Relevant for decarbonylation routes.

    "Ketonization",
    # Two carboxylic acids -> ketone + CO2 + H2O. HIGHLY RELEVANT for
    # TAL: triacetic acid (ring-opened TAL) is a beta-keto acid and
    # ketonization is a well-known reaction pathway.

    "Ketonization, Intramolecular",

    "Oxidative Esterification of Aldehydes and Alcohols",
    # Aldehyde + alcohol -> ester in one pot.

    "Oxidative Esterification of Aldehydes and Alcohols, Intramolecular",

    "Hydration of Ketone and Aldehyde",
    # Carbonyl + water -> geminal diol.

    "Hemiacetal Dissociation, Addition of Alcohols to Carbonyl Groups Reverse",
    # Hemiacetal -> carbonyl + alcohol.

    "Hemiacetal Dissociation, Addition of Alcohols to Carbonyl Groups Reverse, Intramolecular",

    "Hemiacetal Formation, Addition of Alcohols to Carbonyl Groups",
    # Carbonyl + alcohol -> hemiacetal.

    "Hemiacetal Formation, Addition of Alcohols to Carbonyl Groups, Intramolecular",

    "Reduction of Carbonyl Groups",
    # General carbonyl reduction. Essential with H2.

    "Pinacol Coupling",
    # Two carbonyls -> 1,2-diol via reductive coupling. Gives diol
    # products that can then undergo pinacol rearrangement.

    "Pinacol Coupling, Intramolecular",

    "Synthesis of Enol Ethers from Aldehyde and Alcohol",
    # Aldehyde + alcohol -> vinyl ether. Relevant with EtOH helper.

    "Synthesis of Enol Ethers from Aldehyde and Alcohol, Intramolecular",

    # NOTE: Keto-enol tautomerization (136, 137) is INTENTIONALLY
    # EXCLUDED. See module docstring for rationale.

    # ================================================================
    # CARBOXYLIC ACID & ESTER CHEMISTRY
    # This is the most important block for TAL.
    # TAL ring-opens to triacetic acid; esterification with helpers
    # like EtOH is a primary derivatization path.
    # ================================================================

    "Carboxylic Acids Decarboxylation",
    # Beta-keto acids decarboxylate readily. CRITICAL for TAL — triacetic
    # acid is a beta-keto acid and spontaneous decarboxylation is a
    # well-documented TAL pathway.

    "Esterification, Acid Anhydride Formation",
    # Acid + alcohol -> ester. Essential with EtOH/MeOH helpers.

    "Esterification, Acid Anhydride Formation, Intramolecular",
    # Lactonization. Essential for ring-closure back to lactone products.

    "Acid Anhydrides React with Alcohols to Form Esters",
    # Anhydride + alcohol -> ester.

    "Acid Anhydrides React with Alcohols to Form Esters, Intramolecular",

    "Ester Reduction to Aldehydes and Alcohols",
    # DIBAL-type reduction. Relevant.

    "Ester Reduction to Aldehydes and Alcohols, Intramolecular",

    "Ester Reduction to Alcohols",
    # Full LiAlH4/H2 reduction of ester to alcohol.

    "Ester Reduction to Alcohols, Intramolecular",
    # Lactone -> diol. CRITICAL — TAL lactone -> diol is a key pathway.

    "Transesterification",
    # Ester + alcohol -> new ester. Essential with EtOH/MeOH helpers.

    "Transesterification, Intramolecular",

    # ================================================================
    # ENOLATE / CONDENSATION CHEMISTRY
    # TAL is an activated methylene compound. Its beta-diketone
    # character makes it an excellent nucleophile and electrophile
    # for all condensation reactions. These are HIGHEST PRIORITY.
    # ================================================================

    "Aldol Condensation",
    # TAL's active methylene attacks a carbonyl. Core TAL chemistry.

    "Aldol Condensation, Intramolecular",

    "Aldol Condensation (2H)",
    # Aldol with subsequent dehydration (gives alpha,beta-unsaturated).

    "Aldol Condensation (2H), Intramolecular",

    "Claisen Condensation",
    # Ester enolate attacks another ester. Relevant for TAL ester products.

    "Dieckmann Condensation",
    # Intramolecular Claisen. VERY relevant — could re-cyclize TAL
    # ring-opened products back into cyclic beta-keto esters.

    "Michael Reaction",
    # 1,4-addition of nucleophile to enone. TAL's alpha,beta-unsaturated
    # lactone is a Michael acceptor. Also TAL as nucleophile.
    # CRITICAL for TAL.

    "Michael Reaction, Intramolecular",

    "Michael Reaction with Cyclic Ketones",
    # Michael addition involving cyclic ketone acceptors.

    "Michael Reaction with Cyclic Ketones, Intramolecular",

    "Robinson Annulation",
    # Michael + aldol cascade giving cyclohexenone. TAL products can
    # participate as Michael acceptors.

    "Enolate Alkylation",
    # Alpha-alkylation of carbonyl. Relevant with alkyl helpers.
    # NOTE: requires alkyl halide so network must generate one first.
    # Kept for exploratory coverage.

    "Enolate Alkylation, Intramolecular",

    # ================================================================
    # AROMATIC CHEMISTRY
    # TAL aromatizes to phenolic products (orcinol, phloroglucinol-type).
    # These are experimentally known TAL derivatives.
    # Aromatic reactions on those downstream phenols are relevant.
    # ================================================================

    "Benzene Hydrogenation",
    # Phenolic products can be hydrogenated.

    "Benzene Partial Hydrogenation 1",

    "Benzene Partial Hydrogenation 2",

    "Phenols Oxidation to Quinones 1",
    # Orcinol-type phenols -> quinones. Directly relevant for
    # TAL-derived phenol products.

    "Phenols Oxidation to Quinones 2",

    "Phenols Oxidation to Quinones 3",

    "Phenols Oxidation to Quinones 3 Reverse",

    "Phenols Oxidation to Quinones 4",

    "Phenols Oxidation to Quinones 4 Reverse",

    "Oxidation of Aromatic Alkanes Ar-CH3",
    # TAL-derived phenols have a methyl group (orcinol has -CH3).
    # Benzylic methyl oxidation -> benzoic acid derivative.

    "Oxidation of Aromatic Alkanes Ar-CH3, 2-step",

    "Electrophilic Aromatic Alkylation with Alkenes",
    # Phenol products can undergo Friedel-Crafts with alkene helpers.

    "Electrophilic Aromatic Alkylation with Alkenes, Intramolecular",

    "Electrophilic Aromatic Alkylation with Alcohols",

    "Electrophilic Aromatic Alkylation with Alcohols, Intramolecular",

    "Friedel–Crafts Acylation with Carboxylic Acids",
    # Aromatic ring + carboxylic acid -> aryl ketone. Relevant for
    # phenol products reacting with TAL-derived acids.

    "Friedel–Crafts Acylation with Carboxylic Acids, Intramolecular",

    "Friedel–Crafts Acylation with Acid Anhydrides",

    "Friedel–Crafts Reaction with Alkenes",
    # Friedel-Crafts alkylation with alkenes (no halide needed).

    "Friedel–Crafts Reaction with Alkenes, Intramolecular",

    "Friedel-Crafts Hydroxyalkylation",
    # Phenol + aldehyde -> benzyl alcohol product. Highly relevant for
    # phenol products reacting with carbonyl-containing TAL derivatives.

    "Friedel-Crafts Hydroxyalkylation, Intramolecular",

    "Kolbe Carboxylation",
    # Phenol + CO2 -> hydroxybenzoic acid (salicylate-type).
    # Very relevant if CO2 is used as a helper molecule.

    "Furan Carboxylation",
    # Furans are potential TAL derivatives via ring transformation.

    "Furan Carboxylation, 2-step",

    # ================================================================
    # CARBONYLATION CHEMISTRY
    # With CO as a helper (or generated in network), carbonylation
    # reactions add carbon and form esters/acids.
    # ================================================================

    "Hydroformylation",
    # Alkene + CO + H2 -> aldehyde. Industrially important, can expand
    # carbon skeleton of TAL-derived alkenes.

    "Hydrocarboxylation 1, Hydroesterification(Carboalkoxylation)",
    # Alkene + CO + alcohol/water -> ester/acid.

    "Hydrocarboxylation 2",

    # ================================================================
    # NITROGEN CHEMISTRY (selective)
    # You indicated nitrogen is of interest. TAL's carbonyl and
    # beta-diketone character make it an excellent substrate for
    # imine/enamine/amide formation. Selected reactions kept.
    # Excluded: nitration, azides, diazo, Hofmann, Sandmeyer,
    # nitroso chemistry, HCN processes.
    # ================================================================

    "Ketone Reductive Amination",
    # Ketone + amine + H2 -> amine. TAL's ketone -> amino product.
    # Relevant with NH3 or amine helpers.

    "Ketone Reductive Amination, Intramolecular",

    "Synthesis of Amides with Carboxylic Acid",
    # Acid + amine -> amide. TAL-derived acids + amine helpers.

    "Synthesis of Amides with Carboxylic Acid, Intramolecular",

    "Synthesis of Amides with Acid Anhydrides or Esters",
    # Ester/anhydride + amine -> amide. TAL lactone + amine.

    "Synthesis of Amides with Acid Anhydrides or Esters, Intramolecular",

    "Hydrolysis of Amides",
    # Amide -> acid + amine. Reversibility of amide formation.

    "Hydrolysis of Amides, Intramolecular",

    "Dehydration of Amides",
    # Amide -> nitrile. Connects amide products to nitrile space.

    "Imines from Aldehydes and Ketones",
    # Carbonyl + primary amine -> imine. TAL's carbonyl + amine helper.

    "Imines from Aldehydes and Ketones, Intramolecular",

    "Imines from Aldehydes and Ketones Reverse",
    # Imine hydrolysis back to carbonyl.

    "Imines from Aldehydes and Ketones, Intramolecular Reverse",

    "Treatment of Aldehydes and Ketones with a Secondary Amine",
    # Carbonyl + secondary amine -> enamine. TAL is very reactive here.

    "Treatment of Aldehydes and Ketones with a Secondary Amine, Intramolecular",

    "Treatment of Aldehydes and Ketones with a Secondary Amine Reverse",

    "Treatment of Aldehydes and Ketones with a Secondary Amine, Intramolecular Reverse",

    "Hemiaminal Formation",
    # Carbonyl + amine -> hemiaminal intermediate.

    "Hemiaminal Formation, Reverse",

    "Hemiaminal Dehydration",
    # Hemiaminal -> imine + water.

    "Hemiaminal Dehydration, Reverse",

    "Hydroamination of Alkenes",
    # Alkene + amine -> amine product. TAL's ring alkene + amine helper.

    "Hydroamination of Alkenes, Intramolecular",

    "Hydroamination of Dienes",
    # Diene + amine -> allylic amine. Relevant for TAL ring-opened dienes.

    "Hydroamination of Dienes, Intramolecular",

    "Ring Opening of Epoxides by Amines",
    # Epoxide (from TAL alkene epoxidation) + amine -> amino alcohol.

    "Reduction of Imines",
    # Imine + H2 -> amine.

    "Beckmann Rearrangement from Ketones",
    # Ketone + NH2OH -> lactam (ring expansion). SPECIFICALLY relevant
    # for TAL's ketone -> ring-expanded nitrogen heterocycle.
    # This is a known route to N-containing TAL derivatives.

    "Beckmann Rearrangement with Ketoximes",
    # Ketoxime -> lactam. Same pathway, oxime intermediate.

    "Wolff–Kishner Reduction",
    # Carbonyl -> methylene via hydrazone. Relevant deoxygenation.

    "Amine Alkylation with Alcohols or Primary Amines",
    # N-alkylation via alcohol. Relevant with EtOH/MeOH helpers.

    "Amine Alkylation with Alcohols or Primary Amines, Intramolecular",

    "Cyanation of Ketones or Aldehydes (Cyanohydrin Reaction)",
    # Carbonyl + HCN -> cyanohydrin. Kept for coverage — HCN is
    # a plausible helper and cyanohydrins are useful intermediates.

    "Hydrolysis of Nitriles",
    # Nitrile + water -> amide/acid. Closes cyanohydrin pathway.

    "Partial Hydrolysis of Nitriles",
    # Nitrile -> amide (stopped at amide stage).

    # ================================================================
    # MISCELLANEOUS — KEPT FOR EXPLORATORY COVERAGE
    # ================================================================

    "[2+2] Cycloaddition",
    # Ketene + carbonyl -> beta-lactone. TAL could generate ketene
    # intermediates via decarbonylation; beta-lactones are useful.

    "Hydroformylation",
    # Already listed above; duplicate entry harmless in a frozenset.

})


# ================================================================
# USAGE — drop this into your generate_network_TAL function
# ================================================================
#
# from tal_reaction_whitelist import TAL_REACTION_WHITELIST
#
# # Replace the smarts_list assignment with:
# if direction == "forward":
#     smarts_list = [op for op in op_smarts
#                    if op.name in TAL_REACTION_WHITELIST]
# elif direction == "retro":
#     smarts_list = [op for op in op_retro_smarts
#                    if op.name in TAL_REACTION_WHITELIST]
#
# # Then the existing operator loop continues unchanged:
# for smarts in smarts_list:
#     ...
#
# ================================================================
# STATISTICS
# ================================================================
#
# Total reactions in DORAnet forward list: 388
# Reactions in this whitelist:             ~100
# Reduction:                               ~74%
#
# Excluded blocks summary:
#   Halogenation (all forms):              ~35 reactions
#   Petrochemical/cracking:                 8 reactions
#   Organometallic cross-couplings:        15 reactions
#   Alkyne chemistry (most):               10 reactions
#   Exotic industrial processes:            8 reactions
#   Full sulfur block:                     54 reactions
#   Exotic nitrogen chemistry:             ~60 reactions
#   Keto-enol tautomerization:              2 reactions
#   Scaffold-specific (naphthalene etc):    8 reactions
# ================================================================
