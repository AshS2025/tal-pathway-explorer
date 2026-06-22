"""
TAL Bio GLUCOSE -> TAL Pathway Whitelist (auto-generated, STRICT)
=================================================================

Filtered from DORAnet's JN1224MIN ruleset (1224 rules) by
scripts/build_bio_glucose_to_tal_whitelist.py.

A rule is kept iff applying it to one literature intermediate
produces the NEXT literature intermediate (canonical SMILES match).
This is the strict whitelist for validating the glucose -> TAL
biosynthetic pathway end-to-end.

Kept: 80 / 1224 rules.

Literature pathway:
#    1. glucose                        O=CC(O)C(O)C(O)C(O)CO
#    2. glucose-6-phosphate            O=CC(O)C(O)C(O)C(O)COP(=O)(O)O
#    3. fructose-6-phosphate           O=C(CO)C(O)C(O)C(O)COP(=O)(O)O
#    4. fructose-1,6-bisphosphate      O=C(COP(=O)(O)O)C(O)C(O)C(O)COP(=O)(O)O
#    5. dihydroxyacetone phosphate     O=C(CO)COP(=O)(O)O
#    6. glyceraldehyde-3-phosphate     O=CC(O)COP(=O)(O)O
#    7. 1,3-bisphosphoglycerate        O=C(OP(=O)(O)O)C(O)COP(=O)(O)O
#    8. 3-phosphoglycerate             O=C(O)C(O)COP(=O)(O)O
#    9. 2-phosphoglycerate             O=C(O)C(CO)OP(=O)(O)O
#   10. phosphoenolpyruvate            C=C(OP(=O)(O)O)C(=O)O
#   11. pyruvate                       CC(=O)C(=O)O
#   12. acetyl-CoA                     CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(...
#   13. malonyl-CoA                    CC(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@@H](n2cnc3c(N)ncnc...
#   14. acetoacetyl-CoA                CC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)O...
#   15. 3,5-dioxohexanoyl-CoA          CC(=O)CC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=...
#   16. TAL                            Cc1cc(O)cc(=O)o1

DO NOT edit by hand. Re-run scripts/build_bio_glucose_to_tal_whitelist.py
to regenerate.
"""

TAL_BIO_GLUCOSE_TO_TAL_WHITELIST = frozenset({
    "rule0001",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0006",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0007",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0014",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0015",  # glucose -> glucose-6-phosphate; fructose-6-phosphate -> fructose-1,6-bisphosphate
    "rule0016",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0017",  # glucose -> glucose-6-phosphate; fructose-6-phosphate -> fructose-1,6-bisphosphate
    "rule0019",  # 2-phosphoglycerate -> phosphoenolpyruvate
    "rule0023",  # acetyl-CoA -> malonyl-CoA
    "rule0024",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0028",  # glucose-6-phosphate -> fructose-6-phosphate; dihydroxyacetone phosphate -> glyceraldehyde-3-phosphate; 3-phosphoglycerate -> 2-phosphoglycerate
    "rule0029",  # fructose-6-phosphate -> fructose-1,6-bisphosphate; 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0036",  # fructose-6-phosphate -> fructose-1,6-bisphosphate; 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0037",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0046",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0053",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0081",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0085",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0086",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate; pyruvate -> acetyl-CoA
    "rule0087",  # acetoacetyl-CoA -> 3,5-dioxohexanoyl-CoA
    "rule0095",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0096",  # glucose -> glucose-6-phosphate; fructose-6-phosphate -> fructose-1,6-bisphosphate
    "rule0097",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0104",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0115",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate; pyruvate -> acetyl-CoA
    "rule0124",  # fructose-6-phosphate -> fructose-1,6-bisphosphate; 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0125",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0126",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate; acetoacetyl-CoA -> 3,5-dioxohexanoyl-CoA
    "rule0129",  # 3-phosphoglycerate -> 2-phosphoglycerate
    "rule0150",  # glyceraldehyde-3-phosphate -> 1,3-bisphosphoglycerate
    "rule0182",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0215",  # glucose -> glucose-6-phosphate; fructose-6-phosphate -> fructose-1,6-bisphosphate
    "rule0216",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0248",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0250",  # glucose-6-phosphate -> fructose-6-phosphate; dihydroxyacetone phosphate -> glyceraldehyde-3-phosphate
    "rule0256",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0263",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0267",  # glyceraldehyde-3-phosphate -> 1,3-bisphosphoglycerate
    "rule0269",  # pyruvate -> acetyl-CoA
    "rule0274",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0275",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate; pyruvate -> acetyl-CoA
    "rule0306",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0349",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate; pyruvate -> acetyl-CoA
    "rule0350",  # acetoacetyl-CoA -> 3,5-dioxohexanoyl-CoA
    "rule0384",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0394",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0402",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0408",  # phosphoenolpyruvate -> pyruvate
    "rule0431",  # glucose -> glucose-6-phosphate; fructose-6-phosphate -> fructose-1,6-bisphosphate; 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0454",  # 2-phosphoglycerate -> phosphoenolpyruvate
    "rule0461",  # 2-phosphoglycerate -> phosphoenolpyruvate
    "rule0472",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0473",  # glucose -> glucose-6-phosphate; fructose-6-phosphate -> fructose-1,6-bisphosphate
    "rule0503",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate; 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0517",  # glucose -> glucose-6-phosphate; fructose-6-phosphate -> fructose-1,6-bisphosphate; 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0549",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0569",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0584",  # phosphoenolpyruvate -> pyruvate
    "rule0604",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate; 3-phosphoglycerate -> 2-phosphoglycerate
    "rule0605",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate; 3-phosphoglycerate -> 2-phosphoglycerate
    "rule0730",  # acetyl-CoA -> malonyl-CoA
    "rule0731",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule0767",  # fructose-6-phosphate -> fructose-1,6-bisphosphate; 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0768",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0846",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0847",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule0891",  # 3,5-dioxohexanoyl-CoA -> TAL
    "rule1027",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule1028",  # phosphoenolpyruvate -> pyruvate
    "rule1062",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule1086",  # phosphoenolpyruvate -> pyruvate
    "rule1090",  # phosphoenolpyruvate -> pyruvate
    "rule1117",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate; pyruvate -> acetyl-CoA
    "rule1118",  # phosphoenolpyruvate -> pyruvate; malonyl-CoA -> acetoacetyl-CoA
    "rule1134",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule1135",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule1149",  # fructose-1,6-bisphosphate -> dihydroxyacetone phosphate
    "rule1215",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule1221",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
    "rule1222",  # 1,3-bisphosphoglycerate -> 3-phosphoglycerate
})
