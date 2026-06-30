"""
TAL Bio Polyketide-Only Whitelist
=================================

The minimum set of JN1224MIN enzymatic rules needed to describe TAL
biosynthesis from acetyl-CoA via the 2-pyrone synthase (2-PS)
polyketide pathway. Extracted from
tal_bio_glucose_to_tal_whitelist.py — only the four polyketide steps
are kept; glycolysis upstream of acetyl-CoA is removed.

Pathway covered
---------------
  acetyl-CoA  +  CO2  +  ATP   --(acetyl-CoA carboxylase)-->
      malonyl-CoA  +  ADP  +  Pi

  acetyl-CoA  +  malonyl-CoA   --(Claisen condensation 1)-->
      acetoacetyl-CoA  +  CO2  +  CoA

  acetoacetyl-CoA  +  malonyl-CoA   --(Claisen condensation 2)-->
      3,5-dioxohexanoyl-CoA  +  CO2  +  CoA

  3,5-dioxohexanoyl-CoA   --(intramolecular cyclization, 2-PS)-->
      TAL  +  CoA

Use when the search target is acetyl-CoA <- TAL (the CEO-scoped
retro). The rule count (7) is small enough that even multi-substrate
chemistry stays tractable once the cofactor pool is also restricted
to polyketide-relevant cofactors.

DO NOT edit by hand. Re-extract from
tal_bio_glucose_to_tal_whitelist.py if regeneration is needed.
"""

TAL_BIO_POLYKETIDE_WHITELIST = frozenset({
    "rule0023",   # acetyl-CoA -> malonyl-CoA  (carboxylation)
    "rule0730",   # acetyl-CoA -> malonyl-CoA  (variant)
    "rule1118",   # malonyl-CoA -> acetoacetyl-CoA  (Claisen condensation 1)
    "rule0087",   # acetoacetyl-CoA -> 3,5-dioxohexanoyl-CoA  (Claisen 2)
    "rule0126",   # acetoacetyl-CoA -> 3,5-dioxohexanoyl-CoA  (variant)
    "rule0350",   # acetoacetyl-CoA -> 3,5-dioxohexanoyl-CoA  (variant)
    "rule0891",   # 3,5-dioxohexanoyl-CoA -> TAL  (intramolecular cyclization)
})
