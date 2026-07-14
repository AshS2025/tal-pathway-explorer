"""
equilibrator_client.py
======================

Thin wrapper around equilibrator_api's ComponentContribution that lets
us score reaction free energy from SMILES-based reaction strings.

WHY WE NEED A WRAPPER
---------------------
equilibrator_api works on Compound objects that are looked up from a
KEGG-derived database. Our pathway files use SMILES strings that
often don't cleanly map to KEGG IDs (especially CoA-tethered
polyketide intermediates that were never assigned KEGG entries).

The wrapper:
  1. Maintains a small SMILES → KEGG ID map for the common cofactors
     we know are in equilibrator's DB
  2. For SMILES not in the map, falls back to InChI-key lookup
  3. Caches every reaction score by its normalized SMILES string
  4. Returns None gracefully for reactions where any compound is
     unmappable, rather than crashing the pipeline

HOW TO USE
----------
    from equilibrator_client import EquilibratorClient
    client = EquilibratorClient()               # ~20s one-time init
    dg = client.dG_prime("CC(=O)C(=O)O>>CC(O)C(=O)O")
    # dg is a float (kJ/mol) or None if unmappable
"""

from __future__ import annotations

import threading
from typing import Optional

from rdkit import Chem


# --- Hardcoded SMILES → KEGG ID map for compounds we KNOW are in the
# equilibrator database. The RDKit-canonicalized SMILES are the keys;
# KEGG ids are the values. All the common cofactors and small
# molecules our chem + bio pipelines emit.
#
# Adding to this map is the main way to make equilibrator work with
# more of our chemistry. When we hit an unmappable intermediate, look
# it up on KEGG and add an entry here.
_SMILES_TO_KEGG_RAW = {
    # water and inorganics
    "O": "C00001",              # water
    "O=C=O": "C00011",          # CO2
    "[H][H]": "C00282",         # H2
    "O=O": "C00007",            # O2
    "N": "C00014",              # ammonia
    "OP(=O)(O)O": "C00009",     # phosphate

    # currency cofactors
    "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O":
        "C00024",               # acetyl-CoA
    "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O":
        "C00083",               # malonyl-CoA
    "CC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O":
        "C00332",               # acetoacetyl-CoA
    "OC1O[C@@H](COP(=O)(O)OP(=O)(O)OC[C@@H]2OC([n+]3ccc(C(N)=O)cc3)"
    "[C@H](O)[C@H]2O)[C@@H](O)[C@H]1O":
        "C00003",               # NAD+  (approximate SMILES)
}


def _canon(smi: str) -> Optional[str]:
    """Canonicalise a SMILES string. Returns None for invalid input."""
    if not smi:
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


# Canonicalise the map keys once at import time so lookups are cheap.
_SMILES_TO_KEGG = {}
for _raw, _kegg in _SMILES_TO_KEGG_RAW.items():
    _c = _canon(_raw)
    if _c is not None:
        _SMILES_TO_KEGG[_c] = _kegg


class EquilibratorClient:
    """
    Wrapper around equilibrator_api's ComponentContribution.

    Init cost: ~20 seconds (loads the equilibrator SQLite DB into
    memory). One instance can be reused across the whole app session.

    Thread-safe on `.dG_prime` calls via an internal lock, so it can
    be stashed in Streamlit's @st.cache_resource without extra care.
    """

    def __init__(self):
        from equilibrator_api import ComponentContribution
        self._cc = ComponentContribution()
        self._compound_cache = {}      # canon smiles -> Compound or None
        self._reaction_cache = {}      # normalized rxn string -> float or None
        self._lock = threading.Lock()

    # ------------------------------------------------------------------
    # Compound lookup
    # ------------------------------------------------------------------
    # We store BOTH the Compound object and a string "accession" that
    # `parse_reaction_formula` can read. Formats it accepts:
    #     "kegg:C00001"          — namespace prefix + accession
    #     "coco:{compound_id}"   — the internal compound-cache namespace
    # We use the latter as a fallback when only the InChI-key matched.
    def _find_compound(self, smi: str):
        """SMILES → (Compound, accession_string) or (None, None). Cached."""
        canon = _canon(smi)
        if canon is None:
            return None, None
        if canon in self._compound_cache:
            return self._compound_cache[canon]

        cpd = None
        acc = None

        # 1. try the hardcoded SMILES→KEGG map first — cleanest path
        if canon in _SMILES_TO_KEGG:
            kegg = _SMILES_TO_KEGG[canon]
            cpd = self._cc.get_compound(f"kegg:{kegg}")
            if cpd is not None:
                acc = f"kegg:{kegg}"

        # 2. fallback: equilibrator's InChI-key search
        if cpd is None:
            try:
                inchi_key = Chem.InchiToInchiKey(
                    Chem.MolToInchi(Chem.MolFromSmiles(canon))
                )
            except Exception:
                inchi_key = None
            if inchi_key:
                # try RDKit's neutral form and equilibrator's -M form
                for candidate in (inchi_key, inchi_key[:-1] + "M"):
                    results = self._cc.search_compound_by_inchi_key(candidate)
                    if results:
                        cpd = results[0]
                        # coco:{id} is the internal namespace that
                        # parse_reaction_formula recognises for
                        # compound-cache lookups without an external ID.
                        acc = f"coco:{cpd.id}"
                        break

        self._compound_cache[canon] = (cpd, acc)
        return cpd, acc

    # ------------------------------------------------------------------
    # Reaction ΔG'°
    # ------------------------------------------------------------------
    def dG_prime(self, rxn_smiles: str) -> Optional[float]:
        """
        Compute standard transformed reaction Gibbs free energy at pH 7
        for a reaction given in `reactants>>products` SMILES format.
        Multiple reactants/products separated by `.` as usual.

        Returns kJ/mol as a float, or None if any compound in the
        reaction can't be mapped to equilibrator's database.
        """
        if rxn_smiles in self._reaction_cache:
            return self._reaction_cache[rxn_smiles]

        with self._lock:
            # double-check after acquiring the lock in case another
            # thread computed it while we were waiting
            if rxn_smiles in self._reaction_cache:
                return self._reaction_cache[rxn_smiles]

            try:
                lhs, rhs = rxn_smiles.split(">>", 1)
            except ValueError:
                self._reaction_cache[rxn_smiles] = None
                return None

            reactant_smis = [s for s in lhs.split(".") if s.strip()]
            product_smis  = [s for s in rhs.split(".") if s.strip()]

            # Look up every side. If any compound is unmappable we
            # can't score this reaction.
            reactant_pairs = [self._find_compound(s) for s in reactant_smis]
            product_pairs  = [self._find_compound(s) for s in product_smis]
            if any(cpd is None for cpd, _ in reactant_pairs + product_pairs):
                self._reaction_cache[rxn_smiles] = None
                return None

            # Build the equilibrator reaction formula using accessions:
            # "kegg:C00011 + kegg:C00282 = kegg:C00001"
            lhs_str = " + ".join(acc for _, acc in reactant_pairs)
            rhs_str = " + ".join(acc for _, acc in product_pairs)
            formula = f"{lhs_str} = {rhs_str}"
            try:
                rxn = self._cc.parse_reaction_formula(formula)
                # NOTE: we don't gate on rxn.is_balanced(). DORAnet's
                # operators often omit protons or water that
                # equilibrator's atomic-balance check flags as
                # unbalanced. The ΔG'° values still make sense for
                # relative pathway comparison, and refusing to score
                # unbalanced reactions means most bio pathways get no
                # score at all.
                dg = self._cc.standard_dg_prime(rxn)
                # dg is a pint Measurement like `(-45.2 +/- 1.1) kJ/mol`.
                # Extract the nominal value in kJ/mol.
                dg_kJ_per_mol = float(dg.value.magnitude)
                self._reaction_cache[rxn_smiles] = dg_kJ_per_mol
                return dg_kJ_per_mol
            except Exception:
                self._reaction_cache[rxn_smiles] = None
                return None

    # ------------------------------------------------------------------
    # Pathway scoring
    # ------------------------------------------------------------------
    def score_pathway(self, reaction_smiles_list):
        """Given a pathway's list of reaction SMILES strings, return
        (max_dg, avg_dg, coverage_ratio) where:
          - max_dg  : max ΔG'° across scoreable reactions, kJ/mol (or None)
          - avg_dg  : mean ΔG'° across scoreable reactions, kJ/mol (or None)
          - coverage_ratio : fraction of reactions we could score (0..1)
        """
        scores = [self.dG_prime(r) for r in reaction_smiles_list]
        real = [s for s in scores if s is not None]
        coverage = len(real) / max(1, len(scores))
        if not real:
            return None, None, coverage
        return max(real), sum(real) / len(real), coverage
