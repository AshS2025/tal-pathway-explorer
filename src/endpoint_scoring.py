"""
Endpoint "interestingness" scoring for open-exploration mode.

Open exploration expands chemistry forward from a starter without a
specific target. The resulting network contains every reachable
molecule — typically hundreds to thousands. To surface what's worth
looking at, we score each molecule by a transparent heuristic blending
structural complexity, functional-group diversity, and a carbon-count
gate.

The criteria are intentionally simple and defensible so the ranking
can be explained in a report. They are also pluggable — pass different
keyword args to score_interestingness, or add a new term, and the
score updates without touching call sites.

References:
  Bertz, S. H. (1981). "The first general index of molecular complexity."
    J. Am. Chem. Soc. 103, 3599-3601.  (RDKit: GraphDescriptors.BertzCT)
  Ertl, P., & Schuffenhauer, A. (2009). "Estimation of synthetic
    accessibility score..." J. Cheminform. 1, 8.  (not used in V1; add
    via the sascorer contrib when available.)
"""

import dataclasses
from typing import Iterable, Union

from rdkit import Chem
from rdkit.Chem import GraphDescriptors, Lipinski


# Curated set of functional groups for diversity counting. Kept small
# and mutually distinct enough that the count is a meaningful "how many
# functional handles does this molecule offer" signal rather than a
# tally of overlapping SMARTS hits.
_FG_SMARTS = {
    "hydroxyl":   "[OX2H]",
    "carbonyl":   "[CX3]=[OX1]",
    "carboxylic": "[CX3](=O)[OX2H]",
    "ester":      "[#6][CX3](=O)[OX2][#6]",
    "ether":      "[OD2]([#6])[#6]",
    "amine":      "[NX3;H2,H1;!$(NC=O)]",
    "amide":      "[NX3][CX3](=[OX1])",
    "aldehyde":   "[CX3H1](=O)[#6]",
    "nitrile":    "[NX1]#[CX2]",
}
_FG_PATTERNS = {name: Chem.MolFromSmarts(s) for name, s in _FG_SMARTS.items()}


@dataclasses.dataclass(frozen=True)
class EndpointScore:
    """Decomposed interestingness score for a single endpoint molecule."""
    smiles: str
    score: float
    carbons: int
    bertz: float
    n_functional_groups: int
    n_aromatic_rings: int
    has_oxygen: bool
    breakdown: dict


def score_interestingness(
    mol_or_smiles: Union[str, Chem.Mol],
    *,
    carbon_window: tuple = (4, 12),
    require_oxygen: bool = True,
    bertz_saturate: float = 300.0,
    fg_saturate: int = 4,
) -> EndpointScore:
    """
    Score one molecule's interestingness on [0, 1].

    Pipeline (geometric — a zero on any axis zeroes the score):

      1. Carbon window gate     : 0 if C count outside (carbon_window)
      2. Oxygen gate            : 0 if require_oxygen and no O atom
      3. Complexity term        : min(Bertz / bertz_saturate, 1.0)
      4. FG diversity term      : min(n_distinct_FGs / fg_saturate, 1.0)
      5. Aromatic bonus         : 1.0 if any aromatic ring, else 0.6

      score = complexity * fg_diversity * aromatic_bonus

    Why product (not sum): a molecule with high Bertz but no functional
    handles isn't an interesting endpoint, it's an exotic scaffold with
    nothing to do. Multiplying enforces "good on every axis."

    Why those defaults:
      carbon_window (4, 12) — TAL is 6C; common derivatives 4-12C
      bertz_saturate 300    — saturates at "moderate complexity"
                              (TAL ~140, sorbic acid ~80, phloroglucinol ~200)
      fg_saturate 4         — 4 distinct FGs ≈ a richly functionalized
                              small molecule
    """
    if isinstance(mol_or_smiles, str):
        smiles_in = mol_or_smiles
        mol = Chem.MolFromSmiles(smiles_in)
    else:
        mol = mol_or_smiles
        smiles_in = Chem.MolToSmiles(mol) if mol is not None else ""

    if mol is None:
        return EndpointScore(
            smiles=smiles_in, score=0.0, carbons=0, bertz=0.0,
            n_functional_groups=0, n_aromatic_rings=0, has_oxygen=False,
            breakdown={"reason": "unparseable SMILES"},
        )

    smiles = Chem.MolToSmiles(mol)
    carbons = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
    has_oxygen = any(a.GetSymbol() == "O" for a in mol.GetAtoms())

    if not (carbon_window[0] <= carbons <= carbon_window[1]):
        return EndpointScore(
            smiles=smiles, score=0.0, carbons=carbons, bertz=0.0,
            n_functional_groups=0, n_aromatic_rings=0,
            has_oxygen=has_oxygen,
            breakdown={"reason": f"C={carbons} outside {carbon_window}"},
        )

    if require_oxygen and not has_oxygen:
        return EndpointScore(
            smiles=smiles, score=0.0, carbons=carbons, bertz=0.0,
            n_functional_groups=0, n_aromatic_rings=0, has_oxygen=False,
            breakdown={"reason": "no oxygen (require_oxygen=True)"},
        )

    bertz = float(GraphDescriptors.BertzCT(mol))
    n_fgs = sum(
        1 for pat in _FG_PATTERNS.values() if mol.HasSubstructMatch(pat)
    )
    n_aro = int(Lipinski.NumAromaticRings(mol))

    complexity_term = min(bertz / bertz_saturate, 1.0)
    fg_term = min(n_fgs / fg_saturate, 1.0) if fg_saturate > 0 else 1.0
    aromatic_term = 1.0 if n_aro > 0 else 0.6

    score = complexity_term * fg_term * aromatic_term

    return EndpointScore(
        smiles=smiles,
        score=score,
        carbons=carbons,
        bertz=bertz,
        n_functional_groups=n_fgs,
        n_aromatic_rings=n_aro,
        has_oxygen=has_oxygen,
        breakdown={
            "complexity": complexity_term,
            "fg_diversity": fg_term,
            "aromatic_bonus": aromatic_term,
        },
    )


def rank_network_endpoints(
    network,
    *,
    exclude_smiles: Iterable[str] = (),
    top_n: int = 20,
    **scoring_kwargs,
) -> list[EndpointScore]:
    """
    Score every molecule in a DORAnet network and return the top_n by
    interestingness. Skips any SMILES in `exclude_smiles` (starters,
    helpers, cofactors) after canonicalization.
    """
    excluded_canonical = set()
    for smi in exclude_smiles:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            excluded_canonical.add(Chem.MolToSmiles(m))

    scored: list[EndpointScore] = []
    seen: set[str] = set()
    for mol in network.mols:
        # DORAnet MolDatBasicV1 exposes the canonical SMILES on
        # `.smiles` (and the uid is also a SMILES). No `.item` wrapper.
        smi = getattr(mol, "smiles", None) or getattr(mol, "uid", None)
        if not smi:
            continue
        # Re-canonicalize so excluded_canonical comparison is apples-to-apples.
        rd = Chem.MolFromSmiles(smi)
        if rd is None:
            continue
        smi = Chem.MolToSmiles(rd)
        if smi in seen or smi in excluded_canonical:
            continue
        seen.add(smi)
        es = score_interestingness(rd, **scoring_kwargs)
        if es.score > 0:
            scored.append(es)

    scored.sort(key=lambda s: s.score, reverse=True)
    return scored[:top_n]
