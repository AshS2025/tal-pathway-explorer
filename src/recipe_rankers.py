"""
recipe_rankers.py
RecipeRanker implementations for DORAnet priority-queue expansion.

WHAT A RANKER IS
----------------
A RecipeRanker is a callable that scores candidate (operator, reactants)
combinations BEFORE the reaction is actually run. The priority-queue
strategy uses the score to choose which recipes to expand first; combined
with `beam_size`, that yields beam-search pruning, which keeps network
growth tractable at 5+ generations.

CONVENTION (verified against doranet/strategies.py:1017-1028, :790-797):
    Higher rank value = higher priority. DORAnet's RecipeHeap keeps the
    top items by rank and pops them in descending order.

DESIGN PHILOSOPHY: COMPOSABILITY
--------------------------------
Different users will care about different criteria — similarity to a
target, molecular weight fit, atom economy, cost, thermodynamic ΔG, etc.
Each criterion lives in its own small ranker class with its own knobs.
`WeightedCompositeRanker` combines any subset with user-chosen weights,
so the call site fully controls the search direction. Don't hard-code
a single objective here.

PRODUCT- vs. REACTANT-SIDE SCORING
----------------------------------
Forward goal-directed search should rank by predicted *product*
similarity to target, not reactant similarity (small starting fragments
will always look unlike a large target). But the strategy hands the
ranker only (operator, reactants) — products don't exist yet. To get
product-side semantics, `ForwardProductTanimotoRanker` applies the
operator's RDKit reaction itself to predict products, then fingerprints
those predictions. This doubles the operator's RunReactants work per
ranked recipe; that cost is bounded by the heap and is the price of
correct semantics.

Retrosynthesis is fundamentally different: fragments produced by
backward expansion will always have low Tanimoto vs. the parent target,
so global similarity is the wrong signal. Forward rankers in this file
should not be used for retro. (A retro-appropriate ranker — MCS
preservation or catalog matching — is planned but not implemented here.)
"""

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

from doranet import interfaces


# =====================================================================
# Helpers
# =====================================================================

def _morgan_fp(rdkit_mol, radius, n_bits):
    return AllChem.GetMorganFingerprintAsBitVect(
        rdkit_mol, radius, nBits=n_bits
    )


# =====================================================================
# Criterion: Tanimoto similarity of predicted products to a target
# =====================================================================

class ForwardProductTanimotoRanker(interfaces.RecipeRanker):
    """
    Forward-mode ranker. Predicts products by running the recipe's
    operator on its reactants, then returns the MAX Morgan-fingerprint
    Tanimoto similarity between any predicted product and `target_smiles`.

    Why MAX over products: reactions often produce one "real" product and
    one or more byproducts; we want recipes whose main product is close
    to the target, not recipes whose byproducts happen to be small water-
    like molecules with high baseline similarity to nothing.

    Why predict products in the ranker (instead of scoring reactants):
    forward goal-directed search needs the score to reflect what we'd be
    adding to the network if we expanded this recipe. Reactants are
    typically smaller fragments whose similarity to a large target is
    meaningless. (See module docstring.)

    Caching: predicted-product fingerprints are cached by canonical SMILES
    so the same product structure isn't fingerprinted twice across the
    whole expansion.

    Parameters
    ----------
    target_smiles : str
        Target molecule (e.g. TAL).
    radius : int
        Morgan fingerprint radius (default 2 = ECFP4-like).
    n_bits : int
        Folded fingerprint bit length (default 2048).
    """

    def __init__(self, target_smiles, radius=2, n_bits=2048):
        target_mol = Chem.MolFromSmiles(target_smiles)
        if target_mol is None:
            raise ValueError(f"Invalid target SMILES: {target_smiles!r}")
        self.target_smiles = target_smiles
        self.radius = radius
        self.n_bits = n_bits
        self.target_fp = _morgan_fp(target_mol, radius, n_bits)
        self._fp_cache = {}

    def _fp_for_predicted(self, rdkit_mol):
        try:
            Chem.SanitizeMol(rdkit_mol)
        except Exception:
            return None
        smi = Chem.MolToSmiles(rdkit_mol)
        fp = self._fp_cache.get(smi)
        if fp is None:
            fp = _morgan_fp(rdkit_mol, self.radius, self.n_bits)
            self._fp_cache[smi] = fp
        return fp

    def __call__(self, recipe, min_rank=None):
        op = recipe.operator.item
        if not hasattr(op, "rdkitrxn"):
            # Non-RDKit operator type; we can't predict products.
            return 0.0

        reactant_mols = []
        for r in recipe.reactants:
            if not isinstance(r.item, interfaces.MolDatRDKit):
                return 0.0
            reactant_mols.append(r.item.rdkitmol)

        try:
            product_sets = op.rdkitrxn.RunReactants(tuple(reactant_mols))
        except Exception:
            return 0.0

        max_sim = 0.0
        for product_set in product_sets:
            for product in product_set:
                fp = self._fp_for_predicted(product)
                if fp is None:
                    continue
                sim = DataStructs.TanimotoSimilarity(self.target_fp, fp)
                if sim > max_sim:
                    max_sim = sim
        return max_sim

    @property
    def meta_required(self):
        return interfaces.MetaKeyPacket()


# =====================================================================
# Retro: feedstock-proximity ranker (max sim across feedstock pool)
# =====================================================================

class FeedstockProximityRanker(interfaces.RecipeRanker):
    """
    Retro-mode ranker. In retro expansion, an operator's RDKit reaction
    object generates predicted upstream substrates from a downstream
    product. This ranker scores each recipe by the MAX Tanimoto
    similarity between any predicted substrate and any molecule in a
    user-supplied feedstock pool.

    Why this works for retro
    ------------------------
    Forward goal-directed search via Tanimoto-to-target only works
    because forward chemistry GROWS molecules toward the target. Retro
    chemistry FRAGMENTS molecules backward; the resulting upstream
    pieces don't look like the final target. But they DO need to
    eventually look like a feedstock. Pulling the beam toward recipes
    whose predicted substrates resemble glucose (or any other feedstock)
    gives the search a directional signal that the forward
    Tanimoto-to-target heuristic explicitly lacks for retro.

    Pool semantics
    --------------
    Score = max over (predicted_substrate, feedstock) pairs. A recipe
    whose predicted substrates match ANY feedstock highly scores high;
    we don't require matching all feedstocks. Start with a single-entry
    pool (just glucose) to validate the path; later add glycerol,
    acetate, acetyl-CoA, etc.

    Parameters
    ----------
    feedstock_smiles_list : list[str]
        Pool of candidate feedstock SMILES. Empty pool raises.
    radius, n_bits : int
        Morgan fingerprint params (defaults match ForwardProductTanimotoRanker).
    """

    def __init__(self, feedstock_smiles_list, radius=2, n_bits=2048):
        if not feedstock_smiles_list:
            raise ValueError("FeedstockProximityRanker needs at least one feedstock")
        self.radius = radius
        self.n_bits = n_bits
        self.feedstock_fps = []
        for smi in feedstock_smiles_list:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                raise ValueError(f"Invalid feedstock SMILES: {smi!r}")
            self.feedstock_fps.append(_morgan_fp(m, radius, n_bits))
        self._fp_cache = {}

    def _fp_for_predicted(self, rdkit_mol):
        try:
            Chem.SanitizeMol(rdkit_mol)
        except Exception:
            return None
        smi = Chem.MolToSmiles(rdkit_mol)
        fp = self._fp_cache.get(smi)
        if fp is None:
            fp = _morgan_fp(rdkit_mol, self.radius, self.n_bits)
            self._fp_cache[smi] = fp
        return fp

    def __call__(self, recipe, min_rank=None):
        op = recipe.operator.item
        if not hasattr(op, "rdkitrxn"):
            return 0.0
        reactant_mols = []
        for r in recipe.reactants:
            if not isinstance(r.item, interfaces.MolDatRDKit):
                return 0.0
            reactant_mols.append(r.item.rdkitmol)
        try:
            product_sets = op.rdkitrxn.RunReactants(tuple(reactant_mols))
        except Exception:
            return 0.0

        max_sim = 0.0
        for product_set in product_sets:
            for product in product_set:
                fp = self._fp_for_predicted(product)
                if fp is None:
                    continue
                for feed_fp in self.feedstock_fps:
                    sim = DataStructs.TanimotoSimilarity(feed_fp, fp)
                    if sim > max_sim:
                        max_sim = sim
        return max_sim

    @property
    def meta_required(self):
        return interfaces.MetaKeyPacket()


# =====================================================================
# Criterion: Molecular-weight fit to a target MW
# =====================================================================

class ProductMWRanker(interfaces.RecipeRanker):
    """
    Forward-mode ranker. Predicts products by running the operator on
    its reactants, then scores each predicted product by how close its
    molecular weight is to `target_mw`.

    Score per product: max(0, 1 - |MW_product - target_mw| / tolerance).
    A product exactly at `target_mw` scores 1.0; one more than
    `tolerance` Da away scores 0. The ranker returns the MAX over
    predicted products (same rationale as ForwardProductTanimotoRanker).

    Useful as a complementary criterion to Tanimoto — molecules near
    target MW that AREN'T similar (scaffold hops) are interesting in
    their own right.

    Parameters
    ----------
    target_mw : float
        Target molecular weight in Daltons (e.g. TAL = 126.11).
    tolerance : float
        MW window over which the score linearly drops from 1 to 0.
        Default 50 Da.
    """

    def __init__(self, target_mw, tolerance=50.0):
        if tolerance <= 0:
            raise ValueError("tolerance must be positive")
        self.target_mw = float(target_mw)
        self.tolerance = float(tolerance)
        self._mw_cache = {}

    def _mw_for_predicted(self, rdkit_mol):
        try:
            Chem.SanitizeMol(rdkit_mol)
        except Exception:
            return None
        smi = Chem.MolToSmiles(rdkit_mol)
        mw = self._mw_cache.get(smi)
        if mw is None:
            mw = Descriptors.MolWt(rdkit_mol)
            self._mw_cache[smi] = mw
        return mw

    def __call__(self, recipe, min_rank=None):
        op = recipe.operator.item
        if not hasattr(op, "rdkitrxn"):
            return 0.0

        reactant_mols = []
        for r in recipe.reactants:
            if not isinstance(r.item, interfaces.MolDatRDKit):
                return 0.0
            reactant_mols.append(r.item.rdkitmol)

        try:
            product_sets = op.rdkitrxn.RunReactants(tuple(reactant_mols))
        except Exception:
            return 0.0

        max_score = 0.0
        for product_set in product_sets:
            for product in product_set:
                mw = self._mw_for_predicted(product)
                if mw is None:
                    continue
                score = max(0.0, 1.0 - abs(mw - self.target_mw) / self.tolerance)
                if score > max_score:
                    max_score = score
        return max_score

    @property
    def meta_required(self):
        return interfaces.MetaKeyPacket()


# =====================================================================
# Criterion: thermodynamic favorability of the predicted reaction
# =====================================================================

class ProductThermoRanker(interfaces.RecipeRanker):
    """
    Score a recipe by the predicted reaction ΔH (or ΔG, depending on what
    `molecule_thermo_calculator` returns) of running its operator. Soft
    mirror of the existing `Rxn_dH_Filter`: same stoichiometric math, but
    used to steer the beam instead of to hard-reject.

    Why this exists alongside Rxn_dH_Filter
    ---------------------------------------
    Rxn_dH_Filter runs at reaction-time on REAL products — authoritative,
    used for hard rejection. ProductThermoRanker runs at recipe-time on
    PREDICTED products — approximate, used for soft ordering. Together
    they implement the same "same criterion, two layers" pattern as
    Tanimoto: ranker picks the best 50 of 800 candidates per generation;
    filter then rejects any whose real thermo is unacceptable.

    Without this ranker, thermodynamics gates only AFTER you've spent
    beam budget on a recipe. With it, recipes likely to fail the thermo
    filter are deprioritized before they reach expansion.

    Math
    ----
    Replicates `Chem_Rxn_dH_Calculator` (doranet/modules/synthetic/
    generate_network.py:23-90):

        raw_dH = Σ(stoi_p · Hf_predicted_p) − Σ(stoi_r · Hf_r)
        if direction == "forward":  dH = (raw_dH + enthalpy_correction) / n_steps
        if direction == "retro":    dH = (-raw_dH + enthalpy_correction) / n_steps

    Then maps dH to a score in [0, 1]:

        dH <= dH_favorable:            score = 1.0   (great)
        dH >= dH_favorable + penalty:  score = 0.0   (uninteresting)
        in between:                    linear interp

    Multi-product-set handling: an operator may produce several alternate
    product sets from one (operator, reactants) tuple. We score each set
    and return the MAX (best alternative wins), same convention as
    ForwardProductTanimotoRanker.

    Missing data: if any Hf lookup returns None for a predicted product
    or a reactant, that product set scores 0 (uncertain → deprioritized,
    but not dropped from the heap — None would mean "drop").

    Parameters
    ----------
    molecule_thermo_calculator : Callable[[str], float | None]
        Same interface as the existing `Mole_Hf` callback: takes a SMILES
        string and returns its standard formation energy (Hf or Gf —
        the math is the same, only the interpretation differs). Returning
        None signals "no data" and the recipe gets score 0.
    direction : str
        "forward" or "retro". Mirrors Chem_Rxn_dH_Calculator's direction
        handling.
    dH_favorable : float
        Threshold at and below which the score saturates at 1.0. Default
        0 (anything exothermic is fully credited).
    penalty_range : float
        How many units above dH_favorable before the score reaches 0.
        Default 30 (matches the hard-filter range with some headroom).
        Units follow whatever your calculator returns — kcal/mol or
        kJ/mol — same caveat as elsewhere in the pipeline.
    """

    def __init__(
        self,
        molecule_thermo_calculator,
        direction="forward",
        dH_favorable=0.0,
        penalty_range=30.0,
    ):
        if molecule_thermo_calculator is None:
            raise ValueError(
                "ProductThermoRanker requires a thermo calculator; got None. "
                "Provide the same Mole_Hf callable you pass to "
                "Chem_Rxn_dH_Calculator."
            )
        if direction not in ("forward", "retro"):
            raise ValueError(f"direction must be 'forward' or 'retro', got {direction!r}")
        if penalty_range <= 0:
            raise ValueError("penalty_range must be positive")
        self.Mole_Hf = molecule_thermo_calculator
        self.direction = direction
        self.dH_favorable = float(dH_favorable)
        self.penalty_range = float(penalty_range)
        # Hf is expensive (equilibrator calls hit a database/network).
        # Cache results by canonical SMILES so we don't re-query.
        self._hf_cache = {}

    def _hf_for(self, smiles):
        hf = self._hf_cache.get(smiles)
        if hf is None and smiles not in self._hf_cache:
            hf = self.Mole_Hf(smiles)
            self._hf_cache[smiles] = hf
        return hf

    def _score_dH(self, dH):
        if dH <= self.dH_favorable:
            return 1.0
        if dH >= self.dH_favorable + self.penalty_range:
            return 0.0
        return 1.0 - (dH - self.dH_favorable) / self.penalty_range

    def __call__(self, recipe, min_rank=None):
        op = recipe.operator.item
        if not hasattr(op, "rdkitrxn"):
            return 0.0
        op_meta = recipe.operator.meta
        if op_meta is None:
            return 0.0
        try:
            reactants_stoi = op_meta["reactants_stoi"]
            products_stoi = op_meta["products_stoi"]
            n_steps = op_meta["number_of_steps"]
        except KeyError:
            return 0.0
        correction = op_meta.get("enthalpy_correction") or 0

        # Sum reactant Hf contribution once — it doesn't change across
        # alternate product sets from the same recipe.
        rea_sum = 0.0
        for idx, r in enumerate(recipe.reactants):
            if not isinstance(r.item, interfaces.MolDatRDKit):
                return 0.0
            hf = self._hf_for(r.item.smiles)
            if hf is None:
                return 0.0
            rea_sum += hf * reactants_stoi[idx]

        reactant_mols = [r.item.rdkitmol for r in recipe.reactants]
        try:
            product_sets = op.rdkitrxn.RunReactants(tuple(reactant_mols))
        except Exception:
            return 0.0
        if not product_sets:
            return 0.0

        max_score = 0.0
        for product_set in product_sets:
            # Stoichiometry vector is positional — if RunReactants returns
            # a product tuple of unexpected length, skip this set rather
            # than misalign indices.
            if len(product_set) != len(products_stoi):
                continue
            pro_sum = 0.0
            ok = True
            for idx, product in enumerate(product_set):
                try:
                    Chem.SanitizeMol(product)
                except Exception:
                    ok = False
                    break
                hf = self._hf_for(Chem.MolToSmiles(product))
                if hf is None:
                    ok = False
                    break
                pro_sum += hf * products_stoi[idx]
            if not ok:
                continue

            raw = pro_sum - rea_sum
            if self.direction == "forward":
                dH = (raw + correction) / n_steps
            else:  # retro
                dH = (-raw + correction) / n_steps
            score = self._score_dH(dH)
            if score > max_score:
                max_score = score
        return max_score

    @property
    def meta_required(self):
        # We need operator stoichiometry and correction metadata to do
        # the dH math — declare them so the strategy fetches them.
        return interfaces.MetaKeyPacket(
            operator_keys={
                "reactants_stoi",
                "products_stoi",
                "enthalpy_correction",
                "number_of_steps",
            }
        )


# =====================================================================
# Composite: weighted sum of any number of criterion rankers
# =====================================================================

class WeightedCompositeRanker(interfaces.RecipeRanker):
    """
    Combine multiple criterion rankers into a single score via a
    user-defined weighted sum. Each criterion ranker is expected to
    return a non-negative number; the composite returns
    sum(weight_i * score_i) / sum(weight_i).

    This is the seam where user preference enters the search. Build the
    composite at the call site with whatever weights make sense for the
    current task:

        ranker = WeightedCompositeRanker([
            (ForwardProductTanimotoRanker("CC(=O)CC(=O)O"), 1.0),
            (ProductMWRanker(target_mw=126.0, tolerance=80.0), 0.3),
        ])

    DORAnet ships its own `CompositeRecipeRanker` (lexicographic, tuple-
    valued) which is appropriate when one criterion strictly dominates
    another. Use that when you have a hard preference order. Use this
    one when you want to trade criteria off smoothly.

    Parameters
    ----------
    components : list of (RecipeRanker, float) tuples
        (ranker, weight). Weights need not sum to 1 — the composite
        normalizes by the weight total.
    """

    def __init__(self, components):
        if not components:
            raise ValueError("Need at least one (ranker, weight) component")
        weight_sum = sum(w for _, w in components)
        if weight_sum <= 0:
            raise ValueError("Sum of weights must be positive")
        self.components = tuple(components)
        self._weight_sum = float(weight_sum)

    def __call__(self, recipe, min_rank=None):
        total = 0.0
        for ranker, weight in self.components:
            score = ranker(recipe)
            if score is None:
                continue
            total += weight * float(score)
        return total / self._weight_sum

    @property
    def meta_required(self):
        keys = interfaces.MetaKeyPacket()
        for ranker, _ in self.components:
            keys = keys + ranker.meta_required
        return keys
