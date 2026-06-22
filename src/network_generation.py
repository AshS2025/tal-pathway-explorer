"""
network_generation.py
TAL biosynthetic & chemical pathway network generation via DORAnet.

Changes from previous version:
  - Minimum_Carbon_Count_Filter: replaced Python atom loop with native
    rdMolDescriptors C++ descriptors (faster at scale)
  - Forward pipeline reordered: cheap structural filters first,
    thermo calculator last (unchanged in retro — already logical)
  - Default max_molecular_weight tightened 800 → 500 Da
    (TAL = 126 Da; 800 was too permissive for TAL derivative space)
  - Cartesian depth warning added for gen > 3
  - Cleaned up commented-out code; future hooks marked clearly
"""

import doranet as dn
from doranet.modules.synthetic.Reaction_Smarts_Forward import op_smarts
from doranet.modules.synthetic.Reaction_Smarts_Retro import op_retro_smarts
from doranet.modules.synthetic.generate_network import (
    get_smiles_from_file,
    Max_Atoms_Filter,
    Ring_Issues_Filter,
    Enol_filter_forward,
    Enol_filter_retro,
    Check_balance_filter,
    Allowed_Elements_Filter,
    Chem_Rxn_dH_Calculator,
    Rxn_dH_Filter,
    Cross_Reaction_Filter,
    Retro_Not_Aromatic_Filter,
)
# Bio side: cofactor table, ruleset paths, and the canonical SMILES set
# used to identify which molecules in a recipe are cofactors vs. true
# substrates. We import directly so we stay in lockstep with DORAnet's
# definitions; if they update the cofactor list, we pick that up.
from doranet.modules.enzymatic.generate_network import (
    AVAILABLE_RULESETS,
    cofactors_clean_dict,
    cofactors_clean,
    clean_SMILES,
)
import pandas as pd
from datetime import datetime
import time
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from doranet import metadata, interfaces
import collections.abc
import dataclasses
from tal_reaction_whitelist import TAL_REACTION_WHITELIST
from tal_bio_reaction_whitelist import TAL_BIO_REACTION_WHITELIST
from recipe_rankers import ForwardProductTanimotoRanker

# Bio cofactors DORAnet excludes by default — these tend to over-fire
# and produce nonsense intermediates. Keep the same default here.
_DEFAULT_EXCLUDED_COFACTORS = ("CARBONYL_CoF", "AMINO_CoF")


# ============ CUSTOM FILTERS ============

class Molecular_Weight_Filter(metadata.ReactionFilterBase):
    """
    Reject reactions if any product exceeds max molecular weight.
    Uses a single RDKit descriptor lookup — cheapest possible check,
    so this runs first in the pipeline.
    """

    def __init__(self, max_weight=500):
        self.max_weight = max_weight

    def __call__(self, recipe):
        for mol in recipe.products:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"Filter only works with MolDatRDKit, not {type(mol.item)}"
                )
            if Descriptors.MolWt(mol.item.rdkitmol) > self.max_weight:
                return False
        return True

    @property
    def meta_required(self):
        return interfaces.MetaKeyPacket()


class Minimum_Carbon_Count_Filter(metadata.ReactionFilterBase):
    """
    Reject reactions if any product has fewer than min_carbons carbon atoms.

    OPTIMIZATION: replaced Python atom iteration with native RDKit C++
    descriptors. CalcNumAliphaticCarbons + CalcNumAromaticCarbons together
    cover every carbon in the molecule without double-counting, and run
    in compiled C++ rather than interpreted Python — meaningfully faster
    when evaluating large numbers of product nodes.

    Previous implementation (slow):
        sum(1 for atom in mol.item.rdkitmol.GetAtoms() if atom.GetAtomicNum() == 6)

    Current implementation (fast):
        rdMolDescriptors.CalcNumAliphaticCarbons(rdmol)
        + rdMolDescriptors.CalcNumAromaticCarbons(rdmol)
    """

    def __init__(self, min_carbons=0):
        self.min_carbons = min_carbons

    def __call__(self, recipe):
        for mol in recipe.products:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"Filter only works with MolDatRDKit, not {type(mol.item)}"
                )
            rdmol = mol.item.rdkitmol
            # Count all carbons via atomic number — rdMolDescriptors
            # has no CalcNumAliphaticCarbons; iterating atoms is
            # equivalent and avoids the missing-symbol crash.
            carbon_count = sum(
                1 for a in rdmol.GetAtoms() if a.GetAtomicNum() == 6
            )
            if carbon_count < self.min_carbons:
                return False
        return True

    @property
    def meta_required(self):
        return interfaces.MetaKeyPacket()


class Product_Tanimoto_Filter(metadata.ReactionFilterBase):
    """
    Reject reactions where no product is within `min_similarity` Tanimoto
    of the target molecule (Morgan fingerprint, radius 2, 2048 bits).

    This is a HARD filter (accept / reject), distinct from the soft
    ordering done by a ranker. The ranker only decides which recipes get
    expanded first; this filter decides whether a reaction's products
    are even allowed into the network. Use both together: ranker to
    spend the beam budget wisely, filter to keep the network clean.

    Runs on real products (after the operator executes), so it sees
    authoritative structures — unlike the ranker, which works from
    predicted products.

    Parameters
    ----------
    target_smiles : str
        Target molecule.
    min_similarity : float
        Threshold in [0, 1]. A reaction is kept iff at least one of its
        products has Tanimoto similarity >= this value vs. the target.
        Typical values: ~0.45 for scaffold-hopping exploration, ~0.7
        for lead-optimization-style refinement.
    radius : int
        Morgan radius (default 2).
    n_bits : int
        Fingerprint bit length (default 2048).
    """

    def __init__(self, target_smiles, min_similarity, radius=2, n_bits=2048):
        target_mol = Chem.MolFromSmiles(target_smiles)
        if target_mol is None:
            raise ValueError(f"Invalid target SMILES: {target_smiles!r}")
        if not 0.0 <= min_similarity <= 1.0:
            raise ValueError("min_similarity must be in [0, 1]")
        self.target_fp = AllChem.GetMorganFingerprintAsBitVect(
            target_mol, radius, nBits=n_bits
        )
        self.min_similarity = float(min_similarity)
        self.radius = radius
        self.n_bits = n_bits
        self._fp_cache = {}

    def _fp_for(self, mol_data):
        smi = mol_data.smiles
        fp = self._fp_cache.get(smi)
        if fp is None:
            fp = AllChem.GetMorganFingerprintAsBitVect(
                mol_data.rdkitmol, self.radius, nBits=self.n_bits
            )
            self._fp_cache[smi] = fp
        return fp

    def __call__(self, recipe):
        for mol in recipe.products:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                continue
            sim = DataStructs.TanimotoSimilarity(
                self.target_fp, self._fp_for(mol.item)
            )
            if sim >= self.min_similarity:
                return True
        return False

    @property
    def meta_required(self):
        return interfaces.MetaKeyPacket()


# ============ BIO-AWARE RECIPE FILTER ============

@dataclasses.dataclass(frozen=True)
class Bio_Single_Substrate_Filter(interfaces.RecipeFilter):
    """
    Recipe filter that enforces "one non-cofactor substrate per reaction"
    on bio operators ONLY. Chem operators pass through untouched.

    This mirrors DORAnet's bio `Reaction_Type_Filter` (enzymatic
    generate_network.py:134) but is conditional on `is_bio` operator
    meta, so it composes safely with chem operators in the same network.

    Why we need it:
      Without it, two cofactors can react together as a "valid" bio
      recipe (CO2 + H2O on a carboxylase, e.g.) — which inflates the
      network with biochemically meaningless reactions.

    Why it can be opted out (bio_allow_multiple_reactants=True):
      The polyketide route to TAL is 3 × acetyl-CoA → TAL. All three
      reactants are cofactors. With this filter enforcing single-non-
      cofactor we'd reject the entire polyketide condensation. For TAL
      forward expansion you want it OFF.
    """

    cofactor_smiles_set: frozenset

    def __call__(self, recipe) -> bool:
        if recipe.operator.meta is None:
            return True
        if not recipe.operator.meta.get("is_bio", False):
            # Chem operator — no extra constraint.
            return True
        # Bio operator — count non-cofactor reactants.
        non_cofactor = set()
        for mol in recipe.reactants:
            if mol.meta is None:
                continue
            smi = mol.meta.get("SMILES")
            if smi is None:
                continue
            if clean_SMILES(smi) not in self.cofactor_smiles_set:
                non_cofactor.add(clean_SMILES(smi))
        return len(non_cofactor) == 1

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={"is_bio"}, molecule_keys={"SMILES"}
        )


def _load_bio_ops_into_network(
    network,
    engine,
    ruleset,
    direction,
    excluded_cofactors,
    whitelist,
):
    """
    Load bio rules from the TSV, pre-filter through the TAL bio
    whitelist, and register each as a network operator with meta
    populated for BOTH the bio Reactants/Comments columns AND the
    chem-style meta keys (ring_issue, kekulize_flag, allowed_elements,
    etc.) so the existing chem reaction_plan filters treat bio
    operators as no-ops where appropriate.

    Returns (n_loaded, n_in_whitelist).
    """
    if ruleset not in AVAILABLE_RULESETS:
        raise ValueError(
            f"Unknown bio ruleset '{ruleset}'. "
            f"Available: {list(AVAILABLE_RULESETS.keys())}"
        )
    rules = pd.read_csv(AVAILABLE_RULESETS[ruleset], sep="\t")
    n_whitelisted = 0
    n_loaded = 0
    for idx, raw_smarts in enumerate(rules["SMARTS"]):
        name = rules["Name"][idx]
        if name not in whitelist:
            continue
        n_whitelisted += 1
        reactant_types = rules["Reactants"][idx].split(";")
        product_types = rules["Products"][idx].split(";")
        if excluded_cofactors and (
            set(excluded_cofactors) & set(reactant_types)
            or set(excluded_cofactors) & set(product_types)
        ):
            continue

        # ---- direction handling -----------------------------------
        # JN1224MIN's SMARTS are written forward: LHS = substrates,
        # RHS = products. For retro, we want to fire enzymes "in
        # reverse" — i.e., from a target product, find the substrate(s)
        # that would have produced it. We do that by flipping the
        # SMARTS at the >> token AND swapping the Reactants/Products
        # cofactor-slot labels so the meta stays consistent with what
        # the operator now consumes/produces.
        #
        # Note: most enzymatic reactions are operationally reversible
        # under their natural conditions (Haldane), so flipping the
        # SMARTS is biochemically defensible. The exception is irreversible
        # decarboxylations etc — those would need an explicit irreversible
        # flag we don't carry yet. Flagged as a known limitation for now.
        if direction == "retro":
            lhs, rhs = raw_smarts.split(">>", 1)
            smarts_str = f"{rhs}>>{lhs}"
            reactants_meta = rules["Products"][idx]
            products_meta = rules["Reactants"][idx]
        else:
            smarts_str = raw_smarts
            reactants_meta = rules["Reactants"][idx]
            products_meta = rules["Products"][idx]

        n_reactants = len(smarts_str.split(">>")[0].split("."))
        n_products = len(smarts_str.split(">>")[1].split("."))
        meta = {
            "name": name,
            # Bio-specific meta — used by enzyme lookup, also kept so
            # the bio recipe filter can introspect Reactants if needed.
            "Reactants": reactants_meta,
            "Products": products_meta,
            "Comments": rules["Comments"][idx],  # UniProt enzyme IDs
            "SMARTS": smarts_str,
            "is_bio": True,
            # Chem-compatible defaults so the chem reaction_plan
            # filters don't crash on bio recipes:
            "reactants_stoi": (1,) * n_reactants,
            "products_stoi": (1,) * n_products,
            "enthalpy_correction": 0,
            "ring_issue": False,
            "kekulize_flag": False,
            "Retro_Not_Aromatic": False,
            "number_of_steps": 1,
            "allowed_elements": ("All",),
            "Reaction_type": "Enzymatic",
            "Reaction_direction": direction,
        }
        network.add_op(
            engine.op.rdkit(smarts_str, kekulize=False, drop_errors=True),
            meta=meta,
        )
        n_loaded += 1
    return n_loaded, n_whitelisted


# ============ CUSTOM GENERATE FUNCTION ============

# Cartesian strategy becomes impractical beyond this depth on local hardware.
# Priority queue is recommended for gen > 3.
_CARTESIAN_GEN_WARNING_THRESHOLD = 3


def generate_network_tal(
    job_name="default_job",
    starters=False,
    helpers=False,
    gen=2,
    direction="forward",
    molecule_thermo_calculator=None,
    max_rxn_thermo_change=15,
    max_atoms=None,
    max_molecular_weight=500,   # Tightened from 800 — TAL is 126 Da
    allow_multiple_reactants="default",
    targets=None,
    strategy="cartesian",
    min_carbons=0,
    # --- priority-queue / beam-search controls ---
    recipe_ranker=None,
    beam_size=50,
    heap_size=None,
    min_product_similarity=None,
    # --- chem / bio toggles ---
    include_chem=True,
    include_bio=True,
    bio_ruleset="JN1224MIN",
    bio_allow_multiple_reactants=False,
    excluded_cofactors=_DEFAULT_EXCLUDED_COFACTORS,
    # --- future hooks ---
    # cost_calculator=None,     # TODO: add cost-based branch pruning
):
    """
    Generate a TAL biosynthetic/chemical pathway network via DORAnet.

    Parameters
    ----------
    job_name : str
        Label for output files.
    starters : str
        Path to file containing starter molecule SMILES.
    helpers : str
        Path to file containing helper/co-reactant SMILES.
    gen : int
        Number of expansion generations.
        Cartesian strategy will warn if gen > 3.
    direction : str
        "forward" or "retro".
    molecule_thermo_calculator : object
        Thermo calculator passed to Chem_Rxn_dH_Calculator.
    max_rxn_thermo_change : float
        Maximum absolute reaction enthalpy change allowed (units match
        your thermo calculator — confirm kcal/mol vs kJ/mol).
    max_atoms : dict
        Per-element atom count limits, e.g. {"C": 50, "O": 8, "N": 2}.
    max_molecular_weight : float
        Max product molecular weight in Daltons. Default tightened to
        500 Da (was 800). Adjust if your targets are larger.
    allow_multiple_reactants : bool or "default"
        Controls cross-reaction filtering.
    targets : str or list
        Target SMILES for priority_queue ranker and post-expansion check.
    strategy : str
        "cartesian" — exhaustive, exponential scaling.
        "priority_queue" — beam-search style; uses `recipe_ranker` and
        `beam_size` to prune. Required for gen > 3 in practice.
    min_carbons : int
        Minimum carbon count for any product molecule to be retained.
    recipe_ranker : RecipeRanker or None
        Used only when strategy="priority_queue". Decides which
        candidate recipes get expanded first each iteration. If None
        and strategy="priority_queue" and direction="forward" and a
        target is given, a ForwardProductTanimotoRanker is auto-built
        against the first target. For retro, defaults to None (default
        heap order) — global similarity is the wrong signal for retro
        fragments; supply your own ranker (e.g. catalog-matching) if
        you want direction.
    beam_size : int
        Number of top-ranked recipes expanded per iteration in
        priority_queue mode. This is the actual pruning knob. Default
        50; smaller = more aggressive pruning. Ignored in cartesian.
    heap_size : int or None
        Bound on the rank-ordered candidate heap. None = unbounded;
        a number caps memory but may drop low-ranked recipes that
        could matter later. Ignored in cartesian.
    min_product_similarity : float or None
        If set (and target given), adds a Product_Tanimoto_Filter to
        the forward reaction_plan that rejects reactions whose products
        are all below this Tanimoto cutoff to the target. Independent
        of the ranker. Range 0..1; try 0.45 for exploration, 0.7 for
        lead-optimization-style search. Ignored in retro.
    include_chem : bool
        If True, load chem operators (filtered by TAL_REACTION_WHITELIST).
    include_bio : bool
        If True, pre-load all bio cofactors as coreactant helpers and
        load bio operators (filtered by TAL_BIO_REACTION_WHITELIST).
        Bio operators carry a "Comments" meta with semicolon-separated
        UniProt enzyme IDs, so the network records which enzymes are
        compatible with each step.
    bio_ruleset : str
        "JN1224MIN" (1224 rules, default) or "JN3604IMT" (3604 rules).
    bio_allow_multiple_reactants : bool
        Default False — applies the bio Single-Substrate filter, which
        is the regime that keeps non-polyketide bio expansions
        tractable (DORAnet's normal assumption: one variable substrate
        per recipe, cofactors fill the other slots).

        Set to True ONLY when expanding through polyketide chemistry,
        where a single reaction legitimately consumes multiple
        substrates that all happen to be cofactor-like (e.g. the TAL
        biosynthesis route: 1 × acetyl-CoA + 2 × malonyl-CoA → TAL,
        every reactant is a cofactor/CoA-thioester). Without this flag
        the polyketide step gets silently filtered out; with it,
        non-polyketide bio expansions become combinatorially expensive.
    excluded_cofactors : tuple
        Cofactor IDs to exclude from injection. Default mirrors DORAnet's
        ("CARBONYL_CoF", "AMINO_CoF") which over-fire.
    """

    if not starters:
        raise Exception("At least one starter is needed to generate a network")

    # --- depth safety warning ---
    if strategy == "cartesian" and gen > _CARTESIAN_GEN_WARNING_THRESHOLD:
        import warnings
        warnings.warn(
            f"\n[PRUNING WARNING] cartesian strategy with gen={gen} grows "
            f"exponentially and may exhaust local memory.\n"
            f"  → Consider strategy='priority_queue' for gen > "
            f"{_CARTESIAN_GEN_WARNING_THRESHOLD}.\n"
            f"  → Or tighten max_molecular_weight / max_atoms before proceeding.",
            RuntimeWarning,
            stacklevel=2,
        )

    starters = get_smiles_from_file(starters)
    helpers = get_smiles_from_file(helpers)
    targets = get_smiles_from_file(targets)

    print(f"\n{'='*60}")
    print(f"TAL NETWORK GENERATION")
    print(f"{'='*60}")
    print(f"Job name:              {job_name}")
    print(f"Direction:             {direction}")
    print(f"Strategy:              {strategy}")
    print(f"Generations:           {gen}")
    print(f"Max molecular weight:  {max_molecular_weight} Da")
    print(f"Max atoms:             {max_atoms}")
    print(f"Min carbons:           {min_carbons}")
    print(f"Max rxn thermo change: {max_rxn_thermo_change}")
    if strategy == "priority_queue" and targets:
        print(f"Priority queue targets: {targets}")
    print(f"Job started: {datetime.now()}")
    start_time = time.time()

    # ── engine & network setup ──────────────────────────────────────────
    engine = dn.create_engine()
    network = engine.new_network()

    # IMPORTANT pre-loading order so the bundle coreactants filter
    # below treats all helpers (chem) + cofactors (bio) the same way:
    #   1. user-supplied chem helpers
    #   2. bio cofactors (if include_bio)
    #   3. starters
    # Every molecule with index < my_start_i is treated as a coreactant
    # by `engine.filter.bundle.coreactants(tuple(range(my_start_i)))`.
    if helpers:
        for smiles in helpers:
            network.add_mol(engine.mol.rdkit(smiles))

    n_cofactors_loaded = 0
    if include_bio:
        for cof_id, cof_smi in cofactors_clean_dict.items():
            if cof_id in excluded_cofactors:
                continue
            network.add_mol(
                engine.mol.rdkit(cof_smi),
                meta={
                    "SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(cof_smi)),
                },
            )
            n_cofactors_loaded += 1

    my_start_i = -1
    for smiles in starters:
        if my_start_i == -1:
            my_start_i = network.add_mol(
                engine.mol.rdkit(smiles),
                meta={"SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(smiles))},
            )
        else:
            network.add_mol(
                engine.mol.rdkit(smiles),
                meta={"SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(smiles))},
            )

    # ── chem operator loading ──────────────────────────────────────────
    n_chem_ops = 0
    if include_chem:
        if direction == "forward":
            smarts_list = [
                op for op in op_smarts if op.name in TAL_REACTION_WHITELIST
            ]
        elif direction == "retro":
            # DORAnet retro ops are named the same as their forward
            # counterparts, so the TAL chem whitelist filters both
            # sides symmetrically.
            smarts_list = [
                op for op in op_retro_smarts if op.name in TAL_REACTION_WHITELIST
            ]
        else:
            smarts_list = []

        for smarts in smarts_list:
            meta = {
                "name": smarts.name,
                "reactants_stoi": smarts.reactants_stoi,
                "products_stoi": smarts.products_stoi,
                "enthalpy_correction": smarts.enthalpy_correction,
                "ring_issue": smarts.ring_issue,
                "kekulize_flag": smarts.kekulize_flag,
                "Retro_Not_Aromatic": smarts.Retro_Not_Aromatic,
                "number_of_steps": smarts.number_of_steps,
                "allowed_elements": smarts.allowed_elements,
                "Reaction_type": smarts.reaction_type,
                "Reaction_direction": direction,
                "is_bio": False,
            }
            if smarts.kekulize_flag is False:
                network.add_op(
                    engine.op.rdkit(smarts.smarts, drop_errors=True),
                    meta=meta,
                )
            elif smarts.kekulize_flag is True:
                network.add_op(
                    engine.op.rdkit(
                        smarts.smarts, kekulize=True, drop_errors=True
                    ),
                    meta=meta,
                )
            n_chem_ops += 1

    # ── bio operator loading ───────────────────────────────────────────
    n_bio_ops = 0
    if include_bio:
        n_bio_ops, n_bio_whitelisted = _load_bio_ops_into_network(
            network,
            engine,
            ruleset=bio_ruleset,
            direction=direction,
            excluded_cofactors=set(excluded_cofactors),
            whitelist=TAL_BIO_REACTION_WHITELIST,
        )
        print(
            f"Bio: ruleset={bio_ruleset}, cofactors loaded={n_cofactors_loaded}, "
            f"whitelist={n_bio_whitelisted}, ops added={n_bio_ops}"
        )

    print(
        f"Operators loaded: chem={n_chem_ops}, bio={n_bio_ops}, "
        f"total={len(network.ops)}"
    )

    # ── strategy setup ──────────────────────────────────────────────────
    # NOTE: engine.strat.pq is PriorityQueueStrategyBasic (verified in
    # doranet/interfaces.py:1525 and engine.py:200). Its constructor
    # takes (network, num_procs) — the ranker is passed to .expand(),
    # not the constructor. There is no built-in `smiles_distance`
    # ranker; we provide our own in recipe_rankers.py.
    if strategy == "cartesian":
        strat = engine.strat.cartesian(network)

    elif strategy == "priority_queue":
        strat = engine.strat.pq(network)
        if recipe_ranker is None and direction == "forward" and targets:
            target_for_ranker = (
                targets[0] if isinstance(targets, list) else targets
            )
            recipe_ranker = ForwardProductTanimotoRanker(target_for_ranker)
            print(
                f"Priority queue: auto-built ForwardProductTanimotoRanker "
                f"targeting {target_for_ranker}"
            )
        elif recipe_ranker is not None:
            print(f"Priority queue: using user-supplied {type(recipe_ranker).__name__}")
        else:
            print(
                "Priority queue: no ranker (default heap order). "
                "Beam pruning still active via beam_size."
            )
        print(f"  beam_size={beam_size}, heap_size={heap_size}")

    else:
        raise ValueError(
            f"Unknown strategy '{strategy}'. Use 'cartesian' or 'priority_queue'."
        )

    # ── atom limit preprocessing ────────────────────────────────────────
    periodic_table = Chem.GetPeriodicTable()
    if max_atoms is None:
        max_atoms_dict_num = None
    else:
        max_atoms_dict_num = {
            periodic_table.GetAtomicNumber(k): v for k, v in max_atoms.items()
        }

    # ── filter pipelines ────────────────────────────────────────────────
    #
    # PIPELINE ORDER PRINCIPLE:
    #   Cheapest filters first → most candidates rejected before reaching
    #   expensive filters → fewer calls to Chem_Rxn_dH_Calculator overall.
    #
    # FORWARD pipeline:
    #   1. Molecular_Weight_Filter   — single descriptor call, cheapest
    #   2. Allowed_Elements_Filter   — simple atomic number check
    #   3. Max_Atoms_Filter          — atom count iteration
    #   4. Ring_Issues_Filter        — structural check
    #   5. Enol_filter_forward       — structural check
    #   6. Check_balance_filter      — stoichiometry check
    #   7. Minimum_Carbon_Count_Filter — fast via rdMolDescriptors (if used)
    #   8. Chem_Rxn_dH_Calculator    — expensive, correctly last
    #   9. Rxn_dH_Filter             — gates on dH result, correctly last
    #
    # RETRO pipeline:
    #   Order unchanged from original — already logical for retro direction.
    #   Retro_Not_Aromatic_Filter is retro-specific and correctly placed early.

    if direction == "forward":
        reaction_plan = (
            Molecular_Weight_Filter(max_molecular_weight)
            >> Allowed_Elements_Filter()
            >> Max_Atoms_Filter(max_atoms_dict_num)
            >> Ring_Issues_Filter()
            >> Enol_filter_forward()
            >> Check_balance_filter()
        )

        if min_carbons > 0:
            reaction_plan = reaction_plan >> Minimum_Carbon_Count_Filter(min_carbons)

        # Product similarity filter — optional hard cutoff on product
        # Tanimoto vs target. Placed AFTER cheap structural filters and
        # BEFORE the expensive thermo calc.
        if min_product_similarity is not None and targets:
            target_for_filter = (
                targets[0] if isinstance(targets, list) else targets
            )
            reaction_plan = reaction_plan >> Product_Tanimoto_Filter(
                target_for_filter, min_product_similarity
            )

        # TODO: add cost-based branch pruning here once cost_calculator is ready
        # e.g.: reaction_plan = reaction_plan >> Cost_Pruning_Filter(cost_calculator)

        reaction_plan = (
            reaction_plan
            >> Chem_Rxn_dH_Calculator("dH", "forward", molecule_thermo_calculator)
            >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
        )

        recipe_filter = None

    elif direction == "retro":
        reaction_plan = (
            Max_Atoms_Filter(max_atoms_dict_num)
            >> Molecular_Weight_Filter(max_molecular_weight)
            >> Ring_Issues_Filter()
            >> Retro_Not_Aromatic_Filter()
            >> Enol_filter_retro()
            >> Allowed_Elements_Filter()
            >> Check_balance_filter()
        )

        if min_carbons > 0:
            reaction_plan = reaction_plan >> Minimum_Carbon_Count_Filter(min_carbons)

        # TODO: add cost-based branch pruning here once cost_calculator is ready

        reaction_plan = (
            reaction_plan
            >> Chem_Rxn_dH_Calculator("dH", "retro", molecule_thermo_calculator)
            >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
        )

        recipe_filter = Cross_Reaction_Filter(tuple(range(my_start_i)))

    # ── allow_multiple_reactants override ───────────────────────────────
    if allow_multiple_reactants != "default":
        if allow_multiple_reactants is True:
            recipe_filter = None
        elif allow_multiple_reactants is False:
            recipe_filter = Cross_Reaction_Filter(tuple(range(my_start_i)))

    # ── bio-specific recipe filter ──────────────────────────────────────
    # Only stacked when bio is enabled AND we don't want the polyketide-
    # style cofactor+cofactor combos. The filter is conditional on
    # is_bio meta so it leaves chem recipes alone.
    if include_bio and not bio_allow_multiple_reactants:
        bio_filter = Bio_Single_Substrate_Filter(
            cofactor_smiles_set=frozenset(cofactors_clean)
        )
        if recipe_filter is None:
            recipe_filter = bio_filter
        else:
            recipe_filter = recipe_filter & bio_filter

    # ── expansion ───────────────────────────────────────────────────────
    bundle_filter = engine.filter.bundle.coreactants(tuple(range(my_start_i)))
    ini_number = len(network.mols)

    # Cartesian's expand() accepts num_iter directly (it wraps it in a
    # max_iter hook internally). PriorityQueueStrategyBasic's expand()
    # has no num_iter parameter; we must pass the iter cap via
    # global_hooks=[engine.hook.max_iter(gen)], and supply ranker,
    # beam_size, heap_size there too.
    if strategy == "cartesian":
        strat.expand(
            num_iter=gen,
            reaction_plan=reaction_plan,
            bundle_filter=bundle_filter,
            recipe_filter=recipe_filter,
            save_unreactive=False,
        )
    else:
        strat.expand(
            reaction_plan=reaction_plan,
            bundle_filter=bundle_filter,
            recipe_filter=recipe_filter,
            recipe_ranker=recipe_ranker,
            beam_size=beam_size,
            heap_size=heap_size,
            global_hooks=[engine.hook.max_iter(gen)],
            save_unreactive=False,
        )

    # ── target detection ────────────────────────────────────────────────
    if targets is not None:
        print("\nChecking for targets...")
        to_check = set()
        if isinstance(targets, str):
            to_check.add(Chem.MolToSmiles(Chem.MolFromSmiles(targets)))
        else:
            for t in targets:
                to_check.add(Chem.MolToSmiles(Chem.MolFromSmiles(t)))

        targets_found = []
        for mol in network.mols:
            if (
                network.reactivity[network.mols.i(mol.uid)] is True
                and mol.uid in to_check
            ):
                targets_found.append(mol.uid)
                print(f"  ✓ Target found: {mol.uid}")

        if not targets_found:
            print("  ✗ No targets found")

    # ── summary ─────────────────────────────────────────────────────────
    elapsed = (time.time() - start_time) / 60
    print(f"\n{'='*60}")
    print(f"NETWORK GENERATION COMPLETE")
    print(f"{'='*60}")
    print(f"Generations:                  {gen}")
    print(f"Operators:                    {len(network.ops)}")
    print(f"Molecules before expansion:   {ini_number}")
    print(f"Molecules after expansion:    {len(network.mols)}")
    print(f"Reactions:                    {len(network.rxns)}")
    print(f"Time elapsed:                 {elapsed:.2f} minutes")
    print(f"{'='*60}\n")

    network.save_to_file(f"{job_name}_{direction}_saved_network")

    return network








# #WAIT actually im going to use the cosmetic expansion from yesterday to modify instead of starting from scratch 
# #deleting all the unnecessary stuff - no more functional group or alkene preservation filters --> these things can be used in post reaction processing if needed


# import doranet as dn
# from doranet.modules.synthetic.Reaction_Smarts_Forward import op_smarts
# from doranet.modules.synthetic.Reaction_Smarts_Retro import op_retro_smarts
# from doranet.modules.synthetic.generate_network import (
#     get_smiles_from_file,
#     Max_Atoms_Filter,
#     Ring_Issues_Filter,
#     Enol_filter_forward,
#     Enol_filter_retro,
#     Check_balance_filter,
#     Allowed_Elements_Filter,
#     Chem_Rxn_dH_Calculator,
#     Rxn_dH_Filter,
#     Cross_Reaction_Filter,
#     Retro_Not_Aromatic_Filter,
# )
# from datetime import datetime
# import time
# from rdkit import Chem
# from rdkit.Chem import Descriptors
# from doranet import metadata, interfaces
# import collections.abc
# from tal_reaction_whitelist import TAL_REACTION_WHITELIST #THIS IS USED FOR CUSTOM WHITELIST


# # ============ CUSTOM FILTERS ============

# class Molecular_Weight_Filter(metadata.ReactionFilterBase):
#     """Reject reactions if products exceed max molecular weight."""
    
#     def __init__(self, max_weight=400):
#         self.max_weight = max_weight
    
#     def __call__(self, recipe):
#         for mol in recipe.products:
#             if not isinstance(mol.item, interfaces.MolDatRDKit):
#                 raise NotImplementedError(
#                     f"Filter only works with MolDatRDKit, not {type(mol.item)}"
#                 )
            
#             mw = Descriptors.MolWt(mol.item.rdkitmol)
            
#             if mw > self.max_weight:
#                 return False
        
#         return True
    
#     @property
#     def meta_required(self):
#         return interfaces.MetaKeyPacket()


# class Minimum_Carbon_Count_Filter(metadata.ReactionFilterBase):
#     """Reject molecules with fewer than 8 carbons."""
    
#     def __init__(self, min_carbons=0):
#         self.min_carbons = min_carbons
    
#     def __call__(self, recipe):
#         for mol in recipe.products:
#             if not isinstance(mol.item, interfaces.MolDatRDKit):
#                 raise NotImplementedError(
#                     f"Filter only works with MolDatRDKit, not {type(mol.item)}"
#                 )
            
#             # Count carbons
#             carbon_count = sum(
#                 1 for atom in mol.item.rdkitmol.GetAtoms() 
#                 if atom.GetAtomicNum() == 6
#             )
            
#             if carbon_count < self.min_carbons:
#                 return False
        
#         return True
    
#     @property
#     def meta_required(self):
#         return interfaces.MetaKeyPacket()


# # ============ CUSTOM GENERATE FUNCTION ============

# def generate_network_tal(
#     job_name="default_job",
#     starters=False,
#     helpers=False,
#     gen=2, 
#     direction="forward",
#     molecule_thermo_calculator=None,
#     max_rxn_thermo_change=15,
#     max_atoms=None,  # {"C": 50, "O": 8, "N": 2}
#     max_molecular_weight=800,  # NEW: Adjustable MW cap
#     allow_multiple_reactants="default",
#     targets=None,
#     strategy="cartesian",  # NEW: "cartesian" or "priority_queue"
#     #preserve_delta9=True,  # NEW: Preserve Δ9 alkene
#     min_carbons=0,  # NEW: Minimum carbon count
#     #enforce_functional_groups=True,  # NEW: Whitelist functional groups
# ):
#     """
#     Enhanced generate_network for TAL acid derivatives.
    
#     Parameters:
#     -----------
#     max_atoms : dict
#         Atom count limits, e.g., {"C": 50, "O": 8, "N": 2}
    
#     max_molecular_weight : float
#         Max molecular weight in Daltons (default: 800)
    
#     strategy : str
#         "cartesian" (blind, exhaustive) or "priority_queue" (targeted)
#         If "priority_queue", auto-detects target from targets parameter
    
#     min_carbons : int
#         Minimum carbon count in molecules (default: 8)
    
#     """
    
#     if not starters:
#         raise Exception("At least one starter is needed to generate a network")

#     starters = get_smiles_from_file(starters)
#     helpers = get_smiles_from_file(helpers)
#     targets = get_smiles_from_file(targets)

#     print(f"\n{'='*60}")
#     print(f"TAL NETWORK GENERATION")
#     print(f"{'='*60}")
#     print(f"Job name: {job_name}")
#     print(f"Job type: synthetic network expansion {direction}")
#     print(f"Strategy: {strategy}")
#     if strategy == "priority_queue" and targets:
#         print(f"Priority queue targets: {targets}")
#     print(f"Atom limits: {max_atoms}")
#     print(f"Max molecular weight: {max_molecular_weight} Da")
#     #print(f"Preserve Δ9 alkene: {preserve_delta9}")
#     print(f"Min carbons: {min_carbons}")
#     #print(f"Enforce functional groups: {enforce_functional_groups}")
#     print("Job started on:", datetime.now())
#     start_time = time.time()

#     engine = dn.create_engine()
#     network = engine.new_network()

#     if helpers:
#         for smiles in helpers:
#             network.add_mol(engine.mol.rdkit(smiles))

#     my_start_i = -1
#     for smiles in starters:
#         if my_start_i == -1:
#             my_start_i = network.add_mol(engine.mol.rdkit(smiles))
#         else:
#             network.add_mol(engine.mol.rdkit(smiles))

#     if direction == "forward":
#         smarts_list = [op for op in op_smarts if op.name in TAL_REACTION_WHITELIST]
#         #print("SMARTS LIST")
#         #for smarts in smarts_list:
#            # print(smarts.name)
            
#         #TS IS NEW 
#     elif direction == "retro":
#         smarts_list = op_retro_smarts

#     for smarts in smarts_list:
#         if smarts.kekulize_flag is False:
#             network.add_op(
#                 engine.op.rdkit(smarts.smarts, drop_errors=True),
#                 meta={
#                     "name": smarts.name,
#                     "reactants_stoi": smarts.reactants_stoi,
#                     "products_stoi": smarts.products_stoi,
#                     "enthalpy_correction": smarts.enthalpy_correction,
#                     "ring_issue": smarts.ring_issue,
#                     "kekulize_flag": smarts.kekulize_flag,
#                     "Retro_Not_Aromatic": smarts.Retro_Not_Aromatic,
#                     "number_of_steps": smarts.number_of_steps,
#                     "allowed_elements": smarts.allowed_elements,
#                     "Reaction_type": smarts.reaction_type,
#                     "Reaction_direction": direction,
#                 },
#             )
#         if smarts.kekulize_flag is True:
#             network.add_op(
#                 engine.op.rdkit(smarts.smarts, kekulize=True, drop_errors=True),
#                 meta={
#                     "name": smarts.name,
#                     "reactants_stoi": smarts.reactants_stoi,
#                     "products_stoi": smarts.products_stoi,
#                     "enthalpy_correction": smarts.enthalpy_correction,
#                     "ring_issue": smarts.ring_issue,
#                     "kekulize_flag": smarts.kekulize_flag,
#                     "Retro_Not_Aromatic": smarts.Retro_Not_Aromatic,
#                     "number_of_steps": smarts.number_of_steps,
#                     "allowed_elements": smarts.allowed_elements,
#                     "Reaction_type": smarts.reaction_type,
#                     "Reaction_direction": direction,
#                 },
#             )

#     print(f"Number of operators added to network: {len(network.ops)}")

#     # ====== CHOOSE STRATEGY ======
#     if strategy == "cartesian":
#         strat = engine.strat.cartesian(network)
#     elif strategy == "priority_queue":
#         if not targets:
#             raise ValueError(
#                 "priority_queue strategy requires targets parameter to be set"
#             )
#         # Use first target as the ranker target
#         target_for_ranker = targets[0] if isinstance(targets, list) else targets
#         ranker = engine.ranker.smiles_distance(target_for_ranker)
#         strat = engine.strat.priority_queue(network, ranker=ranker)
#         print(f"Using priority queue strategy, targeting: {target_for_ranker}")
#     else:
#         raise ValueError(
#             f"Unknown strategy: {strategy}. Use 'cartesian' or 'priority_queue'"
#         )
#     # =============================

#     periodic_table = Chem.GetPeriodicTable()

#     if max_atoms is None:
#         max_atoms_dict_num = None
#     else:
#         max_atoms_dict_num = dict()
#         for key in max_atoms:
#             max_atoms_dict_num[periodic_table.GetAtomicNumber(key)] = max_atoms[key]

#     # Build reaction plan with custom filters
#     if direction == "forward":
#         reaction_plan = (
#             Max_Atoms_Filter(max_atoms_dict_num)
#             >> Molecular_Weight_Filter(max_molecular_weight)
#             >> Ring_Issues_Filter()
#             >> Enol_filter_forward()
#             >> Check_balance_filter()
#             >> Allowed_Elements_Filter()
#         )
        
#         # Add optional filters
#        #if preserve_delta9:
#           #  reaction_plan = reaction_plan >> Preserve_Delta9_Alkene_Filter()
        
#         if min_carbons > 0:
#             reaction_plan = reaction_plan >> Minimum_Carbon_Count_Filter(min_carbons)
        
#         #  enforce_functional_groups:
#             #reaction_plan = reaction_plan >> Functional_Group_Whitelist_Filter()
        
#         # Add thermodynamics filters at the end
#         reaction_plan = (
#             reaction_plan
#             >> Chem_Rxn_dH_Calculator("dH", "forward", molecule_thermo_calculator)
#             >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
#         )
#         #HOLD UP - is this from original or is this something i added - I DO NOT REMEMBER ADDING THIS
#         #Should thermodynamics be filter or analysis --> actually im 80% sure this is from original so i wont touch
        
#         recipe_filter = None

#     elif direction == "retro":
#         reaction_plan = (
#             Max_Atoms_Filter(max_atoms_dict_num)
#             >> Molecular_Weight_Filter(max_molecular_weight)
#             >> Ring_Issues_Filter()
#             >> Retro_Not_Aromatic_Filter()
#             >> Enol_filter_retro()
#             >> Allowed_Elements_Filter()
#             >> Check_balance_filter()
#         )
        
#         # Add optional filters
#         #if preserve_delta9:
#          #   reaction_plan = reaction_plan >> Preserve_Delta9_Alkene_Filter()"""
        
#         if min_carbons > 0:
#             reaction_plan = reaction_plan >> Minimum_Carbon_Count_Filter(min_carbons)
        
#         #if enforce_functional_groups:
#          #   reaction_plan = reaction_plan >> Functional_Group_Whitelist_Filter()"""
        
#         # Add thermodynamics filters at the end
#         reaction_plan = (
#             reaction_plan
#             >> Chem_Rxn_dH_Calculator("dH", "retro", molecule_thermo_calculator)
#             >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
#         )
        
#         recipe_filter = Cross_Reaction_Filter(tuple(range(my_start_i)))

#     if allow_multiple_reactants != "default":
#         if allow_multiple_reactants is True:
#             recipe_filter = None
#         elif allow_multiple_reactants is False:
#             recipe_filter = Cross_Reaction_Filter(tuple(range(my_start_i)))

#     bundle_filter = engine.filter.bundle.coreactants(tuple(range(my_start_i)))

#     ini_number = len(network.mols)

#     strat.expand(
#         num_iter=gen,
#         reaction_plan=reaction_plan,
#         bundle_filter=bundle_filter,
#         recipe_filter=recipe_filter,
#         save_unreactive=False,
#     )

#     if targets is not None:
#         print("\nChecking for targets...")
#         to_check = set()
#         if isinstance(targets, str):
#             to_check.add(Chem.MolToSmiles(Chem.MolFromSmiles(targets)))
#         else:
#             for i in targets:
#                 to_check.add(Chem.MolToSmiles(Chem.MolFromSmiles(i)))

#         targets_found = []
#         for mol in network.mols:
#             if (
#                 network.reactivity[network.mols.i(mol.uid)] is True
#                 and mol.uid in to_check
#             ):
#                 targets_found.append(mol.uid)
#                 print(f"  ✓ Target found: {mol.uid}")
        
#         if len(targets_found) == 0:
#             print(f"  ✗ No targets found")

#     print("\n" + "="*60)
#     print("NETWORK GENERATION COMPLETE")
#     print("="*60)
#     print(f"Number of generations: {gen}")
#     print(f"Number of operators: {len(network.ops)}")
#     print(f"Number of molecules before expansion: {ini_number}")
#     print(f"Number of molecules after expansion: {len(network.mols)}")
#     print(f"Number of reactions: {len(network.rxns)}")

#     end_time = time.time()
#     elapsed_time = (end_time - start_time) / 60
#     print(f"Time used: {elapsed_time:.2f} minutes")
#     print("="*60 + "\n")

#     network.save_to_file(f"{job_name}_{direction}_saved_network")

#     return network