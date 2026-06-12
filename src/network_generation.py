#WAIT actually im going to use the cosmetic expansion from yesterday to modify instead of starting from scratch 
#deleting all the unnecessary stuff - no more functional group or alkene preservation filters --> these things can be used in post reaction processing if needed


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
from datetime import datetime
import time
from rdkit import Chem
from rdkit.Chem import Descriptors
from doranet import metadata, interfaces
import collections.abc
from tal_reaction_whitelist import TAL_REACTION_WHITELIST #THIS IS USED FOR CUSTOM WHITELIST


# ============ CUSTOM FILTERS ============

class Molecular_Weight_Filter(metadata.ReactionFilterBase):
    """Reject reactions if products exceed max molecular weight."""
    
    def __init__(self, max_weight=400):
        self.max_weight = max_weight
    
    def __call__(self, recipe):
        for mol in recipe.products:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"Filter only works with MolDatRDKit, not {type(mol.item)}"
                )
            
            mw = Descriptors.MolWt(mol.item.rdkitmol)
            
            if mw > self.max_weight:
                return False
        
        return True
    
    @property
    def meta_required(self):
        return interfaces.MetaKeyPacket()


class Minimum_Carbon_Count_Filter(metadata.ReactionFilterBase):
    """Reject molecules with fewer than 8 carbons."""
    
    def __init__(self, min_carbons=0):
        self.min_carbons = min_carbons
    
    def __call__(self, recipe):
        for mol in recipe.products:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"Filter only works with MolDatRDKit, not {type(mol.item)}"
                )
            
            # Count carbons
            carbon_count = sum(
                1 for atom in mol.item.rdkitmol.GetAtoms() 
                if atom.GetAtomicNum() == 6
            )
            
            if carbon_count < self.min_carbons:
                return False
        
        return True
    
    @property
    def meta_required(self):
        return interfaces.MetaKeyPacket()


# ============ CUSTOM GENERATE FUNCTION ============

def generate_network_tal(
    job_name="default_job",
    starters=False,
    helpers=False,
    gen=2, 
    direction="forward",
    molecule_thermo_calculator=None,
    max_rxn_thermo_change=15,
    max_atoms=None,  # {"C": 50, "O": 8, "N": 2}
    max_molecular_weight=800,  # NEW: Adjustable MW cap
    allow_multiple_reactants="default",
    targets=None,
    strategy="cartesian",  # NEW: "cartesian" or "priority_queue"
    #preserve_delta9=True,  # NEW: Preserve Δ9 alkene
    min_carbons=0,  # NEW: Minimum carbon count
    #enforce_functional_groups=True,  # NEW: Whitelist functional groups
):
    """
    Enhanced generate_network for TAL acid derivatives.
    
    Parameters:
    -----------
    max_atoms : dict
        Atom count limits, e.g., {"C": 50, "O": 8, "N": 2}
    
    max_molecular_weight : float
        Max molecular weight in Daltons (default: 800)
    
    strategy : str
        "cartesian" (blind, exhaustive) or "priority_queue" (targeted)
        If "priority_queue", auto-detects target from targets parameter
    
    min_carbons : int
        Minimum carbon count in molecules (default: 8)
    
    """
    
    if not starters:
        raise Exception("At least one starter is needed to generate a network")

    starters = get_smiles_from_file(starters)
    helpers = get_smiles_from_file(helpers)
    targets = get_smiles_from_file(targets)

    print(f"\n{'='*60}")
    print(f"TAL NETWORK GENERATION")
    print(f"{'='*60}")
    print(f"Job name: {job_name}")
    print(f"Job type: synthetic network expansion {direction}")
    print(f"Strategy: {strategy}")
    if strategy == "priority_queue" and targets:
        print(f"Priority queue targets: {targets}")
    print(f"Atom limits: {max_atoms}")
    print(f"Max molecular weight: {max_molecular_weight} Da")
    #print(f"Preserve Δ9 alkene: {preserve_delta9}")
    print(f"Min carbons: {min_carbons}")
    #print(f"Enforce functional groups: {enforce_functional_groups}")
    print("Job started on:", datetime.now())
    start_time = time.time()

    engine = dn.create_engine()
    network = engine.new_network()

    if helpers:
        for smiles in helpers:
            network.add_mol(engine.mol.rdkit(smiles))

    my_start_i = -1
    for smiles in starters:
        if my_start_i == -1:
            my_start_i = network.add_mol(engine.mol.rdkit(smiles))
        else:
            network.add_mol(engine.mol.rdkit(smiles))

    if direction == "forward":
        smarts_list = [op for op in op_smarts if op.name in TAL_REACTION_WHITELIST]
        #print("SMARTS LIST")
        #for smarts in smarts_list:
           # print(smarts.name)
            
        #TS IS NEW 
    elif direction == "retro":
        smarts_list = op_retro_smarts

    for smarts in smarts_list:
        if smarts.kekulize_flag is False:
            network.add_op(
                engine.op.rdkit(smarts.smarts, drop_errors=True),
                meta={
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
                },
            )
        if smarts.kekulize_flag is True:
            network.add_op(
                engine.op.rdkit(smarts.smarts, kekulize=True, drop_errors=True),
                meta={
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
                },
            )

    print(f"Number of operators added to network: {len(network.ops)}")

    # ====== CHOOSE STRATEGY ======
    if strategy == "cartesian":
        strat = engine.strat.cartesian(network)
    elif strategy == "priority_queue":
        if not targets:
            raise ValueError(
                "priority_queue strategy requires targets parameter to be set"
            )
        # Use first target as the ranker target
        target_for_ranker = targets[0] if isinstance(targets, list) else targets
        ranker = engine.ranker.smiles_distance(target_for_ranker)
        strat = engine.strat.priority_queue(network, ranker=ranker)
        print(f"Using priority queue strategy, targeting: {target_for_ranker}")
    else:
        raise ValueError(
            f"Unknown strategy: {strategy}. Use 'cartesian' or 'priority_queue'"
        )
    # =============================

    periodic_table = Chem.GetPeriodicTable()

    if max_atoms is None:
        max_atoms_dict_num = None
    else:
        max_atoms_dict_num = dict()
        for key in max_atoms:
            max_atoms_dict_num[periodic_table.GetAtomicNumber(key)] = max_atoms[key]

    # Build reaction plan with custom filters
    if direction == "forward":
        reaction_plan = (
            Max_Atoms_Filter(max_atoms_dict_num)
            >> Molecular_Weight_Filter(max_molecular_weight)
            >> Ring_Issues_Filter()
            >> Enol_filter_forward()
            >> Check_balance_filter()
            >> Allowed_Elements_Filter()
        )
        
        # Add optional filters
       #if preserve_delta9:
          #  reaction_plan = reaction_plan >> Preserve_Delta9_Alkene_Filter()
        
        if min_carbons > 0:
            reaction_plan = reaction_plan >> Minimum_Carbon_Count_Filter(min_carbons)
        
        #  enforce_functional_groups:
            #reaction_plan = reaction_plan >> Functional_Group_Whitelist_Filter()
        
        # Add thermodynamics filters at the end
        reaction_plan = (
            reaction_plan
            >> Chem_Rxn_dH_Calculator("dH", "forward", molecule_thermo_calculator)
            >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
        )
        #HOLD UP - is this from original or is this something i added - I DO NOT REMEMBER ADDING THIS
        #Should thermodynamics be filter or analysis --> actually im 80% sure this is from original so i wont touch
        
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
        
        # Add optional filters
        #if preserve_delta9:
         #   reaction_plan = reaction_plan >> Preserve_Delta9_Alkene_Filter()"""
        
        if min_carbons > 0:
            reaction_plan = reaction_plan >> Minimum_Carbon_Count_Filter(min_carbons)
        
        #if enforce_functional_groups:
         #   reaction_plan = reaction_plan >> Functional_Group_Whitelist_Filter()"""
        
        # Add thermodynamics filters at the end
        reaction_plan = (
            reaction_plan
            >> Chem_Rxn_dH_Calculator("dH", "retro", molecule_thermo_calculator)
            >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
        )
        
        recipe_filter = Cross_Reaction_Filter(tuple(range(my_start_i)))

    if allow_multiple_reactants != "default":
        if allow_multiple_reactants is True:
            recipe_filter = None
        elif allow_multiple_reactants is False:
            recipe_filter = Cross_Reaction_Filter(tuple(range(my_start_i)))

    bundle_filter = engine.filter.bundle.coreactants(tuple(range(my_start_i)))

    ini_number = len(network.mols)

    strat.expand(
        num_iter=gen,
        reaction_plan=reaction_plan,
        bundle_filter=bundle_filter,
        recipe_filter=recipe_filter,
        save_unreactive=False,
    )

    if targets is not None:
        print("\nChecking for targets...")
        to_check = set()
        if isinstance(targets, str):
            to_check.add(Chem.MolToSmiles(Chem.MolFromSmiles(targets)))
        else:
            for i in targets:
                to_check.add(Chem.MolToSmiles(Chem.MolFromSmiles(i)))

        targets_found = []
        for mol in network.mols:
            if (
                network.reactivity[network.mols.i(mol.uid)] is True
                and mol.uid in to_check
            ):
                targets_found.append(mol.uid)
                print(f"  ✓ Target found: {mol.uid}")
        
        if len(targets_found) == 0:
            print(f"  ✗ No targets found")

    print("\n" + "="*60)
    print("NETWORK GENERATION COMPLETE")
    print("="*60)
    print(f"Number of generations: {gen}")
    print(f"Number of operators: {len(network.ops)}")
    print(f"Number of molecules before expansion: {ini_number}")
    print(f"Number of molecules after expansion: {len(network.mols)}")
    print(f"Number of reactions: {len(network.rxns)}")

    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print(f"Time used: {elapsed_time:.2f} minutes")
    print("="*60 + "\n")

    network.save_to_file(f"{job_name}_{direction}_saved_network")

    return network