"""
build_bio_glucose_to_tal_whitelist.py
Build a STRICT whitelist of JN1224MIN rules that catalyze the
literature glucose -> TAL pathway, step by step.

Difference vs build_bio_whitelist.py
------------------------------------
  - build_bio_whitelist.py: keep a rule if it produces ANY sanitizable
    product on ANY TAL-chemistry seed. Broad. 348 rules pass.
  - this script: keep a rule if applying it to intermediate N produces
    intermediate N+1 EXACTLY (canonical SMILES match). Strict. Expect
    ~15-25 rules.

Goal
----
Pare 1224 rules down to the minimum set needed for the literature
glycolysis + polyketide-biosynthesis pathway. Lets the retro-bio
search validate end-to-end on glucose <- TAL without the cofactor
combinatorial explosion that 348 rules cause.

Method
------
1. Define the literature intermediate sequence (glucose -> ... -> TAL).
2. For each consecutive (current, next) pair:
     - For each rule in JN1224MIN, try putting `current` in any Any
       slot with appropriate cofactors. Check if any product set
       contains `next` (canonical SMILES match).
     - If a rule fires the transformation, record its name AND which
       step it catalyzes.
3. Save:
     - frozenset of unique rule names -> TAL_BIO_GLUCOSE_TO_TAL_WHITELIST
     - dict {rule_name -> step_label} for diagnostics

Run from project root:
    python scripts/build_bio_glucose_to_tal_whitelist.py
"""

import csv
import os
import sys
import time

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")


# --- paths --------------------------------------------------------------
DORANET_BIO = (
    "C:/Users/ashvi/anaconda3/Lib/site-packages/doranet/modules/enzymatic"
)
RULES_TSV = os.path.join(DORANET_BIO, "JN1224MIN_rules.tsv")
COFACTORS_TSV = os.path.join(DORANET_BIO, "all_cofactors.tsv")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_PY = os.path.normpath(
    os.path.join(HERE, "..", "src", "tal_bio_glucose_to_tal_whitelist.py")
)


# --- literature pathway -------------------------------------------------
# Each tuple: (label, SMILES). The sequence walks forward from glucose
# down glycolysis, then via acetyl-CoA carboxylation and polyketide
# condensation to TAL.
#
# Phosphorylated intermediates are written as the deprotonated forms
# with explicit -OPO(=O)(O)O groups. We rely on canonical-SMILES
# round-tripping to match the bio rule outputs.
PATHWAY = [
    ("glucose",                 "OCC(O)C(O)C(O)C(O)C=O"),
    ("glucose-6-phosphate",     "O=CC(O)C(O)C(O)C(O)COP(=O)(O)O"),
    ("fructose-6-phosphate",    "O=C(CO)C(O)C(O)C(O)COP(=O)(O)O"),
    ("fructose-1,6-bisphosphate","O=C(COP(=O)(O)O)C(O)C(O)C(O)COP(=O)(O)O"),
    ("dihydroxyacetone phosphate","OCC(=O)COP(=O)(O)O"),
    ("glyceraldehyde-3-phosphate","O=CC(O)COP(=O)(O)O"),
    ("1,3-bisphosphoglycerate", "O=C(OP(=O)(O)O)C(O)COP(=O)(O)O"),
    ("3-phosphoglycerate",      "O=C(O)C(O)COP(=O)(O)O"),
    ("2-phosphoglycerate",      "O=C(O)C(OP(=O)(O)O)CO"),
    ("phosphoenolpyruvate",     "O=C(O)C(=C)OP(=O)(O)O"),
    ("pyruvate",                "CC(=O)C(=O)O"),
    # Acetyl-CoA and Malonyl-CoA use DORAnet's cofactor SMILES (looked
    # up at runtime to guarantee canonical form matches the cofactor
    # pool).
    ("acetyl-CoA",              None),       # filled from cofactors
    ("malonyl-CoA",             None),       # filled from constant
    # Polyketide intermediates en route to TAL.
    ("acetoacetyl-CoA",         None),       # filled from constant
    ("3,5-dioxohexanoyl-CoA",   None),       # filled from constant
    ("TAL",                     "Cc1cc(O)cc(=O)o1"),
]

# Cofactor IDs to look up in all_cofactors.tsv for the special slots
ACETYL_COA_ID  = "ACETYL-COA"     # DORAnet uses ACETYL-COA (with hyphens)

# Malonyl-CoA is NOT registered as a cofactor in DORAnet's all_cofactors.tsv
# — it's a regular Any-slot substrate. Built by hand by carboxylating the
# acetyl group of acetyl-CoA: -CC(=O)S... becomes OC(=O)CC(=O)S...
CoA_TAIL = (
    "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)
MALONYL_COA_SMILES   = "OC(=O)CC(=O)"        + CoA_TAIL
# 2-PS (2-pyrone synthase) builds the polyketide via two
# decarboxylative Claisen condensations:
#   acetyl-CoA + malonyl-CoA  -> acetoacetyl-CoA + CO2 + CoA
#   acetoacetyl-CoA + malonyl-CoA -> 3,5-dioxohexanoyl-CoA + CO2 + CoA
# Then 3,5-dioxohexanoyl-CoA cyclizes intramolecularly to TAL + CoA.
ACETOACETYL_COA_SMILES        = "CC(=O)CC(=O)"           + CoA_TAIL
DIOXOHEXANOYL_COA_SMILES      = "CC(=O)CC(=O)CC(=O)"     + CoA_TAIL


def load_cofactors():
    cof = {}
    with open(COFACTORS_TSV, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            cid, _name, smi = row[0], row[1], row[2]
            cof[cid] = smi
    return cof


def load_rules():
    with open(RULES_TSV, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield row["Name"], row["Reactants"].split(";"), row["SMARTS"]


def canon(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return Chem.MolToSmiles(m)


def resolve_pathway(cofactors):
    """Replace SMILES=None entries with cofactor lookups; canonicalize all."""
    resolved = []
    for label, smi in PATHWAY:
        if smi is None:
            if label == "acetyl-CoA":
                smi = cofactors.get(ACETYL_COA_ID)
            elif label == "malonyl-CoA":
                smi = MALONYL_COA_SMILES
            elif label == "acetoacetyl-CoA":
                smi = ACETOACETYL_COA_SMILES
            elif label == "3,5-dioxohexanoyl-CoA":
                smi = DIOXOHEXANOYL_COA_SMILES
            if smi is None:
                raise RuntimeError(f"SMILES for {label} not found")
        c = canon(smi)
        if c is None:
            raise RuntimeError(f"SMILES for {label} failed to parse: {smi}")
        resolved.append((label, c))
    return resolved


def build_reactant_tuples(reactant_types, current_mol, cof_mols, helper_mols):
    """
    Yield candidate reactant tuples for testing this rule. `current_mol`
    is the pathway intermediate at this step. cof_mols is the {id: Mol}
    dict for cofactor slots. helper_mols is a small set of common
    fillers (water, the next intermediate's neighbors) we try at extra
    Any slots when there's more than one.
    """
    any_indices = [i for i, t in enumerate(reactant_types) if t == "Any"]
    cof_indices = [i for i, t in enumerate(reactant_types) if t != "Any"]

    # Map cofactor slots to their molecules; abort if any unknown
    base = [None] * len(reactant_types)
    for i in cof_indices:
        cid = reactant_types[i]
        if cid not in cof_mols:
            return
        base[i] = cof_mols[cid]

    if not any_indices:
        yield tuple(base)
        return

    if len(any_indices) == 1:
        i = any_indices[0]
        t = list(base)
        t[i] = current_mol
        yield tuple(t)
        return

    # 2+ Any slots: put current_mol at each Any slot in turn, fill the
    # rest with each helper mol.
    for varying in any_indices:
        for fill in helper_mols + [current_mol]:
            t = list(base)
            for j in any_indices:
                t[j] = current_mol if j == varying else fill
            yield tuple(t)


def rule_produces_target(rxn, reactant_tuples, target_canon):
    """Return True if any reactant tuple produces target_canon."""
    for reactants in reactant_tuples:
        try:
            products = rxn.RunReactants(reactants)
        except Exception:
            continue
        if not products:
            continue
        for product_set in products:
            for prod in product_set:
                try:
                    Chem.SanitizeMol(prod)
                    smi = Chem.MolToSmiles(prod)
                except Exception:
                    continue
                if smi == target_canon:
                    return True
    return False


def main():
    print(f"Loading cofactors from {COFACTORS_TSV}")
    cofactors = load_cofactors()
    print(f"  loaded {len(cofactors)} cofactor SMILES")

    # Resolve pathway (fills cofactor SMILES, canonicalizes everything)
    print("\nResolving literature pathway...")
    pathway = resolve_pathway(cofactors)
    for label, smi in pathway:
        print(f"  {label:30s} {smi[:50]}{'...' if len(smi)>50 else ''}")

    # Pre-parse cofactor + helper molecules
    cof_mols = {}
    for cid, smi in cofactors.items():
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            cof_mols[cid] = m
    helper_mols = [Chem.MolFromSmiles(s) for s in ("O", "[H][H]")]
    helper_mols = [m for m in helper_mols if m is not None]

    print(f"\nLoading rules from {RULES_TSV}")
    rules = list(load_rules())
    print(f"  loaded {len(rules)} rules")

    print("\nTesting each (step) x (rule) pair...")
    t0 = time.perf_counter()

    # rule_name -> list of step_labels it catalyzes
    rule_to_steps: dict[str, list[str]] = {}

    for step_idx in range(len(pathway) - 1):
        cur_label, cur_smi = pathway[step_idx]
        next_label, next_smi = pathway[step_idx + 1]
        step_label = f"{cur_label} -> {next_label}"
        cur_mol = Chem.MolFromSmiles(cur_smi)
        if cur_mol is None:
            print(f"  [skip] {step_label}: bad SMILES")
            continue

        n_match = 0
        for rule_name, reactant_types, smarts in rules:
            rxn = AllChem.ReactionFromSmarts(smarts)
            if rxn is None:
                continue
            tuples = list(
                build_reactant_tuples(
                    reactant_types, cur_mol, cof_mols, helper_mols
                )
            )
            if not tuples:
                continue
            if rule_produces_target(rxn, tuples, next_smi):
                rule_to_steps.setdefault(rule_name, []).append(step_label)
                n_match += 1

        print(f"  step {step_idx+1:2d}: {step_label:55s}  matches={n_match}")

    elapsed = time.perf_counter() - t0
    print(f"\nFinished in {elapsed:.1f}s")
    print(f"  unique rules kept: {len(rule_to_steps)}")

    write_whitelist(rule_to_steps, pathway, len(rules))
    print(f"\nWrote {OUT_PY}")


def write_whitelist(rule_to_steps, pathway, total_rules):
    rule_lines = []
    for name in sorted(rule_to_steps.keys()):
        steps = rule_to_steps[name]
        rule_lines.append(f'    "{name}",  # ' + "; ".join(steps))
    body_lines = "\n".join(rule_lines)

    pathway_lines = "\n".join(
        f"#   {i+1:2d}. {label:30s} {smi[:55]}{'...' if len(smi)>55 else ''}"
        for i, (label, smi) in enumerate(pathway)
    )
    content = f'''"""
TAL Bio GLUCOSE -> TAL Pathway Whitelist (auto-generated, STRICT)
=================================================================

Filtered from DORAnet's JN1224MIN ruleset ({total_rules} rules) by
scripts/build_bio_glucose_to_tal_whitelist.py.

A rule is kept iff applying it to one literature intermediate
produces the NEXT literature intermediate (canonical SMILES match).
This is the strict whitelist for validating the glucose -> TAL
biosynthetic pathway end-to-end.

Kept: {len(rule_to_steps)} / {total_rules} rules.

Literature pathway:
{pathway_lines}

DO NOT edit by hand. Re-run scripts/build_bio_glucose_to_tal_whitelist.py
to regenerate.
"""

TAL_BIO_GLUCOSE_TO_TAL_WHITELIST = frozenset({{
{body_lines}
}})
'''
    with open(OUT_PY, "w", encoding="utf-8") as f:
        f.write(content)


if __name__ == "__main__":
    main()
