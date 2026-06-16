"""
build_bio_whitelist.py
Auto-generate the TAL bio reaction whitelist by prefiltering JN1224MIN's
1224 enzymatic rules against a TAL chemistry seed set.

Why auto-generation:
  Unlike the chem operators (388 rules with human-readable names that map
  to textbook chemistry), the bio rules are SMARTS-encoded named only
  "rule0001", "rule0002", ... A manual whitelist is intractable. Instead
  we ask each rule: "can you fire on at least one TAL-relevant seed?"
  Rules that can't are dropped.

Method:
  1. Load JN1224MIN_rules.tsv and all_cofactors.tsv from DORAnet.
  2. For each rule:
       - Parse SMARTS via RDKit.
       - The rule's Reactants column is a ";"-separated list whose entries
         are either "Any" (a substrate slot) or a cofactor ID.
       - Build candidate reactant tuples by:
           * substituting cofactor SMILES at cofactor slots,
           * trying each seed at each Any slot (with TAL filling the
             other Any slots if there are multiple).
       - Run RDKit RunReactants on each candidate. If ANY combination
         yields at least one sanitizable product, the rule passes.
  3. Write the passing names to src/tal_bio_reaction_whitelist.py as a
     frozenset, mirroring the chem whitelist's interface.

Run from project root:
    python scripts/build_bio_whitelist.py
"""

import csv
import os
import sys
import time
from itertools import product as iproduct

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")  # silence the per-reaction sanitize warnings

# --- paths --------------------------------------------------------------
DORANET_BIO = (
    "C:/Users/ashvi/anaconda3/Lib/site-packages/doranet/modules/enzymatic"
)
RULES_TSV = os.path.join(DORANET_BIO, "JN1224MIN_rules.tsv")
COFACTORS_TSV = os.path.join(DORANET_BIO, "all_cofactors.tsv")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_PY = os.path.normpath(
    os.path.join(HERE, "..", "src", "tal_bio_reaction_whitelist.py")
)

# --- TAL chemistry seed set --------------------------------------------
# Each seed: SMILES + short label. The label is only for the output
# header — what matters is that at least one of these reacts.
SEEDS = [
    ("CC1=CC(O)=CC(=O)O1",       "TAL (target)"),
    ("CC(=O)C(=O)O",             "pyruvate"),
    ("CC(=O)O",                  "acetate"),
    ("CC(=O)CC(=O)O",            "acetoacetate"),
    ("CC(=O)CC(=O)CC(=O)O",      "triacetic acid"),
    ("CC(=O)C=O",                "methylglyoxal"),
    ("CC=O",                     "acetaldehyde"),
    ("CCO",                      "ethanol"),
    ("CO",                       "methanol"),
    ("CC(=O)C",                  "acetone"),
    ("CC(=O)C(=O)C",             "diacetyl"),
    ("OC(=O)CC(=O)O",            "malonate"),
    ("OCC(O)C(=O)O",             "glycerate"),
    ("Cc1cc(O)cc(O)c1",          "orcinol"),
    ("Cc1cc(O)cc(O)c1C(=O)O",    "orsellinic acid"),
    ("O=C1OC=CC(O)=C1",          "4-hydroxy-2-pyrone"),
    ("CCc1ccc(C)oc1=O",          "3-ethyl-6-methyl-2-pyrone"),
]


def load_cofactors():
    """Return dict: cofactor_id -> canonical SMILES."""
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
    """Yield (name, reactants_list, smarts) for each rule."""
    with open(RULES_TSV, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row["Name"]
            reactants = row["Reactants"].split(";")
            smarts = row["SMARTS"]
            yield name, reactants, smarts


def prepare_seed_mols():
    out = []
    for smi, label in SEEDS:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            print(f"  WARNING: seed {label} ({smi}) failed to parse, skipping")
            continue
        out.append((m, label, smi))
    return out


def prepare_cofactor_mols(cofactors):
    """Pre-parse cofactor SMILES. Return dict id -> Mol."""
    out = {}
    for cid, smi in cofactors.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        out[cid] = m
    return out


def candidate_reactant_tuples(reactant_types, seed_mols, cof_mols, tal_mol):
    """
    Given a rule's reactant type list (e.g. ["Any", "NAD_CoF"]),
    yield tuples of RDKit Mols to try as reactants.

    Strategy:
      - Replace cofactor slots with their parsed Mol (skip rule if any
        cofactor is missing from our table).
      - For Any slots: if there's a single Any slot, try each seed.
        If there are multiple, fix TAL at the non-varying Any slots
        and walk one Any slot through all seeds at a time.
    """
    any_indices = [i for i, t in enumerate(reactant_types) if t == "Any"]
    cof_indices = [i for i, t in enumerate(reactant_types) if t != "Any"]

    base = [None] * len(reactant_types)
    for i in cof_indices:
        cid = reactant_types[i]
        if cid not in cof_mols:
            return  # unknown cofactor; cannot test this rule
        base[i] = cof_mols[cid]

    if not any_indices:
        # All slots are cofactors. Yield the all-cofactor tuple once.
        yield tuple(base)
        return

    if len(any_indices) == 1:
        i = any_indices[0]
        for seed_mol, _label, _smi in seed_mols:
            t = list(base)
            t[i] = seed_mol
            yield tuple(t)
        return

    # 2+ Any slots: walk one slot through seeds, fill the rest with TAL.
    for varying in any_indices:
        for seed_mol, _label, _smi in seed_mols:
            t = list(base)
            for j in any_indices:
                t[j] = seed_mol if j == varying else tal_mol
            yield tuple(t)


def rule_fires(rxn, reactant_tuples):
    """Return True if any reactant tuple yields a sanitizable product."""
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
                except Exception:
                    continue
                # Got at least one sanitizable product. Done.
                return True
    return False


def main():
    print(f"Loading cofactors from {COFACTORS_TSV}")
    cofactors = load_cofactors()
    cof_mols = prepare_cofactor_mols(cofactors)
    print(f"  loaded {len(cof_mols)} cofactor Mols")

    print(f"\nPreparing {len(SEEDS)} TAL-chemistry seeds")
    seed_mols = prepare_seed_mols()
    tal_mol = seed_mols[0][0]  # TAL is first
    print(f"  using {len(seed_mols)} seeds")

    print(f"\nLoading rules from {RULES_TSV}")
    rules = list(load_rules())
    print(f"  loaded {len(rules)} rules")

    print("\nPrefiltering...")
    t0 = time.perf_counter()
    kept = []
    unknown_cof_rules = 0
    parse_failures = 0
    for n, (name, reactant_types, smarts) in enumerate(rules):
        if n and n % 100 == 0:
            print(
                f"  [{n}/{len(rules)}] kept so far: {len(kept)} "
                f"({time.perf_counter()-t0:.1f}s)"
            )

        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
        except Exception:
            parse_failures += 1
            continue
        if rxn is None:
            parse_failures += 1
            continue

        tuples = list(
            candidate_reactant_tuples(
                reactant_types, seed_mols, cof_mols, tal_mol
            )
        )
        if not tuples:
            unknown_cof_rules += 1
            continue

        if rule_fires(rxn, tuples):
            kept.append(name)

    elapsed = time.perf_counter() - t0
    print(f"\nDone in {elapsed:.1f}s")
    print(f"  rules total:               {len(rules)}")
    print(f"  parse failures:            {parse_failures}")
    print(f"  skipped (unknown cofactor):{unknown_cof_rules}")
    print(f"  kept (fired on a seed):    {len(kept)}")

    write_whitelist(kept, len(rules))
    print(f"\nWrote {OUT_PY}")


def write_whitelist(kept_names, total):
    seed_lines = "\n".join(
        f"#   {label:30s} {smi}" for smi, label in SEEDS
    )
    body_lines = "\n".join(f'    "{n}",' for n in kept_names)
    content = f'''"""
TAL Bio Reaction Whitelist for DORAnet (auto-generated)
=======================================================

Filtered from DORAnet's JN1224MIN ruleset ({total} rules) by
scripts/build_bio_whitelist.py. A rule is kept iff applying its SMARTS
to at least one TAL-chemistry seed yields a sanitizable product.

Kept: {len(kept_names)} / {total} rules.

Seeds used:
{seed_lines}

Usage in network_generation:
    from tal_bio_reaction_whitelist import TAL_BIO_REACTION_WHITELIST
    # when iterating bio rules:
    if rule_name in TAL_BIO_REACTION_WHITELIST:
        network.add_op(...)

DO NOT edit by hand. Re-run scripts/build_bio_whitelist.py to regenerate.
"""

TAL_BIO_REACTION_WHITELIST = frozenset({{
{body_lines}
}})
'''
    with open(OUT_PY, "w", encoding="utf-8") as f:
        f.write(content)


if __name__ == "__main__":
    main()
