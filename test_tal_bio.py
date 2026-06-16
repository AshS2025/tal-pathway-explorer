"""
test_tal_bio.py
Standalone smoke test for DORAnet's biological (enzymatic) network on TAL.

Goal: confirm that the JN ruleset can reach TAL from a primary-metabolism
seed (pyruvate), with cofactors auto-injected by the bio module.

Why pyruvate as starter (not acetyl-CoA):
  Acetyl-CoA, CoA, ATP, NAD(P)(H), CO2, water, H+, O2 are all in DORAnet's
  cofactor list (modules/enzymatic/all_cofactors.tsv) — they get
  auto-injected as helpers. Pyruvate is a non-cofactor seed that connects
  to acetyl-CoA via pyruvate dehydrogenase rules, which the network can
  then condense into pyranones.

Why allow_multiple_reactants=True:
  TAL is built from 3 × acetyl-CoA, all of which are cofactors. The default
  Reaction_Type_Filter requires exactly ONE non-cofactor substrate per
  reaction — so cofactor+cofactor condensations (the entire polyketide
  machinery) would be rejected. We bypass the filter by setting
  allow_multiple_reactants=True.

Run from project root:
    python test_tal_bio.py
"""

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

from doranet.modules.enzymatic.generate_network import generate_network
from rdkit import Chem


JOB_NAME = "tal_bio_test"
STARTER = "CC(=O)C(=O)O"             # pyruvate
TARGET = "CC1=CC(O)=CC(=O)O1"        # TAL


def banner(text):
    print("\n" + "=" * 70)
    print(text)
    print("=" * 70)


def run(ruleset, gen, max_atoms):
    banner(f"BIO NETWORK | ruleset={ruleset} | gen={gen} | max_atoms={max_atoms}")
    t0 = time.perf_counter()
    network = generate_network(
        job_name=f"{JOB_NAME}_{ruleset}_g{gen}",
        starters=STARTER,
        gen=gen,
        direction="forward",
        max_atoms=max_atoms,
        allow_multiple_reactants=True,
        targets=TARGET,
        ruleset=ruleset,
    )
    elapsed = time.perf_counter() - t0
    print(f"\nElapsed: {elapsed:.1f} s")
    print(f"Molecules: {len(network.mols)}  Reactions: {len(network.rxns)}")

    # Did we reach TAL?
    target_canonical = Chem.MolToSmiles(Chem.MolFromSmiles(TARGET))
    mol_smiles = set()
    for m in network.mols:
        try:
            mol_smiles.add(Chem.MolToSmiles(Chem.MolFromSmiles(m.uid)))
        except Exception:
            pass
    found = target_canonical in mol_smiles
    print(f"\nTAL canonical: {target_canonical}")
    print(f"TAL in network: {found}")
    return network, found


def main():
    # Smaller ruleset first: 1224 rules, faster, sanity check.
    # gen=2 is a tight depth — pyruvate → acetyl-CoA → first condensation.
    # gen=3 to give the polyketide path room to grow.
    # Keep max_atoms generous: TAL has 6 C, 4 O. Intermediates may be bigger.
    run("JN1224MIN", gen=2, max_atoms={"C": 20, "O": 10, "N": 5})
    run("JN1224MIN", gen=3, max_atoms={"C": 20, "O": 10, "N": 5})
    # If JN1224MIN doesn't reach TAL, try the full ruleset.
    # 3604 rules × cartesian expansion will be slower — start with gen=2.
    run("JN3604IMT", gen=2, max_atoms={"C": 20, "O": 10, "N": 5})


if __name__ == "__main__":
    main()
