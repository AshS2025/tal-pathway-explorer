"""
test_rmg_pruning.py
===================

WHAT THIS IS
------------
A before/after comparison. Runs the SAME TAL network expansion twice:
  1. Without RMG  — no thermo data, no thermo pruning
  2. With RMG     — every candidate reaction is screened by its
                    computed ΔH; endothermic reactions get rejected

WHY THIS MATTERS
----------------
Thermo pruning is a HARD filter (binary reject/keep), not a soft
ranking. Reactions whose ΔH exceeds the threshold never enter the
network — so the next generation doesn't have to consider them as
substrates either. This compounds over generations: pruning reduces
the molecule pool, which reduces per-rule trial count, which reduces
the next generation's pool, etc. The hope is that the speedup from a
smaller network outweighs the cost of the RMG queries.

WHAT TO LOOK FOR IN THE OUTPUT
------------------------------
  Run 1 (no thermo):  fast expansion, many reactions, may include
                      stuff that wouldn't actually run in a lab.
  Run 2 (with RMG):   slower for the first few generations (RMG
                      queries cost ~5-100 ms each, cache helps after).
                      Fewer reactions kept. Smaller resulting network.

  Key metric: REACTIONS_PRUNED / REACTIONS_NO_THERMO. A high ratio
  (say >30%) means thermo is doing meaningful work. A low ratio
  (say <5%) means the current threshold is too loose.

UNITS REMINDER
--------------
RMG returns Hf in kJ/mol. The Chem_Rxn_dH_Calculator subtracts to
get ΔH in kJ/mol. The Rxn_dH_Filter rejects reactions whose ΔH is
more positive than max_rxn_thermo_change — which must therefore
ALSO be in kJ/mol. (15 kJ/mol is fairly permissive; you may want
to tighten to 5-10 kJ/mol once you see how it behaves.)
"""

import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
sys.stdout.reconfigure(encoding="utf-8")

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from network_generation import generate_network_tal
from rmg_thermo import RMGThermoClient


STARTER_FILE = "test_rmg_pruning_starter.smi"
HELPER_FILE  = "test_rmg_pruning_helpers.smi"

TAL = "Cc1cc(O)cc(=O)o1"
GEN = 3                              # deeper — where thermo pruning starts paying off
THERMO_THRESHOLD_KJ = 15.0           # kJ/mol


def write_smi(path, lines):
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


# Common kwargs for both runs — only the calculator changes between
# them, so any difference in network size is attributable to thermo
# pruning (not to other random pipeline differences).
def common_kwargs(thermo_calc):
    return dict(
        starters=STARTER_FILE,
        helpers=HELPER_FILE,
        gen=GEN,
        direction="forward",
        molecule_thermo_calculator=thermo_calc,
        max_rxn_thermo_change=THERMO_THRESHOLD_KJ,
        max_atoms={"C": 12, "O": 5, "N": 0, "S": 0},
        max_molecular_weight=300,
        allow_multiple_reactants="default",
        strategy="cartesian",
        min_carbons=0,
        include_chem=True,
        include_bio=False,
    )


def main():
    print("=" * 64)
    print(" RMG thermo pruning — before/after comparison")
    print("=" * 64)
    print(f" starter = TAL, gen = {GEN}, threshold = {THERMO_THRESHOLD_KJ} kJ/mol")

    write_smi(STARTER_FILE, [TAL])
    write_smi(HELPER_FILE,  ["O", "[H][H]"])

    # ----------------- RUN 1: no thermo (baseline) --------------------
    print("\n[1/2] Run WITHOUT thermo (calculator = None)...")
    t0 = time.time()
    network_no_thermo = generate_network_tal(
        job_name="rmg_pruning_off",
        **common_kwargs(thermo_calc=None),
    )
    t_no_thermo = time.time() - t0
    n_mols_no = len(network_no_thermo.mols)
    n_rxns_no = len(network_no_thermo.rxns)
    print(f"  baseline: {n_mols_no} mols, {n_rxns_no} reactions, "
          f"{t_no_thermo:.1f}s")

    # ----------------- RUN 2: RMG thermo pruning ON -------------------
    print("\n[2/2] Run WITH RMG thermo pruning...")
    print("      (spawning RMG server — first query is slow, then it's fast)")

    with RMGThermoClient() as calc:
        t0 = time.time()
        network_with_thermo = generate_network_tal(
            job_name="rmg_pruning_on",
            **common_kwargs(thermo_calc=calc),
        )
        t_with_thermo = time.time() - t0
        n_mols_yes = len(network_with_thermo.mols)
        n_rxns_yes = len(network_with_thermo.rxns)
        print(f"  with thermo: {n_mols_yes} mols, {n_rxns_yes} reactions, "
              f"{t_with_thermo:.1f}s")
        # Peek inside the cache to report how many unique molecules
        # RMG actually had to look up.
        n_unique_queries = len(calc._cache)
        n_hits_real_value = sum(1 for v in calc._cache.values() if v is not None)
        n_hits_no_data = n_unique_queries - n_hits_real_value
        print(f"  RMG queries: {n_unique_queries} unique molecules")
        print(f"    with Hf available:        {n_hits_real_value}")
        print(f"    NO_THERMO (RMG couldn't): {n_hits_no_data}")

    # ----------------- COMPARISON -------------------------------------
    print()
    print("=" * 64)
    print(" Comparison")
    print("=" * 64)
    print(f"  Molecules  : {n_mols_no:5d} (off)  →  {n_mols_yes:5d} (on)   "
          f"diff = {n_mols_no - n_mols_yes:+d}")
    print(f"  Reactions  : {n_rxns_no:5d} (off)  →  {n_rxns_yes:5d} (on)   "
          f"diff = {n_rxns_no - n_rxns_yes:+d}")
    if n_rxns_no > 0:
        pct = 100.0 * (n_rxns_no - n_rxns_yes) / n_rxns_no
        print(f"  Reactions pruned by thermo: {pct:.1f}%")
    print(f"  Runtime    : {t_no_thermo:5.1f}s (off)  →  "
          f"{t_with_thermo:5.1f}s (on)")
    print()
    print(" If 'Reactions pruned by thermo' is >20%, thermo is doing")
    print(" meaningful work and we should turn it on by default for")
    print(" all expansions. If <5%, we might want a tighter threshold")
    print(" (try max_rxn_thermo_change=5.0 kJ/mol) or the threshold is")
    print(" already weeding out the obvious bad reactions and the")
    print(" remaining ones are mostly thermo-favorable already.")


if __name__ == "__main__":
    main()
