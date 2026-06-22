"""
test_rmg_bridge.py
==================

WHAT THIS IS
------------
A "smoke test" — the smallest possible script that confirms the
Python 3 ↔ Python 2.7 RMG bridge actually works on your machine
before you wire it into the full network expansion.

WHY RUN THIS FIRST
------------------
The bridge has a lot of moving parts: a subprocess spawn, an env
path, the database load, stdin/stdout protocol. If any one of them
is broken, you want to catch it here — not 5 minutes into a network
expansion that fails for an unrelated-looking reason.

HOW TO RUN
----------
From the project root (the tal-pathway-explorer folder):

    python test_rmg_bridge.py

EXPECTED OUTPUT (rough)
-----------------------
Hf values (in kJ/mol) for a few well-known molecules. RMG should
give numbers close to literature:
    water        ≈ -241 to -286 kJ/mol  (gas or liquid phase)
    methane      ≈ -74 kJ/mol
    ethanol      ≈ -235 kJ/mol
    TAL          something around -400 to -500 kJ/mol (highly oxidized)

If you see numbers in roughly those ranges, the bridge works. If you
see all "NO_THERMO" results, something is wrong inside RMG (probably
database location). If the script hangs forever, the protocol or pipe
buffering is wrong — Ctrl+C and ping me.
"""

import os
import sys
import time

# Make `src/` importable without installing the package.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from rmg_thermo import RMGThermoClient


# A small set of well-known molecules with rough literature Hf values
# (gas-phase, 298 K). Treat these as ballpark sanity checks, not
# precise validation.
KNOWN = [
    ("water (H2O)",         "O",                  -241.8),
    ("methane (CH4)",       "C",                   -74.6),
    ("methanol (CH3OH)",    "CO",                 -200.9),
    ("ethanol (CH3CH2OH)",  "CCO",                -234.8),
    ("acetic acid",         "CC(=O)O",            -432.0),
    ("acetone",             "CC(=O)C",            -218.5),
    # TAL is large and oxidized — no clean literature reference,
    # but we expect it to be in the -400 to -500 kJ/mol ballpark.
    ("TAL",                 "Cc1cc(O)cc(=O)o1",    None),
    ("sorbic acid",         "CC=CC=CC(=O)O",       None),
]


def main():
    print("=" * 64)
    print(" RMG bridge smoke test")
    print("=" * 64)
    print()
    print("Spawning RMG server (this takes ~10 seconds while the")
    print("ThermoDatabase loads)...")
    print()

    t0 = time.time()
    # Use a `with` block so the server gets cleanly shut down even if
    # something below raises an exception.
    with RMGThermoClient() as calc:
        startup = time.time() - t0
        print(f"  Server READY  (took {startup:.1f}s)")
        print()

        print(f"  {'molecule':24s}  {'SMILES':20s}  "
              f"{'Hf (kJ/mol)':>14s}  {'literature':>12s}")
        print("  " + "-" * 76)

        for label, smiles, lit in KNOWN:
            t_query = time.time()
            value = calc(smiles)
            dt = time.time() - t_query

            if value is None:
                hf_str = "NO_THERMO"
            else:
                hf_str = f"{value:.2f}"

            lit_str = "—" if lit is None else f"{lit:.1f}"

            print(f"  {label:24s}  {smiles:20s}  {hf_str:>14s}  "
                  f"{lit_str:>12s}    ({dt*1000:.0f}ms)")

        print()
        # Cache sanity check — second call should be effectively free.
        t = time.time()
        _ = calc("CCO")
        cache_ms = (time.time() - t) * 1000
        print(f"  Cache hit on ethanol: {cache_ms:.2f}ms "
              f"(should be near zero)")

    print()
    print("Bridge test complete. If you see Hf values close to the")
    print("'literature' column above, the bridge is working and you")
    print("can plug RMGThermoClient into generate_network_tal as")
    print("`molecule_thermo_calculator=calc`.")


if __name__ == "__main__":
    main()
