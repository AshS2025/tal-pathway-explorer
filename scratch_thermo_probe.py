"""
scratch_thermo_probe.py — phase 3
Verify the 1.4 GB cache still loads under the downgraded API (0.6.x), and
that LocalCompoundCache can register a novel SMILES.
"""
import sys, time
sys.stdout.reconfigure(encoding="utf-8")

from equilibrator_api import ComponentContribution, Q_
from rdkit import Chem

t0 = time.perf_counter()
cc = ComponentContribution()
print(f"[{time.perf_counter()-t0:5.1f}s] cc ready (api 0.6.x)")

# --- (A) cache still works on known compounds? ---
print("\n--- (A) cache loads, known InChIKey lookups still work ---")
for name, smi in [("ethanol","CCO"), ("water","O"), ("H2","[H][H]"), ("TAL","CC1=CC(O)=CC(=O)O1")]:
    ik = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
    try:
        result = cc.search_compound_by_inchi_key(ik)
        print(f"  {name:8s} {ik}  -> {result[:1] if result else '[]'}")
    except Exception as e:
        print(f"  {name:8s} ERROR {type(e).__name__}: {e}")

# --- (B) LocalCompoundCache registration of a novel SMILES ---
print("\n--- (B) LocalCompoundCache available? ---")
try:
    from equilibrator_assets.local_compound_cache import LocalCompoundCache
    print("  import OK")
except Exception as e:
    print(f"  import FAIL: {type(e).__name__}: {e}")
    sys.exit(1)

import tempfile, pathlib
tmp = pathlib.Path(tempfile.mkdtemp(prefix="lemnisca_thermo_probe_"))
db_path = tmp / "local.sqlite"
print(f"  local cache dir: {tmp}")

t = time.perf_counter()
lcc = LocalCompoundCache()
# 0.6 API: load_cache(filename) initializes a fresh sqlite or opens existing
try:
    lcc.load_cache(str(db_path))
    print(f"  load_cache OK ({time.perf_counter()-t:.2f}s)")
except Exception as e:
    print(f"  load_cache FAIL: {type(e).__name__}: {e}")
    sys.exit(1)

# Now try to register a novel pyranone
novel = "CCc1c(O)cc(C)oc1=O"  # 3-Et-TAL-OH
print(f"\n  registering novel SMILES: {novel}")
t = time.perf_counter()
try:
    df = lcc.add_compounds([novel])
    print(f"  add_compounds returned ({time.perf_counter()-t:.2f}s):")
    print(df)
except Exception as e:
    print(f"  add_compounds FAIL: {type(e).__name__}: {e}")

# Try to look it up after registration
print("\n  post-registration lookup:")
try:
    cpd = lcc.get_compound(f"smiles:{novel}")
    print(f"  get_compound -> {cpd}")
except Exception as e:
    print(f"  get_compound FAIL: {type(e).__name__}: {e}")

print("\ndone.")
