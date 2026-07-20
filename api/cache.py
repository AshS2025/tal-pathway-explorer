"""
api/cache.py — result cache keyed on canonical input parameters.

Generation (network expansion) and ranking are both slow, so re-running
with identical inputs is pure waste. We hash the parameters that actually
affect the OUTPUT and stash the result under that hash; an identical
request later is served instantly.

"Identical" is stricter than a naive dict comparison, and that's the whole
point of this module:
  * `job_name` is EXCLUDED — it's a per-run UUID, so leaving it in would
    mean no two requests ever collide (the classic caching bug).
  * SMILES are CANONICALIZED with the same call the app already trusts
    (Chem.MolToSmiles(Chem.MolFromSmiles(x)) — see validate_config), so two
    spellings of one molecule collide. This keeps isomeric detail, so cis
    and trans stay DISTINCT — we never merge two real stereoisomers.
  * whitelists are ORDER-INDEPENDENT (sorted + de-duped); order doesn't
    change which operators fire.
  * numeric fields are normalized to float so 200 and 200.0 match.

RULE OF THUMB: every PipelineConfig field except job_name belongs in the
generation key. If you add a new knob to the config, add it here too — the
safe default is "include it," otherwise a changed knob silently returns a
stale result.

Storage is a bounded in-memory dict TODAY. It lives entirely behind this
module (like jobs.py's store), so it can be swapped for Redis / disk later
by reimplementing ResultCache's get/put — the API layer never changes.
It's in-memory, so it resets whenever the server process restarts.
"""
from __future__ import annotations

import hashlib
import json
import threading
from collections import OrderedDict
from typing import Any, Optional

from rdkit import Chem


# --- canonicalization helpers ------------------------------------------------

def _canon_smiles(smi: str) -> str:
    """RDKit-canonical SMILES, or the stripped original if it won't parse
    (validation happens elsewhere; we just want a stable key here)."""
    if not smi:
        return ""
    mol = Chem.MolFromSmiles(smi.strip())
    return Chem.MolToSmiles(mol) if mol is not None else smi.strip()


def _canon_smiles_list(items: Optional[list]) -> Optional[list]:
    """Canonicalize each SMILES, then sort + de-dupe (order/duplication of
    a co-reactant pool doesn't change the network)."""
    if items is None:
        return None
    return sorted({_canon_smiles(s) for s in items if s and s.strip()})


def _canon_names(items: Optional[list]) -> Optional[list]:
    """Sort + de-dupe a whitelist of operator/rule names (order-independent)."""
    if items is None:
        return None
    return sorted({s.strip() for s in items if s and s.strip()})


# Config fields that affect the generation OUTPUT. Deliberately omits
# job_name (per-run UUID) and the SMILES/whitelist fields (handled
# specially above). Numeric fields are cast to float when hashed.
_GEN_SCALAR_FIELDS = (
    "domain", "direction", "generations", "strategy", "beam_size",
    "max_molecular_weight", "max_atoms_c", "max_atoms_o", "max_atoms_n",
    "max_rxn_dh", "enable_rmg", "enable_equilibrator",
    "equilibrator_prune_max_abs_dg", "feasibility_prune_threshold",
)


def _hash(payload: dict) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def generation_key(config: Any) -> str:
    """Stable hash of everything that determines a generation's pathways."""
    payload: dict = {
        "starter": _canon_smiles(config.starter_smiles),
        "target": _canon_smiles(config.target_smiles),
        "helpers": _canon_smiles_list(config.helpers),
        "bio_whitelist": _canon_names(config.bio_whitelist),
        "chem_whitelist": _canon_names(config.chem_whitelist),
    }
    for f in _GEN_SCALAR_FIELDS:
        v = getattr(config, f)
        # normalize numeric types so 3 == 3.0 and 200 == 200.0 (bool stays bool)
        payload[f] = float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else v
    return _hash(payload)


def _canon_weights(w: Optional[dict]) -> Optional[dict]:
    """Weights sorted by key with float-normalized values (None = defaults)."""
    if not w:
        return None
    return {k: float(v) for k, v in sorted(w.items())}


def ranking_key(
    config: Any,
    weights: Optional[dict],
    layer_weights: Optional[dict],
    lemnisca_weights: Optional[dict],
) -> str:
    """Ranking depends on the generated pathways PLUS the weight tiers, so
    the key NESTS the generation key and adds the (canonical) weights.
    Change a weight → new key → re-rank; change nothing → instant hit."""
    return _hash({
        "gen": generation_key(config),
        "weights": _canon_weights(weights),
        "layer_weights": _canon_weights(layer_weights),
        "lemnisca_weights": _canon_weights(lemnisca_weights),
    })


# --- the store ---------------------------------------------------------------

class ResultCache:
    """Bounded, thread-safe LRU cache for generation and ranking results.

    Two independent maps (generation vs ranking) sharing one lock. Payloads
    are plain JSON-serializable dicts, so a future Redis/disk backend just
    reimplements get/put.
    """

    def __init__(self, max_entries: int = 64):
        self._gen: "OrderedDict[str, dict]" = OrderedDict()
        self._rank: "OrderedDict[str, dict]" = OrderedDict()
        self._max = max_entries
        self._lock = threading.Lock()

    def _get(self, store: "OrderedDict[str, dict]", key: str) -> Optional[dict]:
        with self._lock:
            if key not in store:
                return None
            store.move_to_end(key)          # mark most-recently-used
            return store[key]

    def _put(self, store: "OrderedDict[str, dict]", key: str, value: dict) -> None:
        with self._lock:
            store[key] = value
            store.move_to_end(key)
            while len(store) > self._max:
                store.popitem(last=False)   # evict least-recently-used

    # generation: payload = {"pathways", "diagnostics", "job_name"}
    def get_generation(self, key: str) -> Optional[dict]:
        return self._get(self._gen, key)

    def put_generation(self, key: str, payload: dict) -> None:
        self._put(self._gen, key, payload)

    # ranking: payload = {"ranked_pathways", "diagnostics"}
    def get_ranking(self, key: str) -> Optional[dict]:
        return self._get(self._rank, key)

    def put_ranking(self, key: str, payload: dict) -> None:
        self._put(self._rank, key, payload)

    def clear(self) -> None:
        with self._lock:
            self._gen.clear()
            self._rank.clear()


# module-level singleton (mirrors jobs.py's module-level store)
cache = ResultCache()
