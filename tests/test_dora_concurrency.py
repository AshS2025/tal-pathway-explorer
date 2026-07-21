"""The DORA-XGB client talks to its subprocess over a single one-in/one-out
pipe. The API runs jobs on a thread pool, so concurrent callers must be
serialized or their replies cross. This proves the lock does that — without
needing the real dora_xgb env: we subclass the client, skip the subprocess,
and make the pipe round-trip (_query) detect any overlap."""
import threading
import time

from dora_xgb_client import DoraXGBClient


class _Probe(DoraXGBClient):
    """A DoraXGBClient with the real locking logic but a fake pipe. _query
    records how many threads are inside it at once; if the lock in
    feasibility() works, that count never exceeds 1."""

    def __init__(self):
        self._closed = False
        self._proc = None
        self._cache = {}
        self._lock = threading.Lock()
        self._active = 0
        self.max_active = 0

    def _query(self, rxn_smiles):
        self._active += 1
        self.max_active = max(self.max_active, self._active)
        time.sleep(0.005)          # widen the window so a broken lock is caught
        self._active -= 1
        return float(len(rxn_smiles))   # deterministic stand-in "score"


def test_feasibility_serializes_concurrent_callers():
    client = _Probe()
    rxns = [f"A{i}>>B{i}" for i in range(20)]
    results: dict = {}

    def worker(r):
        results[r] = client.feasibility(r)

    threads = [threading.Thread(target=worker, args=(r,)) for r in rxns]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    # never two threads inside the pipe round-trip at once
    assert client.max_active == 1
    # every caller got ITS OWN correct answer (no crossed replies)
    assert all(results[r] == float(len(r)) for r in rxns)


def test_cache_hit_skips_the_pipe():
    client = _Probe()
    client.feasibility("A>>B")
    calls_before = client.max_active
    # second identical call should come from cache, not re-enter _query
    hits = {"n": 0}
    orig = client._query

    def counting_query(r):
        hits["n"] += 1
        return orig(r)

    client._query = counting_query
    client.feasibility("A>>B")      # cached -> _query must NOT run
    assert hits["n"] == 0
    assert calls_before == 1
