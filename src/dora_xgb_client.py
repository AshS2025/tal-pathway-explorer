"""
dora_xgb_client.py
==================

Python-3.13 bridge to the DORA-XGB enzymatic-reaction feasibility model,
which is trapped in its own `dora_xgb` conda env (legacy deps: xgboost
1.6.2, mordred, old numpy). This mirrors RMGThermoClient: we spawn a
long-running server process inside that env and talk to it over
stdin/stdout pipes.

Usage
-----
    from dora_xgb_client import DoraXGBClient
    client = DoraXGBClient()                       # spawns server, ~5-10s
    s = client.feasibility("CC=O.OC(=O)...>>...")  # 0-1 or None
    client.close()                                  # or use a `with` block

The score is a probability in [0, 1]: higher = more likely to be an
enzymatically feasible reaction. None means DORA-XGB couldn't score it.
"""
from __future__ import annotations

import os
import subprocess
import threading
from typing import Optional


DEFAULT_DORA_PYTHON = r"C:\Users\ashvi\anaconda3\envs\dora_xgb\python.exe"
DEFAULT_SERVER_SCRIPT = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "scripts", "dora_xgb_server.py",
    )
)


class DoraXGBClient:
    """Talks to a long-running DORA-XGB server for reaction feasibility.

    `feasibility(rxn_smiles)` returns a float in [0, 1] or None. Results
    are cached, so scoring the same reaction twice is free.
    """

    def __init__(
        self,
        dora_python: str = DEFAULT_DORA_PYTHON,
        server_script: str = DEFAULT_SERVER_SCRIPT,
    ):
        self._closed = True
        self._proc = None
        self._cache: dict[str, Optional[float]] = {}
        # The subprocess speaks a strict one-line-in / one-line-out protocol
        # over a single pipe, so only one thread may hold that conversation
        # at a time — otherwise concurrent callers cross their replies. The
        # API runs jobs on a thread pool, so this lock is load-bearing.
        self._lock = threading.Lock()

        if not os.path.isfile(dora_python):
            raise FileNotFoundError(
                f"dora_xgb env python not found at: {dora_python}\n"
                "Create it with:  conda create -n dora_xgb python=3.9 -y  "
                "then  pip install DORA-XGB"
            )
        if not os.path.isfile(server_script):
            raise FileNotFoundError(f"Server script not found at: {server_script}")

        self._closed = False
        self._proc = subprocess.Popen(
            [dora_python, server_script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

        # Wait for the READY handshake (model load takes several seconds).
        first = self._proc.stdout.readline().strip()
        if first != "READY":
            stderr_blob = ""
            try:
                stderr_blob = self._proc.stderr.read()
            except Exception:
                pass
            self.close()
            raise RuntimeError(
                f"DORA-XGB server didn't start cleanly. First stdout line: "
                f"{first!r}\nstderr:\n{stderr_blob}"
            )

    def feasibility(self, rxn_smiles: str) -> Optional[float]:
        """Feasibility score in [0, 1] for a 'reactants>>products' reaction,
        or None if it can't be scored."""
        if self._closed:
            raise RuntimeError("DoraXGBClient is already closed.")
        # Hold the lock across the cache check AND the pipe round-trip: two
        # threads must never interleave writes/reads on the shared pipe, and
        # the cache read/write is guarded for free.
        with self._lock:
            if rxn_smiles in self._cache:
                return self._cache[rxn_smiles]
            value = self._query(rxn_smiles)
            self._cache[rxn_smiles] = value
            return value

    def _query(self, rxn_smiles: str) -> Optional[float]:
        """One write+read round-trip on the subprocess pipe. The caller MUST
        hold self._lock — the protocol is strictly one line in, one line out,
        so overlapping calls would cross their replies."""
        try:
            self._proc.stdin.write(rxn_smiles + "\n")
            self._proc.stdin.flush()
        except (BrokenPipeError, ValueError):
            raise RuntimeError("DORA-XGB server pipe is broken — process may have crashed.")

        reply = self._proc.stdout.readline().strip()
        if not reply:
            raise RuntimeError("DORA-XGB server gave no reply — process may have died.")

        if reply == "NO_SCORE":
            return None
        try:
            return float(reply)
        except ValueError:
            return None

    # convenience: make the instance callable, like RMGThermoClient
    def __call__(self, rxn_smiles: str) -> Optional[float]:
        return self.feasibility(rxn_smiles)

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        if self._proc is None:
            return
        try:
            if self._proc.stdin and not self._proc.stdin.closed:
                self._proc.stdin.close()
            self._proc.wait(timeout=5.0)
        except Exception:
            try:
                self._proc.kill()
            except Exception:
                pass

    def __enter__(self) -> "DoraXGBClient":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass
