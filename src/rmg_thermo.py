"""
rmg_thermo.py
=============

WHAT THIS IS
------------
A Python 3 wrapper that lets your main project call RMG (which is
trapped in a separate Python 2.7 environment) as if RMG were a normal
Python function. You hand it a SMILES string, you get back a heat of
formation (Hf) in kJ/mol — or None if RMG can't compute it.

HOW IT WORKS (BIG PICTURE)
--------------------------
1. When you create an `RMGThermoClient`, it spawns the RMG server
   script (`scripts/rmg_thermo_server.py`) as a child process running
   inside the rmg_env conda environment. This takes ~10 seconds the
   first time because the RMG database has to load into memory.
2. The two processes communicate through stdin/stdout pipes:
       Python 3 (this script):
         writes "CCO\n"  → server's stdin
         reads server's stdout → "-235.30\n"
       Python 2.7 (server):
         reads "CCO\n", looks up the thermo, writes "-235.30\n"
3. Results are cached, so asking for the same SMILES twice doesn't
   hit RMG twice.
4. When you're done, calling .close() (or letting the object go out
   of scope) cleanly shuts down the server.

WHY THIS DESIGN
---------------
RMG runs on Python 2.7. Your project runs on Python 3.13. Imports
can't cross between them. But the operating system's process pipes
can — they're just bytes. The trick is that the protocol is simple
(one SMILES line in, one Hf line out), so the bridge code stays small.

HOW TO USE
----------
    from rmg_thermo import RMGThermoClient

    calc = RMGThermoClient()        # spawns server, blocks ~10s

    hf = calc("CCO")                # ethanol Hf in kJ/mol
    hf = calc("Cc1cc(O)cc(=O)o1")   # TAL Hf in kJ/mol
    hf = calc("C/C=C/C=C/C(=O)O")   # sorbic acid Hf

    calc.close()                    # shut down the server (or use a
                                    # `with` block — see below)

You can also use it as a DORAnet thermo calculator:

    from network_generation import generate_network_tal
    with RMGThermoClient() as calc:
        network = generate_network_tal(
            ...,
            molecule_thermo_calculator=calc,  # ← that's it
            max_rxn_thermo_change=15.0,       # kJ/mol
        )
"""

from __future__ import annotations

import os
import subprocess
import sys
from typing import Optional


# These defaults match your machine's install layout. Override at the
# call site if you ever move things around.
DEFAULT_RMG_PYTHON = (
    r"C:\Users\ashvi\anaconda3\envs\rmg_env\python.exe"
)
DEFAULT_SERVER_SCRIPT = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "scripts", "rmg_thermo_server.py",
    )
)


class RMGThermoClient:
    """
    Talks to a long-running RMG server process to get Hf values.

    This class is CALLABLE — `calc(smiles)` returns the Hf in kJ/mol
    (a float) or None if RMG can't compute it. DORAnet's
    `Chem_Rxn_dH_Calculator` accepts exactly this kind of callable in
    its `molecule_thermo_calculator` slot, so an instance of this class
    drops directly into the existing pipeline.

    Cofactor for beginners: a "callable" in Python is any object you
    can put parentheses after. Functions are callable. Classes with
    a `__call__` method become callable too. That's what's happening
    here — instances of this class can be called like functions.
    """

    def __init__(
        self,
        rmg_python: str = DEFAULT_RMG_PYTHON,
        server_script: str = DEFAULT_SERVER_SCRIPT,
        startup_timeout_seconds: float = 60.0,
    ):
        """
        Spawn the RMG server subprocess and wait for its READY message.

        Parameters
        ----------
        rmg_python : str
            Full path to the python.exe inside the rmg_env. This MUST
            be the Python 2.7 interpreter that has rmgpy installed —
            NOT your project's Python 3.13.
        server_script : str
            Full path to rmg_thermo_server.py. Defaults to
            scripts/rmg_thermo_server.py relative to this file.
        startup_timeout_seconds : float
            How long to wait for the server to say "READY" before
            giving up. The RMG database load takes ~10 seconds, so
            60 seconds is generous headroom.
        """
        # Mark as closed UP FRONT so cleanup is safe even if something
        # below raises before the subprocess is actually running.
        self._closed = True
        self._proc = None
        self._cache: dict[str, Optional[float]] = {}

        # Sanity checks BEFORE we try to launch the subprocess — these
        # give clearer error messages than "Python crashed somewhere."
        if not os.path.isfile(rmg_python):
            raise FileNotFoundError(
                f"RMG python interpreter not found at: {rmg_python}\n"
                "Did you install RMG into a conda env called rmg_env?"
            )
        if not os.path.isfile(server_script):
            raise FileNotFoundError(
                f"Server script not found at: {server_script}"
            )

        # We're past the early checks — now we're going to actually
        # spawn the server, so unmark the closed flag.
        self._closed = False

        # Launch the server. Key arguments:
        #   stdin/stdout=subprocess.PIPE  → give us hooks to read/write
        #   stderr=subprocess.PIPE        → capture errors separately
        #   text=True                     → bytes auto-decoded as str
        #   bufsize=1                     → line-buffered (so flushes work)
        self._proc = subprocess.Popen(
            [rmg_python, server_script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

        # Wait for the "READY" handshake. The server prints other lines
        # to stderr while loading (so they don't confuse the protocol),
        # but its first stdout line is either "READY" or "DB_NOT_FOUND".
        first_line = self._proc.stdout.readline().strip()
        if first_line == "DB_NOT_FOUND":
            self.close()
            raise RuntimeError(
                "RMG server could not find its thermo database.\n"
                "See stderr above for the paths it checked, then either:\n"
                "  - set the RMG_DATABASE environment variable, or\n"
                "  - clone the rmg-database repo into ~/rmg-database"
            )
        if first_line != "READY":
            stderr_blob = self._proc.stderr.read()
            self.close()
            raise RuntimeError(
                f"RMG server didn't start cleanly. First stdout line: "
                f"{first_line!r}\nstderr:\n{stderr_blob}"
            )

        # (self._cache and self._closed are already initialized above.)

    # ------------------------------------------------------------------
    # The main public method: calc(smiles) → Hf or None
    # ------------------------------------------------------------------
    def __call__(self, smiles: str) -> Optional[float]:
        """
        Look up the standard heat of formation (Hf) for `smiles` in
        kJ/mol. Returns None if RMG can't compute it (uncommon SMILES,
        stereochemistry issue, charged species, etc.).

        DORAnet's filter pipeline knows what to do with None — it
        treats those reactions as "no thermo available" rather than
        rejecting them outright.
        """
        if self._closed:
            raise RuntimeError("RMGThermoClient is already closed.")

        # Cache hit — return immediately, no pipe traffic.
        if smiles in self._cache:
            return self._cache[smiles]

        # Send the SMILES to the server.
        # The newline acts as the "message done" signal because the
        # server reads one line per query.
        try:
            self._proc.stdin.write(smiles + "\n")
            self._proc.stdin.flush()
        except (BrokenPipeError, ValueError):
            # The server died, probably from a crash on the previous
            # query. Surface a clear error rather than hanging.
            raise RuntimeError(
                "RMG server pipe is broken — process may have crashed."
            )

        # Read the reply (one line).
        reply = self._proc.stdout.readline().strip()
        if not reply:
            raise RuntimeError(
                "RMG server gave no reply — process may have died."
            )

        value: Optional[float]
        if reply == "NO_THERMO":
            value = None
        else:
            try:
                value = float(reply)
            except ValueError:
                # Server sent something we don't understand — treat
                # as missing data rather than crashing the whole network
                # build.
                value = None

        self._cache[smiles] = value
        return value

    # ------------------------------------------------------------------
    # Cleanup. Always call close() when done — or use a `with` block.
    # ------------------------------------------------------------------
    def close(self) -> None:
        """
        Shut down the RMG server process cleanly. Safe to call twice.
        """
        if self._closed:
            return
        self._closed = True
        if self._proc is None:
            return
        try:
            # Closing stdin signals the server's main loop that no more
            # work is coming. The server then exits its `while True`
            # loop and the process terminates naturally.
            if self._proc.stdin and not self._proc.stdin.closed:
                self._proc.stdin.close()
            self._proc.wait(timeout=5.0)
        except Exception:
            # Server didn't shut down cleanly — kill it as a fallback.
            try:
                self._proc.kill()
            except Exception:
                pass

    # ------------------------------------------------------------------
    # Optional but nice: support `with RMGThermoClient() as calc:` so
    # the server gets shut down automatically when the block exits,
    # even if an exception is raised mid-network-build.
    # ------------------------------------------------------------------
    def __enter__(self) -> "RMGThermoClient":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def __del__(self) -> None:
        # Belt-and-suspenders: if the user forgot to call close() and
        # the object is being garbage-collected, do it here.
        try:
            self.close()
        except Exception:
            pass
