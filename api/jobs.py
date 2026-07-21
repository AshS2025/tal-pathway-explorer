"""
api/jobs.py — tiny in-memory job store + background execution.

A "run" is one generate(+rank) job. Runs live in a plain dict; the slow
pipeline work runs on a small thread pool so HTTP handlers never block
(clients poll for status instead). Both the dict and the thread pool are
deliberately simple and live ONLY behind this module — they can be
swapped for Redis / a task queue later without changing the API.
"""
from __future__ import annotations

import glob
import os
import threading
import traceback
import uuid
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

# Pipeline work is CPU-heavy (and single-process via the in-process Pool
# shim), so keep the pool small. Bump later or swap for a task queue.
_executor = ThreadPoolExecutor(max_workers=2)
_lock = threading.Lock()
_runs: "OrderedDict[str, Run]" = OrderedDict()

# Cap on retained runs. Each Run holds its full result payload, so an
# uncapped dict is a slow memory leak on a long-lived server. Past this
# cap the OLDEST TERMINAL runs are evicted (and their disk files deleted);
# in-flight runs are never evicted. Bump freely — it's a memory/history
# trade-off, not a correctness knob.
_MAX_RUNS = 200
_TERMINAL_STATUSES = {"generated", "ranked", "error"}


@dataclass
class Run:
    """One run's state. status flows:
    pending → generating → generated → ranking → ranked  (or → error)."""
    id: str
    status: str = "pending"
    config: Any = None                        # PipelineConfig (needed to rank)
    pathways: Optional[list] = None           # unranked generation result
    ranked_pathways: Optional[list] = None    # ranking result
    diagnostics: dict = field(default_factory=dict)
    error: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "status": self.status,
            "pathways": self.pathways,
            "ranked_pathways": self.ranked_pathways,
            "diagnostics": self.diagnostics,
            "error": self.error,
        }


def create_run() -> Run:
    run = Run(id=uuid.uuid4().hex)
    with _lock:
        _runs[run.id] = run
        _evict_locked()
    return run


def _evict_locked() -> None:
    """Drop oldest TERMINAL runs until back under the cap, deleting their
    on-disk artefacts too. Must hold _lock. In-flight runs
    (pending/generating/ranking) are skipped so we never yank a run out
    from under an active worker or a polling client."""
    overflow = len(_runs) - _MAX_RUNS
    if overflow <= 0:
        return
    for rid, run in list(_runs.items()):        # oldest first
        if overflow <= 0:
            break
        if run.status in _TERMINAL_STATUSES:
            del _runs[rid]
            _delete_run_files_locked(run)
            overflow -= 1


def _delete_run_files_locked(run: "Run") -> None:
    """Delete a run's artefact files — but ONLY if no other live run shares
    the same job_name. Cache hits reuse an earlier run's job_name/files, so
    two runs can point at one set of files; deleting them out from under a
    survivor would break its rank/graph. Must hold _lock."""
    job = getattr(run.config, "job_name", None)
    if not job:
        return
    if any(getattr(r.config, "job_name", None) == job for r in _runs.values()):
        return                                  # still referenced — keep files
    purge_job_files(job)


def purge_job_files(job_name: str) -> None:
    """Remove every '{job_name}_*' artefact (pathways, network, graph html,
    reaxys, ...). Safe: job_names are fixed-length UUIDs, so one is never a
    prefix of another."""
    for path in glob.glob(f"{job_name}_*"):
        try:
            os.remove(path)
        except OSError:
            pass


def sweep_orphan_api_artifacts() -> None:
    """Delete leftover API run artefacts from previous server sessions.
    Safe to call at startup: the in-memory run store is empty then, so
    every 'api_*' file on disk is an orphan. Prevents disk growth from
    accumulating across restarts."""
    for path in glob.glob("api_*"):
        try:
            os.remove(path)
        except OSError:                         # dirs / locked files: skip
            pass


def get_run(run_id: str) -> Optional[Run]:
    with _lock:
        return _runs.get(run_id)


def set_status(run: Run, status: str) -> None:
    with _lock:
        run.status = status


def run_in_background(run: Run, worker: Callable[["Run"], None]) -> None:
    """Execute `worker(run)` on the thread pool. The worker mutates `run`
    (its result fields + status). If it raises, the run is marked errored
    with the exception message (and the traceback is logged server-side)."""
    def _task() -> None:
        try:
            worker(run)
        except Exception as e:            # noqa: BLE001 - report any failure
            traceback.print_exc()
            with _lock:
                run.status = "error"
                run.error = f"{type(e).__name__}: {e}"

    _executor.submit(_task)
