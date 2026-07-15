"""
api/jobs.py — tiny in-memory job store + background execution.

A "run" is one generate(+rank) job. Runs live in a plain dict; the slow
pipeline work runs on a small thread pool so HTTP handlers never block
(clients poll for status instead). Both the dict and the thread pool are
deliberately simple and live ONLY behind this module — they can be
swapped for Redis / a task queue later without changing the API.
"""
from __future__ import annotations

import threading
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any, Callable, Optional

# Pipeline work is CPU-heavy (and single-process via the in-process Pool
# shim), so keep the pool small. Bump later or swap for a task queue.
_executor = ThreadPoolExecutor(max_workers=2)
_lock = threading.Lock()
_runs: "dict[str, Run]" = {}


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
    return run


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
