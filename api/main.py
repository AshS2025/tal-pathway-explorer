"""
api/main.py — FastAPI backend ("the waiter") over src/pipeline.py.

Run from the PROJECT ROOT:
    uvicorn api.main:app --port 8000
Interactive docs once running: http://localhost:8000/docs

Endpoints:
    GET  /health               liveness
    POST /runs                 start a generation job         -> {run_id, status}
    GET  /runs/{run_id}         poll a run's full state
    POST /runs/{run_id}/rank    start a ranking job on the run -> {run_id, status}

Generation and ranking are slow, so they run in the background (see
api/jobs.py) and the client polls GET /runs/{id} for status + results.
"""
from __future__ import annotations

import sys
from pathlib import Path

# --- UTF-8 stdout: DORAnet prints unicode arrows; Windows cp1252 crashes. ---
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except (AttributeError, ValueError):
    pass

# --- make src/ importable ---
_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / "src"))

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

# --- Windows multiprocessing safety ---------------------------------------
# Replace multiprocessing.Pool with an in-process shim BEFORE doranet is
# imported, so DORAnet's ranking never spawns subprocesses that re-import
# this module (the same fix the Streamlit app uses).
import multiprocessing as _mp


class _InProcessAsyncResult:
    def __init__(self, value): self._value = value
    def get(self, timeout=None): return self._value
    def wait(self, timeout=None): pass
    def ready(self): return True
    def successful(self): return True


class _InProcessPool:
    def __init__(self, *a, **k): pass
    def __enter__(self): return self
    def __exit__(self, *a): pass
    def map(self, f, it): return [f(x) for x in it]
    def imap(self, f, it):
        for x in it:
            yield f(x)
    def imap_unordered(self, f, it): return self.imap(f, it)
    def apply(self, f, args=(), kwds=None): return f(*args, **(kwds or {}))
    def apply_async(self, func, args=(), kwds=None, callback=None,
                    error_callback=None):
        try:
            r = func(*args, **(kwds or {}))
            if callback:
                callback(r)
            return _InProcessAsyncResult(r)
        except Exception as e:
            if error_callback:
                error_callback(e)
            raise
    def starmap(self, f, it): return [f(*a) for a in it]
    def close(self): pass
    def join(self): pass
    def terminate(self): pass


_mp.Pool = _InProcessPool
from doranet.modules.post_processing import post_processing as _dpp
_dpp.Pool = _InProcessPool
# --------------------------------------------------------------------------

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from pipeline import (
    PipelineConfig, run_pipeline, rank_pathways,
    validate_config, cleanup_job_files,
)
from . import jobs
from .schemas import GenerateRequest, RankRequest

app = FastAPI(title="TAL Pathway Explorer API", version="0.1.0")

# DEV: let the React dev server (a different origin/port) call us.
# Tighten allow_origins before deploying anywhere public.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/runs")
def start_run(req: GenerateRequest):
    """Start a generation job. The run's UUID doubles as its pipeline
    job_name, so concurrent runs never clobber each other's files."""
    run = jobs.create_run()
    config = PipelineConfig.from_dict(
        {**req.model_dump(), "job_name": f"api_{run.id}"}
    )
    err = validate_config(config)
    if err:
        raise HTTPException(status_code=422, detail=err)
    run.config = config
    jobs.set_status(run, "generating")

    def _worker(r: jobs.Run) -> None:
        cleanup_job_files(config.job_name)
        result = run_pipeline(config)
        if not result.ok:
            r.status, r.error = "error", result.error
            return
        r.pathways = result.to_dict()["ranked_pathways"]
        r.diagnostics = result.diagnostics
        r.status = "generated"

    jobs.run_in_background(run, _worker)
    return {"run_id": run.id, "status": run.status}


@app.get("/runs/{run_id}")
def get_run(run_id: str):
    """Poll a run's full state: status, generated pathways, ranked
    pathways, and any error."""
    run = jobs.get_run(run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="run not found")
    return run.to_dict()


@app.post("/runs/{run_id}/rank")
def start_rank(run_id: str, req: RankRequest):
    """Start a ranking job on an already-generated run, using the three
    weight tiers (all optional)."""
    run = jobs.get_run(run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="run not found")
    if run.status not in ("generated", "ranked"):
        raise HTTPException(
            status_code=409,
            detail=f"run not ready to rank (status={run.status})",
        )
    config = run.config
    jobs.set_status(run, "ranking")

    def _worker(r: jobs.Run) -> None:
        result = rank_pathways(
            config,
            weights=req.weights,
            layer_weights=req.layer_weights,
            lemnisca_weights=req.lemnisca_weights,
        )
        if not result.ok:
            r.status, r.error = "error", result.error
            return
        r.ranked_pathways = result.to_dict()["ranked_pathways"]
        r.diagnostics = {**r.diagnostics, **result.diagnostics}
        r.status = "ranked"

    jobs.run_in_background(run, _worker)
    return {"run_id": run.id, "status": run.status}
