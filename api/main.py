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
from fastapi.responses import HTMLResponse

from pipeline import (
    PipelineConfig, run_pipeline, rank_pathways,
    validate_config, cleanup_job_files,
)
from pathway_tools import load_pathways_from_file
from visualize_pathways import visualize_pathways
from enzyme_info import annotate_pathways, enzyme_ids_for_rule
from pathway_scoring import is_bio_op
import uniprot_client
from . import jobs
from .cache import cache, generation_key, ranking_key
from .schemas import GenerateRequest, RankRequest

app = FastAPI(title="TAL Pathway Explorer API", version="0.1.0")

# Shown to the user when a job hits its wall-clock deadline (see jobs.py).
_GEN_TIMEOUT_MSG = (
    "This search took too long and was stopped. The search space is probably "
    "too large — try lowering Generations, Beam size, or Max MW, or narrowing "
    "the whitelist, then run again."
)
_RANK_TIMEOUT_MSG = (
    "Ranking took too long and was stopped. Ranking is single-threaded and "
    "slow on large pathway sets — try tightening the search so fewer pathways "
    "are generated, then rank again."
)

# DORA-XGB feasibility model lives in a separate conda env; spawning it
# takes a few seconds, so we keep one long-running client and reuse it
# across ranking jobs. Lazily created on first use.
_dora_client = None


def _get_dora_client():
    """Return a shared DoraXGBClient, or None if its env isn't set up."""
    global _dora_client
    if _dora_client is None:
        from dora_xgb_client import DoraXGBClient
        _dora_client = DoraXGBClient()
    return _dora_client

# DEV: let the React dev server (a different origin/port) call us.
# Tighten allow_origins before deploying anywhere public.
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Reclaim disk from previous sessions' leftover run files on startup. The
# run store is empty at this point, so every api_* artefact is an orphan.
jobs.sweep_orphan_api_artifacts()


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
    if not req.enable_dora:
        config.feasibility_prune_threshold = None   # prune off unless enabled
    err = validate_config(config)
    if err:
        raise HTTPException(status_code=422, detail=err)

    # ---- generation cache ----
    # If an identical generation already ran (same inputs, canonicalized —
    # see cache.py) AND its pathway file is still on disk, adopt that run's
    # job_name so the graph/rank endpoints resolve against the existing
    # artifacts, and serve its result instantly (no expansion).
    gkey = generation_key(config)
    hit = cache.get_generation(gkey)
    if hit and Path(f"{hit['job_name']}_pathways.txt").exists():
        config.job_name = hit["job_name"]     # point at the cached artifacts
        run.config = config
        run.pathways = hit["pathways"]
        run.diagnostics = {**hit["diagnostics"], "gen_cache_hit": True}
        jobs.set_status(run, "generated")
        return {"run_id": run.id, "status": run.status}

    run.config = config
    jobs.set_status(run, "generating")

    def _worker(r: jobs.Run) -> None:
        cleanup_job_files(config.job_name)
        dora = None
        if req.enable_dora:
            try:
                dora = _get_dora_client()
            except Exception as e:  # env missing → generate without the prune
                print(f"DORA-XGB unavailable, generating without feasibility prune: {e}")
        result = run_pipeline(config, dora_client=dora)
        if not result.ok:
            jobs.complete(r, "error", error=result.error)
            return
        pathways = result.to_dict()["ranked_pathways"]
        annotate_pathways(pathways)   # per-step enzyme counts (bio steps)
        # Cache the (valid) result even if we finished past the deadline — a
        # future identical request then returns instantly and won't time out.
        cache.put_generation(gkey, {
            "pathways": pathways,
            "diagnostics": result.diagnostics,
            "job_name": config.job_name,
        })
        jobs.complete(r, "generated", pathways=pathways, diagnostics=result.diagnostics)

    jobs.run_in_background(run, _worker, timeout=jobs.GENERATION_TIMEOUT_S,
                           timeout_message=_GEN_TIMEOUT_MSG)
    return {"run_id": run.id, "status": run.status}


@app.get("/runs/{run_id}/graph", response_class=HTMLResponse)
def get_graph(run_id: str, view: str = "all"):
    """Render the interactive pathway graph as self-contained HTML.
    view="all" shows every pathway; view="top5" shows the 5 shortest.
    The React Graph tab embeds this in an <iframe>."""
    run = jobs.get_run(run_id)
    if run is None:
        raise HTTPException(status_code=404, detail="run not found")
    if run.status not in ("generated", "ranking", "ranked"):
        raise HTTPException(status_code=409, detail="no pathways to graph yet")
    cfg = run.config

    if view == "top5":
        raw = load_pathways_from_file(cfg.job_name)
        order = sorted(range(1, len(raw) + 1), key=lambda i: raw[i - 1].num_steps)
        pathway_filter = order[:5]
    else:
        pathway_filter = "all"

    op_labels = _bio_op_enzyme_labels(run.ranked_pathways or run.pathways)

    path = visualize_pathways(
        job_name=cfg.job_name,
        starter_smiles=cfg.starter_smiles,
        target_smiles=cfg.target_smiles,
        starter_label="starter",
        target_label="target",
        helpers=cfg.helpers,
        pathway_filter=pathway_filter,
        output_html=f"{cfg.job_name}_graph_{view}.html",
        op_labels=op_labels,
        top_n_threshold=10 ** 9,   # we control the view; disable auto top-N
    )
    with open(path, encoding="utf-8") as f:
        return HTMLResponse(content=f.read())


def _bio_op_enzyme_labels(pathways) -> dict:
    """Map each distinct bio rule in these pathways to a representative
    enzyme name (the first that resolves), for the graph's edge tooltips.
    Cached + degrades gracefully — if UniProt is unreachable, a rule just
    keeps its rule-id label with no enzyme name."""
    rules = {
        name
        for p in (pathways or [])
        for name in p.get("reaction_names", [])
        if is_bio_op(name)
    }
    labels: dict = {}
    for rule in rules:
        ids = enzyme_ids_for_rule(rule)
        if not ids:
            continue
        recs = uniprot_client.resolve(ids[:5])   # one name is enough for a label
        name = next((r["protein_name"] for r in recs if r.get("protein_name")), None)
        if name:
            labels[rule] = name
    return labels


@app.get("/rules/{rule_name}/enzymes")
def rule_enzymes(rule_name: str, limit: int = 25):
    """Resolve the enzymes for one bio rule to human-readable metadata
    (protein name, EC, reaction, gene, organism) via UniProt. On-demand and
    capped, because a rule can list thousands of accessions — the frontend
    calls this when the user expands a bio step."""
    ids = enzyme_ids_for_rule(rule_name)
    limit = max(0, min(limit, 200))
    enzymes = uniprot_client.resolve(ids[:limit])
    return {
        "rule": rule_name,
        "total": len(ids),                 # enzymes annotated on the rule
        "shown": len(enzymes),             # of those, how many resolved in UniProt
        "truncated": len(ids) > limit,
        "enzymes": enzymes,
    }


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

    # ---- ranking cache ----
    # Keyed on the generation key PLUS the weight tiers, so re-ranking the
    # same pathways with the same weights is instant; changing a weight is
    # a miss and re-ranks (without re-expanding).
    rkey = ranking_key(config, req.weights, req.layer_weights, req.lemnisca_weights)
    hit = cache.get_ranking(rkey)
    if hit:
        run.ranked_pathways = hit["ranked_pathways"]
        run.diagnostics = {**run.diagnostics, **hit["diagnostics"], "rank_cache_hit": True}
        jobs.set_status(run, "ranked")
        return {"run_id": run.id, "status": run.status}

    jobs.set_status(run, "ranking")

    def _worker(r: jobs.Run) -> None:
        result = rank_pathways(
            config,
            weights=req.weights,
            layer_weights=req.layer_weights,
            lemnisca_weights=req.lemnisca_weights,
        )
        if not result.ok:
            jobs.complete(r, "error", error=result.error)
            return
        ranked = result.to_dict()["ranked_pathways"]
        annotate_pathways(ranked)     # per-step enzyme counts (bio steps)
        cache.put_ranking(rkey, {
            "ranked_pathways": ranked,
            "diagnostics": result.diagnostics,
        })
        jobs.complete(r, "ranked", ranked_pathways=ranked,
                      diagnostics={**r.diagnostics, **result.diagnostics})

    jobs.run_in_background(run, _worker, timeout=jobs.RANKING_TIMEOUT_S,
                           timeout_message=_RANK_TIMEOUT_MSG)
    return {"run_id": run.id, "status": run.status}
