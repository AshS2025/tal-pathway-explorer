import { useEffect, useState } from "react";
import {
  startRun,
  getRun,
  startRank,
  graphUrl,
  type GenerateRequest,
  type RankRequest,
  type RunState,
} from "./api";
import InputsForm from "./components/InputsForm";
import PathwaysView from "./components/PathwaysView";
import RankSidebar from "./components/RankSidebar";
import "./App.css";

type Tab = "inputs" | "pathways" | "graph";

export default function App() {
  const [runId, setRunId] = useState<string | null>(null);
  const [run, setRun] = useState<RunState | null>(null);
  const [tab, setTab] = useState<Tab>("inputs");
  const [err, setErr] = useState<string | null>(null);
  const [graphView, setGraphView] = useState<"all" | "top5">("all");

  const status = run?.status;
  const transient = status === "generating" || status === "ranking";

  // Poll the run while it's still working (generating or ranking).
  useEffect(() => {
    if (!runId || !transient) return;
    const t = setInterval(async () => {
      try {
        setRun(await getRun(runId));
      } catch (e) {
        setErr(String(e));
      }
    }, 2500);
    return () => clearInterval(t);
  }, [runId, transient]);

  async function handleRun(req: GenerateRequest) {
    setErr(null);
    try {
      const { run_id } = await startRun(req);
      setRunId(run_id);
      setRun({
        id: run_id,
        status: "generating",
        pathways: null,
        ranked_pathways: null,
        diagnostics: {},
        error: null,
      });
      setTab("pathways");
    } catch (e) {
      setErr(String(e));
    }
  }

  async function handleRank(req: RankRequest) {
    if (!runId) return;
    setErr(null);
    try {
      await startRank(runId, req);
      setRun((r) => (r ? { ...r, status: "ranking" } : r));
      setTab("pathways");
    } catch (e) {
      setErr(String(e));
    }
  }

  const pathways = run?.ranked_pathways ?? run?.pathways ?? [];
  const isRanked = !!run?.ranked_pathways;

  return (
    <div className="app">
      <RankSidebar
        onRank={handleRank}
        disabled={!run || status === "generating"}
        busy={status === "ranking"}
      />

      <main className="main">
        <h1>⚗️ TAL Pathway Explorer</h1>
        {err && <div className="error">{err}</div>}

        <nav className="tabs">
          <button className={tab === "inputs" ? "active" : ""} onClick={() => setTab("inputs")}>
            Inputs
          </button>
          {run && (
            <button
              className={tab === "pathways" ? "active" : ""}
              onClick={() => setTab("pathways")}
            >
              Pathways
            </button>
          )}
          {run && (
            <button className={tab === "graph" ? "active" : ""} onClick={() => setTab("graph")}>
              Graph
            </button>
          )}
        </nav>

        <section className="content">
          {tab === "inputs" && <InputsForm onRun={handleRun} busy={status === "generating"} />}

          {tab === "pathways" && run && (
            <>
              {status === "generating" && <p className="muted">Generating pathways…</p>}
              {status === "error" && <div className="error">{run.error}</div>}
              {status === "ranking" && (
                <p className="muted">Ranking in progress… (showing unranked meanwhile)</p>
              )}
              {Number(run.diagnostics?.feasibility_pruned) > 0 && (
                <p className="muted">
                  DORA-XGB pruned {String(run.diagnostics.feasibility_pruned)} infeasible
                  bio pathway(s) during generation.
                </p>
              )}
              {(run.diagnostics?.gen_cache_hit || run.diagnostics?.rank_cache_hit) && (
                <p className="muted">
                  ⚡ Served from cache — identical inputs
                  {run.diagnostics?.rank_cache_hit ? " and weights" : ""}, no
                  re-{run.diagnostics?.rank_cache_hit ? "ranking" : "expansion"} needed.
                </p>
              )}
              {(status === "generated" || status === "ranking" || status === "ranked") && (
                <PathwaysView pathways={pathways} ranked={isRanked} />
              )}
            </>
          )}

          {tab === "graph" && run && (
            <div className="graph">
              <div className="graph-toggle">
                <button
                  className={graphView === "all" ? "active" : ""}
                  onClick={() => setGraphView("all")}
                >
                  All pathways
                </button>
                <button
                  className={graphView === "top5" ? "active" : ""}
                  onClick={() => setGraphView("top5")}
                >
                  Top 5
                </button>
              </div>
              {status === "generated" || status === "ranked" ? (
                <iframe
                  className="graph-frame"
                  title="pathway graph"
                  src={graphUrl(run.id, graphView)}
                />
              ) : (
                <p className="muted">Graph appears once generation finishes.</p>
              )}
            </div>
          )}
        </section>
      </main>
    </div>
  );
}
