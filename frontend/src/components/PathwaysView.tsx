import { useState } from "react";
import { fetchRuleEnzymes, type Pathway, type RuleEnzymes } from "../api";

const fmt = (x: number | null | undefined) =>
  x === null || x === undefined ? "—" : x.toFixed(2);

const comp = (p: Pathway, key: string) => {
  const v = p.lemnisca_components?.[key];
  return v === undefined ? "—" : v.toFixed(2);
};

const trunc = (s: string, n = 42) => (s.length <= n ? s : s.slice(0, n - 3) + "...");

// Expandable enzyme list for one bio step. Lazily fetches from UniProt (via
// our backend) the first time it's opened, so we only pay for what's viewed.
function EnzymeList({ rule, count }: { rule: string; count: number }) {
  const [open, setOpen] = useState(false);
  const [data, setData] = useState<RuleEnzymes | null>(null);
  const [loading, setLoading] = useState(false);
  const [err, setErr] = useState<string | null>(null);

  async function toggle() {
    const next = !open;
    setOpen(next);
    if (next && !data && !loading) {
      setLoading(true);
      setErr(null);
      try {
        setData(await fetchRuleEnzymes(rule));
      } catch (e) {
        setErr(e instanceof Error ? e.message : String(e));
      } finally {
        setLoading(false);
      }
    }
  }

  return (
    <>
      <button className="enz" type="button" onClick={toggle}>
        {" · "}
        {count} enzyme{count === 1 ? "" : "s"} {open ? "▾" : "▸"}
      </button>
      {open && (
        <div className="enz-panel">
          {loading && <span className="muted">Loading from UniProt…</span>}
          {err && <span className="error">{err}</span>}
          {data && (
            <>
              {data.shown < data.total && (
                <p className="muted">
                  Showing {data.shown} of {data.total} annotated enzymes
                  {data.shown === 0 ? " (none resolved in UniProt)" : ""}.
                </p>
              )}
              {data.enzymes.length > 0 && (
                <table className="enz-table">
                  <thead>
                    <tr>
                      <th>Enzyme</th>
                      <th>EC</th>
                      <th>Reaction</th>
                      <th>Gene</th>
                      <th>Organism</th>
                    </tr>
                  </thead>
                  <tbody>
                    {data.enzymes.map((e) => (
                      <tr key={e.accession}>
                        <td>{e.protein_name || e.accession}</td>
                        <td>{e.ec.join(", ") || "—"}</td>
                        <td>
                          {e.reactions[0] ?? "—"}
                          {e.reactions_truncated ? " …" : ""}
                        </td>
                        <td>{e.gene || "—"}</td>
                        <td>
                          <i>{e.organism || "—"}</i>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              )}
            </>
          )}
        </div>
      )}
    </>
  );
}

function Steps({ p }: { p: Pathway }) {
  return (
    <ol className="steps">
      {p.reaction_smiles.map((smi, i) => {
        const [lhs, rhs] = smi.split(">>");
        const dh = p.reaction_enthalpies[i];
        const enz = p.reaction_enzymes?.[i];   // null = chem step; 0 = bio, no enzyme
        return (
          <li key={i}>
            <span className="op">{p.reaction_names[i]}</span>
            {dh !== null && dh !== undefined ? (
              <span className="dh"> ΔH={dh.toFixed(1)}</span>
            ) : null}
            {enz === null || enz === undefined ? null : enz === 0 ? (
              <span className="enz-none"> · no enzyme (possibly spontaneous)</span>
            ) : (
              <EnzymeList rule={p.reaction_names[i]} count={enz} />
            )}
            <div className="rxn">
              {(lhs ?? "").split(".").map(trunc).join(" + ")} <b>→</b>{" "}
              {(rhs ?? "").split(".").map(trunc).join(" + ")}
            </div>
          </li>
        );
      })}
    </ol>
  );
}

interface Props {
  pathways: Pathway[];
  ranked: boolean;
}

export default function PathwaysView({ pathways, ranked }: Props) {
  if (!pathways.length) return <p className="muted">No pathways.</p>;

  return (
    <div>
      {ranked ? (
        <>
          <p className="muted">
            Ranked by the final blended score = geomean(DORAnet chemistry,
            Lemnisca viability). A catastrophic intermediate gates a route to 0.
          </p>
          <table className="grid">
            <thead>
              <tr>
                <th>Rank</th>
                <th>Final</th>
                <th>DORAnet</th>
                <th>Lemnisca</th>
                <th>Stability</th>
                <th>Diversity</th>
                <th>Steps</th>
              </tr>
            </thead>
            <tbody>
              {pathways.map((p) => (
                <tr key={p.rank}>
                  <td>{p.rank}</td>
                  <td>{fmt(p.blended_score)}</td>
                  <td>{comp(p, "doranet")}</td>
                  <td>{fmt(p.lemnisca_score)}</td>
                  <td>{comp(p, "stability")}</td>
                  <td>{comp(p, "diversity")}</td>
                  <td>{p.reaction_smiles.length}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </>
      ) : (
        <p className="muted">
          Unranked — sorted by step count. Set weights in the sidebar and click
          Rank to score them.
        </p>
      )}

      <h3>Details</h3>
      {pathways.map((p, idx) => (
        <details key={p.rank} open={idx === 0}>
          <summary>
            {ranked
              ? `Rank ${p.rank} — final ${fmt(p.blended_score)} — ${p.reaction_smiles.length} steps`
              : `#${p.rank} — ${p.reaction_smiles.length} steps`}
          </summary>
          <Steps p={p} />
        </details>
      ))}
    </div>
  );
}
