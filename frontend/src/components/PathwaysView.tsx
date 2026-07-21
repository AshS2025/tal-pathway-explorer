import type { Pathway } from "../api";

const fmt = (x: number | null | undefined) =>
  x === null || x === undefined ? "—" : x.toFixed(2);

const comp = (p: Pathway, key: string) => {
  const v = p.lemnisca_components?.[key];
  return v === undefined ? "—" : v.toFixed(2);
};

const trunc = (s: string, n = 42) => (s.length <= n ? s : s.slice(0, n - 3) + "...");

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
              <span className="enz"> · {enz} enzyme{enz === 1 ? "" : "s"}</span>
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
