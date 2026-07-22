import { useEffect, useState } from "react";
import { fetchRuleEnzymes, type Pathway, type RuleEnzymes } from "../api";

const fmt = (x: number | null | undefined) =>
  x === null || x === undefined ? "—" : x.toFixed(2);

const comp = (p: Pathway, key: string) => {
  const v = p.lemnisca_components?.[key];
  return v === undefined ? "—" : v.toFixed(2);
};

// Minimum distinct enzymes to build the route (shared multifunctional
// enzymes counted once). Falls back to the bio-step count for older data.
const enzymeCount = (p: Pathway) =>
  p.min_enzymes ??
  (p.reaction_enzymes ?? []).filter((e) => e !== null).length;

const trunc = (s: string, n = 42) => (s.length <= n ? s : s.slice(0, n - 3) + "...");

// Enzyme metadata table for one bio rule. Fetches from UniProt (via our
// backend) when mounted — and it's only mounted once its step is expanded,
// so we pay the lookup only for steps the user actually opens.
function EnzymeTable({ rule }: { rule: string }) {
  const [data, setData] = useState<RuleEnzymes | null>(null);
  const [loading, setLoading] = useState(true);
  const [err, setErr] = useState<string | null>(null);

  useEffect(() => {
    let alive = true;
    setLoading(true);
    setErr(null);
    fetchRuleEnzymes(rule)
      .then((d) => alive && setData(d))
      .catch((e) => alive && setErr(e instanceof Error ? e.message : String(e)))
      .finally(() => alive && setLoading(false));
    return () => {
      alive = false;
    };
  }, [rule]);

  if (loading) return <p className="muted">Loading enzymes from UniProt…</p>;
  if (err) return <p className="error">{err}</p>;
  if (!data) return null;

  return (
    <>
      {data.enzymes.length > 0 && (
        <table className="enz-table">
          <thead>
            <tr>
              <th>Enzyme</th>
              <th>EC</th>
              <th>MW (kDa)</th>
              <th>Reaction</th>
              <th>Gene</th>
              <th>Organism</th>
            </tr>
          </thead>
          <tbody>
            {data.enzymes.map((e) => {
              const uniprotUrl = `https://www.uniprot.org/uniprotkb/${e.accession}/entry`;
              return (
                <tr key={e.accession}>
                  <td>
                    {e.deleted ? (
                      <>
                        {e.accession}{" "}
                        <span className="muted">(deleted in UniProt)</span>
                      </>
                    ) : (
                      e.protein_name || e.accession
                    )}
                  </td>
                  <td>
                    {e.ec.length ? e.ec.join(", ") : <span className="muted">no EC</span>}
                  </td>
                  <td>
                    {e.mass ? (e.mass / 1000).toFixed(1) : <span className="muted">—</span>}
                  </td>
                  <td>
                    {e.reaction_count === 0 ? (
                      <span className="muted">no reaction</span>
                    ) : e.reaction_count <= 3 ? (
                      e.reactions.map((r, i) => <div key={i}>{r}</div>)
                    ) : (
                      <a href={uniprotUrl} target="_blank" rel="noreferrer">
                        {e.reaction_count} reactions — view in UniProt ↗
                      </a>
                    )}
                  </td>
                  <td>{e.gene || "—"}</td>
                  <td>
                    <i>{e.organism || "—"}</i>
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      )}
      {data.total > data.shown && (
        <p className="muted">
          Showing the first {data.shown} of {data.total} annotated enzymes.{" "}
          <a href={data.uniprot_url} target="_blank" rel="noreferrer">
            View all in UniProt ↗
          </a>
        </p>
      )}
      {data.shown === 0 && (
        <p className="muted">No enzymes resolved in UniProt for this step.</p>
      )}
    </>
  );
}

// One step in a pathway — itself a dropdown. Collapsed shows the operator +
// enzyme badge; expanded shows the reaction and (for bio steps) the UniProt
// enzyme table, fetched lazily on first open.
function StepItem({
  name,
  smi,
  dh,
  enz,
}: {
  name: string;
  smi: string;
  dh: number | null;
  enz: number | null | undefined;
}) {
  const [open, setOpen] = useState(false);
  const [lhs, rhs] = smi.split(">>");
  const isChem = enz === null || enz === undefined;
  const isBioNoEnzyme = enz === 0;

  return (
    <details
      className="step"
      onToggle={(e) => setOpen((e.currentTarget as HTMLDetailsElement).open)}
    >
      <summary>
        <span className="op">{name}</span>
        {dh !== null && dh !== undefined ? (
          <span className="dh"> ΔH={dh.toFixed(1)}</span>
        ) : null}
        {isChem ? null : isBioNoEnzyme ? (
          <span className="enz-none"> · no enzyme (possibly spontaneous)</span>
        ) : (
          <span className="enz"> · {enz} enzyme{enz === 1 ? "" : "s"}</span>
        )}
      </summary>
      <div className="step-body">
        <div className="rxn">
          {(lhs ?? "").split(".").map(trunc).join(" + ")} <b>→</b>{" "}
          {(rhs ?? "").split(".").map(trunc).join(" + ")}
        </div>
        {open &&
          (isChem ? (
            <p className="muted">Chemical operator — no enzyme data.</p>
          ) : isBioNoEnzyme ? (
            <p className="enz-none">No known enzyme; this step may be spontaneous.</p>
          ) : (
            <EnzymeTable rule={name} />
          ))}
      </div>
    </details>
  );
}

function Steps({ p }: { p: Pathway }) {
  return (
    <div className="steps">
      {p.reaction_smiles.map((smi, i) => (
        <StepItem
          key={i}
          name={p.reaction_names[i]}
          smi={smi}
          dh={p.reaction_enthalpies[i]}
          enz={p.reaction_enzymes?.[i]}
        />
      ))}
    </div>
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
                <th>Enz. load</th>
                <th>Steps</th>
                <th>Min. enzymes</th>
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
                  <td>{comp(p, "enzyme_load")}</td>
                  <td>{p.reaction_smiles.length}</td>
                  <td>{enzymeCount(p)}</td>
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
