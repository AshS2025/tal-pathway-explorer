import { useState } from "react";
import type { RankRequest } from "../api";

function Weight({
  label,
  value,
  onChange,
}: {
  label: string;
  value: number;
  onChange: (v: number) => void;
}) {
  return (
    <label className="weight">
      <span>{label}</span>
      <input
        type="number"
        min={0}
        max={10}
        value={value}
        onChange={(e) => onChange(Number(e.target.value))}
      />
    </label>
  );
}

interface Props {
  onRank: (req: RankRequest) => void;
  disabled: boolean; // no run generated yet
  busy: boolean; // a rank is in progress
}

export default function RankSidebar({ onRank, disabled, busy }: Props) {
  // tier 2 — layer blend
  const [wDoranet, setWDoranet] = useState(2);
  const [wLemnisca, setWLemnisca] = useState(1);
  // tier 1 — lemnisca components
  const [wStability, setWStability] = useState(1);
  const [wDiversity, setWDiversity] = useState(1);
  const [wFeasibility, setWFeasibility] = useState(1);
  const [enableDora, setEnableDora] = useState(false);
  // tier 0 — DORAnet internals
  const [wSteps, setWSteps] = useState(4);
  const [wThermo, setWThermo] = useState(2);
  const [wByprod, setWByprod] = useState(2);
  const [wAtom, setWAtom] = useState(1);

  function submit() {
    const lemnisca: Record<string, number> = {
      stability: wStability,
      diversity: wDiversity,
    };
    if (enableDora) lemnisca.feasibility = wFeasibility;
    const req: RankRequest = {
      layer_weights: { doranet: wDoranet, lemnisca: wLemnisca },
      lemnisca_weights: lemnisca,
      weights: {
        number_of_steps: wSteps,
        reaction_thermo: wThermo,
        by_product_number: wByprod,
        atom_economy: wAtom,
      },
      enable_dora: enableDora,
    };
    onRank(req);
  }

  return (
    <aside className="sidebar">
      <h2>🏆 Ranking</h2>
      {disabled ? (
        <p className="muted">Run a search first, then rank the pathways here.</p>
      ) : (
        <>
          <h4>Layer blend</h4>
          <Weight label="DORAnet (chemistry)" value={wDoranet} onChange={setWDoranet} />
          <Weight label="Lemnisca (viability)" value={wLemnisca} onChange={setWLemnisca} />

          <h4>Lemnisca components</h4>
          <Weight label="Stability" value={wStability} onChange={setWStability} />
          <Weight label="Diversity" value={wDiversity} onChange={setWDiversity} />
          <label className="weight">
            <span>Feasibility (DORA-XGB)</span>
            <input
              type="checkbox"
              checked={enableDora}
              onChange={(e) => setEnableDora(e.target.checked)}
            />
          </label>
          {enableDora && (
            <Weight label="Feasibility weight" value={wFeasibility} onChange={setWFeasibility} />
          )}

          <details className="advanced">
            <summary>DORAnet internal weights</summary>
            <Weight label="Steps" value={wSteps} onChange={setWSteps} />
            <Weight label="Thermo" value={wThermo} onChange={setWThermo} />
            <Weight label="By-products" value={wByprod} onChange={setWByprod} />
            <Weight label="Atom economy" value={wAtom} onChange={setWAtom} />
          </details>

          <button className="primary" onClick={submit} disabled={busy}>
            {busy ? "Ranking…" : "Rank pathways"}
          </button>
        </>
      )}
    </aside>
  );
}
