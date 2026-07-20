import { useState } from "react";
import type { GenerateRequest } from "../api";

// Splits a textarea into a clean list of lines (SMILES / rule names).
function lines(text: string): string[] {
  return text
    .split("\n")
    .map((l) => l.trim())
    .filter(Boolean);
}

interface Props {
  onRun: (req: GenerateRequest) => void;
  busy: boolean;
}

export default function InputsForm({ onRun, busy }: Props) {
  const [starter, setStarter] = useState("Cc1cc(O)cc(=O)o1"); // TAL
  const [target, setTarget] = useState("CC=CC=CC(=O)O"); // sorbic acid
  const [domain, setDomain] = useState("chem");
  const [direction, setDirection] = useState("bidirectional");
  const [generations, setGenerations] = useState(3);

  // advanced
  const [strategy, setStrategy] = useState("priority_queue");
  const [beamSize, setBeamSize] = useState(1000);
  const [maxMW, setMaxMW] = useState(200);
  const [maxC, setMaxC] = useState(10);
  const [maxO, setMaxO] = useState(5);
  const [maxN, setMaxN] = useState(2);
  const [maxDH, setMaxDH] = useState(15);
  const [helpers, setHelpers] = useState("O\n[H][H]");
  const [chemWhitelist, setChemWhitelist] = useState("");
  const [bioWhitelist, setBioWhitelist] = useState("");
  const [enableDora, setEnableDora] = useState(false);
  const [feasThreshold, setFeasThreshold] = useState(0.5);

  const includeChem = domain === "chem" || domain === "both";
  const includeBio = domain === "bio" || domain === "both";

  function submit() {
    const req: GenerateRequest = {
      starter_smiles: starter.trim(),
      target_smiles: target.trim(),
      domain,
      direction,
      generations,
      strategy,
      beam_size: beamSize,
      max_molecular_weight: maxMW,
      max_atoms_c: maxC,
      max_atoms_o: maxO,
      max_atoms_n: maxN,
      max_rxn_dh: maxDH,
      helpers: lines(helpers).length ? lines(helpers) : ["O", "[H][H]"],
      chem_whitelist: includeChem && chemWhitelist.trim() ? lines(chemWhitelist) : null,
      bio_whitelist: includeBio && bioWhitelist.trim() ? lines(bioWhitelist) : null,
      enable_dora: enableDora,
      feasibility_prune_threshold: feasThreshold,
    };
    onRun(req);
  }

  return (
    <div className="form">
      <label>
        Starter SMILES
        <input value={starter} onChange={(e) => setStarter(e.target.value)} />
      </label>
      <label>
        Target SMILES
        <input value={target} onChange={(e) => setTarget(e.target.value)} />
      </label>

      <div className="row">
        <label>
          Domain
          <select value={domain} onChange={(e) => setDomain(e.target.value)}>
            <option value="chem">Chem</option>
            <option value="bio">Bio</option>
            <option value="both">Both</option>
          </select>
        </label>
        <label>
          Direction
          <select value={direction} onChange={(e) => setDirection(e.target.value)}>
            <option value="bidirectional">Bidirectional</option>
            <option value="forward">Forward</option>
            <option value="retro">Retro</option>
          </select>
        </label>
        <label>
          Generations
          <input
            type="number"
            min={1}
            max={6}
            value={generations}
            onChange={(e) => setGenerations(Number(e.target.value))}
          />
        </label>
      </div>

      <details className="advanced">
        <summary>Advanced</summary>

        <div className="row">
          <label>
            Strategy
            <select value={strategy} onChange={(e) => setStrategy(e.target.value)}>
              <option value="priority_queue">Priority queue</option>
              <option value="cartesian">Cartesian</option>
            </select>
          </label>
          <label>
            Beam size
            <input
              type="number"
              value={beamSize}
              onChange={(e) => setBeamSize(Number(e.target.value))}
            />
          </label>
        </div>

        <div className="row">
          <label>
            Max MW
            <input type="number" value={maxMW} onChange={(e) => setMaxMW(Number(e.target.value))} />
          </label>
          <label>
            Max C
            <input type="number" value={maxC} onChange={(e) => setMaxC(Number(e.target.value))} />
          </label>
          <label>
            Max O
            <input type="number" value={maxO} onChange={(e) => setMaxO(Number(e.target.value))} />
          </label>
          <label>
            Max N
            <input type="number" value={maxN} onChange={(e) => setMaxN(Number(e.target.value))} />
          </label>
          <label>
            Max |ΔH|
            <input type="number" value={maxDH} onChange={(e) => setMaxDH(Number(e.target.value))} />
          </label>
        </div>

        <label>
          Helper molecules (one SMILES per line)
          <textarea rows={2} value={helpers} onChange={(e) => setHelpers(e.target.value)} />
        </label>

        {includeChem && (
          <label>
            Chem reaction whitelist (blank = built-in TAL default)
            <textarea
              rows={3}
              value={chemWhitelist}
              onChange={(e) => setChemWhitelist(e.target.value)}
              placeholder="One operator name per line"
            />
          </label>
        )}
        {includeBio && (
          <label>
            Bio rule whitelist (blank = built-in default)
            <textarea
              rows={3}
              value={bioWhitelist}
              onChange={(e) => setBioWhitelist(e.target.value)}
              placeholder="rule1118&#10;rule0087&#10;rule0891"
            />
          </label>
        )}

        <label className="checkbox-row">
          <input
            type="checkbox"
            checked={enableDora}
            onChange={(e) => setEnableDora(e.target.checked)}
          />
          Prune infeasible bio pathways with DORA-XGB
        </label>
        {enableDora && (
          <label>
            Feasibility threshold — drop pathways with any bio step below this (0–1)
            <input
              type="number"
              min={0}
              max={1}
              step={0.05}
              value={feasThreshold}
              onChange={(e) => setFeasThreshold(Number(e.target.value))}
            />
          </label>
        )}
      </details>

      <button className="primary" onClick={submit} disabled={busy}>
        {busy ? "Running…" : "Run"}
      </button>
    </div>
  );
}
