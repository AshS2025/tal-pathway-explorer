// src/api.ts — thin, typed HTTP client for the TAL Pathway Explorer API.
// Every function here is one "order to the waiter" (an HTTP request to FastAPI).

const BASE = "http://localhost:8000";

// ---- request shapes (mirror api/schemas.py) ----
export interface GenerateRequest {
  starter_smiles: string;
  target_smiles: string;
  domain: string; // "chem" | "bio" | "both"
  direction: string; // "forward" | "retro" | "bidirectional"
  generations: number;
  strategy?: string;
  beam_size?: number;
  max_molecular_weight?: number;
  max_atoms_c?: number;
  max_atoms_o?: number;
  max_atoms_n?: number;
  max_rxn_dh?: number;
  helpers?: string[];
  bio_whitelist?: string[] | null;
  chem_whitelist?: string[] | null;
  // DORA-XGB feasibility prune (bio only), applied during generation
  enable_dora?: boolean;
  feasibility_prune_threshold?: number;
}

export interface RankRequest {
  weights?: Record<string, number> | null; // tier-0 DORAnet internals
  layer_weights?: Record<string, number> | null; // tier-2 DORAnet vs Lemnisca
  lemnisca_weights?: Record<string, number> | null; // tier-1 stability/diversity
}

// ---- response shapes (mirror RankedPathway / Run.to_dict) ----
// NOTE: num_steps is a Python @property, so it's NOT in the JSON — derive it
// on the client as reaction_smiles.length.
export interface Pathway {
  rank: number;
  final_score: number | null; // DORAnet raw composite (null = unranked)
  atomic_economy: number;
  pathway_byproduct_count: number;
  intermediate_byproducts: Record<string, number>;
  reaction_smiles: string[];
  reaction_names: string[];
  reaction_enthalpies: (number | null)[];
  // per-step enzyme count for BIO steps (null for chem steps; 0 = no known
  // enzyme, possibly spontaneous). Parallel to reaction_names.
  reaction_enzymes?: (number | null)[];
  equilibrator_max_dg: number | null;
  equilibrator_avg_dg: number | null;
  equilibrator_coverage: number;
  lemnisca_score: number | null; // Lemnisca sub-score
  blended_score: number | null; // final ranking key
  lemnisca_components: Record<string, number>;
}

export type RunStatus =
  | "pending"
  | "generating"
  | "generated"
  | "ranking"
  | "ranked"
  | "error";

export interface RunState {
  id: string;
  status: RunStatus;
  pathways: Pathway[] | null; // unranked (generation result)
  ranked_pathways: Pathway[] | null; // ranking result
  diagnostics: Record<string, unknown>;
  error: string | null;
}

async function jsonOrThrow(res: Response) {
  if (!res.ok) {
    let detail = res.statusText;
    try {
      detail = (await res.json()).detail ?? detail;
    } catch {
      /* ignore */
    }
    throw new Error(detail);
  }
  return res.json();
}

export async function startRun(
  req: GenerateRequest,
): Promise<{ run_id: string; status: RunStatus }> {
  const res = await fetch(`${BASE}/runs`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(req),
  });
  return jsonOrThrow(res);
}

export async function getRun(id: string): Promise<RunState> {
  return jsonOrThrow(await fetch(`${BASE}/runs/${id}`));
}

export async function startRank(
  id: string,
  req: RankRequest,
): Promise<{ run_id: string; status: RunStatus }> {
  const res = await fetch(`${BASE}/runs/${id}/rank`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(req),
  });
  return jsonOrThrow(res);
}

// The interactive graph HTML lives at this URL (backend endpoint added
// separately); the Graph tab embeds it in an <iframe>.
export function graphUrl(id: string, view: "all" | "top5"): string {
  return `${BASE}/runs/${id}/graph?view=${view}`;
}
