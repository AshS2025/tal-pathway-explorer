"""
api/schemas.py — request models = the API's input contract.

These mirror the fields of pipeline.PipelineConfig (generation) and the
three weight tiers (ranking). Responses are the pipeline's own
`to_dict()` output, so they are NOT re-typed here — the dataclasses stay
the single source of truth.
"""
from __future__ import annotations

from typing import Optional

from pydantic import BaseModel, Field


class GenerateRequest(BaseModel):
    """Inputs for the generation phase (POST /runs)."""
    starter_smiles: str
    target_smiles: str
    domain: str = "chem"                 # "chem" | "bio" | "both"
    direction: str = "bidirectional"     # "forward" | "retro" | "bidirectional"
    generations: int = 3
    strategy: str = "priority_queue"     # "priority_queue" | "cartesian"
    beam_size: int = 1000
    max_molecular_weight: float = 200.0
    max_atoms_c: int = 10
    max_atoms_o: int = 5
    max_atoms_n: int = 2
    max_rxn_dh: float = 15.0
    helpers: list[str] = Field(default_factory=lambda: ["O", "[H][H]"])
    bio_whitelist: Optional[list[str]] = None
    chem_whitelist: Optional[list[str]] = None
    enable_rmg: bool = False
    enable_equilibrator: bool = False
    equilibrator_prune_max_abs_dg: float = 100.0
    # DORA-XGB feasibility prune (bio only): drop pathways whose weakest
    # bio step is below this feasibility. Only applied when enable_dora.
    enable_dora: bool = False
    feasibility_prune_threshold: float = 0.5


class RankRequest(BaseModel):
    """Inputs for the ranking phase (POST /runs/{id}/rank). Each weight
    dict is optional; omitted → the pipeline's built-in defaults."""
    weights: Optional[dict] = None            # tier-0: DORAnet internals
    layer_weights: Optional[dict] = None      # tier-2: DORAnet vs Lemnisca
    lemnisca_weights: Optional[dict] = None   # tier-1: stability, diversity
