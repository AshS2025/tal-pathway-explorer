"""
pipeline.py
===========

UI-agnostic orchestration for the TAL pathway-search pipeline.

This module is the ONE place the whole pipeline lives — validation,
direction handling, bio-whitelist parsing, expansion,
network merge, pathway trace, DORAnet ranking, equilibrator
decoration + pruning. Both `streamlit_app.py` (today's frontend) and
any future frontend (FastAPI/React, CLI, notebook) import this same
module and get identical behaviour.

The pipeline is expressed via two dataclasses:

    PipelineConfig — inputs (SMILES, direction, limits, whitelists...)
    PipelineResult — outputs (ranked pathways, timing, error, paths)

Both are dataclass instances but serializable to plain dicts via
`config.to_dict()` / `result.to_dict()` so a REST layer can accept /
emit JSON without extra glue.

Client lifecycle (RMG / equilibrator) is NOT owned here — callers pass
already-constructed clients in as arguments. That lets Streamlit cache
them via `@st.cache_resource`, a FastAPI backend keep them as module-
level singletons, and a CLI create-and-close them per invocation.
"""

from __future__ import annotations

import dataclasses
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, List, Optional

from rdkit import Chem

import csv

import doranet
from network_generation import generate_network_tal
from pathway_tools import load_pathways_from_file, parse_reaction_string
from pathway_scoring import (
    RankedPathway, generate_base_rankings, decorate_with_equilibrator,
)
from recipe_rankers import (
    FeedstockProximityRanker, ForwardProductTanimotoRanker,
)
from doranet.modules.post_processing.post_processing import (
    pathway_finder, pretreat_networks,
)
# Native DORAnet enzymatic expansion. We call this DIRECTLY for bio
# instead of going through generate_network_tal's engine path: native
# treats cofactors as operator SLOTS (bounded ~100 reactions), whereas
# the wrapper's engine path with multi-reactant enabled lets all ~44
# cofactors free-combine and the network explodes to millions of
# reactions. The native path is what produced the working Jun-30 data.
from doranet.modules.enzymatic.generate_network import (
    generate_network as _native_bio_generate_network,
    AVAILABLE_RULESETS as _BIO_AVAILABLE_RULESETS,
)


DEFAULT_BIO_WHITELIST = ("rule1118", "rule0087", "rule0891")

# Malonyl-CoA — the polyketide chain-extender unit consumed by the
# Claisen condensation rules (rule1118, rule0087). It is NOT in DORAnet's
# default cofactor table (only acetyl-CoA is), so bio expansion produces
# ZERO reactions unless malonyl-CoA is supplied some other way. We inject
# it as a freely-available HELPER whenever bio operators are enabled — the
# same thing the working tal_centered_combined recipe did. Without this,
# no bio pathway ever connects.
MALONYL_COA = (
    "OC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)"
    "OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)(O)O"
)


# ============================================================
# Config / result dataclasses
# ============================================================
@dataclass
class PipelineConfig:
    """Everything the pipeline needs to run, expressed frontend-agnostic.

    Any UI (Streamlit, React, CLI) that can populate this dataclass can
    drive the pipeline. Keep fields JSON-serializable — no callables,
    no Numpy arrays, no runtime objects.
    """
    # required
    starter_smiles: str
    target_smiles: str
    # domain: which operator sets fire — "chem", "bio", or "both".
    # This is independent of search direction (below).
    domain: str = "chem"
    # direction of search: "forward" (expand from starter toward target),
    # "retro" (expand back from target toward starter), or "bidirectional"
    # (both — meet-in-the-middle). NO starter/target swap is applied; the
    # molecules are always used exactly as entered.
    direction: str = "bidirectional"
    # expansion depth per side of the bidirectional search
    generations: int = 3
    # "cartesian" (exhaustive) or "priority_queue" (target-guided beam)
    strategy: str = "priority_queue"
    beam_size: int = 1000
    # atom/mass/thermo caps applied during expansion
    max_molecular_weight: float = 200.0
    max_atoms_c: int = 10
    max_atoms_o: int = 5
    max_atoms_n: int = 2
    max_rxn_dh: float = 15.0
    # operator whitelists — if None, each falls back to its built-in TAL
    # default inside network_generation. When supplied (e.g. from a UI
    # upload), the list FULLY REPLACES the default.
    bio_whitelist: Optional[List[str]] = None
    chem_whitelist: Optional[List[str]] = None
    # thermo toggles
    enable_rmg: bool = False
    enable_equilibrator: bool = False
    # equilibrator post-hoc pruning threshold (kJ/mol); pathways with
    # any bio step |ΔG'°| exceeding this get dropped
    equilibrator_prune_max_abs_dg: float = 100.0
    # file prefix for temp artefacts written to disk during a run
    job_name: str = "pipeline_job"

    # -------- helpers ---------------------------------------------
    @property
    def max_atoms(self) -> dict:
        """Element-symbol → cap dict, in the shape generate_network_tal expects."""
        return {"C": self.max_atoms_c,
                "O": self.max_atoms_o,
                "N": self.max_atoms_n}

    @property
    def include_chem(self) -> bool:
        return self.domain in ("chem", "both")

    @property
    def include_bio(self) -> bool:
        return self.domain in ("bio", "both")

    @property
    def directions(self) -> list:
        """Which expansion directions to run. 'bidirectional' expands
        from both ends (meet-in-the-middle); the others run one side."""
        if self.direction == "bidirectional":
            return ["forward", "retro"]
        return [self.direction]

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "PipelineConfig":
        # tolerate unknown keys (e.g. UI-only fields) by filtering to known ones
        known = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in known})


@dataclass
class PipelineResult:
    """Structured output. Serializable for JSON responses."""
    ok: bool
    error: Optional[str] = None
    elapsed_seconds: float = 0.0
    n_pathways: int = 0
    ranked_pathways: List[RankedPathway] = field(default_factory=list)
    pathway_file_path: Optional[str] = None
    diagnostics: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "ok": self.ok,
            "error": self.error,
            "elapsed_seconds": self.elapsed_seconds,
            "n_pathways": self.n_pathways,
            "pathway_file_path": self.pathway_file_path,
            "diagnostics": self.diagnostics,
            "ranked_pathways": [dataclasses.asdict(p) for p in self.ranked_pathways],
        }


# ============================================================
# Validation
# ============================================================
def validate_config(config: PipelineConfig) -> Optional[str]:
    """
    Return None if the config is runnable, or an error message string
    explaining why not. Kept SEPARATE from run_pipeline() so a UI can
    show live validation errors as the user types.
    """
    if not config.starter_smiles or not config.starter_smiles.strip():
        return "Starter SMILES is required."
    if not config.target_smiles or not config.target_smiles.strip():
        return "Target SMILES is required."

    starter_mol = Chem.MolFromSmiles(config.starter_smiles.strip())
    target_mol  = Chem.MolFromSmiles(config.target_smiles.strip())
    if starter_mol is None:
        return f"Invalid starter SMILES: `{config.starter_smiles}`"
    if target_mol is None:
        return f"Invalid target SMILES: `{config.target_smiles}`"

    if Chem.MolToSmiles(starter_mol) == Chem.MolToSmiles(target_mol):
        return (
            "Starter and target are the same molecule — there's no "
            "pathway to find. Change the target, or leave it empty for "
            "open exploration."
        )
    if config.domain not in ("chem", "bio", "both"):
        return f"Domain must be 'chem', 'bio', or 'both', got {config.domain!r}."
    if config.direction not in ("forward", "retro", "bidirectional"):
        return (
            f"Direction must be 'forward', 'retro', or 'bidirectional', "
            f"got {config.direction!r}."
        )
    if config.strategy not in ("cartesian", "priority_queue"):
        return (
            f"Strategy must be 'cartesian' or 'priority_queue', got "
            f"{config.strategy!r}."
        )
    if config.generations < 1 or config.generations > 8:
        return "Generations per side must be between 1 and 8."
    return None


# ============================================================
# The pipeline
# ============================================================
def run_pipeline(
    config: PipelineConfig,
    *,
    thermo_calc: Optional[Callable[[str], Optional[float]]] = None,
) -> PipelineResult:
    """
    GENERATION phase: expand → merge → trace → parse. Returns the found
    pathways UNRANKED (sorted by step count).

    Ranking is deliberately NOT done here. DORAnet's pathway_ranking is
    slow (minutes, single-threaded under Windows/Streamlit) and its
    criterion weights are something the user tunes AFTER seeing the
    pathways — so it's a separate, opt-in step: call `rank_pathways()`
    with the chosen weights once these pathways are on screen.

    Parameters
    ----------
    config : PipelineConfig
        Frontend-agnostic settings.
    thermo_calc : callable or None
        A SMILES → Hf(kJ/mol) callable (e.g. an RMG client). When
        supplied and config.enable_rmg is True, gets passed through as
        DORAnet's molecule_thermo_calculator during expansion, so
        every reaction gets a real ΔH stored on the network.

    Returns
    -------
    PipelineResult
        `.ranked_pathways` holds the UNRANKED pathways (final_score is
        None); `.diagnostics["ranked"]` is False.
    """
    err = validate_config(config)
    if err:
        return PipelineResult(ok=False, error=err)

    t0 = time.time()

    # No starter/target swap. Domain (chem/bio) and direction
    # (forward/retro/bidirectional) are independent, and the molecules
    # are used exactly as entered: `starter` is always what the user
    # starts from, `target` is always what they want to reach.
    starter = config.starter_smiles.strip()
    target  = config.target_smiles.strip()

    limits = {
        "max_mw": float(config.max_molecular_weight),
        "max_atoms": config.max_atoms,
        "max_dh": float(config.max_rxn_dh),
    }
    effective_thermo = thermo_calc if config.enable_rmg else None

    # ---- expand + trace ----
    try:
        _expand_and_trace(
            starter=starter,
            target=target,
            job_name=config.job_name,
            gen_count=config.generations,
            strategy=config.strategy,
            beam_size=config.beam_size,
            limits=limits,
            include_chem=config.include_chem,
            include_bio=config.include_bio,
            directions=config.directions,
            bio_whitelist=config.bio_whitelist,
            chem_whitelist=config.chem_whitelist,
            thermo_calc=effective_thermo,
        )
    except Exception as e:
        return PipelineResult(
            ok=False,
            error=f"{type(e).__name__}: {e}",
            elapsed_seconds=time.time() - t0,
        )

    pathway_file = f"{config.job_name}_pathways.txt"
    if not os.path.exists(pathway_file):
        return PipelineResult(
            ok=True,
            error=None,
            elapsed_seconds=time.time() - t0,
            n_pathways=0,
            diagnostics={
                "reason": "no_pathways_found",
                "note": "Network built successfully but no pathways connect "
                        "starter to target within the search depth.",
            },
        )

    # ---- parse (UNRANKED) ----
    # Generation stops here — see the function docstring for why ranking
    # is a separate step. Return the raw pathways sorted by step count so
    # the UI can display them instantly.
    pathways = _parse_unranked_pathways(config.job_name)
    return PipelineResult(
        ok=True,
        error=None,
        elapsed_seconds=time.time() - t0,
        n_pathways=len(pathways),
        ranked_pathways=pathways,
        pathway_file_path=pathway_file,
        diagnostics={
            "ranked": False,
            "domain": config.domain,
            "directions": config.directions,
        },
    )


def rank_pathways(
    config: PipelineConfig,
    *,
    weights: Optional[dict] = None,
    thermo_calc: Optional[Callable[[str], Optional[float]]] = None,
    equilibrator_client: Optional[Any] = None,
) -> PipelineResult:
    """
    RANKING phase: run DORAnet's pathway_ranking (with adjustable
    criterion `weights`) on the pathways produced by a prior
    `run_pipeline()` call, then optionally decorate/prune with
    equilibrator. Requires `{config.job_name}_pathways.txt` to exist.

    Parameters
    ----------
    weights : dict or None
        DORAnet criterion weights (number_of_steps, reaction_thermo,
        by_product_number, atom_economy, ...). None uses DEFAULT_WEIGHTS.
    thermo_calc : callable or None
        RMG-style SMILES → Hf calculator, used only when
        config.enable_rmg is True.
    equilibrator_client : EquilibratorClient or None
        When supplied and config.enable_equilibrator is True, ranked
        pathways get per-bio-step ΔG'° values and are pruned by
        config.equilibrator_prune_max_abs_dg.
    """
    t0 = time.time()
    pathway_file = f"{config.job_name}_pathways.txt"
    if not os.path.exists(pathway_file):
        return PipelineResult(
            ok=False,
            error="No pathways to rank — generate pathways first.",
        )

    starter = config.starter_smiles.strip()
    target  = config.target_smiles.strip()
    effective_thermo = thermo_calc if config.enable_rmg else None

    try:
        ranked = generate_base_rankings(
            starter=starter,
            target=target,
            helpers=["O", "[H][H]"],
            job_name=config.job_name,
            weights=weights,
            molecule_thermo_calculator=effective_thermo,
        )
    except Exception as e:
        # Ranking is optional. Fall back to raw pathways sorted by length.
        ranked = _fallback_ranked(config.job_name)
        return PipelineResult(
            ok=True,
            error=f"ranking failed ({type(e).__name__}: {e}); unranked list returned",
            elapsed_seconds=time.time() - t0,
            n_pathways=len(ranked),
            ranked_pathways=ranked,
            pathway_file_path=pathway_file,
            diagnostics={"ranked": False},
        )

    # ---- equilibrator decoration + optional pruning ----
    n_before_eq = len(ranked)
    n_pruned = 0
    if config.enable_equilibrator and equilibrator_client is not None:
        ranked = decorate_with_equilibrator(
            ranked,
            equilibrator_client,
            max_abs_dg_threshold=float(config.equilibrator_prune_max_abs_dg),
        )
        n_pruned = n_before_eq - len(ranked)

    return PipelineResult(
        ok=True,
        error=None,
        elapsed_seconds=time.time() - t0,
        n_pathways=len(ranked),
        ranked_pathways=ranked,
        pathway_file_path=pathway_file,
        diagnostics={
            "ranked": True,
            "equilibrator_pruned": n_pruned,
        },
    )


# ============================================================
# Internals — expansion + trace + fallback
# ============================================================
def _write_smi(path: str, lines: List[str]) -> None:
    with open(path, "w") as f:
        for s in lines:
            f.write(s + "\n")


def _register_bio_ruleset(whitelist: Optional[List[str]], job_name: str) -> str:
    """
    Build a filtered copy of DORAnet's JN1224MIN rule table containing
    ONLY the whitelisted rules, register it under a ruleset name, and
    return that name.

    DORAnet's native enzymatic.generate_network selects operators by
    ruleset name (not by an in-memory whitelist), so restricting to our
    3 polyketide rules means writing a filtered TSV and registering it.
    Keeping the whitelist tight is also what keeps expansion bounded —
    the Claisen-2 *variants* (rule0126, rule0350) explode, so they must
    stay out of the list.
    """
    wl = set(whitelist) if whitelist else set(DEFAULT_BIO_WHITELIST)
    src = (Path(doranet.__file__).parent / "modules" / "enzymatic"
           / "JN1224MIN_rules.tsv")
    out = Path(f"{job_name}_bio_ruleset.tsv").resolve()
    with open(src, encoding="utf-8") as fin, \
         open(out, "w", encoding="utf-8", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames,
                                delimiter="\t")
        writer.writeheader()
        for row in reader:
            if row["Name"] in wl:
                writer.writerow(row)
    ruleset_name = f"PIPELINE_BIO_{job_name}"
    _BIO_AVAILABLE_RULESETS[ruleset_name] = out
    return ruleset_name


def _expand_and_trace(
    *,
    starter: str,
    target: str,
    job_name: str,
    gen_count: int,
    strategy: str,
    beam_size: int,
    limits: dict,
    include_chem: bool,
    include_bio: bool,
    directions: List[str],
    bio_whitelist: Optional[List[str]],
    chem_whitelist: Optional[List[str]],
    thermo_calc,
) -> None:
    """
    Run expansion + pretreat + pathway_finder. Writes the pathway file
    as a side effect. Returns nothing — callers check for the file's
    existence.

    Domain and direction are fully independent:
      - `include_chem` / `include_bio` select the operator set(s); either
        or both may be enabled (both = a mixed chem+bio network).
      - `directions` is a list drawn from {"forward", "retro"}. Each
        enabled domain is expanded in each requested direction, and all
        resulting networks are merged before the pathway trace. A
        two-direction list is the meet-in-the-middle search.

    No starter/target swap happens anywhere: forward always expands from
    `starter`, retro always expands back from `target`.
    """
    # Helpers = freely-available co-reactants that don't count as pathway
    # steps. When bio operators are on, malonyl-CoA MUST be a helper: the
    # Claisen rules consume it, it's not a DORAnet cofactor, and it's what
    # test_tal_centered_combined.py passed to pretreat/pathway_finder.
    # Without it, no bio pathway ever connects (pathway_finder rejects any
    # route that consumes an un-suppliable molecule).
    helpers = ["O", "[H][H]"]
    if include_bio:
        helpers = helpers + [MALONYL_COA]

    starter_file = f"{job_name}_starter.smi"
    helper_file  = f"{job_name}_helpers.smi"        # chem expansion: O + H2
    target_file  = f"{job_name}_target.smi"
    _write_smi(starter_file, [starter])
    _write_smi(helper_file, ["O", "[H][H]"])
    _write_smi(target_file, [target])

    do_forward = "forward" in directions
    do_retro   = "retro" in directions

    networks = []

    if include_chem:
        chem_kwargs = dict(
            include_chem=True,
            include_bio=False,
            max_atoms=limits["max_atoms"],
            max_molecular_weight=limits["max_mw"],
            max_rxn_thermo_change=limits["max_dh"],
            molecule_thermo_calculator=thermo_calc,
            chem_whitelist=chem_whitelist,
        )
        # retro gets a smaller beam — feedstock-proximity is a coarser
        # guide than product-Tanimoto, so we spend fewer expansions there.
        beam_retro = max(50, beam_size // 5)
        if do_forward:
            if strategy == "priority_queue":
                networks.append(generate_network_tal(
                    job_name=f"{job_name}_chem_fwd",
                    starters=starter_file, helpers=helper_file,
                    gen=gen_count, direction="forward",
                    strategy="priority_queue",
                    targets=target,
                    recipe_ranker=ForwardProductTanimotoRanker(target),
                    beam_size=int(beam_size),
                    **chem_kwargs,
                ))
            else:
                networks.append(generate_network_tal(
                    job_name=f"{job_name}_chem_fwd",
                    starters=starter_file, helpers=helper_file,
                    gen=gen_count, direction="forward",
                    strategy="cartesian", **chem_kwargs,
                ))
        if do_retro:
            if strategy == "priority_queue":
                networks.append(generate_network_tal(
                    job_name=f"{job_name}_chem_retro",
                    starters=target_file, helpers=helper_file,
                    gen=gen_count, direction="retro",
                    strategy="priority_queue",
                    targets=starter,
                    recipe_ranker=FeedstockProximityRanker([starter]),
                    beam_size=int(beam_retro),
                    **chem_kwargs,
                ))
            else:
                networks.append(generate_network_tal(
                    job_name=f"{job_name}_chem_retro",
                    starters=target_file, helpers=helper_file,
                    gen=gen_count, direction="retro",
                    strategy="cartesian", **chem_kwargs,
                ))

    if include_bio:
        # Bio expansion goes through DORAnet's NATIVE enzymatic
        # generate_network (see import note). The JN1224MIN operators are
        # substrate->product, so they fire FORWARD; native retro yields 0
        # reactions (it doesn't flip the SMARTS). We therefore always
        # expand bio forward, regardless of `directions` — it's the only
        # productive bio direction.
        #
        # Starters = the user's starter PLUS malonyl-CoA, the polyketide
        # chain-extender the Claisen rules consume. Malonyl-CoA is not a
        # DORAnet cofactor, so it must be supplied as a co-starter; native
        # then treats the real cofactors (CoA, CO2, ...) as operator slots
        # and the expansion stays bounded (~100 reactions for gen=3).
        ruleset_name = _register_bio_ruleset(bio_whitelist, job_name)
        bio_starter_file = f"{job_name}_bio_starters.smi"
        _write_smi(bio_starter_file, [starter, MALONYL_COA])
        networks.append(_native_bio_generate_network(
            job_name=f"{job_name}_bio_fwd",
            starters=bio_starter_file,
            gen=gen_count,
            direction="forward",
            allow_multiple_reactants=True,
            targets=target,
            ruleset=ruleset_name,
            max_rxn_thermo_change=limits["max_dh"],
        ))

    if not networks:
        raise RuntimeError(
            "No networks generated. Need at least one domain "
            "(include_chem or include_bio) AND at least one direction "
            "(forward or retro)."
        )

    # Longest traceable path = sum of the per-side depths we actually ran.
    # Bidirectional (two sides) reaches gen_count*2; a single side reaches
    # gen_count.
    reach = gen_count * len(directions)
    pretreat_networks(
        networks=networks,
        starters=[starter],
        helpers=helpers,
        total_generations=reach,
        job_name=job_name,
    )
    pathway_finder(
        starters=[starter],
        helpers=helpers,
        target=[target],
        search_depth=reach,
        max_num_rxns=reach + 3,
        job_name=job_name,
    )


def _parse_unranked_pathways(job_name: str) -> List[RankedPathway]:
    """Parse the raw pathway file into RankedPathway objects sorted by
    step count, WITHOUT running DORAnet's scorer. `final_score` is set to
    None to signal "not yet ranked"; score-derived fields (atom economy,
    by-products) stay at defaults until rank_pathways() fills them in.

    Used by the generation phase so the UI can show pathways immediately.
    """
    try:
        raw = load_pathways_from_file(job_name)
    except FileNotFoundError:
        return []
    raw_sorted = sorted(raw, key=lambda p: p.num_steps)
    out = []
    for idx, p in enumerate(raw_sorted, 1):
        smiles_list, names_list, dh_list = [], [], []
        for rxn_str in p.reactions:
            parsed = parse_reaction_string(rxn_str)
            smiles_list.append(
                ".".join(parsed["reactants"]) + ">>" + ".".join(parsed["products"])
            )
            names_list.append(parsed["op_name"])
            dh_list.append(parsed["dH"])
        out.append(RankedPathway(
            rank=idx,
            final_score=None,               # None == unranked
            atomic_economy=0.0,
            pathway_byproduct_count=0,
            intermediate_byproducts={},
            reaction_smiles=smiles_list,
            reaction_names=names_list,
            reaction_enthalpies=dh_list,
        ))
    return out


def _fallback_ranked(job_name: str) -> List[RankedPathway]:
    """When DORAnet's ranker fails, build a minimum-viable ranked list
    from the raw pathway file, ordered by step count. Fields the ranker
    normally computes (final_score, atomic_economy, byproducts) get
    zero/empty defaults."""
    try:
        raw = load_pathways_from_file(job_name)
    except FileNotFoundError:
        return []
    raw_sorted = sorted(raw, key=lambda p: p.num_steps)
    ranked = []
    for idx, p in enumerate(raw_sorted, 1):
        # extract each step's SMILES / op / dH from the pathway_tools objs
        smiles_list, names_list, dh_list = [], [], []
        for rxn_str in p.reactions:
            parsed = parse_reaction_string(rxn_str)
            smiles_list.append(
                ".".join(parsed["reactants"]) + ">>" + ".".join(parsed["products"])
            )
            names_list.append(parsed["op_name"])
            dh_list.append(parsed["dH"])
        ranked.append(RankedPathway(
            rank=idx,
            final_score=float(len(raw_sorted) - idx + 1),  # crude ordering
            atomic_economy=0.0,
            pathway_byproduct_count=0,
            intermediate_byproducts={},
            reaction_smiles=smiles_list,
            reaction_names=names_list,
            reaction_enthalpies=dh_list,
        ))
    return ranked


# ============================================================
# House-keeping (side-effect files)
# ============================================================
def cleanup_job_files(job_name: str) -> None:
    """Remove per-run artefacts. Frontends should call this before Run
    so stale results from a previous invocation don't linger."""
    prefixes = [
        f"{job_name}_pathways.txt",
        f"{job_name}_ranked_pathways.txt",
        f"{job_name}_network_pretreated.json",
        f"{job_name}_reaxys_batch_query.txt",
        f"{job_name}_reaxys_batch_result.csv",
    ]
    for p in prefixes:
        try:
            os.remove(p)
        except FileNotFoundError:
            pass
