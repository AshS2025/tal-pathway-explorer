"""
enzyme_info.py — enzyme annotations for bio (enzymatic) reaction rules.

Each JN1224MIN rule lists, in the ruleset's `Comments` column, the UniProt
protein accessions known to catalyze that transformation. Some rules list
none — e.g. spontaneous cyclizations/lactonizations (like the TAL
ring-closure that offloads CoA) that need no dedicated enzyme, or gaps in
the source curation. Either way "no enzyme" is meaningful to a biologist, so
we surface it rather than hide it.

This module loads that per-rule enzyme table once (cached) and exposes the
count per pathway step. Chem steps have no enzyme concept and return None.
"""
from __future__ import annotations

import csv
import functools
from pathlib import Path
from typing import List, Optional

import doranet

from pathway_scoring import is_bio_op


def _ruleset_path() -> Path:
    return (Path(doranet.__file__).parent / "modules" / "enzymatic"
            / "JN1224MIN_rules.tsv")


@functools.lru_cache(maxsize=1)
def _rule_enzymes() -> dict:
    """rule name -> list of UniProt accessions (possibly empty)."""
    out: dict = {}
    with open(_ruleset_path(), encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            ids = [x.strip() for x in (row.get("Comments") or "").split(";")
                   if x.strip()]
            out[row["Name"]] = ids
    return out


def enzyme_ids_for_rule(rule_name: str) -> List[str]:
    """UniProt accessions catalyzing this rule (empty list if none/unknown)."""
    return _rule_enzymes().get(rule_name, [])


def enzyme_count_for_step(op_name: str) -> Optional[int]:
    """Enzyme count for one pathway step, or None if the step is CHEMICAL
    (no enzyme concept). A bio step returns its count; 0 means "no known
    enzyme" (possibly spontaneous)."""
    if not is_bio_op(op_name):
        return None
    return len(enzyme_ids_for_rule(op_name))


def minimum_enzyme_count(reaction_names: list) -> int:
    """Minimum distinct enzymes needed to build a pathway's bio steps: one
    per DISTINCT bio rule that has a known enzyme.

    Repeating the SAME rule (the same reaction) reuses the same enzyme, so
    it collapses to one. We deliberately do NOT collapse *different* rules
    that merely share a candidate enzyme — whether a promiscuous enzyme
    truly catalyzes both reactions needs biological judgment we don't have,
    so we stay conservative. Chem steps and no-known-enzyme bio steps (e.g.
    the possibly-spontaneous rule0891) don't count.
    """
    distinct = {n for n in reaction_names
                if is_bio_op(n) and enzyme_ids_for_rule(n)}
    return len(distinct)


def annotate_pathways(pathways: list) -> None:
    """Mutate each pathway dict in place, adding:
      * `reaction_enzymes`: list parallel to `reaction_names` — None for chem
        steps, an int candidate-count for bio steps (0 = no known enzyme).
      * `min_enzymes`: minimum distinct enzymes to build the whole route
        (shared multifunctional enzymes counted once)."""
    for p in pathways:
        names = p.get("reaction_names", [])
        p["reaction_enzymes"] = [enzyme_count_for_step(n) for n in names]
        p["min_enzymes"] = minimum_enzyme_count(names)
