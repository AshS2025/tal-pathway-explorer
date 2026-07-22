"""
uniprot_client.py — resolve UniProt accessions to enzyme metadata.

The bio rules (JN1224MIN) list UniProt accessions per rule; this turns an
accession like `P42765` into the human-readable fields a biologist wants:
protein name, EC number(s), the reaction equation(s) the enzyme catalyzes,
gene name, and source organism. All of it comes straight from UniProt's
public REST API (no key needed) — the authoritative source, so we never
rely on the model's memory for a mapping.

`resolve(accessions)` batches the lookups into few HTTP calls, caches
results, and degrades gracefully: if UniProt is unreachable it returns
whatever is cached rather than raising. The parsing (`parse_tsv`) is a pure
function so it can be tested without network.

Honest caveats:
  * Many enzymes are multifunctional — they list several EC numbers and
    reactions, and we can't tell from the accession alone which one matches
    the specific rule. We surface the enzyme's reactions (capped) and let
    the biologist judge.
  * ~1% of the ids in the ruleset aren't UniProt accessions (gene names,
    GenBank ids); those simply won't resolve and are omitted.

Storage is an in-memory dict behind this module (like the other caches),
swappable for a disk/Redis cache later — UniProt records are stable, so a
persistent cache would be a reasonable upgrade.
"""
from __future__ import annotations

import re
import threading
import urllib.parse
import urllib.request
from typing import Dict, List, Optional

_BASE = "https://rest.uniprot.org/uniprotkb/search"
_WEB = "https://www.uniprot.org/uniprotkb"
_FIELDS = ("accession,protein_name,ec,gene_names,organism_name,mass,"
           "cc_catalytic_activity")
# We only show reactions inline when an enzyme has few of them (a
# multifunctional enzyme with many reactions is ambiguous — we can't tell
# which matches the rule — so the UI links out to UniProt instead). We still
# report the true total via `reaction_count`.
_MAX_INLINE_REACTIONS = 3
_BATCH = 80                 # accessions per HTTP request (keeps the URL sane)
_TIMEOUT = 25               # seconds per request

_cache: Dict[str, dict] = {}
_lock = threading.Lock()

_REACTION_RE = re.compile(r"Reaction=(.*?);")

# Official UniProt accession syntax. ~1% of the ids in the ruleset are gene
# names / GenBank ids (e.g. "NCED52", "ADAMTS2") that aren't accessions —
# and a single invalid token like `accession:NCED52` makes UniProt reject
# the WHOLE batch query with HTTP 400. Filter to real accessions first.
_ACCESSION_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)


def _is_accession(s: str) -> bool:
    return bool(_ACCESSION_RE.match(s or ""))


def _primary_name(protein_names: str) -> str:
    """UniProt lists the recommended name first, then alternates in parens —
    take the part before the first ' ('."""
    return (protein_names or "").split(" (")[0].strip()


def _scientific_organism(organism: str) -> str:
    """Drop the '(Human)' common-name suffix, keep the scientific name."""
    return (organism or "").split(" (")[0].strip()


def parse_tsv(tsv: str) -> List[dict]:
    """Parse a UniProt TSV response (with our _FIELDS columns) into records.
    Pure function — no network — so it's unit-testable against a fixture."""
    lines = tsv.splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    idx = {name: i for i, name in enumerate(header)}
    out: List[dict] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")

        def get(col: str) -> str:
            i = idx.get(col)
            return cols[i] if i is not None and i < len(cols) else ""

        # A demerged/removed UniProt entry comes back with the literal name
        # "deleted" and no other data. Keep it visible (the accession is
        # still meaningful to a biologist) but flag it rather than showing
        # "deleted" as if it were a protein name.
        raw_name = get("Protein names").strip()
        deleted = raw_name.lower() == "deleted"

        mass_raw = get("Mass").replace(",", "").strip()   # Daltons
        mass = int(mass_raw) if mass_raw.isdigit() else None

        ec = [e.strip() for e in get("EC number").split(";") if e.strip()]
        reactions = [r.strip() for r in _REACTION_RE.findall(get("Catalytic activity"))]
        # de-dupe reactions while preserving order
        seen, uniq = set(), []
        for r in reactions:
            if r not in seen:
                seen.add(r)
                uniq.append(r)
        out.append({
            "accession": get("Entry"),
            "protein_name": "" if deleted else _primary_name(raw_name),
            "deleted": deleted,
            "ec": ec,
            "mass": mass,                                      # Daltons (None if unknown)
            "gene": get("Gene Names").split(" ")[0].strip(),   # first gene symbol
            "organism": _scientific_organism(get("Organism")),
            # inline reactions only when there are few; `reaction_count` is
            # the true total so the UI can link out when there are more.
            "reactions": uniq[:_MAX_INLINE_REACTIONS],
            "reaction_count": len(uniq),
        })
    return out


def _fetch_tsv(accessions: List[str]) -> str:
    """One HTTP call to UniProt for a batch of accessions. Raises on failure."""
    query = " OR ".join(f"accession:{a}" for a in accessions)
    url = _BASE + "?" + urllib.parse.urlencode({
        "query": query, "fields": _FIELDS, "format": "tsv", "size": 500,
    })
    req = urllib.request.Request(url, headers={"User-Agent": "TAL-Pathway-Explorer"})
    with urllib.request.urlopen(req, timeout=_TIMEOUT) as r:
        return r.read().decode("utf-8")


def uniprot_search_url(accessions: List[str], cap: int = 100) -> str:
    """A uniprot.org search URL covering these accessions — used for the
    'view all in UniProt' link when we only show the first few inline."""
    valid = [a for a in dict.fromkeys(accessions) if _is_accession(a)][:cap]
    if not valid:
        return _WEB
    query = " OR ".join(f"(accession:{a})" for a in valid)
    return _WEB + "?" + urllib.parse.urlencode({"query": query})


def resolve(accessions: List[str], limit: Optional[int] = None) -> List[dict]:
    """Resolve accessions to enzyme records, in the given order. Cached and
    batched into few HTTP calls. Non-accession ids (gene names etc.) are
    dropped first — they'd 400 the batch — then the list is capped to
    `limit`. Unresolvable ids (or a failed fetch) are omitted; this never
    raises on a network error."""
    # de-dupe, keep order, and drop non-accession ids (they'd 400 the batch)
    want = list(dict.fromkeys(a for a in accessions if _is_accession(a)))
    if limit is not None:
        want = want[:limit]
    with _lock:
        missing = [a for a in want if a not in _cache]

    for i in range(0, len(missing), _BATCH):
        batch = missing[i:i + _BATCH]
        try:
            records = parse_tsv(_fetch_tsv(batch))
        except Exception as e:               # network/SSL/parse — degrade, don't crash
            print(f"UniProt lookup failed for a batch of {len(batch)}: {e}")
            records = []
        by_acc = {r["accession"]: r for r in records}
        with _lock:
            for a in batch:                  # cache misses (None) too, so we don't refetch
                _cache[a] = by_acc.get(a)

    with _lock:
        return [_cache[a] for a in want if _cache.get(a)]


def clear_cache() -> None:
    with _lock:
        _cache.clear()
