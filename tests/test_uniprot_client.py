"""UniProt response parsing — tested against a real captured TSV fixture, so
no network is needed and the assertions can't drift from the real format.
resolve()'s caching/batching is tested with a stubbed fetch."""
from conftest import FIXTURES

import uniprot_client


def _fixture_tsv() -> str:
    return (FIXTURES / "uniprot_sample.tsv").read_text(encoding="utf-8")


def test_parses_all_records():
    recs = uniprot_client.parse_tsv(_fixture_tsv())
    assert {r["accession"] for r in recs} == {"P42765", "Q64428", "D2I940"}


def test_primary_protein_name_stripped_of_alternates():
    recs = {r["accession"]: r for r in uniprot_client.parse_tsv(_fixture_tsv())}
    # P42765's full name has many parenthetical alternates; keep the primary
    assert recs["P42765"]["protein_name"] == "3-ketoacyl-CoA thiolase, mitochondrial"
    assert recs["D2I940"]["protein_name"] == "Acetyl-CoA acetyltransferase"


def test_ec_numbers_split():
    recs = {r["accession"]: r for r in uniprot_client.parse_tsv(_fixture_tsv())}
    assert "2.3.1.9" in recs["P42765"]["ec"]
    assert "2.3.1.16" in recs["P42765"]["ec"]


def test_human_readable_reaction_extracted():
    recs = {r["accession"]: r for r in uniprot_client.parse_tsv(_fixture_tsv())}
    # the reaction equation is pulled out of the verbose catalytic-activity blob
    assert recs["P42765"]["reactions"][0] == "an acyl-CoA + acetyl-CoA = a 3-oxoacyl-CoA + CoA"
    assert recs["D2I940"]["reactions"] == ["succinyl-CoA + acetyl-CoA = 3-oxoadipyl-CoA + CoA"]
    assert recs["D2I940"]["reaction_count"] == 1


def test_gene_and_organism():
    recs = {r["accession"]: r for r in uniprot_client.parse_tsv(_fixture_tsv())}
    assert recs["P42765"]["gene"] == "ACAA2"
    assert recs["P42765"]["organism"] == "Homo sapiens"        # common name dropped
    assert recs["Q64428"]["organism"] == "Rattus norvegicus"


def test_blank_ec_column_is_ok():
    # D2I940 has an empty EC column but an EC lives inside its reaction text;
    # the parser must not choke on the blank column.
    recs = {r["accession"]: r for r in uniprot_client.parse_tsv(_fixture_tsv())}
    assert recs["D2I940"]["ec"] == []
    assert recs["D2I940"]["reactions"]                          # still got the reaction


def test_multifunctional_enzyme_reports_true_count_and_caps_inline():
    recs = {r["accession"]: r for r in uniprot_client.parse_tsv(_fixture_tsv())}
    # Q64428 is multifunctional with dozens of reactions: inline list is
    # capped, but the true total is reported so the UI can link out.
    assert len(recs["Q64428"]["reactions"]) <= uniprot_client._MAX_INLINE_REACTIONS
    assert recs["Q64428"]["reaction_count"] > 3


def test_deleted_entry_is_flagged_not_named():
    tsv = (
        "Entry\tProtein names\tEC number\tGene Names\tOrganism\tCatalytic activity\n"
        "A0A072ZQE7\tdeleted\t\t\t\t\n"
    )
    rec = uniprot_client.parse_tsv(tsv)[0]
    assert rec["accession"] == "A0A072ZQE7"
    assert rec["deleted"] is True
    assert rec["protein_name"] == ""          # not the literal "deleted"
    assert rec["ec"] == [] and rec["reaction_count"] == 0


def test_empty_input():
    assert uniprot_client.parse_tsv("") == []


def test_resolve_caches_and_batches(monkeypatch):
    calls = {"n": 0}

    def fake_fetch(accessions):
        calls["n"] += 1
        return _fixture_tsv()   # returns the 3 fixture records regardless

    uniprot_client.clear_cache()
    monkeypatch.setattr(uniprot_client, "_fetch_tsv", fake_fetch)

    recs = uniprot_client.resolve(["P42765", "Q64428", "D2I940"])
    assert {r["accession"] for r in recs} == {"P42765", "Q64428", "D2I940"}
    assert calls["n"] == 1                     # one batched call

    # second call for the same accessions -> served from cache, no new fetch
    uniprot_client.resolve(["P42765"])
    assert calls["n"] == 1


def test_is_accession_accepts_real_rejects_gene_names():
    assert uniprot_client._is_accession("P42765")
    assert uniprot_client._is_accession("A0A023HHB5")
    # gene names / GenBank ids from the ruleset's messy tail
    assert not uniprot_client._is_accession("NCED52")
    assert not uniprot_client._is_accession("ADAMTS2")
    assert not uniprot_client._is_accession("")


def test_resolve_drops_non_accessions_before_query(monkeypatch):
    """A non-accession token like NCED52 would 400 the whole batch, so it
    must be filtered out before the query is built."""
    seen = {}

    def capture(accessions):
        seen["batch"] = list(accessions)
        return _fixture_tsv()

    uniprot_client.clear_cache()
    monkeypatch.setattr(uniprot_client, "_fetch_tsv", capture)
    uniprot_client.resolve(["P42765", "NCED52", "Q64428"])
    assert "NCED52" not in seen["batch"]
    assert "P42765" in seen["batch"] and "Q64428" in seen["batch"]


def test_resolve_survives_network_failure(monkeypatch):
    def boom(accessions):
        raise OSError("network down")

    uniprot_client.clear_cache()
    monkeypatch.setattr(uniprot_client, "_fetch_tsv", boom)
    # must not raise — just returns nothing resolvable
    assert uniprot_client.resolve(["P42765"]) == []
