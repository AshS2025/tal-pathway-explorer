"""Per-rule enzyme annotations sourced from JN1224MIN's Comments column."""
import enzyme_info


def test_known_rule_enzyme_counts():
    # counts confirmed against the shipped JN1224MIN ruleset
    assert enzyme_info.enzyme_count_for_step("rule0087") == 7
    assert enzyme_info.enzyme_count_for_step("rule1118") == 29


def test_enzymeless_rule_is_zero_not_none():
    # rule0891 (TAL ring-closure / CoA offload) has no annotated enzyme
    assert enzyme_info.enzyme_count_for_step("rule0891") == 0


def test_chem_step_returns_none():
    # chem operators have human names, not rule IDs -> no enzyme concept
    assert enzyme_info.enzyme_count_for_step("Dehydration of Alcohol") is None
    assert enzyme_info.enzyme_count_for_step("Aldol Condensation") is None


def test_enzyme_ids_are_uniprot_like():
    ids = enzyme_info.enzyme_ids_for_rule("rule0087")
    assert len(ids) == 7
    assert all(isinstance(x, str) and x for x in ids)


def test_annotate_pathways_adds_parallel_list():
    pathways = [{
        "reaction_names": ["rule0087", "Dehydration of Alcohol", "rule0891"],
    }]
    enzyme_info.annotate_pathways(pathways)
    assert pathways[0]["reaction_enzymes"] == [7, None, 0]


def test_annotate_handles_missing_names():
    pathways = [{}]
    enzyme_info.annotate_pathways(pathways)
    assert pathways[0]["reaction_enzymes"] == []
