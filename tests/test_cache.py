"""api/cache.py — canonical key logic + the bounded LRU store.
Verifies the caching correctness properties we designed for: job_name never
affects the key, equivalent SMILES collide, real parameter changes miss."""
import pytest

from pipeline import PipelineConfig
from api.cache import generation_key, ranking_key, ResultCache

BASE = dict(
    starter_smiles="CC=CC=CC(=O)O", target_smiles="Cc1cc(O)cc(=O)o1",
    domain="chem", direction="bidirectional", generations=3,
    max_molecular_weight=200.0, helpers=["O", "[H][H]"], job_name="api_abc",
)


def cfg(**over):
    return PipelineConfig.from_dict({**BASE, **over})


# ------------------------------------------------------- generation key
def test_job_name_does_not_affect_key():
    assert generation_key(cfg()) == generation_key(cfg(job_name="api_zzz_different"))


@pytest.mark.parametrize("over", [
    {"max_molecular_weight": 300.0},
    {"max_rxn_dh": 25.0},
    {"generations": 4},
    {"feasibility_prune_threshold": 0.7},
    {"domain": "bio"},
    {"direction": "forward"},
])
def test_changing_any_expansion_param_changes_key(over):
    assert generation_key(cfg()) != generation_key(cfg(**over))


def test_equivalent_smiles_spelling_collides():
    assert generation_key(cfg(target_smiles="Cc1cc(O)cc(=O)o1")) == \
        generation_key(cfg(target_smiles="c1(C)cc(O)cc(=O)o1"))


def test_stereochemistry_stays_distinct():
    # cis vs trans are different molecules -> must NOT collide
    assert generation_key(cfg(starter_smiles="C/C=C/C=C/C(=O)O")) != \
        generation_key(cfg(starter_smiles="C/C=C\\C=C/C(=O)O"))


def test_int_vs_float_numeric_normalized():
    assert generation_key(cfg(max_molecular_weight=200)) == \
        generation_key(cfg(max_molecular_weight=200.0))


def test_whitelist_order_independent_but_content_matters():
    assert generation_key(cfg(bio_whitelist=["r1", "r2"])) == \
        generation_key(cfg(bio_whitelist=["r2", "r1"]))
    assert generation_key(cfg(bio_whitelist=["r1", "r2"])) != \
        generation_key(cfg(bio_whitelist=["r1", "r3"]))


# ----------------------------------------------------------- ranking key
def test_ranking_key_nests_generation():
    c = cfg()
    assert ranking_key(c, None, None, None) != ranking_key(cfg(max_molecular_weight=300.0), None, None, None)


def test_ranking_key_reacts_to_weights():
    c = cfg()
    base = ranking_key(c, None, None, None)
    assert base == ranking_key(c, None, None, None)                       # stable
    assert base != ranking_key(c, None, None, {"stability": 2.0})         # weight change -> miss


def test_ranking_key_weight_order_independent():
    c = cfg()
    assert ranking_key(c, None, None, {"stability": 1.0, "diversity": 2.0}) == \
        ranking_key(c, None, None, {"diversity": 2.0, "stability": 1.0})


# ------------------------------------------------------------ ResultCache
def test_store_get_and_miss():
    rc = ResultCache(max_entries=3)
    rc.put_generation("a", {"pathways": [1], "job_name": "j"})
    assert rc.get_generation("a")["pathways"] == [1]
    assert rc.get_generation("missing") is None


def test_lru_eviction():
    rc = ResultCache(max_entries=3)
    for k in ("a", "b", "c"):
        rc.put_generation(k, {"x": k})
    rc.get_generation("a")               # touch 'a' -> most recently used
    rc.put_generation("d", {"x": "d"})   # over cap -> evict LRU ('b')
    assert rc.get_generation("b") is None
    assert rc.get_generation("a") is not None
    assert rc.get_generation("d") is not None


def test_generation_and_ranking_stores_independent():
    rc = ResultCache()
    rc.put_generation("k", {"x": 1})
    assert rc.get_ranking("k") is None
