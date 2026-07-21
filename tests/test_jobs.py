"""Run-store eviction: the fix for the memory + disk leaks.

Runs are capped so the in-memory store can't grow forever, and an evicted
run's files are deleted — but never when another live run still shares that
job_name (cache hits reuse an earlier run's files). In-flight runs are never
evicted. jobs.py is stdlib-only, so these run without doranet/FastAPI."""
import time
import types

import pytest

from api import jobs


def _wait(pred, timeout=3.0):
    end = time.time() + timeout
    while time.time() < end:
        if pred():
            return True
        time.sleep(0.01)
    return False


@pytest.fixture(autouse=True)
def _isolate_job_store():
    """Each test gets a clean store and its own cap; restore afterwards."""
    saved_max = jobs._MAX_RUNS
    jobs._runs.clear()
    yield
    jobs._runs.clear()
    jobs._MAX_RUNS = saved_max


def _terminal(job_name, status="generated"):
    r = jobs.create_run()
    r.config = types.SimpleNamespace(job_name=job_name)
    r.status = status
    return r


def test_cap_evicts_oldest_terminal_runs():
    jobs._MAX_RUNS = 3
    ids = [_terminal(f"job{i}").id for i in range(5)]
    assert len(jobs._runs) == 3
    assert ids[0] not in jobs._runs and ids[1] not in jobs._runs   # oldest gone
    assert ids[4] in jobs._runs                                    # newest kept


def test_inflight_run_is_never_evicted():
    jobs._MAX_RUNS = 2
    busy = jobs.create_run()
    busy.config = types.SimpleNamespace(job_name="busy")
    busy.status = "generating"                 # in-flight, and the oldest
    for i in range(4):
        _terminal(f"j{i}")
    assert busy.id in jobs._runs               # survived despite being oldest
    assert busy.status == "generating"


def test_evicted_run_files_are_deleted(tmp_path):
    jobs._MAX_RUNS = 1
    f = tmp_path / "api_del_pathways.txt"
    f.write_text("x")
    r = jobs.create_run()
    r.config = types.SimpleNamespace(job_name=str(tmp_path / "api_del"))
    r.status = "generated"
    jobs.create_run()                          # over cap -> evicts r -> deletes files
    assert not f.exists()


def test_shared_job_name_files_are_preserved(tmp_path):
    jobs._MAX_RUNS = 2
    f = tmp_path / "api_share_pathways.txt"
    f.write_text("x")
    job = str(tmp_path / "api_share")
    _terminal(job)                             # r1 uses job
    _terminal(job)                             # r2 ALSO uses job (cache-hit clone)
    jobs.create_run()                          # evicts r1, but r2 still shares job
    assert f.exists()                          # files kept for the survivor


def test_job_times_out_with_friendly_message():
    run = jobs.create_run()

    def slow_worker(r):
        time.sleep(0.6)                     # longer than the timeout below
        jobs.complete(r, "generated", pathways=[1])

    jobs.run_in_background(run, slow_worker, timeout=0.05,
                           timeout_message="too long, tighten params")
    assert _wait(lambda: run.status == "error"), "should have timed out"
    assert run.timed_out is True
    assert "tighten params" in run.error

    # the slow worker finishes AFTER the deadline — it must NOT resurrect the run
    time.sleep(0.7)                         # let slow_worker (0.6s) run to completion
    assert run.status == "error" and "tighten params" in run.error
    assert run.pathways is None             # late result was discarded by complete()


def test_fast_job_completes_normally():
    run = jobs.create_run()

    def quick_worker(r):
        jobs.complete(r, "generated", pathways=[1, 2])

    jobs.run_in_background(run, quick_worker, timeout=5)
    assert _wait(lambda: run.status == "generated")
    assert run.timed_out is False
    assert run.pathways == [1, 2]


def test_worker_exception_marks_error():
    run = jobs.create_run()

    def boom(r):
        raise ValueError("kaboom")

    jobs.run_in_background(run, boom, timeout=5)
    assert _wait(lambda: run.status == "error")
    assert "kaboom" in run.error
    assert run.timed_out is False


def test_purge_job_files_only_matches_exact_prefix(tmp_path):
    keep = tmp_path / "api_abcd_pathways.txt"      # different job, must survive
    drop = tmp_path / "api_abc_pathways.txt"
    keep.write_text("k")
    drop.write_text("d")
    jobs.purge_job_files(str(tmp_path / "api_abc"))
    assert not drop.exists()
    assert keep.exists()                       # 'api_abc' is not a prefix of 'api_abcd_'
