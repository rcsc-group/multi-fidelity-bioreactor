"""Tests for scripts/run_key.py: one key per run, an index of keys, and the duplicate check in submit_slurm.

A run key is the sha256 of (physics params, binary identity, parent run's key). Two runs with the same key give
the same results, so the second one is not submitted.
"""
import json
import os
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from scripts import run_key as rk
from scripts.simulate import submit_slurm, DuplicateRun
from tests.conftest import CANONICAL_PARAMS, PROJECT_ROOT

P = {**CANONICAL_PARAMS, "run_id": "a", "fidelity": 5}


def _bin(tmp_path, name="BioReactor-mpi", content=b"binary v1"):
    b = tmp_path / "bin" / name
    b.parent.mkdir(parents=True, exist_ok=True)
    b.write_bytes(content)
    return str(b)


# -- the key ---------------------------------------------------------------------------------------------------

def test_key_ignores_run_id_and_underscore_keys():
    k = rk.run_key(P, "sha256:x", None)
    assert rk.run_key(P | {"run_id": "b", "_source_run": "s", "_canonical_run_dir": "/d"}, "sha256:x", None) == k


def test_key_changes_with_physics():
    assert rk.run_key(P | {"omega_b": P["omega_b"] * 1.01}, "sha256:x", None) != rk.run_key(P, "sha256:x", None)


def test_key_ignores_float_noise_and_int_float_spelling():
    noisy = P | {"omega_b": P["omega_b"] * (1 + 1e-15), "fidelity": 5.0}
    assert rk.run_key(noisy, "sha256:x", None) == rk.run_key(P, "sha256:x", None)


def test_key_changes_with_binary_and_parent():
    k = rk.run_key(P, "sha256:x", None)
    assert rk.run_key(P, "sha256:y", None) != k
    assert rk.run_key(P, "sha256:x", "k_parent") != k


# -- binary identity -------------------------------------------------------------------------------------------

def test_binary_identity_is_the_file_hash(tmp_path):
    a, b = _bin(tmp_path, "a", b"one"), _bin(tmp_path, "b", b"one")
    assert rk.binary_identity(a) == rk.binary_identity(b) != rk.binary_identity(_bin(tmp_path, "c", b"two"))
    assert rk.binary_identity(a).startswith("sha256:")


def test_binary_identity_unknown_if_missing_or_changed_after_the_run(tmp_path):
    assert rk.binary_identity(str(tmp_path / "gone")) is None
    assert rk.binary_identity(None) is None
    b = _bin(tmp_path)
    assert rk.binary_identity(b, not_after=time.time() + 60) is not None
    assert rk.binary_identity(b, not_after=os.path.getmtime(b) - 60) is None   # rebuilt in place after the run


def test_default_binary_matches_the_templates():
    for template in ("slurm_mpi_template.sh", "slurm_template.sh"):
        path = rk.default_binary(PROJECT_ROOT / "config" / template, PROJECT_ROOT)
        text = (PROJECT_ROOT / "config" / template).read_text()
        assert Path(path).name in text


# -- the index (backfill of old runs) --------------------------------------------------------------------------

def _run(root, rid, params, results=True, mtime=None):
    d = root / rid
    d.mkdir(parents=True)
    (d / "params.json").write_text(json.dumps({**params, "run_id": rid}))
    if results:
        (d / "results.json").write_text("{}")
    if mtime is not None:
        os.utime(d / "params.json", (mtime, mtime))


def test_build_index_keys_fresh_runs_and_restarts_with_a_known_parent(tmp_path):
    root, b = tmp_path / "runs", _bin(tmp_path)
    later = time.time() + 3600
    _run(root, "cold", P | {"_binary": b}, mtime=later)
    _run(root, "cold_twin", P | {"_binary": b}, mtime=later)
    _run(root, "warm", P | {"_binary": b, "t_checkpoint": 10.0, "_source_run": "cold"}, mtime=later)
    _run(root, "warm_orphan", P | {"_binary": b, "t_checkpoint": 10.0}, mtime=later)
    _run(root, "old_binary", P | {"_binary": str(tmp_path / "gone")}, mtime=later)
    idx = {e["run_id"]: e for e in rk.build_index(root)}
    assert idx["cold"]["key"] == idx["cold_twin"]["key"] and idx["cold"]["dedupe"]
    assert idx["warm"]["parent_key"] == idx["cold"]["key"] and idx["warm"]["dedupe"]
    assert not idx["warm_orphan"]["dedupe"]          # restart from an unknown dump
    assert not idx["old_binary"]["dedupe"]           # binary deleted: its identity is lost
    assert {e["run_id"] for e in rk.read_index(root)} == set(idx)


def test_build_index_skips_dirs_without_params(tmp_path):
    root = tmp_path / "runs"
    (root / "empty").mkdir(parents=True)
    assert rk.build_index(root) == []


# -- submit_slurm --------------------------------------------------------------------------------------------

def _submit(tmp_path, params, **kw):
    ok = MagicMock(stdout="Submitted batch job 42\n", returncode=0)
    with patch("subprocess.run", return_value=ok) as run:
        jid = submit_slurm(params, project_root=PROJECT_ROOT, runs_root=tmp_path / "runs",
                           mpi_scratch_root=tmp_path / "mpi", **kw)
    assert run.called
    return jid


def _refused(tmp_path, params, exc=DuplicateRun, **kw):
    with patch("subprocess.run", side_effect=AssertionError("sbatch must not be called")):
        with pytest.raises(exc) as e:
            submit_slurm(params, project_root=PROJECT_ROOT, runs_root=tmp_path / "runs",
                         mpi_scratch_root=tmp_path / "mpi", **kw)
    return e.value


def test_submit_writes_the_key_and_indexes_the_run(tmp_path):
    p = P | {"_binary": _bin(tmp_path)}
    _submit(tmp_path, p)
    written = json.loads((tmp_path / "runs" / "a" / "params.json").read_text())
    assert written["_run_key"] == rk.run_key(p, rk.binary_identity(p["_binary"]), None)
    assert [e["run_id"] for e in rk.read_index(tmp_path / "runs")] == ["a"]


def test_submit_refuses_a_finished_duplicate_and_names_it(tmp_path):
    p = P | {"_binary": _bin(tmp_path)}
    _submit(tmp_path, p)
    (tmp_path / "runs" / "a" / "results.json").write_text("{}")
    e = _refused(tmp_path, p | {"run_id": "b"})
    assert e.run_dir == tmp_path / "runs" / "a"
    assert not (tmp_path / "runs" / "b").exists()


def test_submit_refuses_a_duplicate_still_in_flight(tmp_path):
    p = P | {"_binary": _bin(tmp_path)}
    _submit(tmp_path, p)
    with patch("scripts.run_key.job_alive", return_value=True) as alive:
        e = _refused(tmp_path, p | {"run_id": "b"})
    alive.assert_called_with("42")
    assert not e.finished


def test_a_failed_duplicate_does_not_block_a_retry(tmp_path):
    """A failed run (no results.json, job gone) and its renamed copies (.failed-<job>, indexed by the backfill
    without a job id) must not block the next submission of the same point."""
    p = P | {"_binary": _bin(tmp_path)}
    _submit(tmp_path, p)
    _run(tmp_path / "runs", "a.failed-41", p, results=False)
    rk.append_index(tmp_path / "runs", rk.entry("a.failed-41", p, rk.binary_identity(p["_binary"]), None, None))
    with patch("scripts.run_key.job_alive", return_value=False):
        assert _submit(tmp_path, p | {"run_id": "b"}) == "42"


def test_submit_allows_a_requested_repeat(tmp_path):
    p = P | {"_binary": _bin(tmp_path)}
    _submit(tmp_path, p)
    assert _submit(tmp_path, p | {"run_id": "b"}, allow_duplicate=True) == "42"


def test_submit_allows_a_retry_under_the_same_run_id(tmp_path):
    p = P | {"_binary": _bin(tmp_path)}
    _submit(tmp_path, p)                            # failed: no results.json
    assert _submit(tmp_path, p) == "42"


def test_submit_refuses_to_overwrite_finished_results_with_other_physics(tmp_path):
    p = P | {"_binary": _bin(tmp_path)}
    _submit(tmp_path, p)
    (tmp_path / "runs" / "a" / "results.json").write_text("{}")
    _refused(tmp_path, p | {"omega_b": 9.0}, exc=ValueError)


def test_submit_refuses_an_unknown_binary(tmp_path):
    e = _refused(tmp_path, P | {"_binary": str(tmp_path / "gone")}, exc=ValueError)
    assert "binary" in str(e)
    assert not (tmp_path / "runs" / "a").exists()


def test_submit_with_unknown_binary_allowed_is_never_a_duplicate(tmp_path):
    p = P | {"_binary": str(tmp_path / "gone")}
    _submit(tmp_path, p, allow_unknown_binary=True)
    assert _submit(tmp_path, p | {"run_id": "b"}, allow_unknown_binary=True) == "42"


def test_restart_key_follows_the_dump_s_run(tmp_path):
    b = _bin(tmp_path)
    _submit(tmp_path, P | {"_binary": b})
    dump = tmp_path / "runs" / "a" / "checkpoint.dump"
    dump.write_text("dump")
    (dump.parent / "checkpoint.phase").write_text('{"t": 1.0, "omega_b": 1.0}')
    r = P | {"run_id": "r", "_binary": b, "t_checkpoint": 1.0}
    with patch("scripts.simulate._check_restart_forcing"):
        _submit(tmp_path, r, checkpoint=str(dump))
    entry = {e["run_id"]: e for e in rk.read_index(tmp_path / "runs")}["r"]
    assert entry["parent_run"] == "a"
    assert entry["parent_key"] == {e["run_id"]: e for e in rk.read_index(tmp_path / "runs")}["a"]["key"]


def test_cli_reports_a_duplicate_without_submitting(tmp_path):
    """`simulate.py PARAMS --slurm` (what `make submit` runs) stops at a duplicate: exit 3, the existing run named."""
    import subprocess
    import sys
    p = P | {"_binary": _bin(tmp_path)}
    _submit(tmp_path, p)
    (tmp_path / "runs" / "a" / "results.json").write_text("{}")
    f = tmp_path / "b.json"
    f.write_text(json.dumps(p | {"run_id": "b"}))
    r = subprocess.run([sys.executable, str(PROJECT_ROOT / "scripts" / "simulate.py"), str(f), "--slurm",
                        "--runs-root", str(tmp_path / "runs")], capture_output=True, text=True,
                       env={**os.environ, "PATH": str(tmp_path)})          # no sbatch on PATH
    assert r.returncode == 3, r.stderr
    assert "same run key as a" in r.stdout + r.stderr
