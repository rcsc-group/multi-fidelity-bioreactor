"""Tests for simulate.py — local run + SLURM job submission.

simulate.run_local(params, project_root) runs the binary directly.
simulate.submit_slurm(params, project_root, **sbatch_kwargs) submits via sbatch
  and returns a job_id string.
simulate.wait_for_result(run_dir, timeout, poll) blocks until results.json appears.
simulate.run_trial(params, backend, ...) → float  black-box: submit + wait + return kLa.
"""
import json
import math
import subprocess
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from scripts.simulate import (
    run_local, submit_slurm, wait_for_result, run_trial, get_job_elapsed_seconds,
)
from tests.conftest import CANONICAL_PARAMS, PROJECT_ROOT

SHORT_PARAMS = {**CANONICAL_PARAMS, "run_id": "sim_test", "fidelity": 3}


# ── run_local ─────────────────────────────────────────────────────────────────

@pytest.mark.medium
def test_run_local_creates_run_dir(tmp_path):
    """run_local creates runs/{run_id}/ under project_root."""
    run_local(SHORT_PARAMS, project_root=PROJECT_ROOT, runs_root=tmp_path)
    assert (tmp_path / SHORT_PARAMS["run_id"]).is_dir()


@pytest.mark.medium
def test_run_local_writes_params_json(tmp_path):
    """run_local writes params.json into the run directory."""
    run_local(SHORT_PARAMS, project_root=PROJECT_ROOT, runs_root=tmp_path)
    p = tmp_path / SHORT_PARAMS["run_id"] / "params.json"
    assert p.exists()
    loaded = json.loads(p.read_text())
    assert loaded["omega_b"] == SHORT_PARAMS["omega_b"]


@pytest.mark.medium
def test_run_local_returns_run_dir(tmp_path):
    """run_local returns the Path of the run directory."""
    run_dir = run_local(SHORT_PARAMS, project_root=PROJECT_ROOT, runs_root=tmp_path)
    assert isinstance(run_dir, Path)
    assert run_dir == tmp_path / SHORT_PARAMS["run_id"]


# ── submit_slurm ──────────────────────────────────────────────────────────────

def test_submit_slurm_calls_sbatch(tmp_path):
    """submit_slurm calls sbatch with the correct --export=PARAMS=... argument."""
    fake_result = MagicMock()
    fake_result.stdout = "Submitted batch job 123456\n"
    fake_result.returncode = 0

    with patch("subprocess.run", return_value=fake_result) as mock_run:
        job_id = submit_slurm(SHORT_PARAMS, project_root=PROJECT_ROOT, runs_root=tmp_path,
                               mpi_scratch_root=tmp_path / "mpi_scratch")

    assert mock_run.called
    cmd = mock_run.call_args[0][0]
    assert "sbatch" in cmd[0]
    assert any("PARAMS=" in a for a in cmd), f"No PARAMS= in sbatch args: {cmd}"
    assert job_id == "123456"


def test_submit_slurm_passes_walltime(tmp_path):
    """submit_slurm forwards walltime to sbatch as --time."""
    fake_result = MagicMock()
    fake_result.stdout = "Submitted batch job 42\n"
    fake_result.returncode = 0

    with patch("subprocess.run", return_value=fake_result) as mock_run:
        submit_slurm(SHORT_PARAMS, project_root=PROJECT_ROOT, runs_root=tmp_path,
                     walltime="02:30:00", mpi_scratch_root=tmp_path / "mpi_scratch")

    cmd = mock_run.call_args[0][0]
    assert any("02:30:00" in a for a in cmd), f"walltime not in sbatch args: {cmd}"


# ── wait_for_result ───────────────────────────────────────────────────────────

def test_wait_for_result_returns_when_file_appears(tmp_path):
    """wait_for_result returns results dict once results.json is written."""
    run_dir = tmp_path / "run_wait"
    run_dir.mkdir()
    expected = {"kLa_10": 0.005, "kLa_25": 0.008, "kLa_50": 0.012}

    # write results.json after a brief delay (simulates job completion)
    def _write_late():
        time.sleep(0.2)
        (run_dir / "results.json").write_text(json.dumps(expected))

    import threading
    t = threading.Thread(target=_write_late)
    t.start()

    result = wait_for_result(run_dir, timeout=5, poll=0.05)
    t.join()

    assert result["kLa_25"] == pytest.approx(0.008)


def test_wait_for_result_raises_on_timeout(tmp_path):
    """wait_for_result raises TimeoutError if results.json never appears."""
    run_dir = tmp_path / "run_timeout"
    run_dir.mkdir()

    with pytest.raises(TimeoutError):
        wait_for_result(run_dir, timeout=0.3, poll=0.05)


# ── get_job_elapsed_seconds ────────────────────────────────────────────────────

@pytest.mark.parametrize("sacct_elapsed,expected_seconds", [
    ("00:05:30", 5 * 60 + 30),
    ("01:02:03", 1 * 3600 + 2 * 60 + 3),
    ("1-02:03:04", 1 * 86400 + 2 * 3600 + 3 * 60 + 4),
    ("3-00:00:00", 3 * 86400),
])
def test_get_job_elapsed_seconds_parses_slurm_format(sacct_elapsed, expected_seconds):
    """SLURM's Elapsed format is [D-]HH:MM:SS; must parse both with and
    without a day component."""
    fake_result = MagicMock()
    fake_result.stdout = f"{sacct_elapsed}\n"
    fake_result.returncode = 0

    with patch("subprocess.run", return_value=fake_result) as mock_run:
        seconds = get_job_elapsed_seconds("123456")

    assert mock_run.called
    cmd = mock_run.call_args[0][0]
    assert "sacct" in cmd[0]
    assert "123456" in cmd
    assert seconds == expected_seconds


def test_get_job_elapsed_seconds_returns_none_on_empty_output():
    """If sacct returns nothing (job not yet in accounting), return None
    rather than raising -- callers treat a missing wall time as optional."""
    fake_result = MagicMock()
    fake_result.stdout = "\n"
    fake_result.returncode = 0

    with patch("subprocess.run", return_value=fake_result):
        assert get_job_elapsed_seconds("999999") is None


# ── run_trial ─────────────────────────────────────────────────────────────────

class TestRunTrial:
    _RESULTS = {"kLa_10": 0.081, "kLa_25": 0.083, "kLa_50": 0.060}

    def test_run_trial_local_returns_scalar(self, tmp_path):
        """run_trial(backend='local') returns the kLa_25 scalar as a float."""
        def _fake_run_local(params, project_root, runs_root, **kw):
            run_dir = runs_root / params["run_id"]
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "results.json").write_text(json.dumps(self._RESULTS))
            return run_dir

        with patch("scripts.simulate.run_local", side_effect=_fake_run_local):
            result = run_trial(
                SHORT_PARAMS, backend="local",
                project_root=PROJECT_ROOT, runs_root=tmp_path,
            )

        assert isinstance(result, float)
        assert result == pytest.approx(0.083)

    def test_run_trial_slurm_returns_scalar(self, tmp_path):
        """run_trial(backend='slurm') submits, waits, and returns the kLa_25 scalar."""
        run_dir = tmp_path / SHORT_PARAMS["run_id"]
        run_dir.mkdir(parents=True)

        fake_sbatch = MagicMock()
        fake_sbatch.stdout = "Submitted batch job 99\n"
        fake_sbatch.returncode = 0

        def _write_results(*a, **kw):
            (run_dir / "results.json").write_text(json.dumps(self._RESULTS))
            return fake_sbatch

        with patch("scripts.simulate.submit_slurm", return_value="99") as mock_sub, \
             patch("scripts.simulate.wait_for_result", return_value=self._RESULTS) as mock_wait:
            result = run_trial(
                SHORT_PARAMS, backend="slurm",
                project_root=PROJECT_ROOT, runs_root=tmp_path,
            )

        mock_sub.assert_called_once()
        mock_wait.assert_called_once()
        assert result == pytest.approx(0.083)

    def test_run_trial_returns_nan_on_missing_key(self, tmp_path):
        """run_trial returns NaN when the requested kla_key is absent — does not raise."""
        def _fake_run_local(params, project_root, runs_root, **kw):
            run_dir = runs_root / params["run_id"]
            run_dir.mkdir(parents=True, exist_ok=True)
            (run_dir / "results.json").write_text(json.dumps({"kLa_10": 0.08}))
            return run_dir

        with patch("scripts.simulate.run_local", side_effect=_fake_run_local):
            result = run_trial(
                SHORT_PARAMS, backend="local", kla_key="kLa_25",
                project_root=PROJECT_ROOT, runs_root=tmp_path,
            )

        assert math.isnan(result)

    def test_run_trial_slurm_timeout_raises(self, tmp_path):
        """run_trial raises TimeoutError when results.json never appears."""
        with patch("scripts.simulate.submit_slurm", return_value="77"), \
             patch("scripts.simulate.wait_for_result",
                   side_effect=TimeoutError("timed out")):
            with pytest.raises(TimeoutError):
                run_trial(
                    SHORT_PARAMS, backend="slurm",
                    project_root=PROJECT_ROOT, runs_root=tmp_path,
                )
