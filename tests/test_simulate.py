"""Tests for simulate.py — local run + SLURM job submission.

simulate.run_local(params, project_root) runs the binary directly.
simulate.submit_slurm(params, project_root, **sbatch_kwargs) submits via sbatch
  and returns a job_id string.
simulate.wait_for_result(run_dir, timeout, poll) blocks until results.json appears.
"""
import json
import math
import subprocess
import time
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from scripts.simulate import run_local, submit_slurm, wait_for_result
from tests.conftest import CANONICAL_PARAMS, PROJECT_ROOT

SHORT_PARAMS = {**CANONICAL_PARAMS, "run_id": "sim_test", "fidelity": 3}


# ── run_local ─────────────────────────────────────────────────────────────────

def test_run_local_creates_run_dir(tmp_path):
    """run_local creates runs/{run_id}/ under project_root."""
    run_local(SHORT_PARAMS, project_root=PROJECT_ROOT, runs_root=tmp_path)
    assert (tmp_path / SHORT_PARAMS["run_id"]).is_dir()


def test_run_local_writes_params_json(tmp_path):
    """run_local writes params.json into the run directory."""
    run_local(SHORT_PARAMS, project_root=PROJECT_ROOT, runs_root=tmp_path)
    p = tmp_path / SHORT_PARAMS["run_id"] / "params.json"
    assert p.exists()
    loaded = json.loads(p.read_text())
    assert loaded["omega_b"] == SHORT_PARAMS["omega_b"]


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
        job_id = submit_slurm(SHORT_PARAMS, project_root=PROJECT_ROOT, runs_root=tmp_path)

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
                     walltime="02:30:00")

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
