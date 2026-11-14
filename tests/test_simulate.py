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
    binary_reads_label,
    source_forcing_keys,
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


# ── a restart must say which forcing it continues (pitfall "restart-ramp") ────
# Every restart ramps the forcing from its source to params.json. Binaries built since 2026-11-13 read the
# source from the dump's label (checkpoint.phase) and refuse *_prev keys; older binaries read *_prev and, with
# none, ramp from zero rocking (fig9_l9_rpm37.5_ext1, 2026-11-13). Refused before anything is written.

LABEL = ('{"t": 10.0, "cycles": 4.0, "U_bio": 0.1, "T_bio": 1.0, "omega_b": 3.4, "theta_max": [7, 0, 0], '
         '"Th": 0.0, "Th_d": 0.4}')


def _restart(**extra):
    return {**SHORT_PARAMS, "run_id": "restart_test", "t_checkpoint": 10.0, "t_end": 5.0,
            "restart_continue": 1, "n_mix_cycles": 0, **extra}


def _prev(p):
    return {"omega_b_prev": p["omega_b"], "theta_max_prev": list(p["theta_max"]),
            "phi_angular_prev": list(p["phi_angular"]), "amplitude_h_prev": list(p["amplitude_h"]),
            "phi_horizontal_prev": list(p["phi_horizontal"]), "omega_h_prev": p["omega_h"]}


def _binary(tmp_path, reads_label):
    b = tmp_path / ("new_bin" if reads_label else "old_bin")
    b.write_bytes(b"\x7fELF..." + (b"restart refused" if reads_label else b"") + b"...")
    return str(b)


def _dump(tmp_path, label=LABEL):
    d = tmp_path / "ck"
    d.mkdir(exist_ok=True)
    (d / "checkpoint.dump").write_bytes(b"x")
    if label is not None:
        (d / "checkpoint.phase").write_text(label)
    return str(d / "checkpoint.dump")


def _submit(tmp_path, params, checkpoint=None):
    ok = MagicMock(stdout="Submitted batch job 42\n", returncode=0)
    with patch("subprocess.run", return_value=ok):
        return submit_slurm(params, project_root=PROJECT_ROOT, runs_root=tmp_path / "runs",
                            mpi_scratch_root=tmp_path / "mpi", checkpoint=checkpoint)


def _refused(tmp_path, params, checkpoint=None, match="restart-ramp"):
    with patch("subprocess.run", side_effect=AssertionError("sbatch must not be called")):
        with pytest.raises(ValueError, match=match):
            submit_slurm(params, project_root=PROJECT_ROOT, runs_root=tmp_path / "runs",
                         mpi_scratch_root=tmp_path / "mpi", checkpoint=checkpoint)
    assert not (tmp_path / "runs" / params["run_id"]).exists()


def test_binary_reads_label_from_its_marker(tmp_path):
    assert binary_reads_label(_binary(tmp_path, True))
    assert not binary_reads_label(_binary(tmp_path, False))
    assert binary_reads_label(None)                       # the binary built from this src/


def test_old_binary_restart_without_prev_is_refused(tmp_path):
    _refused(tmp_path, _restart(_binary=_binary(tmp_path, False)))


def test_old_binary_restart_with_prev_is_submitted(tmp_path):
    p = _restart(_binary=_binary(tmp_path, False))
    assert _submit(tmp_path, p | _prev(p)) == "42"


def test_label_binary_restart_with_label_is_submitted(tmp_path):
    assert _submit(tmp_path, _restart(_binary=_binary(tmp_path, True)), _dump(tmp_path)) == "42"


def test_label_binary_refuses_prev_next_to_a_label(tmp_path):
    p = _restart(_binary=_binary(tmp_path, True))
    _refused(tmp_path, p | _prev(p), _dump(tmp_path), match="two sources")


def test_label_binary_refuses_a_dump_without_label(tmp_path):
    _refused(tmp_path, _restart(_binary=_binary(tmp_path, True)), _dump(tmp_path, label=None),
             match="checkpoint.phase")


def test_label_binary_takes_prev_with_an_old_label(tmp_path):
    p = _restart(_binary=_binary(tmp_path, True))
    old = '{"t": 10.0, "t_phase_offset": 0.0, "T_per_st": 2.5}'
    assert _submit(tmp_path, p | _prev(p), _dump(tmp_path, old)) == "42"
    _refused(tmp_path, p | {"run_id": "restart_test2"}, _dump(tmp_path, old), match="old label")


def test_mpi_staging_copies_the_label_with_the_dump(tmp_path):
    ck = _dump(tmp_path)
    ok = MagicMock(stdout="Submitted batch job 44\n", returncode=0)
    with patch("subprocess.run", return_value=ok):
        submit_slurm(_restart(_binary=_binary(tmp_path, True)), project_root=PROJECT_ROOT,
                     runs_root=tmp_path / "runs", mpi_scratch_root=tmp_path / "mpi", checkpoint=ck,
                     template=PROJECT_ROOT / "config" / "slurm_mpi_template.sh")
    assert (tmp_path / "mpi" / "restart_test" / "checkpoint.phase").read_text() == LABEL


def test_source_forcing_keys_follow_the_binary(tmp_path):
    p = dict(SHORT_PARAMS)
    assert source_forcing_keys(p, _binary(tmp_path, True)) == {}
    assert source_forcing_keys(p, _binary(tmp_path, False)) == _prev(p)


def test_replay_restart_needs_no_prev(tmp_path):
    """A replay run solves only scalars on a recorded flow; the forcing does not drive it."""
    ok = MagicMock(stdout="Submitted batch job 43\n", returncode=0)
    with patch("subprocess.run", return_value=ok):
        assert submit_slurm(_restart(rec_dir="/x"), project_root=PROJECT_ROOT, runs_root=tmp_path,
                            mpi_scratch_root=tmp_path / "mpi") == "43"


def test_submit_slurm_serial_restart_stages_the_label(tmp_path):
    """Non-MPI restart: the solver runs in runs/<id>/ and reads checkpoint.phase there, so the dump's label must
    be copied next to params.json (the MPI branch copies it to scratch). Review finding 2026-11-14."""
    src = tmp_path / "src_run"
    src.mkdir()
    (src / "checkpoint.dump").write_bytes(b"\0" * 16)
    (src / "checkpoint.phase").write_text('{"t": 5.0, "omega_b": 3.0}')
    fake = MagicMock(stdout="Submitted batch job 1\n", returncode=0)
    p = {**SHORT_PARAMS, "run_id": "serial_restart", "t_checkpoint": 5.0}
    with patch("subprocess.run", return_value=fake):
        submit_slurm(p, project_root=PROJECT_ROOT, runs_root=tmp_path / "runs", checkpoint=str(src / "checkpoint.dump"),
                     template=PROJECT_ROOT / "config" / "slurm_template.sh", mpi_scratch_root=tmp_path / "mpi")
    staged = tmp_path / "runs" / "serial_restart" / "checkpoint.phase"
    assert staged.is_file() and staged.read_text() == (src / "checkpoint.phase").read_text()
