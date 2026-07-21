"""Tests for scripts/loop.py — wall-clock time tracking through the BO loop.

Wall-clock time isn't otherwise available in results.json (postprocess.py
only knows about simulated quantities, not the SLURM job's real-world
duration), so it must be attached at the point where loop.py still has the
job_id in scope, then persisted into ExperimentData so a future acquisition
function can trade off quality against cost.
"""
from pathlib import Path
from unittest.mock import patch

from scripts.loop import _submit_and_wait, _append_to_ed
from f3dasm import ExperimentData
import pandas as pd


_CFG = {"walltime": "01:00:00", "job_timeout": 10}


def test_submit_and_wait_attaches_wall_time_s(tmp_path):
    """_submit_and_wait adds a wall_time_s key sourced from sacct, alongside
    the simulation's own results.json contents."""
    params = {"run_id": "abc123"}

    with patch("scripts.simulate.submit_slurm", return_value="42") as mock_submit, \
         patch("scripts.simulate.wait_for_result",
               return_value={"kLa_25": 0.01}) as mock_wait, \
         patch("scripts.simulate.get_job_elapsed_seconds", return_value=321.0) as mock_elapsed, \
         patch("scripts.loop._PROJECT_ROOT", tmp_path):
        results = _submit_and_wait(params, _CFG)

    mock_elapsed.assert_called_once_with("42")
    assert results["kLa_25"] == 0.01
    assert results["wall_time_s"] == 321.0


def test_submit_and_wait_tolerates_missing_wall_time(tmp_path):
    """If sacct has no record yet, wall_time_s is None rather than missing
    entirely or raising -- the run's real results must not be discarded."""
    params = {"run_id": "abc123"}

    with patch("scripts.simulate.submit_slurm", return_value="42"), \
         patch("scripts.simulate.wait_for_result", return_value={"kLa_25": 0.01}), \
         patch("scripts.simulate.get_job_elapsed_seconds", return_value=None), \
         patch("scripts.loop._PROJECT_ROOT", tmp_path):
        results = _submit_and_wait(params, _CFG)

    assert results["wall_time_s"] is None


def test_append_to_ed_stores_wall_time_s_column(tmp_path):
    """The ExperimentData output row includes wall_time_s so it's available
    to a future cost-aware acquisition function."""
    params = {
        "run_id": "r1", "omega_b": 2.0, "n_harmonics": 1,
        "theta_max": [5.0, 0.0, 0.0], "phi_angular": [0.0, 0.0, 0.0],
        "omega_h": 0.0, "amplitude_h": [0.0, 0.0, 0.0],
        "phi_horizontal": [0.0, 0.0, 0.0],
        "geometry": {"a": 0.25, "b": 0.071, "n": 4.0},
        "fill_level": 0.5,
    }
    results = {"kLa_25": 0.01, "wall_time_s": 1234.0}

    ed = ExperimentData(
        input_data=pd.DataFrame(),
        output_data=pd.DataFrame(),
    )
    ed = _append_to_ed(ed, params, results, fidelity=7, phase="bo")

    _, out_df = ed.to_pandas()
    assert "wall_time_s" in out_df.columns
    assert out_df["wall_time_s"].iloc[0] == 1234.0
