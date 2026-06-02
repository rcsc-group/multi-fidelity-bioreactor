"""Unit tests for scripts/sample.py.

All tests run without submitting SLURM jobs or executing simulations.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parents[1]))

from scripts.sample import _row_to_params, run_sampling

# ── fixtures ──────────────────────────────────────────────────────────────────

N_MAX = 3  # matches config/param_space.yaml

def _make_row(**overrides) -> pd.Series:
    """Minimal valid ExperimentData input row."""
    data = {
        "omega_b":         3.93,
        "n_harmonics":     1.0,
        "theta_max_0":     7.0,
        "theta_max_1":     0.0,
        "theta_max_2":     0.0,
        "phi_angular_0":   0.0,
        "phi_angular_1":   0.5,
        "phi_angular_2":   1.0,
        "omega_h":         0.0,
        "amplitude_h_0":   0.0,
        "amplitude_h_1":   0.0,
        "amplitude_h_2":   0.0,
        "phi_horizontal_0": 0.0,
        "phi_horizontal_1": 0.0,
        "phi_horizontal_2": 0.0,
        "geometry_a":      0.25,
        "geometry_b":      0.071,
        "geometry_n":      8.0,
        "fill_level":      0.5,
    }
    data.update(overrides)
    return pd.Series(data)


def _minimal_cfg(tmp_path: Path, strategy: str = "latin",
                 n_samples: int = 3, submit: bool = False) -> dict:
    return {
        "experiment_dir": str(tmp_path / "exp"),
        "fidelity": 3,
        "t_buffer": 5.0,
        "walltime": "00:10:00",
        "submit": submit,
        "sampling": {
            "strategy": strategy,
            "n_samples": n_samples,
            "seed": 42,
        },
    }


# ── _row_to_params ────────────────────────────────────────────────────────────

class TestRowToParams:
    def test_returns_dict_with_required_keys(self):
        row = _make_row()
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        for key in ("run_id", "fidelity", "omega_b", "n_harmonics",
                    "theta_max", "phi_angular", "phi_horizontal",
                    "omega_h", "amplitude_h", "geometry", "fill_level", "t_end"):
            assert key in p, f"Missing key: {key}"

    def test_run_id_is_8_char_hex(self):
        row = _make_row()
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert len(p["run_id"]) == 8
        int(p["run_id"], 16)  # raises ValueError if not hex

    def test_run_id_unique_across_calls(self):
        row = _make_row()
        ids = {_row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)["run_id"]
               for _ in range(20)}
        assert len(ids) == 20, "run_ids should be unique"

    def test_fidelity_passed_through(self):
        row = _make_row()
        p = _row_to_params(row, fidelity=5, t_buffer=5.0, n_max=N_MAX)
        assert p["fidelity"] == 5

    def test_omega_b_passed_through(self):
        row = _make_row(omega_b=2.71)
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert p["omega_b"] == pytest.approx(2.71)

    def test_fill_level_passed_through(self):
        row = _make_row(fill_level=0.4)
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert p["fill_level"] == pytest.approx(0.4)

    def test_theta_max_is_list_of_n_max(self):
        row = _make_row(theta_max_0=5.0, theta_max_1=2.0, theta_max_2=1.0)
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert isinstance(p["theta_max"], list)
        assert len(p["theta_max"]) == N_MAX
        assert p["theta_max"][0] == pytest.approx(5.0)
        assert p["theta_max"][1] == pytest.approx(2.0)

    def test_phi_angular_index0_always_zero(self):
        # phi_angular_0 in the row should be ignored; always set to 0.0
        row = _make_row(phi_angular_0=3.14)
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert p["phi_angular"][0] == 0.0, "phi_angular[0] must be 0 (time-origin reference)"

    def test_phi_angular_higher_indices_passed(self):
        row = _make_row(phi_angular_1=1.57, phi_angular_2=3.14)
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert p["phi_angular"][1] == pytest.approx(1.57)
        assert p["phi_angular"][2] == pytest.approx(3.14)

    def test_phi_angular_length_is_n_max(self):
        row = _make_row()
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert len(p["phi_angular"]) == N_MAX

    def test_n_harmonics_at_least_1(self):
        # Even if row value rounds to 0, output must be >= 1
        row = _make_row(n_harmonics=0.3)
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert p["n_harmonics"] >= 1

    def test_n_harmonics_rounded_to_int(self):
        row = _make_row(n_harmonics=2.6)
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert p["n_harmonics"] == 3

    def test_geometry_is_nested_dict(self):
        row = _make_row(geometry_a=0.20, geometry_b=0.06, geometry_n=4.0)
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert isinstance(p["geometry"], dict)
        assert set(p["geometry"].keys()) == {"a", "b", "n"}
        assert p["geometry"]["a"] == pytest.approx(0.20)
        assert p["geometry"]["b"] == pytest.approx(0.06)
        assert p["geometry"]["n"] == pytest.approx(4.0)

    def test_t_end_is_positive(self):
        row = _make_row()
        p = _row_to_params(row, fidelity=3, t_buffer=5.0, n_max=N_MAX)
        assert p["t_end"] > 0

    def test_t_end_increases_with_t_buffer(self):
        row = _make_row()
        p_small = _row_to_params(row, fidelity=3, t_buffer=5.0,  n_max=N_MAX)
        p_large = _row_to_params(row, fidelity=3, t_buffer=50.0, n_max=N_MAX)
        assert p_large["t_end"] > p_small["t_end"]


# ── run_sampling ──────────────────────────────────────────────────────────────

class TestRunSampling:
    def test_unknown_strategy_raises(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, strategy="bogus")
        with pytest.raises(ValueError, match="bogus"):
            run_sampling(cfg)

    def test_all_known_strategies_accepted(self, tmp_path):
        for strategy in ("latin", "random", "sobol"):
            cfg = _minimal_cfg(tmp_path / strategy, strategy=strategy, n_samples=2)
            run_sampling(cfg)  # must not raise

    def test_no_submit_writes_params_json(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, n_samples=2, submit=False)
        run_sampling(cfg)
        runs_root = Path(__file__).parents[1] / "runs"
        written = list(runs_root.glob("*/params.json"))
        assert len(written) >= 2

    def test_correct_number_of_runs_written(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, n_samples=4, submit=False)
        runs_root = Path(__file__).parents[1] / "runs"
        before = set(runs_root.glob("*/params.json"))
        run_sampling(cfg)
        after = set(runs_root.glob("*/params.json"))
        new_runs = after - before
        assert len(new_runs) == 4

    def test_params_json_has_required_keys(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, n_samples=1, submit=False)
        runs_root = Path(__file__).parents[1] / "runs"
        before = set(runs_root.glob("*/params.json"))
        run_sampling(cfg)
        after = set(runs_root.glob("*/params.json"))
        new_file = list(after - before)[0]
        params = json.loads(new_file.read_text())
        for key in ("run_id", "fidelity", "omega_b", "theta_max",
                    "geometry", "fill_level", "t_end"):
            assert key in params, f"Missing key: {key}"

    def test_unique_run_ids(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, n_samples=5, submit=False)
        runs_root = Path(__file__).parents[1] / "runs"
        before = set(runs_root.glob("*/params.json"))
        run_sampling(cfg)
        after = set(runs_root.glob("*/params.json"))
        new_files = after - before
        ids = [json.loads(f.read_text())["run_id"] for f in new_files]
        assert len(ids) == len(set(ids)), "All run_ids must be unique"

    def test_experiment_data_store_created(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, n_samples=2, submit=False)
        run_sampling(cfg)
        exp_dir = Path(cfg["experiment_dir"])
        # f3dasm stores at least an input CSV
        assert any(exp_dir.iterdir()), "Experiment dir should not be empty"

    def test_submit_false_does_not_call_sbatch(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, n_samples=2, submit=False)
        with patch("scripts.simulate.submit_slurm") as mock_sbatch:
            run_sampling(cfg)
            mock_sbatch.assert_not_called()

    def test_submit_true_calls_sbatch_n_times(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, n_samples=3, submit=True)
        with patch("scripts.simulate.submit_slurm", return_value="12345") as mock_sbatch:
            run_sampling(cfg)
            assert mock_sbatch.call_count == 3

    def test_params_within_param_space_bounds(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, n_samples=10, submit=False)
        runs_root = Path(__file__).parents[1] / "runs"
        before = set(runs_root.glob("*/params.json"))
        run_sampling(cfg)
        after = set(runs_root.glob("*/params.json"))
        for path in after - before:
            p = json.loads(path.read_text())
            assert 1.57 <= p["omega_b"] <= 6.28, f"omega_b out of bounds: {p['omega_b']}"
            assert 0.3  <= p["fill_level"] <= 0.7, f"fill_level out of bounds: {p['fill_level']}"
            assert 0.15 <= p["geometry"]["a"] <= 0.35
            assert 0.05 <= p["geometry"]["b"] <= 0.15

    def test_phi_angular_index0_always_zero_in_output(self, tmp_path):
        cfg = _minimal_cfg(tmp_path, n_samples=5, submit=False)
        runs_root = Path(__file__).parents[1] / "runs"
        before = set(runs_root.glob("*/params.json"))
        run_sampling(cfg)
        after = set(runs_root.glob("*/params.json"))
        for path in after - before:
            p = json.loads(path.read_text())
            assert p["phi_angular"][0] == 0.0, \
                f"phi_angular[0] must be 0.0, got {p['phi_angular'][0]}"
