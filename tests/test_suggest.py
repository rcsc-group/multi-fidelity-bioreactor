"""Tests for suggest.py — acquisition function and candidate generation."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from f3dasm import ExperimentData

from scripts.suggest import _compute_y_best, main as suggest_main

N_MAX = 3


def _make_experiment_data(tmp_dir: Path, n_lf: int = 8, n_hf: int = 3) -> Path:
    rng = np.random.default_rng(0)
    rows_in, rows_out = [], []
    for fidelity, n in ((5, n_lf), (7, n_hf)):
        for i in range(n):
            rows_in.append({
                "omega_b": rng.uniform(1.57, 6.28),
                "n_harmonics": 1,
                "theta_max_0": rng.uniform(2, 7), "theta_max_1": 0.0, "theta_max_2": 0.0,
                "phi_angular_0": 0.0, "phi_angular_1": 0.0, "phi_angular_2": 0.0,
                "omega_h": 0.0,
                "amplitude_h_0": 0.0, "amplitude_h_1": 0.0, "amplitude_h_2": 0.0,
                "phi_horizontal_0": 0.0, "phi_horizontal_1": 0.0, "phi_horizontal_2": 0.0,
                "geometry_a": 0.25, "geometry_b": 0.071, "geometry_n": rng.uniform(2, 8),
                "fill_level": rng.uniform(0.3, 0.7),
                "fidelity": fidelity,
            })
            rows_out.append({
                "kLa_10": rng.uniform(0.002, 0.015),
                "kLa_25": rng.uniform(0.003, 0.02),
                "kLa_50": rng.uniform(0.002, 0.018),
            })
    ed = ExperimentData(
        input_data=pd.DataFrame(rows_in),
        output_data=pd.DataFrame(rows_out),
    )
    exp_dir = tmp_dir / "experiment"
    ed.store(str(exp_dir))
    return exp_dir


def test_candidates_within_bounds(tmp_path):
    """All suggested parameters are within param_space.yaml bounds."""
    exp_dir = _make_experiment_data(tmp_path)
    param_space = Path("config/param_space.yaml")

    candidate = suggest_main(str(exp_dir), str(param_space))

    assert 1.57 <= candidate["omega_b"] <= 6.28
    assert 0.3 <= candidate["fill_level"] <= 0.7


def test_phi_angular_zero_always_fixed(tmp_path):
    """phi_angular[0] == 0.0 in every suggested candidate."""
    exp_dir = _make_experiment_data(tmp_path)
    param_space = Path("config/param_space.yaml")

    candidate = suggest_main(str(exp_dir), str(param_space))

    assert candidate["phi_angular"][0] == 0.0


def test_vectors_zero_padded_to_nmax(tmp_path):
    """theta_max, phi_angular, amplitude_h, phi_horizontal all have length N_max."""
    exp_dir = _make_experiment_data(tmp_path)
    param_space = Path("config/param_space.yaml")

    candidate = suggest_main(str(exp_dir), str(param_space))

    for vec in ("theta_max", "phi_angular", "amplitude_h", "phi_horizontal"):
        assert len(candidate[vec]) == N_MAX, f"{vec} must have length {N_MAX}"


class _FakeModel:
    """Returns a fixed HF-scale prediction regardless of input, for testing
    that y_best is computed via the surrogate (not raw LF observations) when
    no HF data exists yet."""

    def __init__(self, mean_value: float):
        self._mean_value = mean_value

    def predict(self, X: np.ndarray):
        n = X.shape[0]
        return np.full(n, self._mean_value), np.full(n, 0.01)


def test_y_best_uses_hf_observations_when_available():
    """With HF rows present, y_best is the raw observed HF max (unchanged
    behaviour) — the surrogate should not be consulted at all here."""
    inp_df = pd.DataFrame({
        "fidelity": [5, 5, 7, 7],
        "omega_b": [2.0, 3.0, 2.5, 3.5],
    })
    out_df = pd.DataFrame({"kLa_25": [0.001, 0.002, 0.010, 0.015]})
    fake_model = _FakeModel(mean_value=999.0)  # would be very wrong if used

    y_best = _compute_y_best(
        inp_df, out_df, fake_model, kla_key="kLa_25",
        hf_fidelity=7, lf_fidelity=5, feature_cols=["omega_b"],
    )

    assert y_best == pytest.approx(0.015)


def test_y_best_does_not_mix_lf_scale_when_no_hf_observations():
    """Root cause of the bug: with zero HF observations, y_best must NOT
    fall back to a raw LF-observed value, since LF and HF are expected to
    have a systematic bias (that's the entire premise of the multi-fidelity
    correction). It must instead use the surrogate's own HF-scale prediction
    evaluated at the LF-observed inputs."""
    inp_df = pd.DataFrame({
        "fidelity": [5, 5, 5],
        "omega_b": [2.0, 3.0, 4.0],
    })
    # Raw LF observations are all much larger than the surrogate's HF-scale
    # estimate -- if the bug were present, y_best would equal 0.5 (the raw
    # LF max), not the model-predicted value.
    out_df = pd.DataFrame({"kLa_25": [0.3, 0.4, 0.5]})
    fake_model = _FakeModel(mean_value=0.02)

    y_best = _compute_y_best(
        inp_df, out_df, fake_model, kla_key="kLa_25",
        hf_fidelity=7, lf_fidelity=5, feature_cols=["omega_b"],
    )

    assert y_best == pytest.approx(0.02)
    assert y_best != pytest.approx(0.5)
