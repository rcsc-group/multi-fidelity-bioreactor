"""Tests for train_surrogate.py — MF-BML surrogate training."""
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from f3dasm import ExperimentData
from f3dasm.design import Domain

from scripts.train_surrogate import main as train_main

N_FEATURES = 18  # flat feature vector (20 input cols minus phi_angular_0 and fidelity)

INPUT_COLS = [
    "omega_b", "n_harmonics",
    "theta_max_0", "theta_max_1", "theta_max_2",
    "phi_angular_0", "phi_angular_1", "phi_angular_2",
    "omega_h",
    "amplitude_h_0", "amplitude_h_1", "amplitude_h_2",
    "phi_horizontal_0", "phi_horizontal_1", "phi_horizontal_2",
    "geometry_a", "geometry_b", "geometry_n",
    "fill_level", "fidelity",
]


def _make_experiment_data(tmp_dir: Path, n_lf: int = 8, n_hf: int = 3) -> Path:
    rng = np.random.default_rng(42)
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
                "fill_level": 0.5,
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


def test_surrogate_trains_on_synthetic_data(tmp_path):
    """Given synthetic LF+HF ExperimentData, model.pkl is written without error."""
    exp_dir = _make_experiment_data(tmp_path)
    model_path = tmp_path / "model.pkl"

    train_main(str(exp_dir), str(model_path))

    assert model_path.exists()


def test_surrogate_prediction_returns_mean_and_variance(tmp_path):
    """Loaded model.pkl predicts (mean, variance) tuple for a new input."""
    exp_dir = _make_experiment_data(tmp_path)
    model_path = tmp_path / "model.pkl"
    train_main(str(exp_dir), str(model_path))

    with open(model_path, "rb") as f:
        model = pickle.load(f)

    x_new = np.zeros((1, N_FEATURES))
    mean, var = model.predict(x_new)
    assert mean.shape == (1,)
    assert var.shape == (1,)
    assert var[0] >= 0.0
