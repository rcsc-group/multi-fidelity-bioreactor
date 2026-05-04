"""Tests for train_surrogate.py — MF-BML surrogate training."""
import csv
import pickle
from pathlib import Path

import numpy as np
import pytest

from scripts.train_surrogate import main as train_main

N_FEATURES = 14  # flat feature vector length (see param_space.yaml)


def _write_synthetic_results(path: Path, n_lf: int = 8, n_hf: int = 3):
    rng = np.random.default_rng(42)
    rows = []
    for fidelity, n in ((7, n_lf), (9, n_hf)):
        for i in range(n):
            row = {
                "run_id": f"f{fidelity}_{i}",
                "fidelity": fidelity,
                "omega_b": rng.uniform(1.57, 6.28),
                "n_harmonics": 1,
                "theta_max_0": rng.uniform(2, 7), "theta_max_1": 0.0, "theta_max_2": 0.0,
                "phi_angular_0": 0.0, "phi_angular_1": 0.0, "phi_angular_2": 0.0,
                "omega_h": 0.0,
                "amplitude_h_0": 0.0, "amplitude_h_1": 0.0, "amplitude_h_2": 0.0,
                "phi_horizontal_0": 0.0, "phi_horizontal_1": 0.0, "phi_horizontal_2": 0.0,
                "geometry_a": 0.25, "geometry_b": 0.071, "geometry_tilt": 0.0,
                "fill_level": 0.5,
                "kLa_25": rng.uniform(0.003, 0.02),
            }
            rows.append(row)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def test_surrogate_trains_on_synthetic_data(tmp_path):
    """Given synthetic LF+HF results.csv, model.pkl is written without error."""
    results_csv = tmp_path / "results.csv"
    _write_synthetic_results(results_csv)
    model_path = tmp_path / "model.pkl"

    train_main(str(results_csv), str(model_path))

    assert model_path.exists()


def test_surrogate_prediction_returns_mean_and_variance(tmp_path):
    """Loaded model.pkl predicts (mean, variance) tuple for a new input."""
    results_csv = tmp_path / "results.csv"
    _write_synthetic_results(results_csv)
    model_path = tmp_path / "model.pkl"
    train_main(str(results_csv), str(model_path))

    with open(model_path, "rb") as f:
        model = pickle.load(f)

    x_new = np.zeros((1, N_FEATURES))
    mean, var = model.predict(x_new)
    assert mean.shape == (1,)
    assert var.shape == (1,)
    assert var[0] >= 0.0
