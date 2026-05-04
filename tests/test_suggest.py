"""Tests for suggest.py — acquisition function and candidate generation."""
import csv
import json
import math
from pathlib import Path

import pytest

from scripts.suggest import main as suggest_main

N_MAX = 3


def _write_minimal_results_csv(path: Path, n_lf: int = 5):
    """Write a minimal results.csv with synthetic LF runs."""
    rows = []
    for i in range(n_lf):
        rows.append({
            "run_id": f"lf_{i}",
            "fidelity": 7,
            "omega_b": 2.0 + i * 0.5,
            "n_harmonics": 1,
            "theta_max_0": 5.0, "theta_max_1": 0.0, "theta_max_2": 0.0,
            "phi_angular_0": 0.0, "phi_angular_1": 0.0, "phi_angular_2": 0.0,
            "omega_h": 0.0,
            "amplitude_h_0": 0.0, "amplitude_h_1": 0.0, "amplitude_h_2": 0.0,
            "phi_horizontal_0": 0.0, "phi_horizontal_1": 0.0, "phi_horizontal_2": 0.0,
            "geometry_a": 0.25, "geometry_b": 0.071, "geometry_tilt": 0.0,
            "fill_level": 0.5,
            "kLa_25": 0.005 + i * 0.001,
        })
    fieldnames = rows[0].keys()
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_candidates_within_bounds(tmp_path):
    """All suggested parameters are within param_space.yaml bounds."""
    results_csv = tmp_path / "results.csv"
    _write_minimal_results_csv(results_csv)
    param_space = Path("config/param_space.yaml")

    candidate = suggest_main(str(results_csv), str(param_space))

    assert 1.57 <= candidate["omega_b"] <= 6.28
    assert 0.3 <= candidate["fill_level"] <= 0.7


def test_phi_angular_zero_always_fixed(tmp_path):
    """phi_angular[0] == 0.0 in every suggested candidate."""
    results_csv = tmp_path / "results.csv"
    _write_minimal_results_csv(results_csv)
    param_space = Path("config/param_space.yaml")

    candidate = suggest_main(str(results_csv), str(param_space))

    assert candidate["phi_angular"][0] == 0.0


def test_vectors_zero_padded_to_nmax(tmp_path):
    """theta_max, phi_angular, amplitude_h, phi_horizontal all have length N_max."""
    results_csv = tmp_path / "results.csv"
    _write_minimal_results_csv(results_csv)
    param_space = Path("config/param_space.yaml")

    candidate = suggest_main(str(results_csv), str(param_space))

    for vec in ("theta_max", "phi_angular", "amplitude_h", "phi_horizontal"):
        assert len(candidate[vec]) == N_MAX, f"{vec} must have length {N_MAX}"
