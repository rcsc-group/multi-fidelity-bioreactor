"""Tests for postprocess.py — kLa extraction from tr_oxy.dat."""
import json
import math
import tempfile
from pathlib import Path

import numpy as np
import pytest

# Import will fail until scripts/postprocess.py exists — that's the RED state.
from scripts.postprocess import main as postprocess_main


def _synthetic_tr_oxy(kla: float, n_points: int = 50) -> np.ndarray:
    """Return synthetic tr_oxy data: columns [time, C*_w_oxy].
    C*(t) = 1 - exp(-kla * t) is the exact first-order solution.
    """
    t = np.linspace(0.1, 10.0, n_points)
    c = 1.0 - np.exp(-kla * t)
    return np.column_stack([t, c])


def test_kla_extracted_from_synthetic_data(tmp_path):
    """Given a synthetic tr_oxy.dat with known kLa=0.01, extraction is within 5%."""
    true_kla = 0.01
    data = _synthetic_tr_oxy(true_kla)
    run_dir = tmp_path / "run_test"
    run_dir.mkdir()
    np.savetxt(run_dir / "tr_oxy.dat", data, header="time C_w_oxy")

    postprocess_main(str(run_dir))

    results = json.loads((run_dir / "results.json").read_text())
    assert math.isfinite(results["kLa_25"])
    assert abs(results["kLa_25"] - true_kla) / true_kla < 0.05


def test_kla_reported_at_saturation_levels(tmp_path):
    """results.json contains keys kLa_10, kLa_25, kLa_50."""
    data = _synthetic_tr_oxy(0.01)
    run_dir = tmp_path / "run_test"
    run_dir.mkdir()
    np.savetxt(run_dir / "tr_oxy.dat", data)

    postprocess_main(str(run_dir))

    results = json.loads((run_dir / "results.json").read_text())
    for key in ("kLa_10", "kLa_25", "kLa_50"):
        assert key in results


def test_kla_returns_nan_on_insufficient_data(tmp_path):
    """tr_oxy.dat with fewer than 5 points → kLa fields are NaN, not a crash."""
    data = _synthetic_tr_oxy(0.01, n_points=3)
    run_dir = tmp_path / "run_test"
    run_dir.mkdir()
    np.savetxt(run_dir / "tr_oxy.dat", data)

    postprocess_main(str(run_dir))

    results = json.loads((run_dir / "results.json").read_text())
    assert math.isnan(results["kLa_25"])
