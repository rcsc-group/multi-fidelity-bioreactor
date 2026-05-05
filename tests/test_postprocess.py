"""Tests for postprocess.py — kLa extraction from BioReactor output files.

Real file formats (space-separated, header line starts with 'i'):
  tr_oxy.dat:          i t oxy_liq_sum oxy_liq_sum2 c_liq_sum ... (12 cols)
  vol_frac_interf.dat: i t f_liq_sum f_liq_interf posY_max posY_min (6 cols)

C* (dimensionless dissolved O2) = oxy_liq_sum / f_liq_sum_mean
First-order kinetic model: dC*/dt = kLa (1 - C*)
"""
import json
import math
import tempfile
from pathlib import Path

import numpy as np
import pytest

from scripts.postprocess import main as postprocess_main

F_LIQ_MEAN = 0.3571   # typical f_liq_sum for fill_level=0.5, geometry as in CANONICAL_PARAMS


def _write_tr_oxy(run_dir: Path, t: np.ndarray, oxy_liq_sum: np.ndarray) -> None:
    """Write a realistic 12-column tr_oxy.dat (unused columns set to 0)."""
    n = len(t)
    zeros = np.zeros(n)
    # columns: i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 ... (8 more zeros)
    data = np.column_stack([
        np.arange(n), t, oxy_liq_sum, zeros, zeros, zeros,
        zeros, zeros, zeros, zeros, zeros, zeros,
    ])
    header = "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 c1_liq_sum c1_liq_sum2 c2_liq_sum c2_liq_sum2 c3_liq_sum c3_liq_sum2"
    np.savetxt(run_dir / "tr_oxy.dat", data, header=header, comments="")


def _write_vol_frac(run_dir: Path, t: np.ndarray, f_liq_sum: float = F_LIQ_MEAN) -> None:
    """Write a 6-column vol_frac_interf.dat with constant f_liq_sum."""
    n = len(t)
    data = np.column_stack([
        np.arange(n), t,
        np.full(n, f_liq_sum),  # f_liq_sum
        np.ones(n) * 0.1,       # f_liq_interf (unused by postprocess)
        np.ones(n) * 0.05,      # posY_max
        np.ones(n) * -0.05,     # posY_min
    ])
    header = "i t f_liq_sum f_liq_interf posY_max posY_min"
    np.savetxt(run_dir / "vol_frac_interf.dat", data, header=header, comments="")


def _synthetic_run(run_dir: Path, kla: float, n_points: int = 50,
                   f_liq: float = F_LIQ_MEAN) -> None:
    """Write synthetic run files with exact first-order kinetics: C*(t) = 1 - exp(-kla*t).

    Time span covers up to 90% saturation so all three thresholds (10/25/50%) are reachable.
    """
    t_max = -math.log(0.10) / kla     # time to reach 90% saturation
    t = np.linspace(0.1, t_max, n_points)
    # C*(t) = oxy_liq_sum / f_liq  → oxy_liq_sum = f_liq * C*(t)
    c_star = 1.0 - np.exp(-kla * t)
    oxy_liq_sum = f_liq * c_star
    _write_tr_oxy(run_dir, t, oxy_liq_sum)
    _write_vol_frac(run_dir, t, f_liq)


def test_kla_extracted_from_synthetic_data(tmp_path):
    """Synthetic run with known kLa=0.01 — extraction within 5%."""
    true_kla = 0.01
    run_dir = tmp_path / "run_test"
    run_dir.mkdir()
    _synthetic_run(run_dir, true_kla)

    postprocess_main(str(run_dir))

    results = json.loads((run_dir / "results.json").read_text())
    assert math.isfinite(results["kLa_25"]), "kLa_25 is not finite"
    assert results["kLa_25"] > 0, "kLa_25 is not positive"
    assert abs(results["kLa_25"] - true_kla) / true_kla < 0.05, (
        f"kLa_25={results['kLa_25']:.4f} too far from true {true_kla}"
    )


def test_kla_reported_at_saturation_levels(tmp_path):
    """results.json must contain kLa_10, kLa_25, kLa_50."""
    run_dir = tmp_path / "run_test"
    run_dir.mkdir()
    _synthetic_run(run_dir, 0.01)

    postprocess_main(str(run_dir))

    results = json.loads((run_dir / "results.json").read_text())
    for key in ("kLa_10", "kLa_25", "kLa_50"):
        assert key in results, f"Missing key {key}"
        assert math.isfinite(results[key]) or math.isnan(results[key]), \
            f"{key} is neither finite nor NaN"


def test_kla_returns_nan_on_insufficient_data(tmp_path):
    """Fewer than 5 rows → NaN output, no crash."""
    run_dir = tmp_path / "run_test"
    run_dir.mkdir()
    _synthetic_run(run_dir, 0.01, n_points=3)

    postprocess_main(str(run_dir))

    results = json.loads((run_dir / "results.json").read_text())
    assert math.isnan(results["kLa_25"]), "Expected NaN for insufficient data"
