"""Tests for the τ₉₈ (98th-percentile shear stress) KPI pipeline.

shear_stress.dat columns: i t tau_98 tau_mean tau_max
All tau values are non-dimensional: [μ_nd × U_bio / L] = [μ₁ × 1/T_bio].
postprocess._compute_tau98_kpis() converts to Pa using τ_Pa = τ_nd × μ_w / T_bio.

These tests cover:
 - Scalar extraction from synthetic shear_stress.dat
 - QSS window correctly excludes the ramp period (t < 3 × T_per_nd)
 - NaN return and no crash when shear_stress.dat is absent
 - Dimensional conversion formula: τ_Pa = τ_nd × μ_w / T_bio
"""
import json
import math
from pathlib import Path

import numpy as np
import pytest

from scripts.postprocess import main as postprocess_main, _t_scales, _compute_tau98_kpis
from tests.conftest import CANONICAL_PARAMS

# re-use synthetic kLa infrastructure from test_postprocess.py
F_LIQ_MEAN = 0.3571

MU_W = 1.0e-3   # Pa·s, water viscosity — matches BioReactor.c mu_w constant


# ── file writers ──────────────────────────────────────────────────────────────

def _write_shear_stress(
    run_dir: Path,
    t: np.ndarray,
    tau_98: np.ndarray,
    tau_mean: np.ndarray | None = None,
    tau_max: np.ndarray | None = None,
) -> None:
    """Write a 5-column shear_stress.dat.

    Columns: i t tau_98 tau_mean tau_max
    If tau_mean/tau_max are not provided, defaults to 0.8×tau_98 / 1.2×tau_98.
    """
    n = len(t)
    if tau_mean is None:
        tau_mean = 0.8 * tau_98
    if tau_max is None:
        tau_max = 1.2 * tau_98
    data = np.column_stack([np.arange(n), t, tau_98, tau_mean, tau_max])
    header = "i t tau_98 tau_mean tau_max"
    np.savetxt(run_dir / "shear_stress.dat", data, header=header, comments="")


def _write_minimal_kla_files(run_dir: Path, t: np.ndarray) -> None:
    """Write minimal tr_oxy.dat and vol_frac_interf.dat with zero oxygen (no kLa)."""
    n = len(t)
    zeros = np.zeros(n)
    tr_data = np.column_stack([np.arange(n), t] + [zeros] * 10)
    tr_header = "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 c1 c1s c2 c2s c3 c3s"
    np.savetxt(run_dir / "tr_oxy.dat", tr_data, header=tr_header, comments="")

    vf_data = np.column_stack([
        np.arange(n), t,
        np.full(n, F_LIQ_MEAN),
        np.full(n, 0.1), np.full(n, 0.05), np.full(n, -0.05),
    ])
    vf_header = "i t f_liq_sum f_liq_interf posY_max posY_min"
    np.savetxt(run_dir / "vol_frac_interf.dat", vf_data, header=vf_header, comments="")


# ── tests ─────────────────────────────────────────────────────────────────────

def test_tau98_extraction_known_value(tmp_path):
    """Constant tau_98_nd in QSS window → tau_98_qss matches expected Pa value ±1%.

    tau_98_qss = median(tau_98_nd_in_QSS_window) × tau_scale
    where tau_scale = μ_w / T_bio.
    """
    run_dir = tmp_path / "run_tau"
    run_dir.mkdir()

    T_bio, T_per_nd = _t_scales(CANONICAL_PARAMS)
    t_ramp = 3.0 * T_per_nd
    # QSS window: t_ramp to t_end (no injection in this test)
    t = np.linspace(0.1, t_ramp + 5.0, 80)

    tau_nd = 0.10  # constant non-dimensional shear stress in QSS window
    tau_98 = np.where(t > t_ramp, tau_nd, tau_nd * 10)  # spike during ramp (should be excluded)
    _write_shear_stress(run_dir, t, tau_98)
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)

    tau_scale = MU_W / T_bio   # Pa
    expected_pa = tau_nd * tau_scale
    assert math.isfinite(result["tau_98_qss"]), "tau_98_qss is not finite"
    assert abs(result["tau_98_qss"] - expected_pa) / expected_pa < 0.01, (
        f"tau_98_qss={result['tau_98_qss']:.4e} Pa, expected {expected_pa:.4e} Pa"
    )


def test_tau98_qss_excludes_ramp(tmp_path):
    """Spike during ramp (t < 3 × T_per_nd) must not inflate tau_98_qss.

    tau_98_max should capture the spike; tau_98_qss should reflect only the
    post-ramp quasi-steady value.
    """
    run_dir = tmp_path / "run_tau_ramp"
    run_dir.mkdir()

    T_bio, T_per_nd = _t_scales(CANONICAL_PARAMS)
    t_ramp = 3.0 * T_per_nd
    t = np.linspace(0.1, t_ramp + 5.0, 80)

    tau_normal = 0.10
    tau_spike  = tau_normal * 1000.0  # 1000× spike during ramp
    tau_98 = np.where(t < t_ramp, tau_spike, tau_normal)
    _write_shear_stress(run_dir, t, tau_98,
                        tau_max=np.where(t < t_ramp, tau_spike * 1.05, tau_normal * 1.05))
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)
    tau_scale = MU_W / T_bio

    assert math.isfinite(result["tau_98_qss"])
    # QSS value should be near tau_normal (within 2×), not inflated by the spike
    assert result["tau_98_qss"] < tau_normal * tau_scale * 2.0, (
        f"tau_98_qss={result['tau_98_qss']:.4e} Pa inflated by ramp spike"
    )
    # tau_98_max should reflect the spike
    assert result["tau_98_max"] > tau_normal * tau_scale * 100.0, (
        f"tau_98_max={result['tau_98_max']:.4e} Pa did not capture ramp spike"
    )


def test_tau98_nan_on_missing_file(tmp_path):
    """No shear_stress.dat → all tau_98_* keys are NaN; no crash.

    Old runs produced before this feature was added must postprocess gracefully.
    """
    run_dir = tmp_path / "run_no_tau"
    run_dir.mkdir()

    T_bio, T_per_nd = _t_scales(CANONICAL_PARAMS)
    t = np.linspace(0.1, 10.0, 50)
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))
    # intentionally do NOT write shear_stress.dat

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)
    assert math.isnan(result["tau_98_qss"]),      "Expected NaN for missing shear_stress.dat"
    assert math.isnan(result["tau_98_mean_qss"]), "Expected NaN for missing shear_stress.dat"
    assert math.isnan(result["tau_98_max"]),       "Expected NaN for missing shear_stress.dat"


def test_tau98_dimensional_conversion(tmp_path):
    """Known tau_nd → tau_Pa = tau_nd × μ_w / T_bio within 1%.

    tau_scale = μ_w × U_bio / L = μ_w / T_bio (since U_bio / L = 1/T_bio).
    This identity must hold exactly for the constraint to be in meaningful Pa units.
    """
    run_dir = tmp_path / "run_tau_dim"
    run_dir.mkdir()

    T_bio, T_per_nd = _t_scales(CANONICAL_PARAMS)
    t_ramp = 3.0 * T_per_nd
    t = np.linspace(t_ramp + 0.1, t_ramp + 5.0, 40)  # all post-ramp, no injection

    tau_nd = 0.25  # arbitrary known value
    _write_shear_stress(run_dir, t, np.full(len(t), tau_nd))
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)

    expected_pa = tau_nd * MU_W / T_bio
    assert abs(result["tau_98_qss"] - expected_pa) / expected_pa < 0.01, (
        f"Dimensional conversion error: got {result['tau_98_qss']:.4e} Pa, "
        f"expected {expected_pa:.4e} Pa  (tau_nd={tau_nd}, T_bio={T_bio:.3f} s)"
    )
