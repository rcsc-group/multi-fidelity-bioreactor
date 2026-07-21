"""Tests for the shear stress percentile KPI pipeline.

shear_stress.dat columns: i t tau_95 tau_98 tau_100 tau_mean
BioReactor.c sets rho1=1, mu1=1/Re_w, so tau_nd = tau_dim / (rho_w * U_bio^2)
-- not tau_dim * T_bio / mu_w. postprocess._compute_tau98_kpis() converts to
Pa using tau_Pa = tau_nd * rho_w * U_bio^2, with U_bio = geometry.a / T_bio.
"""
import json
import math
from pathlib import Path

import numpy as np

from scripts.postprocess import _t_scales, _compute_tau98_kpis
from tests.conftest import CANONICAL_PARAMS

F_LIQ_MEAN = 0.3571
RHO_W = 1.0e3   # kg/m^3, water density — matches postprocess._compute_tau98_kpis


def _tau_scale(params: dict) -> float:
    """tau_Pa = tau_nd * RHO_W * U_bio^2, U_bio = geometry.a / T_bio."""
    T_bio, _ = _t_scales(params)
    L = params.get("geometry", {}).get("a", 0.25)
    U_bio = L / T_bio
    return RHO_W * U_bio**2


# ── file writers ──────────────────────────────────────────────────────────────

def _write_shear_stress(
    run_dir: Path,
    t: np.ndarray,
    tau_98: np.ndarray,
    tau_95: np.ndarray | None = None,
    tau_100: np.ndarray | None = None,
    tau_mean: np.ndarray | None = None,
) -> None:
    """Write a 6-column shear_stress.dat: i t tau_95 tau_98 tau_100 tau_mean."""
    n = len(t)
    if tau_95  is None: tau_95  = 0.95 * tau_98
    if tau_100 is None: tau_100 = 1.05 * tau_98
    if tau_mean is None: tau_mean = 0.8 * tau_98
    data = np.column_stack([np.arange(n), t, tau_95, tau_98, tau_100, tau_mean])
    header = "i t tau_95 tau_98 tau_100 tau_mean"
    np.savetxt(run_dir / "shear_stress.dat", data, header=header, comments="")


def _write_minimal_kla_files(run_dir: Path, t: np.ndarray) -> None:
    """Write minimal tr_oxy.dat and vol_frac_interf.dat with zero oxygen."""
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
    """Constant tau_98_nd in QSS window → tau_98_qss matches expected Pa value ±1%."""
    run_dir = tmp_path / "run_tau"
    run_dir.mkdir()

    T_bio, T_per_nd = _t_scales(CANONICAL_PARAMS)
    t_ramp = 3.0 * T_per_nd
    t = np.linspace(0.1, t_ramp + 5.0, 80)

    tau_nd = 0.10
    tau_98 = np.where(t > t_ramp, tau_nd, tau_nd * 10)
    _write_shear_stress(run_dir, t, tau_98)
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)

    expected_pa = tau_nd * _tau_scale(CANONICAL_PARAMS)
    assert math.isfinite(result["tau_98_qss"]), "tau_98_qss is not finite"
    assert abs(result["tau_98_qss"] - expected_pa) / expected_pa < 0.01, (
        f"tau_98_qss={result['tau_98_qss']:.4e} Pa, expected {expected_pa:.4e} Pa"
    )


def test_tau98_qss_excludes_ramp(tmp_path):
    """Spike during ramp must not inflate tau_98_qss OR tau_98_max.

    tau_98_max is deliberately QSS-window-restricted (see
    _compute_tau98_kpis' docstring: it matches Kim et al.'s "max over one
    period in the quasi-steady regime", explicitly excluding the startup
    transient) -- a spike entirely inside the excluded ramp must not leak
    into either metric.
    """
    run_dir = tmp_path / "run_tau_ramp"
    run_dir.mkdir()

    T_bio, T_per_nd = _t_scales(CANONICAL_PARAMS)
    t_ramp = 3.0 * T_per_nd
    t = np.linspace(0.1, t_ramp + 5.0, 80)

    tau_normal = 0.10
    tau_spike  = tau_normal * 1000.0
    tau_98  = np.where(t < t_ramp, tau_spike, tau_normal)
    tau_100 = np.where(t < t_ramp, tau_spike * 1.05, tau_normal * 1.05)
    _write_shear_stress(run_dir, t, tau_98, tau_100=tau_100)
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)
    tau_scale = _tau_scale(CANONICAL_PARAMS)

    assert math.isfinite(result["tau_98_qss"])
    assert result["tau_98_qss"] < tau_normal * tau_scale * 2.0, (
        f"tau_98_qss={result['tau_98_qss']:.4e} Pa inflated by ramp spike"
    )
    assert result["tau_98_max"] < tau_normal * tau_scale * 2.0, (
        f"tau_98_max={result['tau_98_max']:.4e} Pa leaked the ramp spike "
        f"despite being QSS-restricted"
    )


def test_tau98_max_captures_qss_window_spike(tmp_path):
    """A transient spike INSIDE the QSS window is captured by tau_98_max,
    while tau_98_qss (median) stays robust to it -- confirms _qss_max is a
    real per-window max, not a metric that's collapsed to the same median."""
    run_dir = tmp_path / "run_tau_qss_spike"
    run_dir.mkdir()

    T_bio, T_per_nd = _t_scales(CANONICAL_PARAMS)
    t_ramp = 3.0 * T_per_nd
    t = np.linspace(0.1, t_ramp + 10.0, 160)

    tau_normal = 0.10
    tau_spike  = tau_normal * 50.0
    # Single spike well after the ramp, deep inside the QSS window.
    spike_idx = int(np.argmin(np.abs(t - (t_ramp + 5.0))))
    tau_98 = np.full_like(t, tau_normal)
    tau_98[spike_idx] = tau_spike
    _write_shear_stress(run_dir, t, tau_98)
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)
    tau_scale = _tau_scale(CANONICAL_PARAMS)

    assert result["tau_98_qss"] < tau_normal * tau_scale * 2.0, (
        f"tau_98_qss={result['tau_98_qss']:.4e} Pa inflated by a single spike"
    )
    assert result["tau_98_max"] > tau_normal * tau_scale * 10.0, (
        f"tau_98_max={result['tau_98_max']:.4e} Pa did not capture the "
        f"in-QSS-window spike"
    )


def test_tau98_nan_on_missing_file(tmp_path):
    """No shear_stress.dat → all tau_* keys are NaN; no crash."""
    run_dir = tmp_path / "run_no_tau"
    run_dir.mkdir()

    t = np.linspace(0.1, 10.0, 50)
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)
    for key in ("tau_95_qss", "tau_98_qss", "tau_100_qss",
                "tau_95_max", "tau_98_max", "tau_100_max"):
        assert math.isnan(result[key]), f"Expected NaN for {key} when file missing"


def test_tau98_dimensional_conversion(tmp_path):
    """Known tau_nd → tau_Pa = tau_nd * rho_w * U_bio^2 within 1%."""
    run_dir = tmp_path / "run_tau_dim"
    run_dir.mkdir()

    T_bio, T_per_nd = _t_scales(CANONICAL_PARAMS)
    t_ramp = 3.0 * T_per_nd
    t = np.linspace(t_ramp + 0.1, t_ramp + 5.0, 40)

    tau_nd = 0.25
    _write_shear_stress(run_dir, t, np.full(len(t), tau_nd))
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)

    expected_pa = tau_nd * _tau_scale(CANONICAL_PARAMS)
    assert abs(result["tau_98_qss"] - expected_pa) / expected_pa < 0.01, (
        f"Dimensional conversion error: got {result['tau_98_qss']:.4e} Pa, "
        f"expected {expected_pa:.4e} Pa"
    )


def test_tau95_lt_tau98_lt_tau100(tmp_path):
    """For a realistic distribution, tau_95_qss ≤ tau_98_qss ≤ tau_100_qss."""
    run_dir = tmp_path / "run_tau_order"
    run_dir.mkdir()

    T_bio, T_per_nd = _t_scales(CANONICAL_PARAMS)
    t_ramp = 3.0 * T_per_nd
    t = np.linspace(t_ramp + 0.1, t_ramp + 5.0, 40)

    n = len(t)
    tau_95  = np.full(n, 0.08)
    tau_98  = np.full(n, 0.10)
    tau_100 = np.full(n, 0.15)
    _write_shear_stress(run_dir, t, tau_98, tau_95=tau_95, tau_100=tau_100)
    _write_minimal_kla_files(run_dir, t)
    (run_dir / "params.json").write_text(json.dumps(CANONICAL_PARAMS))

    result = _compute_tau98_kpis(run_dir, CANONICAL_PARAMS)

    assert result["tau_95_qss"] <= result["tau_98_qss"], (
        f"tau_95_qss={result['tau_95_qss']:.4e} > tau_98_qss={result['tau_98_qss']:.4e}"
    )
    assert result["tau_98_qss"] <= result["tau_100_qss"], (
        f"tau_98_qss={result['tau_98_qss']:.4e} > tau_100_qss={result['tau_100_qss']:.4e}"
    )
