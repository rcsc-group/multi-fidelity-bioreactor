"""Unit tests for health_report.py — three new numerical health KPIs.

All tests use synthetic in-memory data so no binary is required.

Columns referenced (0-indexed):
  normf  : 0=i, 1=t, 7=ux_rms, 9=ux_max, 11=uy_rms, 13=uy_max
  vol_frac: 0=i, 1=t, 2=f_liq_sum, 3=f_liq_interf
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.health_report import (
    _parse_logstats,
    kpi_cfl,
    kpi_kinetic_energy,
    kpi_pressure_residual,
)

# ── shared synthetic-data helpers ─────────────────────────────────────────────

_N     = 200
_T_MAX = (_N - 1) * 0.1   # 19.9
_T_RAMP = 2.0
_FIDELITY = 6
_DX = 1.0 / (2**_FIDELITY)   # 1/64 ≈ 0.015625


def _make_normf(n=_N, ux_rms=0.1, uy_rms=0.08, ux_max=0.5, uy_max=0.3) -> np.ndarray:
    """Quasi-steady normf array with constant velocities."""
    data = np.zeros((n, 14))
    for k in range(n):
        data[k, 0] = k * 24   # step index (approx 24 steps per t+=0.1)
        data[k, 1] = k * 0.1  # t
        data[k, 7] = ux_rms
        data[k, 9] = ux_max
        data[k, 11] = uy_rms
        data[k, 13] = uy_max
    return data


def _make_logstats(n=_N, dt=0.004) -> list[tuple[int, float, float]]:
    """Synthetic logstats matching the step indices produced by _make_normf."""
    return [(k * 24, k * 0.1, dt) for k in range(n)]


def _write_logstats(path: Path, rows: list[tuple[int, float, float]]) -> None:
    """Write synthetic logstats.dat in BioReactor format."""
    lines = [
        f"i: {i} t: {t} dt: {dt} #Cells: 4096 Wall clock time (s): 0 CPU time (s): 0\n"
        for i, t, dt in rows
    ]
    path.write_text("".join(lines))


# ── CFL KPI ───────────────────────────────────────────────────────────────────

class TestKpiCFL:
    def test_healthy_cfl_ok(self):
        # CFL = 0.004 * 0.5 / (1/64) = 0.128 — well within threshold 0.6
        logstats = _make_logstats(dt=0.004)
        normf = _make_normf(ux_max=0.5, uy_max=0.3)
        cfl_max, status = kpi_cfl(logstats, normf, fidelity=_FIDELITY)
        assert status == "OK"
        assert cfl_max == pytest.approx(0.128, rel=1e-3)

    def test_high_cfl_fails(self):
        # CFL = 0.02 * 0.5 / (1/64) = 0.64 — exceeds threshold 0.6
        logstats = _make_logstats(dt=0.02)
        normf = _make_normf(ux_max=0.5)
        cfl_max, status = kpi_cfl(logstats, normf, fidelity=_FIDELITY)
        assert status == "FAIL"
        assert cfl_max > 0.6

    def test_fidelity_7_tighter_dx(self):
        # Same dt=0.004, U=0.5, but fidelity=7 → dx=1/128 → CFL=0.256 — still OK
        logstats = _make_logstats(dt=0.004)
        normf = _make_normf(ux_max=0.5)
        cfl_max, status = kpi_cfl(logstats, normf, fidelity=7)
        assert status == "OK"
        assert cfl_max == pytest.approx(0.256, rel=1e-3)

    def test_parse_logstats_roundtrip(self, tmp_path):
        rows = _make_logstats(n=5, dt=0.0042)
        _write_logstats(tmp_path / "logstats.dat", rows)
        parsed = _parse_logstats(tmp_path)
        assert len(parsed) == 5
        assert parsed[0] == (0, pytest.approx(0.0), pytest.approx(0.0042))
        assert parsed[1][2] == pytest.approx(0.0042)


# ── KE quasi-steady KPI ───────────────────────────────────────────────────────

class TestKpiKineticEnergy:
    def test_quasi_steady_ok(self):
        normf = _make_normf()   # constant velocities → ratio ≈ 1.0
        ratio, status = kpi_kinetic_energy(normf, t_ramp=_T_RAMP)
        assert status == "OK"
        assert 0.5 < ratio < 2.0

    def test_decaying_ke_fails(self):
        normf = _make_normf()
        normf[110:, 7] = 0.001   # second half: near-zero ux_rms (flow stagnates)
        normf[110:, 11] = 0.001
        ratio, status = kpi_kinetic_energy(normf, t_ramp=_T_RAMP)
        assert status == "FAIL"
        assert ratio < 0.5

    def test_growing_ke_fails(self):
        normf = _make_normf()
        normf[110:, 7] = 5.0   # second half: 50× larger (blowup)
        normf[110:, 11] = 5.0
        ratio, status = kpi_kinetic_energy(normf, t_ramp=_T_RAMP)
        assert status == "FAIL"
        assert ratio > 2.0

    def test_insufficient_data_skips(self):
        normf = _make_normf(n=3)   # only 3 rows, all before t_ramp
        _, status = kpi_kinetic_energy(normf, t_ramp=_T_RAMP)
        assert status == "SKIP"


# ── Pressure residual KPI ─────────────────────────────────────────────────────

class TestKpiPressureResidual:
    def test_healthy_residuals_ok(self, tmp_path):
        (tmp_path / "pressure_diag.dat").write_text(
            "i t mgp_resa mgu_resa mgp_i mgu_i\n"
            "0 0.0 1e-7 1e-7 4 3\n"
            "24 0.1 2e-7 1.5e-7 5 4\n"
        )
        resa_max, status = kpi_pressure_residual(tmp_path)
        assert status == "OK"
        assert resa_max == pytest.approx(2e-7, rel=1e-3)

    def test_high_residuals_fail(self, tmp_path):
        (tmp_path / "pressure_diag.dat").write_text(
            "i t mgp_resa mgu_resa mgp_i mgu_i\n"
            "0 0.0 5e-4 3e-4 15 12\n"
        )
        resa_max, status = kpi_pressure_residual(tmp_path)
        assert status == "FAIL"
        assert resa_max > 1e-4

    def test_missing_file_skips(self, tmp_path):
        resa_max, status = kpi_pressure_residual(tmp_path)
        assert status == "SKIP"
        assert np.isnan(resa_max)

    def test_threshold_boundary(self, tmp_path):
        # Exactly at threshold (1e-4) should FAIL (tight check: < not <=)
        (tmp_path / "pressure_diag.dat").write_text(
            "i t mgp_resa mgu_resa mgp_i mgu_i\n"
            "0 0.0 1e-4 1e-5 8 6\n"
        )
        _, status = kpi_pressure_residual(tmp_path)
        assert status == "FAIL"
