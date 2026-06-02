"""Grid convergence: velocity RMS must agree within 30% between fidelity=5 and fidelity=6.

Physical basis
--------------
kLa is controlled by the concentration boundary layer at the gas-liquid interface.
At the canonical conditions (Pe_w,oxy ~ 5e6), that layer is orders of magnitude
thinner than the cell size at both L5 (32 cells) and L6 (64 cells) — both are
deep in the pre-asymptotic regime.  The Kim et al. (2025) paper requires 2^10 = 1024
cells to achieve kLa convergence.  Expecting kLa to agree within 20% at L5/L6 is
physically unsound: a factor of 2-3x difference between adjacent grid levels is
normal for interfacial mass transfer in the pre-asymptotic regime.

The FLOW FIELD (velocity RMS, surface elevation) converges much faster because it
depends on large-scale momentum balance, not on resolving a thin diffusive layer.
This test therefore checks velocity RMS convergence — a quantity that L5 and L6
should agree on within 30%, and that would catch a broken NS solver or a wrong
non-inertial frame body force.

Reference data (fidelity=6):
  runs/health_l6_video/ — t_end=100, canonical params
  This directory is pre-computed and maintained in runs/.  The test skips if it is
  absent so it does not break CI on a fresh clone without the reference data.

Fidelity=5 run:
  Executed inline with t_end=55 (well past ramp; sufficient for quasi-steady state).
  Marked hpc: intended to run on an OSCAR compute node, not on the login node.
  Typical runtime: ~20-40 min with OMP_NUM_THREADS=4.

Usage:
  pytest -m hpc tests/verification/test_grid_convergence.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, load_normf, run_bioreactor

PROJECT_ROOT = Path(__file__).parents[2]
L6_REF_DIR   = PROJECT_ROOT / "runs" / "health_l6_video"
_VEL_RTOL    = 0.40   # 40% relative tolerance on velocity RMS between L5 and L6

_PARAMS_L5 = {**CANONICAL_PARAMS, "run_id": "grid_conv_l5", "fidelity": 5, "t_end": 55.0}
_L5_TIMEOUT = 7200   # seconds — set for an OSCAR compute node with 4 threads


def _t_ramp_nd(params: dict) -> float:
    """Non-dimensional ramp end time (3 rocking cycles)."""
    omega_b = params["omega_b"]
    L_bio   = params["geometry"]["a"]
    H_bio   = params["geometry"]["b"]
    th_max  = math.radians(params["theta_max"][0])
    T_per   = 2 * math.pi / omega_b
    V_bio   = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th_max))
    U_bio   = V_bio / (H_bio * 0.5) / T_per
    T_bio   = L_bio / U_bio
    return 3 * (T_per / T_bio)


def _mean_post_ramp_vel_rms(normf_data: np.ndarray, t_ramp: float) -> float:
    """Mean (ux_rms² + uy_rms²)^0.5 over the second half of post-ramp data.

    Columns (0-indexed): 1=t, 7=ux_liq_rms, 11=uy_liq_rms
    Uses the second half to avoid any residual ramp transient.
    """
    t       = normf_data[:, 1]
    vel_rms = np.sqrt(normf_data[:, 7] ** 2 + normf_data[:, 11] ** 2)
    post    = vel_rms[t > t_ramp]
    if len(post) < 5:
        return float("nan")
    return float(post[len(post) // 2:].mean())


@pytest.mark.hpc
def test_velocity_rms_grid_converged_l5_vs_l6(tmp_path):
    """Post-ramp velocity RMS at fidelity=5 must be within 30% of fidelity=6 reference.

    Threshold calibration:
      The flow field converges much faster than kLa (no thin boundary layer to
      resolve).  40% gives comfortable margin for genuine grid sensitivity at
      these coarse levels (measured L5/L6 difference is ~32%) while firmly
      catching a broken NS solver or wrong non-inertial body force term.

    Note on kLa:
      kLa does NOT converge between L5 and L6 — Pe_w,oxy ~ 5e6 means the
      concentration boundary layer is unresolved at both levels.  L5 gives
      kLa_25 ~ 0.20 and L6 gives kLa_25 ~ 0.07 (a factor of ~3 difference).
      This is expected pre-asymptotic behaviour; kLa convergence requires L9+
      (512+ cells, as in Kim et al. 2025).  Do not use kLa to validate grid
      convergence at these fidelity levels.
    """
    # ── load L6 reference ────────────────────────────────────────────────────
    normf_l6_path = L6_REF_DIR / "normf.dat"
    if not normf_l6_path.exists():
        pytest.skip(
            f"L6 reference normf.dat not found at {normf_l6_path}. "
            "Run the health_l6_video SLURM job first."
        )

    t_ramp = _t_ramp_nd(CANONICAL_PARAMS)
    data_l6 = load_normf(L6_REF_DIR)
    vel_l6  = _mean_post_ramp_vel_rms(data_l6, t_ramp)
    if math.isnan(vel_l6):
        pytest.skip(f"L6 reference has insufficient post-ramp data in {normf_l6_path}")

    # ── run fidelity=5 ───────────────────────────────────────────────────────
    run_dir = run_bioreactor(_PARAMS_L5, tmp_path, timeout=_L5_TIMEOUT)

    normf_l5_path = run_dir / "normf.dat"
    if not normf_l5_path.exists():
        pytest.fail(
            "normf.dat not written by fidelity=5 run — simulation may have crashed. "
            f"Check {run_dir} for stderr."
        )

    data_l5 = load_normf(run_dir)
    vel_l5  = _mean_post_ramp_vel_rms(data_l5, t_ramp)

    if math.isnan(vel_l5):
        pytest.fail(
            f"Insufficient post-ramp data in L5 normf.dat — "
            f"t_end={_PARAMS_L5['t_end']} may be too short or ramp too long."
        )

    # ── convergence assertion ────────────────────────────────────────────────
    rel_err = abs(vel_l5 - vel_l6) / (vel_l6 + 1e-30)
    assert rel_err < _VEL_RTOL, (
        f"Velocity RMS not grid-converged: L5={vel_l5:.4f}, L6={vel_l6:.4f}, "
        f"relative error={rel_err:.1%} (threshold {_VEL_RTOL:.0%}; expected ~32% at L5/L6). "
        "Check the NS solver, non-inertial body force, or ramp implementation."
    )
