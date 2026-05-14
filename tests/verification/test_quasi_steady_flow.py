"""Physical verification: post-ramp velocity RMS reaches quasi-steady state.

After the soft-start ramp (N_RAMP_CYCLES = 3 rocking cycles), the driving is
fully applied and the flow should settle into a quasi-periodic driven regime.
The liquid velocity RMS (ux_rms² + uy_rms²)^0.5 must not grow monotonically
through the remainder of the run.

Failure modes caught:
  - Numerical instability (velocity blowing up after ramp)
  - Ramp not working (sudden forcing causes runaway kinetic energy)
  - Non-inertial pseudo-force term with wrong sign (continuous acceleration)
"""
import math

import numpy as np
import pytest

from tests.conftest import CANONICAL_PARAMS, load_normf, run_bioreactor


# ── helpers ───────────────────────────────────────────────────────────────────

def _t_ramp_nd(params: dict) -> float:
    """Non-dimensional ramp end time (3 rocking cycles).

    Replicates BioReactor.c:  t_change_st = N_RAMP_CYCLES * T_per_st
    """
    omega_b = params["omega_b"]
    L_bio   = params["geometry"]["a"]
    H_bio   = params["geometry"]["b"]
    th_max  = math.radians(params["theta_max"][0])

    T_per  = 2 * math.pi / omega_b
    V_bio  = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th_max))
    U_bio  = V_bio / (H_bio * 0.5) / T_per
    T_bio  = L_bio / U_bio
    T_per_st = T_per / T_bio
    return 3 * T_per_st   # N_RAMP_CYCLES = 3


# ── test ──────────────────────────────────────────────────────────────────────

@pytest.mark.medium
def test_velocity_rms_quasi_steady_after_ramp(tmp_path):
    """Post-ramp liquid velocity RMS must stay quasi-steady (ratio in [0.7, 1.5]).

    Physical basis: a linearly ramped sinusoidal forcing produces a driven
    quasi-periodic flow.  Once the ramp is complete the RMS kinetic energy
    oscillates around a stationary mean.  A ratio outside [0.7, 1.5] indicates
    numerical blowup, continuous acceleration, or a decaying/non-driven flow.

    Threshold calibrated from L6 reference run (health_l6_video): ratio = 1.12.
    [0.7, 1.5] gives ~35% margin around the measured value on both sides.

    Columns in normf.dat (0-indexed):
        1: t, 7: ux_liq_rms, 11: uy_liq_rms
    """
    params  = {**CANONICAL_PARAMS, "run_id": "ke_steady", "t_end": 30.0}
    run_dir = run_bioreactor(params, tmp_path)

    data = load_normf(run_dir)
    assert len(data) >= 20, "Too few output rows — simulation may not have run"

    t       = data[:, 1]
    vel_rms = np.sqrt(data[:, 7] ** 2 + data[:, 11] ** 2)

    t_ramp = _t_ramp_nd(params)
    mask   = t > t_ramp
    if mask.sum() < 10:
        pytest.skip("Not enough post-ramp data (increase t_end or reduce omega_b)")

    post_ramp = vel_rms[mask]
    mid       = len(post_ramp) // 2
    ratio     = post_ramp[mid:].mean() / (post_ramp[:mid].mean() + 1e-30)

    assert 0.7 <= ratio <= 1.5, (
        f"Velocity RMS ratio {ratio:.2f} outside [0.7, 1.5] — "
        f"expected quasi-steady driven flow post-ramp"
    )
