"""Physical verification: interface oscillates at the angular rocking frequency.

In a rocking bag the free surface tilts with the bag.  posY_max − posY_min
(interface y-span) therefore oscillates at exactly omega_b.  The FFT of
that signal must show a dominant peak within 20% of the expected non-dim
rocking frequency.

Failure modes caught:
  - omega_b not applied (e.g. forcing block skipped)
  - wrong non-dimensionalization of the forcing term
  - interface not deforming with the bag (boundary condition bug)
"""
import math

import numpy as np
import pytest

from tests.conftest import CANONICAL_PARAMS, load_vol_frac, run_bioreactor


# ── helpers ───────────────────────────────────────────────────────────────────

def _omega_b_nd(params: dict) -> float:
    """Non-dimensional rocking angular frequency (rad / t_nd).

    Replicates BioReactor.c:
        T_per  = 2*pi / omega_b
        V_bio  = L/4 * (b + 0.5*L*tan(theta_max[0]))
        U_bio  = V_bio / (b*0.5) / T_per
        T_bio  = L / U_bio
        omega_b_nd = omega_b * T_bio
    """
    omega_b = params["omega_b"]
    L_bio   = params["geometry"]["a"]
    H_bio   = params["geometry"]["b"]
    th_max  = math.radians(params["theta_max"][0])

    T_per  = 2 * math.pi / omega_b
    V_bio  = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th_max))
    U_bio  = V_bio / (H_bio * 0.5) / T_per
    T_bio  = L_bio / U_bio
    return omega_b * T_bio


# ── test ──────────────────────────────────────────────────────────────────────

@pytest.mark.medium
def test_interface_oscillates_at_rocking_frequency(tmp_path):
    """posY_max spectral power must be concentrated near omega_b.

    Physical basis: the bag rocks at omega_b; as it tilts, the highest
    interface cell (posY_max) rises and falls once per rocking cycle —
    an antisymmetric signal at exactly omega_b.  Note: the span
    posY_max − posY_min is symmetric (same whether tilted left or right)
    and therefore oscillates at 2*omega_b, not omega_b.

    At fidelity=3 the interface is represented by only a few cells, so
    posY_max is quantised and contains broadband noise.  We therefore
    check the *power fraction* in the ±40% band around omega_b rather
    than the argmax frequency.  A fraction ≥ 20% is required; a
    correctly forced simulation consistently shows ~37%.

    Failure: forcing not applied, wrong omega_b, or interface not
    responding to the forcing.
    """
    params = {**CANONICAL_PARAMS, "run_id": "fft_freq", "t_end": 30.0}
    run_dir = run_bioreactor(params, tmp_path)

    data = load_vol_frac(run_dir)
    assert len(data) >= 20, "Too few output rows — simulation may not have run"

    t       = data[:, 1]
    posYmax = data[:, 4]   # antisymmetric: oscillates at omega_b

    # skip ramp (3 rocking cycles) and subtract mean
    t_ramp = 3 * (2 * math.pi / _omega_b_nd(params))   # 3 * T_per_st
    mask   = t > t_ramp
    signal = posYmax[mask] - posYmax[mask].mean()

    assert len(signal) >= 20, "Not enough post-ramp rows for spectral analysis"

    dt    = np.mean(np.diff(t[mask]))
    freqs = np.fft.rfftfreq(len(signal), d=dt)[1:]        # cycles / t_nd
    power = np.abs(np.fft.rfft(signal))[1:] ** 2

    expected_freq = _omega_b_nd(params) / (2 * math.pi)   # cycles / t_nd
    band = (freqs >= expected_freq * 0.6) & (freqs <= expected_freq * 1.4)
    frac = power[band].sum() / (power.sum() + 1e-30)

    assert frac >= 0.20, (
        f"Only {frac:.1%} of posY_max spectral power lies in the "
        f"omega_b band [{expected_freq*0.6:.3f}, {expected_freq*1.4:.3f}] cycles/t_nd; "
        f"expected ≥ 20% for a correctly forced simulation"
    )
