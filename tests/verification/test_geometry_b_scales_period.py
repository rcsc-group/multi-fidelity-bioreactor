"""Physical verification: changing geometry.b shifts the rocking period
exactly as the closed-form non-dimensionalization predicts.

T_per_st = T_per / T_bio depends on H_bio, which depends on geometry.b.
This is a DIFFERENTIAL test (matching the style of
test_geometry_a_wired.py): rather than trusting one absolute formula, it
compares the MEASURED period ratio between two different geometry.b
values against the ratio predicted by the documented, independently-
reasoned formula (H_bio = 2*L_bio*Ly, full height = 2x the half-height
semi-axis geometry.b/L_bio).

A scale bug in H_bio (like the 2026-08-03 missing factor of 2) shifts
T_per_st itself, but NOT necessarily in a way a single-point check would
catch if the bug were, say, "off by a constant additive amount" rather
than multiplicative — checking the RATIO across two distinct geometry.b
values is a stronger, more general guardrail than checking one value in
isolation, and doesn't depend on getting T_per_st's absolute value
"right" for an unrelated reason.

Failure modes caught:
  - Any future error in how geometry.b enters H_bio/V_bio/U_bio/T_bio.
  - A sign or scale error that happens to look right for the specific
    (a=0.25, b=0.071) case this project has used throughout, but is
    wrong in general.
"""
import math
import sys
import pathlib

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, load_vol_frac, run_bioreactor


def _t_per_st_theory(params: dict) -> float:
    """Closed-form non-dim rocking period, full height = 2x geometry.b."""
    omega_b = params["omega_b"]
    L_bio   = params["geometry"]["a"]
    H_bio   = 2 * params["geometry"]["b"]  # full height (BioReactor.c:295)
    th_max  = math.radians(params["theta_max"][0])

    T_per  = 2 * math.pi / omega_b
    V_bio  = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th_max))
    U_bio  = V_bio / (H_bio * 0.5) / T_per
    T_bio  = L_bio / U_bio
    return T_per / T_bio


def _measured_period(run_dir, params) -> float:
    """Dominant FFT period of posY_max - posY_min (interface span), which
    oscillates at 2*omega_b (see test_forcing_frequency.py docstring) —
    so the rocking period itself is 2x the measured span period.
    """
    data = load_vol_frac(run_dir)
    t = data[:, 1]
    span = data[:, 4] - data[:, 5]
    t_per_st_guess = _t_per_st_theory(params)
    t_ramp = 3 * t_per_st_guess
    mask = t > t_ramp
    signal = span[mask] - span[mask].mean()
    dt = np.mean(np.diff(t[mask]))
    freqs = np.fft.rfftfreq(len(signal), d=dt)[1:]
    power = np.abs(np.fft.rfft(signal))[1:] ** 2
    dominant_freq = freqs[np.argmax(power)]
    return 2 / dominant_freq  # span oscillates at 2*omega_b -> period/2


@pytest.mark.medium
def test_doubling_geometry_b_shifts_period_as_theory_predicts(tmp_path):
    """Ratio of measured periods (b=0.10 / b=0.05) must match the ratio
    predicted by the closed-form T_per_st formula, within 15%.
    """
    params_small = {
        **CANONICAL_PARAMS, "run_id": "geo_b_small",
        "geometry": {"a": 0.25, "b": 0.05, "n": 8.0},
        "t_end": 15.0,
    }
    params_large = {
        **CANONICAL_PARAMS, "run_id": "geo_b_large",
        "geometry": {"a": 0.25, "b": 0.10, "n": 8.0},
        "t_end": 15.0,
    }

    run_small = run_bioreactor(params_small, tmp_path)
    run_large = run_bioreactor(params_large, tmp_path)

    measured_small = _measured_period(run_small, params_small)
    measured_large = _measured_period(run_large, params_large)
    measured_ratio = measured_large / measured_small

    theory_ratio = _t_per_st_theory(params_large) / _t_per_st_theory(params_small)

    rel_err = abs(measured_ratio - theory_ratio) / theory_ratio
    assert rel_err < 0.15, (
        f"Measured period ratio (b=0.10/b=0.05) = {measured_ratio:.3f}, "
        f"theory predicts {theory_ratio:.3f} — off by {rel_err:.1%}. "
        "A scale error in H_bio/geometry.b's contribution to T_per_st "
        "would show up here even if it looked fine at a single geometry."
    )
