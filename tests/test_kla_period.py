"""The whole-period kLa must not depend on where in the rocking cycle C* crosses.

Synthetic oxygen curve with a known mean transfer rate k and an in-cycle
modulation at 2 f_b, as measured in the L9/L10 runs (diag_kla_phase_lock.py,
2026-10-01):

    d ln(1 - C*)/dt = -k (1 + a cos(4 pi t / T))
 => ln(1 - C*) = -k t - (k a T / 4 pi) sin(4 pi t / T)

A fit over one whole period recovers k to a few percent at every phase; the
5-sample fit does not. That contrast is the bug being fixed.
"""
import math

import numpy as np
import pytest

from scripts.postprocess import _kla_5pt_at_threshold, _kla_period_at_threshold

T, DT, K, A = 0.6073, 0.02, 0.05, 0.3


def curve(shift):
    t = np.arange(0.0, 40.0, DT)
    lnr = -K * t - (K * A * T / (4 * math.pi)) * np.sin(4 * math.pi * (t + shift) / T)
    return t, 1.0 - np.exp(lnr)


@pytest.mark.parametrize("shift", np.linspace(0, T, 9, endpoint=False))
def test_one_period_fit_recovers_the_mean_rate_at_every_phase(shift):
    t, c = curve(shift)
    assert _kla_period_at_threshold(t, c, 0.25, T) == pytest.approx(K, rel=0.03)


def test_five_sample_fit_is_phase_dependent_on_the_same_curves():
    # negative control: the old estimator misses k by > 15% at some phase
    errs = [abs(_kla_5pt_at_threshold(*curve(s), 0.25) / K - 1)
            for s in np.linspace(0, T, 9, endpoint=False)]
    assert max(errs) > 0.15


def test_window_running_past_the_data_returns_nan():
    t, c = curve(0.0)
    keep = c < 0.26                     # series ends right after the crossing
    assert math.isnan(_kla_period_at_threshold(t[keep], c[keep], 0.25, T))


def test_window_reaching_back_before_release_returns_nan():
    t, c = curve(0.0)
    c = np.where(t < 5.0, 0.0, c)       # nothing before release at t = 5
    k = int(np.argmax(c >= 0.10))
    assert t[k] - 5.0 < T / 2           # the 10% crossing is within half a period of release
    assert math.isnan(_kla_period_at_threshold(t, c, 0.10, T))


def test_threshold_never_reached_returns_nan():
    t, c = curve(0.0)
    assert math.isnan(_kla_period_at_threshold(t, c * 0.2, 0.25, T))
