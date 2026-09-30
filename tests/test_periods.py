"""The end-of-run period boundary must not depend on float round-off.

The solver stops at a period boundary, and chain.py / sweep.py predict that
boundary to schedule the next segment. The old rule, floor(t/T) + 1, turns an
exact multiple k*T into k OR k+1 periods depending on the last bit of t/T.
On CI that made a 10-period segment run 11 while a 20-period run ran 20
(diary 2026-09-30). The rule is now: the smallest boundary >= t, with a
tolerance of 1e-5 period.
"""
import math

import pytest

from scripts.periods import next_period_boundary, PERIOD_TOL

T = 0.607329  # 7 deg, 37.5 rpm non-dimensional period (any value works)


@pytest.mark.parametrize("k", [1, 2, 10, 20, 80, 162])
@pytest.mark.parametrize("jitter", [-1e-9, 0.0, 1e-9, -5e-7, 5e-7])
def test_exact_multiple_is_k_periods_whatever_the_roundoff(k, jitter):
    n, t = next_period_boundary(k * T + jitter, T)
    assert n == k
    assert t == pytest.approx(k * T, abs=1e-12)


@pytest.mark.parametrize("frac", [0.01, 0.3, 0.5, 0.99])
def test_between_boundaries_rounds_up(frac):
    n, _ = next_period_boundary((10 + frac) * T, T)
    assert n == 11


def test_rounded_to_six_decimals_like_the_submit_scripts():
    for k in range(1, 200):
        n, _ = next_period_boundary(round(k * T, 6), T)
        assert n == k, k


def test_never_zero_periods():
    assert next_period_boundary(0.0, T)[0] == 1


def test_tolerance_is_far_below_any_requested_duration():
    assert PERIOD_TOL <= 1e-4 and PERIOD_TOL * T > 1e-6  # > round(...,6) error
    n, _ = next_period_boundary((10 + 10 * PERIOD_TOL) * T, T)
    assert n == 11
