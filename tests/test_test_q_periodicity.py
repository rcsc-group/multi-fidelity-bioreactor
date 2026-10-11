"""Test Q (periodicity per rpm): tests of the pure metric on synthetic fields."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.test_q_periodicity import (classify, floor_error, periodicity_error,  # noqa: E402
                                        bar_from_floor)

T = 1.0
FPP = 100
NP_ = 4


def _fields(drift=0.0, n=6):
    rng = np.random.default_rng(0)
    A = rng.normal(size=(2, n, n))
    B = rng.normal(size=(2, n, n)) + 3.0
    nt = NP_ * FPP + 1
    t = np.arange(nt) * T / FPP
    s = np.sin(2 * np.pi * t / T)[:, None, None, None]
    U = A[None] * s + B[None]
    U = U * (1.0 + drift * t / T)[:, None, None, None]
    mask = np.ones((nt, n, n), bool)
    return t, U, mask


def test_exactly_periodic_is_periodic():
    t, U, m = _fields()
    ph, e, e_all = periodicity_error(t, U, m, T)
    assert len(ph) == (NP_ - 1) * FPP + 1
    assert np.median(e) < 1e-2
    assert classify(np.median(e), 0.01)[0] == "PERIODIC"


def test_periodic_with_offset_period_uses_interpolation():
    t, U, m = _fields()
    # period estimate off the frame grid by a third of a frame: error stays small
    ph, e, _ = periodicity_error(t, U, m, T * (1 + 0.33 / FPP))
    assert np.median(e) < 2e-2


def test_drift_is_aperiodic():
    t, U, m = _fields(drift=0.20)
    _, e, _ = periodicity_error(t, U, m, T)
    assert np.median(e) > 0.10
    assert classify(np.median(e), 0.0)[0] == "APERIODIC"


def test_mask_restricts_cells():
    t, U, m = _fields()
    U2 = U.copy()
    U2[:, :, 0, 0] += np.linspace(0, 50, len(t))[:, None]   # junk in one cell
    m2 = m.copy()
    m2[:, 0, 0] = False
    _, e_masked, e_all = periodicity_error(t, U2, m2, T)
    assert np.median(e_masked) < 1e-2
    assert np.median(e_all) > np.median(e_masked)


def test_floor_order_of_magnitude():
    t, U, m = _fields()
    fl = floor_error(t, U, m)
    # pure sinusoid about offset B: error ~ (2 pi/100)^2/2 * |A sin|/|U| < 0.2%
    ref = (2 * np.pi / FPP) ** 2 / 2
    assert 0.02 * ref < np.median(fl) < 1.5 * ref


def test_floor_bar_rule():
    assert bar_from_floor(0.004) == (0.01, False)
    assert bar_from_floor(0.02) == (0.04, True)
    assert classify(0.03, 0.02)[0] == "PERIODIC"        # bar 4%
    assert classify(0.05, 0.02)[0] == "QUASI-PERIODIC"
    assert classify(0.05, 0.004)[0] == "QUASI-PERIODIC"
    assert classify(0.2, 0.004)[0] == "APERIODIC"
