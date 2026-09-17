"""The cross-level L2 norm must compare solutions, not regridding artifacts.

Kim's Fig. A.16(c) differences each resolution against a finer reference. Two
choices decide whether the number means anything, and both are silent when
wrong: the direction of the regridding, and the mask.

Coarsening the reference DOWN is exact for uniform grids at a power-of-two
ratio -- the coarse cell's value is the mean of the fine cells inside it.
Interpolating the coarse field UP instead would measure how smooth the
interpolant is, and would report a small error for a coarse field that is
merely smooth.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.plot_figA16_c import block_mean, l2_error  # noqa: E402


def _field(n, ux, uy, liquid=True):
    one = np.ones((n, n))
    return (0.0, one * ux, one * uy, one * (1.0 if liquid else 0.0), one)


def test_block_mean_is_the_exact_coarse_cell_average():
    a = np.arange(16, dtype=float).reshape(4, 4)
    got = block_mean(a, 2)
    assert got.shape == (2, 2)
    assert got[0, 0] == pytest.approx((0 + 1 + 4 + 5) / 4)
    assert got[1, 1] == pytest.approx((10 + 11 + 14 + 15) / 4)


def test_identical_fields_have_zero_error():
    ref = _field(8, 2.0, 1.0)
    coarse = _field(4, 2.0, 1.0)
    ex, ey = l2_error(coarse, ref, factor=2)
    assert ex == pytest.approx(0.0)
    assert ey == pytest.approx(0.0)


def test_error_is_normalised_by_the_reference_mean():
    """e = rms(difference) / <u_ref>. A coarse field 10% high against a
    reference of 2.0 must read 0.10, not 0.20."""
    ref = _field(8, 2.0, 1.0)
    coarse = _field(4, 2.2, 1.0)
    ex, _ = l2_error(coarse, ref, factor=2)
    assert ex == pytest.approx(0.10, rel=1e-6)


def test_only_liquid_inside_the_bag_enters_the_norm():
    """A cell that is dry on either grid carries no fluid. Including it drags
    the norm toward whatever the solver left in the void."""
    n = 4
    one = np.ones((n, n))
    f = np.ones((n, n))
    f[0, :] = 0.0                                # one dry row on the coarse grid
    coarse = (0.0, one * 2.0, one * 1.0, f, one)
    big = np.ones((8, 8))
    ref = (0.0, big * 2.0, big * 1.0, big, big)
    ref[1][0:2, :] = 99.0                        # garbage under the dry row
    ex, _ = l2_error(coarse, ref, factor=2)
    assert ex == pytest.approx(0.0)
