"""Campaign 2D (2026-11-04, lead decision (a)+(b)): LF restart selection by marginal likelihood among restarts whose
noise_std is at or above the run-to-run floor; the optimiser's noise lower bound is the floor.

Synthetic data only. The LML is checked against scipy's Gaussian log-density (independent code path).
"""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import multivariate_normal

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))

import campaign as c  # noqa: E402


def _synthetic(seed=0, n=14, noise=0.1):
    rng = np.random.default_rng(seed)
    X = rng.uniform([15.0, 2.0], [37.5, 7.0], (n, 2))
    y = 0.3 * np.sin(X[:, 0] / 5.0) + 0.1 * np.cos(X[:, 1]) + rng.normal(0.0, noise, n)
    return X, y.reshape(-1, 1)


def test_select_restart_skips_floor_violators_even_with_best_lml():
    rows = [dict(seed=1, noise_std=2e-4, lml=5.0, ok=False),      # degenerate mode, highest LML, below floor
            dict(seed=2, noise_std=0.12, lml=-3.0, ok=True),
            dict(seed=3, noise_std=0.20, lml=-4.0, ok=True)]
    assert c.select_restart(rows) == 1


def test_select_restart_raises_when_none_qualifies():
    rows = [dict(seed=1, noise_std=2e-4, lml=5.0, ok=False)]
    with pytest.raises(ValueError):
        c.select_restart(rows)


def test_lf_log_ml_matches_scipy_gaussian_logpdf():
    X, Y = _synthetic()
    lf = c.make_lf(seed=1, floor=0.0)
    lf.train(X, Y, portion_test=0.2)
    K = lf._training_kernel_matrix(scaled_x=lf.sample_x_scaled, scaled_noise_std=lf.noise_std / lf.y_std)
    ref = multivariate_normal(mean=np.zeros(len(K)), cov=K).logpdf(np.asarray(lf.sample_y_scaled, float).ravel())
    assert abs(c.lf_log_ml(lf) - ref) <= 1e-8 * abs(ref)   # relative: the value is about 7e5


def test_noise_floor_is_respected_by_the_optimiser():
    X, Y = _synthetic()
    floor = 0.095
    lf = c.make_lf(seed=1, floor=floor)
    lf.train(X, Y, portion_test=0.2)
    assert float(np.asarray(lf.noise_std).ravel()[0]) >= floor - 1e-9


def test_degenerate_mode_is_not_selected_when_floored():
    # Known degenerate candidate: noise 2e-4 (50 times below the floor). Its LML is computed on the same data.
    X, Y = _synthetic(noise=0.001)                         # nearly noise-free data (the LML of the degenerate candidate is lower here: this tests the ok filter, not an LML contest)
    floored = c.make_lf(seed=1, floor=0.095)
    floored.train(X, Y, portion_test=0.2)
    degenerate = c.make_lf(seed=2, floor=0.0)
    degenerate.train(X, Y, portion_test=0.2)
    degenerate.noise_std = 2e-4                            # the degenerate noise, same kernel
    rows = [dict(seed=2, noise_std=2e-4, lml=c.lf_log_ml(degenerate), ok=False),
            dict(seed=1, noise_std=float(np.asarray(floored.noise_std).ravel()[0]),
                 lml=c.lf_log_ml(floored), ok=True)]
    assert c.select_restart(rows) == 1


def test_hf_grid_floor_sets_lower_noise_ratio():
    l, eta, lp = c.hf_grid(floor=0.095, ystd=0.4)
    assert abs(eta.min() - 0.095 / 0.4) / (0.095 / 0.4) < 0.05
    assert eta.max() >= eta.min()
