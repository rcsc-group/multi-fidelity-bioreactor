"""Guards on the chi -> dtmix pipeline that produces the Kim Fig. 9 replica.

Why this test exists
--------------------
2026-09-15. The mixing-time pipeline (`scripts.postprocess._compute_mixing_metrics`)
has existed for weeks but had never produced a number anyone compared against
Kim -- the Fig 13 sweeps all ran with `t_end` shorter than `t_mix`, so the
tracer was never actually released in any run on disk (runs/a34fc4d4's
tr_oxy.dat is all zeros). Before spending the ~500 wall-clock hours the L7
Fig 9 sweep needs, the pipeline gets pinned against cases whose answer is
known analytically.

The injection invariant
-----------------------
At injection the tracer is 1 over the top half of the liquid and 0 over the
bottom half, so over the liquid

    <c> = 0.5,  <c^2> = 0.5,  sigma^2_max = <c^2> - <c>^2 = 0.25

*exactly*, independent of resolution, rpm and bag geometry. That makes 0.25
a hard reference the code's own `sigma2_max` can be checked against.

The failure mode this protects against
--------------------------------------
`sigma2_max` is taken from the FIRST tr_oxy.dat row with non-zero tracer,
rather than from the analytic 0.25. If that row is not the injection instant
-- coarse output cadence, or an injection that did not fill exactly the top
half of the liquid -- the normalisation reference is wrong.

[CORRECTED 2026-09-15] My first version of this test asserted such a run
yields a too-SHORT dtmix. That is wrong, and the test caught me: it reported
36.44 s against a true 32.09 s, i.e. too LONG. Working it through, for a pure
exponential chi = 1 - exp(-t/tau) the two errors cancel *exactly* -- a first
sample late by delta starts the clock at t_inj+delta but also lowers the
threshold by exp(-delta/tau), and dtmix_code = tau*ln(20) independent of
delta. So dtmix is not merely an unreliable detector of a bad sigma2_max,
it is provably blind to it in the exponential limit.

That is the real reason to assert the invariant DIRECTLY: sigma2_max must
equal 0.25, and no amount of looking at dtmix will tell you when it doesn't.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.postprocess import _compute_mixing_metrics  # noqa: E402

# 32.5 rpm, 7 deg -- the condition with the most reference data.
PARAMS = {
    "omega_b": 2 * math.pi * 32.5 / 60.0,
    "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
    "theta_max": [7.0, 0.0, 0.0],
}
V_LIQ = 0.28          # liquid volume; any positive constant works
TAU_ND = 4.0          # non-dimensional mixing time constant


def _write_run(tmp_path: Path, t: np.ndarray, chi_true: np.ndarray) -> Path:
    """Synthesise a run whose tracer follows a prescribed chi(t).

    Inverts the pipeline's own definitions: tracer mass is conserved so
    <c> = 0.5 throughout, and sigma^2 = 0.25*(1 - chi) gives
    <c^2> = sigma^2 + 0.25. tr_oxy.dat stores volume-weighted sums, i.e.
    <c>*V_liq and <c^2>*V_liq.
    """
    c_mean = np.full_like(t, 0.5)
    c_mean[chi_true < 0] = 0.0                    # pre-injection rows
    sigma2 = 0.25 * (1.0 - np.clip(chi_true, 0.0, 1.0))
    c_mean2 = sigma2 + c_mean**2
    c_mean2[chi_true < 0] = 0.0

    run = tmp_path
    rows = []
    for k, (tk, cm, cm2) in enumerate(zip(t, c_mean, c_mean2)):
        # cols: i t oxy_sum oxy_sum2 c_sum c_sum2 c1_sum c1_sum2 c2_sum c2_sum2 c3_sum c3_sum2
        r = [k, tk, 0, 0, 0, 0, 0, 0, cm * V_LIQ, cm2 * V_LIQ, 0, 0]
        rows.append(" ".join(f"{v:.17g}" for v in r))
    (run / "tr_oxy.dat").write_text(
        "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 c1_liq_sum "
        "c1_liq_sum2 c2_liq_sum c2_liq_sum2 c3_liq_sum c3_liq_sum2\n"
        + "\n".join(rows) + "\n")
    (run / "vol_frac_interf.dat").write_text(
        "i t f_liq_sum f_liq_interf posY_max posY_min\n"
        + "\n".join(f"{k} {tk:.17g} {V_LIQ:.17g} 0 0 0" for k, tk in enumerate(t))
        + "\n")
    return run


def _t_bio() -> float:
    from scripts.postprocess import _t_scales
    return _t_scales(PARAMS)[0]


def test_exponential_mixing_recovers_analytic_dtmix(tmp_path):
    """chi(t) = 1 - exp(-t/tau) crosses each threshold at tau*ln(1/(1-chi))."""
    t_inject = 10.0
    t = np.concatenate([np.linspace(0, t_inject, 50, endpoint=False),
                        t_inject + np.linspace(0, 60.0, 4000)])
    chi_true = np.where(t < t_inject, -1.0,
                        1.0 - np.exp(-(t - t_inject) / TAU_ND))
    run = _write_run(tmp_path, t, chi_true)

    got = _compute_mixing_metrics(run, PARAMS)
    T_bio = _t_bio()
    for thr in (0.50, 0.75, 0.95):
        expected = TAU_ND * math.log(1.0 / (1.0 - thr)) * T_bio
        assert got[f"dtmix_{thr:.2f}"] == pytest.approx(expected, rel=2e-3), (
            f"dtmix_{thr:.2f}: got {got[f'dtmix_{thr:.2f}']:.4f} s, "
            f"analytic {expected:.4f} s")


def test_dtmix_thresholds_are_ordered(tmp_path):
    """Reaching 95% mixing can never be faster than reaching 50%."""
    t_inject = 5.0
    t = np.concatenate([np.linspace(0, t_inject, 20, endpoint=False),
                        t_inject + np.linspace(0, 60.0, 3000)])
    chi_true = np.where(t < t_inject, -1.0,
                        1.0 - np.exp(-(t - t_inject) / TAU_ND))
    got = _compute_mixing_metrics(_write_run(tmp_path, t, chi_true), PARAMS)
    assert got["dtmix_0.50"] <= got["dtmix_0.75"] <= got["dtmix_0.95"]


def test_unreached_threshold_is_nan_not_the_last_sample(tmp_path):
    """A run that stops before 95% mixing must report NaN, not its end time."""
    t_inject = 5.0
    t = np.concatenate([np.linspace(0, t_inject, 20, endpoint=False),
                        t_inject + np.linspace(0, 3.0, 400)])   # only ~chi=0.53
    chi_true = np.where(t < t_inject, -1.0,
                        1.0 - np.exp(-(t - t_inject) / TAU_ND))
    got = _compute_mixing_metrics(_write_run(tmp_path, t, chi_true), PARAMS)
    assert not math.isnan(got["dtmix_0.50"])
    assert math.isnan(got["dtmix_0.95"]), (
        "a threshold never reached within the run must be NaN -- otherwise a "
        "too-short run silently reports a finite, flatteringly small dtmix")


def test_sigma2_max_is_reported(tmp_path):
    """A clean injection must report sigma^2_max = 0.25, the analytic value.

    Reported so it lands in results.json and can be checked per run, instead
    of being an invisible internal normalisation constant.
    """
    t_inject = 10.0
    t = np.concatenate([np.linspace(0, t_inject, 50, endpoint=False),
                        t_inject + np.linspace(0, 60.0, 4000)])
    chi_true = np.where(t < t_inject, -1.0,
                        1.0 - np.exp(-(t - t_inject) / TAU_ND))
    got = _compute_mixing_metrics(_write_run(tmp_path, t, chi_true), PARAMS)
    assert got["sigma2_max"] == pytest.approx(0.25, rel=1e-6), (
        "at injection the tracer is 1 over the top half of the liquid and 0 "
        "over the bottom, so sigma^2_max = 0.25 exactly")


def test_bad_sigma2_max_invalidates_dtmix(tmp_path):
    """sigma^2_max far from 0.25 means the normalisation reference is wrong,
    so dtmix must be refused -- dtmix itself cannot reveal this (see module
    docstring: the errors cancel exactly for exponential mixing)."""
    t = np.linspace(0.0, 60.0, 4000)
    chi_true = 1.0 - 0.70 * np.exp(-t / TAU_ND)     # first sample already chi=0.30
    got = _compute_mixing_metrics(_write_run(tmp_path, t, chi_true), PARAMS)

    assert got["sigma2_max"] == pytest.approx(0.175, rel=1e-6), (
        "the measured sigma^2_max must be surfaced even when it is wrong")
    for thr in (0.50, 0.75, 0.95):
        assert math.isnan(got[f"dtmix_{thr:.2f}"]), (
            f"sigma^2_max = {got['sigma2_max']:.4f} is 30% below the analytic "
            f"0.25, so chi is normalised against the wrong reference and "
            f"dtmix_{thr:.2f} is not a mixing time. Return NaN.")
