"""A sweep point that ran out of clock must extend itself, not sit as a hole.

`submit_fig9.py` sizes t_end as MARGIN x Kim's own dtmix_0.95, on the premise
that a coarser grid mixes faster. The premise holds early and fails late: at
37.5 rpm we cross chi=0.50 in 0.71x Kim's time, chi=0.75 in 1.01x, and never
reach 0.95 inside the window at all. Numerical diffusion buys the bulk
homogenisation and nothing in the tail.

So a point can finish cleanly, pass every invariant, and still report
dtmix_0.95 = NaN. The distinction that matters is WHY it is NaN:

  * sigma^2_max off 0.25  -> the normalisation is wrong; more compute cannot
                             fix it and must not be spent.
  * sigma^2_max sound     -> the physics is fine and the clock ran out; a
                             continuation segment finishes the measurement.

These tests pin that distinction and the extrapolation that sizes the
continuation, because both failures are silent: the first burns a day of
compute to produce the same NaN, the second stops one sample short of the
threshold and looks identical to a point that was never submitted.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(PROJECT_ROOT))
from scripts.autoextend import (  # noqa: E402
    chain_tip, extension_cycles, needs_extension,
)


# ── which points deserve more compute ────────────────────────────────────────

def test_sound_run_that_ran_out_of_clock_is_extended():
    assert needs_extension({"dtmix_0.95": math.nan, "sigma2_max": 0.2423})


def test_completed_point_is_not_extended():
    assert not needs_extension({"dtmix_0.95": 184.2, "sigma2_max": 0.2423})


def test_broken_normalisation_is_not_extended():
    """sigma^2_max away from 0.25 means chi is measured against the wrong
    reference. The tail will never be reached because the target is wrong, so
    a continuation would burn a full segment to return the same NaN."""
    assert not needs_extension({"dtmix_0.95": math.nan, "sigma2_max": 0.11})


def test_missing_key_is_not_extended():
    """No dtmix key at all means postprocess bailed before the mixing block
    (no tr_oxy.dat, too few rows). That is a broken run, not a short one."""
    assert not needs_extension({"sigma2_max": 0.2423})


# ── how much more compute ────────────────────────────────────────────────────

def _truncated_exponential(tau_s: float, t_stop_s: float, n: int = 400):
    """chi(t) = 1 - exp(-t/tau), sampled up to t_stop and no further."""
    t = np.linspace(0.0, t_stop_s, n)
    return t, 1.0 - np.exp(-t / tau_s)


def test_extrapolation_recovers_the_analytic_remainder():
    """For a single exponential the remaining time to 0.95 is known in closed
    form: tau * ln(u_end / 0.05). The estimator must land on it."""
    tau, t_stop = 12.0, 24.0                      # stops at chi = 0.865
    t, chi = _truncated_exponential(tau, t_stop)
    analytic = tau * math.log((1 - chi[-1]) / 0.05)
    got = extension_cycles(t, chi, T_per=1.846, margin=1.0)
    assert got * 1.846 == pytest.approx(analytic, rel=0.05)


def test_margin_is_applied_on_top():
    tau, t_stop = 12.0, 24.0
    t, chi = _truncated_exponential(tau, t_stop)
    bare = extension_cycles(t, chi, T_per=1.846, margin=1.0)
    padded = extension_cycles(t, chi, T_per=1.846, margin=1.5)
    assert padded == pytest.approx(1.5 * bare, rel=1e-6)


def test_target_already_reached_asks_for_nothing():
    t, chi = _truncated_exponential(12.0, 60.0)   # chi = 0.993
    assert extension_cycles(t, chi, T_per=1.846) == 0.0


def test_a_slowing_tail_is_extrapolated_from_the_tail_not_the_whole_curve():
    """Mixing is not a single exponential -- the measured tail relaxes slower
    than the bulk. Fitting the whole curve therefore UNDER-estimates the
    remainder and re-submits a segment that falls short again. The fit must
    use the late window only."""
    t = np.linspace(0.0, 32.0, 500)
    chi = np.where(t < 18.0,
                   1 - np.exp(-t / 8.0),
                   1 - np.exp(-18.0 / 8.0) * np.exp(-(t - 18.0) / 25.0))
    slow_tau = 25.0
    analytic = slow_tau * math.log((1 - chi[-1]) / 0.05)
    got = extension_cycles(t, chi, T_per=1.846, margin=1.0)
    assert got * 1.846 == pytest.approx(analytic, rel=0.10)


def test_a_stalled_tail_falls_back_instead_of_asking_for_infinity():
    """A chi that has flattened gives tau -> inf and an unbounded request.
    Cap it: an unbounded ask is a bug report, not a job."""
    t = np.linspace(0.0, 40.0, 400)
    chi = np.full_like(t, 0.60)
    got = extension_cycles(t, chi, T_per=1.846, max_cycles=300.0)
    assert got == 300.0


# ── where the next segment attaches ──────────────────────────────────────────

def test_chain_tip_is_the_base_run_when_nothing_was_extended(tmp_path):
    (tmp_path / "fig9_l9_rpm30").mkdir()
    assert chain_tip(tmp_path, "fig9_l9_rpm30") == "fig9_l9_rpm30"


def test_chain_tip_follows_the_deepest_extension(tmp_path):
    for name in ("fig9_l9_rpm30", "fig9_l9_rpm30_ext1", "fig9_l9_rpm30_ext2"):
        (tmp_path / name).mkdir()
    assert chain_tip(tmp_path, "fig9_l9_rpm30") == "fig9_l9_rpm30_ext2"


def test_chain_tip_is_not_confused_by_a_different_point(tmp_path):
    """fig9_l9_rpm3 is a prefix of fig9_l9_rpm30 as a plain string. Matching on
    the prefix alone would attach rpm30's continuation to rpm3's chain."""
    for name in ("fig9_l9_rpm3", "fig9_l9_rpm30", "fig9_l9_rpm30_ext1"):
        (tmp_path / name).mkdir()
    assert chain_tip(tmp_path, "fig9_l9_rpm3") == "fig9_l9_rpm3"


def test_results_are_read_from_the_chain_tip(tmp_path):
    """plot_fig9 reads the base run_id. Once a point has been extended, the
    finished measurement lives in the LAST segment -- the base run's own
    results.json still says NaN forever."""
    from scripts.autoextend import tip_results

    base = tmp_path / "fig9_l9_rpm30"
    base.mkdir()
    (base / "results.json").write_text(json.dumps({"dtmix_0.95": math.nan}))
    ext = tmp_path / "fig9_l9_rpm30_ext1"
    ext.mkdir()
    (ext / "results.json").write_text(json.dumps({"dtmix_0.95": 217.1}))

    assert tip_results(tmp_path, "fig9_l9_rpm30")["dtmix_0.95"] == 217.1


# ── units ────────────────────────────────────────────────────────────────────

def test_period_used_for_cycle_counts_is_dimensional():
    """`_t_scales` returns (T_bio [s], T_per_nd [-]) -- its SECOND element is a
    period in non-dimensional units, not seconds. Sizing a segment with it
    inflated the cycle count by T_bio (2.63x at 37.5 rpm) and would have bought
    a 2.6x longer job than the extrapolation asked for."""
    import math as _m
    from scripts.autoextend import period_seconds
    from scripts.postprocess import _t_scales

    params = {"omega_b": 37.5 * 2 * _m.pi / 60.0,
              "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0},
              "theta_max": [7.0, 0.0, 0.0]}
    assert period_seconds(params) == pytest.approx(60.0 / 37.5)
    assert period_seconds(params) != pytest.approx(_t_scales(params)[1], rel=0.1)


def test_untested_rpm_buys_a_bigger_walltime_cushion():
    """min_per_cycle falls back to the max over the rpms measured at that
    (level, ntasks). With a single entry in the table that fallback is not
    conservative at all, and cost per cycle RISES as rpm drops -- 2.9x from
    32.5 to 17.5 rpm at L10/32. A fallback cost therefore needs a cushion
    wide enough to cover that spread, or the extension times out and the
    segment is lost."""
    from scripts.autoextend import walltime_safety

    assert walltime_safety("measured") < walltime_safety("measured_other_rpm_max")
    assert walltime_safety("measured_other_rpm_max") >= 2.9
