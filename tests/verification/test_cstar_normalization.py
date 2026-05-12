"""Verify that C* (dissolved O2 saturation) is correctly normalized to [0, 1].

Root cause of the bug (inherited from upstream rcsc-group/BioReactor):
  1. f_liq[] = (1-cs[])*f[]: cs=1 inside the bag, so this is zero in all bulk
     liquid cells. Only bag-wall cut cells contribute → f_liq_sum ≈ 0.002 instead
     of the true liquid volume ≈ 0.284. Dividing by this tiny number blows C* past 1.
  2. oxy_liq[] = oxy[]*f[]: oxy[] is a Henry's-law mixed quantity; at interface cells
     it contains gas-phase oxygen, inflating the integral.
  3. postprocess.py omits the c_oxy_alpha = 1/30 denominator factor.

After the fix:
  - f_liq[] = f[]  (true liquid volume fraction)
  - oxy_liq[] = f*alpha*oxy/(f*alpha + (1-f))  (liquid-phase-only oxygen)
  - postprocess.py: c_star = oxy_liq_sum / (C_OXY_ALPHA * f_mean)
"""
import sys
import pathlib
import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, run_bioreactor
from scripts.postprocess import _compute_c_star

# t_mix ≈ 48.65 non-dim for canonical params; run just past it to get oxygen data
PARAMS_OXY = {**CANONICAL_PARAMS, "run_id": "cstar_norm", "t_end": 60.0}


@pytest.mark.medium
def test_c_star_bounded_between_0_and_1(tmp_path):
    """C* must stay in [0, 1] at all time steps after oxygen transfer begins.

    Before the fix: f_liq_sum ≈ 0.002 (only wall cut cells) and oxy_liq includes
    gas-side oxygen → C* > 5 at first non-zero row, identical for all thresholds.
    After the fix: C* rises smoothly from 0 toward 1.
    """
    run_dir = run_bioreactor(PARAMS_OXY, tmp_path, timeout=120)
    t, c_star = _compute_c_star(run_dir)

    oxy_rows = c_star[c_star > 0]
    assert len(oxy_rows) > 0, "No oxygen transfer detected — t_end may be < t_mix"
    assert float(c_star.max()) <= 1.0, (
        f"C* exceeds 1.0 (max={c_star.max():.3f}): "
        "f_liq_sum or oxy_liq normalization is wrong"
    )
    assert float(c_star.min()) >= 0.0, f"C* is negative (min={c_star.min():.3f})"


@pytest.mark.medium
def test_kla_values_differ_across_saturation_levels(tmp_path):
    """kLa_10, kLa_25, kLa_50 must not all be identical.

    Before the fix: all three thresholds are crossed at the same row (C* jumps
    above 1 at the first non-zero step) → identical kLa from the same window.
    After the fix: C* rises gradually, thresholds are crossed at different times.
    """
    from scripts.postprocess import main as postprocess_main

    run_dir = run_bioreactor(PARAMS_OXY, tmp_path, timeout=120)
    results = postprocess_main(str(run_dir))

    kla_10 = results["kLa_10"]
    kla_25 = results["kLa_25"]
    kla_50 = results["kLa_50"]

    finite = [v for v in (kla_10, kla_25, kla_50) if v == v]  # exclude NaN
    assert len(finite) >= 2, (
        f"Fewer than 2 finite kLa values — run may be too short: {results}"
    )
    assert len(set(round(v, 6) for v in finite)) > 1, (
        f"kLa_10={kla_10:.4f}, kLa_25={kla_25:.4f}, kLa_50={kla_50:.4f} are all identical: "
        "C* normalization is wrong — all thresholds crossed at the same time step"
    )
