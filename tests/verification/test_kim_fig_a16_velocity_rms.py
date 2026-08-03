"""Physical verification: u'_x,rms and u'_y,rms match Kim et al.'s
published Fig. A.16(a)/(b) at their exact baseline condition.

This is the strongest guardrail in this file: an external, physically-
grounded check against the published paper's reported numbers, not
against our own code's internal formulas (which is what let the
2026-08-03 bug hide since the params.json pipeline's inception -- every
internal ratio/ramp/frequency check stayed self-consistent because they
all derived from the same wrong H_bio formula AND the same wrong
geometry.b default, and self-consistency checks can't detect a constant
error shared by both the code and the test replicating it).

Baseline condition (Main.tex Sec. 4, `theta_b,max=7 deg`, `f_b=32.5rpm`),
run from a fresh cold start out to t/T_p=31 (Fig. A.16's caption window
is t/T_p in [29,31], "over two cycles"). Published values, read directly
off Fig. A.16(a)/(b): u'_x,rms peak ~0.8, u'_y,rms peak ~0.21 -- visible
at every resolution Kim tested, including their coarsest (n_L=2^5), so
this is not a resolution-dependent target.

Failure modes caught:
  - Any future H_bio/U_bio/geometry.b scale bug that distorts the actual
    simulated physics (this is what happened 2026-08-03: with the old
    H_bio formula AND old geometry.b=0.071 default, u'_y,rms matched at
    ~0.22 but u'_x,rms sat at ~0.39 -- half of Kim's value -- since this
    project's params.json pipeline began, because nothing checked
    against an external reference. Fixing only the H_bio formula wasn't
    enough either: it moved u'_x,rms to ~0.49 while pushing u'_y,rms to
    ~0.28, still wrong in both directions, until geometry.b's default
    was corrected too -- see diary.md 2026-08-03 for both confirmation
    runs).
  - Any change (advection scheme, BC, forcing term, embedded solid setup)
    that quietly stops reproducing the actual published physics, even if
    every internal self-consistency check still passes.
"""
import math
import sys
import pathlib

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, load_normf, run_bioreactor

_KIM_UX_RMS_PEAK = 0.80
_KIM_UY_RMS_PEAK = 0.21
_REL_TOL = 0.25  # generous: reading a peak off a published figure by eye


def _t_per_st(params: dict) -> float:
    omega_b = params["omega_b"]
    L_bio   = params["geometry"]["a"]
    H_bio   = 2 * params["geometry"]["b"]  # full height (BioReactor.c:295)
    th_max  = math.radians(params["theta_max"][0])

    T_per  = 2 * math.pi / omega_b
    V_bio  = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th_max))
    U_bio  = V_bio / (H_bio * 0.5) / T_per
    T_bio  = L_bio / U_bio
    return T_per / T_bio


@pytest.mark.hpc
def test_ux_uy_rms_match_kim_baseline(tmp_path):
    """Peak u'_x,rms/u'_y,rms over t/T_p in [29,31] must be within 25% of
    Kim et al.'s published Fig. A.16(a)/(b) values (~0.8, ~0.21).
    """
    omega_b_32_5rpm = 2 * math.pi * 32.5 / 60
    params = {
        **CANONICAL_PARAMS,
        "run_id": "kim_fig_a16_baseline",
        "fidelity": 6,
        "omega_b": omega_b_32_5rpm,
        "theta_max": [7.0, 0.0, 0.0],
        "n_mix_cycles": 0,
        "t_end": 19.0,
    }
    run_dir = run_bioreactor(params, tmp_path, timeout=1800)

    data = load_normf(run_dir)
    assert len(data) >= 20, "Too few output rows — simulation may not have run"

    t = data[:, 1]
    ux_rms = data[:, 7]
    uy_rms = data[:, 11]

    t_per_st = _t_per_st(params)
    t_tp = t / t_per_st
    mask = (t_tp >= 29.0) & (t_tp <= 31.0)
    assert mask.sum() >= 10, (
        f"Only {mask.sum()} samples in t/T_p=[29,31] window — run didn't "
        f"reach the target window (max t/T_p={t_tp.max():.2f}, t_end too "
        f"short for this omega_b/T_per_st combination)"
    )

    ux_peak = ux_rms[mask].max()
    uy_peak = uy_rms[mask].max()

    ux_err = abs(ux_peak - _KIM_UX_RMS_PEAK) / _KIM_UX_RMS_PEAK
    uy_err = abs(uy_peak - _KIM_UY_RMS_PEAK) / _KIM_UY_RMS_PEAK

    assert ux_err < _REL_TOL, (
        f"u'_x,rms peak = {ux_peak:.3f}, Kim et al. report ~{_KIM_UX_RMS_PEAK} "
        f"(off by {ux_err:.1%}, tolerance {_REL_TOL:.0%}). Check both "
        "H_bio's formula (BioReactor.c:295) and geometry.b's default "
        "value -- the 2026-08-03 bug required fixing both together."
    )
    assert uy_err < _REL_TOL, (
        f"u'_y,rms peak = {uy_peak:.3f}, Kim et al. report ~{_KIM_UY_RMS_PEAK} "
        f"(off by {uy_err:.1%}, tolerance {_REL_TOL:.0%})."
    )
