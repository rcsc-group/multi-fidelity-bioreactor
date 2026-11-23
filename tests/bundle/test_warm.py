"""Guards of the warm starts, read from the feature bundle (bars: tools/bundle/bars.yaml).

S3, S15 cross-level warm start: l5 (L5, to 4T) restarted at L6 (xl6, with the pre-refine order of
        CROSS_LEVEL_WARMSTART) must match the cold L6 run over [6T, 8T] (test_cross_level_warmstart bars).
S4      cross-case warm start: condition A's state at 4T, restarted at condition B (caseB_warm), must match a
        cold run at B (caseB_cold) over each run's last 3 periods, both ending at 14 periods of B (x8 bar,
        max(5 %, 2 sd); the warm start settles ~10 periods after the switch: lean diary 2026-11-14).
S8      early tracer release after a warm start: caseB_warm releases its tracer after 1 transition cycle.
"""
from pathlib import Path

import pytest
import numpy as np

from tests.bundle import _series as s
from tests.bundle.conftest import case
from tools.bundle import cases as bc


def rd(bundle, name) -> Path:
    return Path(case(bundle, name)["dir"])


def test_s3_cross_level_warm_start_matches_cold(bundle, bars):
    r = s.compare(bars["xl6"], rd(bundle, "cold"), rd(bundle, "xl6"))
    print("\nxl6:", {m: f"{d:.2e}" for m, (_, _, d, _) in r.items()})
    assert not s.failures(r), "\n".join(s.failures(r))


def test_s4_cross_case_warm_start_matches_cold(bundle, bars):
    r = s.compare(bars["caseB_warm"], rd(bundle, "caseB_cold"), rd(bundle, "caseB_warm"), rpm=bc.RPM_B)
    print("\ncaseB_warm:", {m: f"{d:.2e}" for m, (_, _, d, _) in r.items()})
    assert not s.failures(r), "\n".join(s.failures(r))


def test_s8_tracer_released_after_the_transition(bundle):
    w = rd(bundle, "caseB_warm")
    t, c2 = s.column(w, "tr_oxy.dat", "c2_liq_sum")
    assert c2[0] == 0, "tracer present at the start of a cross-case restart (it must be wiped and re-seeded)"
    assert c2[-1] > 0, "tracer never released after the transition cycle"
    t_rel = t[np.argmax(c2 > 0)]
    T_b = bc.period(bc.RPM_B)
    print(f"\ncaseB_warm: tracer released {(t_rel - t[0]) / T_b:.2f} periods after the restart")


# S16 seam velocity (diary 2026-11-14: mutant no_rescale was NOT caught by S4, which measures the settled state).
# A warm start must carry the PHYSICAL velocity across the restart. The solver's velocity unit is L / T_bio, so
# physical speed ~ v / T_bio. Compared: the liquid rms speed sqrt(ux_rms^2 + uy_rms^2) (normf.dat) on the parent's
# last row at or before the dump vs the child's first row, each divided by its own T_bio (parent: checkpoint.phase
# label; child: restart_forcing.json). Bar 10 %: the parent's row-to-row change near the dump is ~5 % at L6
# (head8: 0.750 -> 0.787 over one output step); without the rescale the jump is |T_A/T_B - 1| = 23 % (32.5 -> 25 rpm,
# seen in bundle mut_rescale2: child 0.799 nondimensional = parent 0.787, i.e. not rescaled).
SEAM = ["caseB_warm", "caseT_warm", "xl6", "cont"]


def _speed(rd, row):
    cols, a = s.table(rd, "normf.dat")
    r = a[row]
    return float(np.hypot(r[cols.index("ux_liq_rms")], r[cols.index("uy_liq_rms")])), float(r[cols.index("t")])


@pytest.mark.parametrize("name", SEAM)
def test_s16_warm_start_keeps_the_physical_speed_at_the_seam(bundle, name):
    import json
    child = rd(bundle, name)
    parent = rd(bundle, bc.CASES[name].parent)
    lab = json.loads((parent / "checkpoint.phase").read_text())
    new = json.loads((child / "restart_forcing.json").read_text())
    cols, a = s.table(parent, "normf.dat")
    t = a[:, cols.index("t")]
    i = int(np.nonzero(t <= lab["t"] + 1e-9)[0][-1])
    vp, tp = _speed(parent, i)
    vc, tc = _speed(child, 0)
    up, uc = vp / lab["T_bio"], vc / new["T_bio"]
    d = abs(uc - up) / up
    print(f"\n{name}: parent {vp:.4f} (t {tp}), child {vc:.4f} (t {tc}); physical rel diff {d:.3f}")
    assert d <= 0.10, f"physical liquid speed jumps by {d:.1%} at the restart (bar 10 %)"


# S17 cross-theta warm start (diary 2026-11-22, pre-registered): half (theta 7, 4T) restarted at theta 2, same rate,
# must match a cold theta-2 run over each run's last 3 periods (the S4 bar, unchanged). S13 and S16 cover the seam.
def test_s17_cross_theta_warm_start_matches_cold(bundle, bars):
    r = s.compare(bars["caseT_warm"], rd(bundle, "caseT_cold"), rd(bundle, "caseT_warm"), theta=bc.THETA_T)
    print("\ncaseT_warm:", {m: f"{d:.2e}" for m, (_, _, d, _) in r.items()})
    assert not s.failures(r), "\n".join(s.failures(r))
