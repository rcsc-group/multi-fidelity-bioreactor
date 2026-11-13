"""Guards of the warm starts, read from the feature bundle (bars: tools/bundle/bars.yaml).

S3, S15 cross-level warm start: l5 (L5, 2T) restarted at L6 (xl6, with the pre-refine order of
        CROSS_LEVEL_WARMSTART) must match the cold L6 run over [2T, 3T] (test_cross_level_warmstart bars).
S4      cross-case warm start: condition A's state at 2T, restarted at condition B (caseB_warm), must match a
        cold run at B (caseB_cold) over each run's last period (x8 bar, 5 %).
S8      early tracer release after a warm start: caseB_warm releases its tracer after 1 transition cycle.
"""
from pathlib import Path

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
