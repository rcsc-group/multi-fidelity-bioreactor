"""Guards of flux record/replay, read from the feature bundle (bars: tools/bundle/bars.yaml).

S1  exact replay: rec recorded the flux of [T, 2T]; replay solves only the scalars over that period and
    must give the coupled run's scalar sums (F1 of tests/test_replay_gate.py, relative 1e-6).
S2  looped replay: the one recorded period, played twice over [T, 3T]. No pass bar (the physics bars live
    in the diary); the guard checks it ran the whole window and its scalars are finite, and prints the
    error against cold for the ledger. Replay builds write no final dump: the end is the last logged time.
"""
from pathlib import Path

import numpy as np

from tests.bundle import _series as s
from tests.bundle.conftest import case
from tools.bundle import cases as bc

T = bc.period(bc.RPM_A)


def rd(bundle, name) -> Path:
    return Path(case(bundle, name)["dir"])


def test_s1_replay_reproduces_the_coupled_scalars(bundle, bars):
    ref, rep = rd(bundle, "rec"), rd(bundle, "replay")
    r = s.compare(bars["replay"], ref, rep)
    print("\nreplay:", {m: f"{d:.2e}" for m, (_, _, d, _) in r.items()})
    assert not s.failures(r), "\n".join(s.failures(r))
    t, _ = s.column(ref, "tr_oxy.dat", "oxy_liq_sum")
    per_T = np.sum((t > T) & (t <= 2 * T))
    _, n = s.series_max("oxy_liq_sum", ref, rep, s.window([1, 2]))
    assert n >= 0.9 * per_T, f"replay shares {n} of rec's {per_T} logged times in [T, 2T]"
    _, c2 = s.column(rep, "tr_oxy.dat", "c2_liq_sum")
    assert np.max(np.abs(c2)) > 0, "trivial comparison: no tracer in the replayed window"


def test_s2_looped_replay_runs_the_window_with_finite_scalars(bundle):
    loop = rd(bundle, "loop")
    t, _ = s.column(loop, "tr_oxy.dat", "oxy_liq_sum")
    dt = float(np.median(np.diff(t)))
    assert 3 * T - 2 * dt <= t[-1] <= 3 * T + s.EPS, f"loop's last logged time {t[-1]:.4f}, window ends at 3T={3 * T:.4f}"
    for col in ("oxy_liq_sum", "c2_liq_sum"):
        _, v = s.column(loop, "tr_oxy.dat", col)
        assert np.all(np.isfinite(v)), f"loop: non-finite {col}"
        d, n = s.series_max(col, rd(bundle, "cold"), loop, s.window([1, 3]))
        print(f"\nloop vs cold, {col}: max rel diff {d:.3e} over {n} times (ledger, no bar)")
