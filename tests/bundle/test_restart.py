"""Guards of the restart mechanisms, read from the feature bundle (bars: tools/bundle/bars.yaml).

S5  same-condition chaining: a run stopped at 4T and continued (cont), and a run stopped early by its
    wall-clock limit and resumed (kill -> kill_resume), must reproduce the unbroken run (cold).
S9  cs/fs are recomputed after every restore: the solver prints "fs sanity [post-solid]: OK".
S10 tracer and oxygen carry across a restart: the restored fields equal the dumped ones.
S12 a run stops on the rocking-period boundary that scripts/periods.py predicts.
"""
import math
from pathlib import Path

import pytest

from tests.bundle import _series as s
from tests.bundle.conftest import case
from tools.bundle import cases as bc

T = bc.period(bc.RPM_A)
RESTARTED = ["cont", "kill_resume", "replay", "loop", "xl6", "caseB_warm"]


def rd(bundle, name) -> Path:
    return Path(case(bundle, name)["dir"])


def _rows_by_time(path: Path, t_max: float) -> dict:
    """{time token: whole line} for the rows with t <= t_max (printed text, so equality is bitwise)."""
    out = {}
    for ln in path.read_text().splitlines()[1:]:
        tok = ln.split()
        if tok and float(tok[1]) <= t_max + s.EPS:
            out[tok[1]] = ln
    return out


def test_s5_stopping_early_leaves_the_history_unchanged(bundle, bars):
    bar = bars["half"]
    _, t1 = s.window(bar["window"])
    cold, half = rd(bundle, "cold"), rd(bundle, "half")
    for f in bar["files"]:
        a, b = _rows_by_time(cold / f, t1), _rows_by_time(half / f, t1)
        common = sorted(set(a) & set(b), key=float)
        assert len(common) >= 0.9 * len(a), f"{f}: half logged {len(b)} rows of cold's {len(a)} before {t1}"
        diff = [t for t in common if a[t] != b[t]]
        assert not diff, f"{f}: {len(diff)} rows differ, first at t={diff[0]}:\n {a[diff[0]]}\n {b[diff[0]]}"


@pytest.mark.parametrize("name", ["cont", "kill_resume"])
def test_s5_continuation_reproduces_the_unbroken_run(bundle, bars, name):
    r = s.compare(bars[name], rd(bundle, "cold"), rd(bundle, name))
    print(f"\n{name}:", {m: f"{d:.2e}" for m, (_, _, d, _) in r.items()})
    assert not s.failures(r), "\n".join(s.failures(r))


def test_s5_early_dump_writes_dump_and_phase(bundle):
    k = rd(bundle, "kill")
    t = s.dump_time(k / "checkpoint.dump")
    lo, hi = bc.N_RAMP * T, (bc.END - 1) * T   # after the ramp, before kill_resume's window [7T, 8T]
    assert lo < t < hi, f"early dump at t={t:.4f}, expected inside ({lo:.4f}, {hi:.4f})"
    assert (k / "checkpoint.phase").is_file(), "early dump wrote no phase sidecar"
    assert "checkpointing at t=" in (k / "stderr.txt").read_text(), "the wall-clock limit did not fire"


@pytest.mark.parametrize("name", RESTARTED)
def test_s9_fs_is_sane_after_every_restore(bundle, name):
    err = (rd(bundle, name) / "stderr.txt").read_text()
    assert "fs sanity [post-solid]: OK" in err, f"{name}: no fs sanity line after the restore"
    assert "ERROR [" not in err


# kill -> kill_resume is not here: the early-dump path writes no pre_dump diagnostic (S5 covers it)
@pytest.mark.parametrize("parent,child", [("half", "cont")])
def test_s10_restore_gives_back_the_dumped_fields(bundle, parent, child):
    """restart_diagnostic_pre_dump.txt (written just before the dump) vs ..._post_restore.txt (just after the
    restore, before init changes anything): every field statistic must be equal to round-off."""
    def stats(p):
        return {k: float(v) for k, v in (ln.split() for ln in p.read_text().splitlines()[1:] if ln.strip())}
    a = stats(rd(bundle, parent) / "restart_diagnostic_pre_dump.txt")
    b = stats(rd(bundle, child) / "restart_diagnostic_post_restore.txt")
    assert set(a) == set(b), f"different fields: {sorted(set(a) ^ set(b))}"
    bad = [f"{k}: {a[k]!r} vs {b[k]!r}" for k in a if s.rel(a[k], b[k]) > 1e-12 and abs(a[k] - b[k]) > 1e-300]
    assert not bad, "restored fields differ from the dumped ones:\n" + "\n".join(bad)


def test_s10_continuation_keeps_the_tracer(bundle):
    _, c2 = s.column(rd(bundle, "cont"), "tr_oxy.dat", "c2_liq_sum")
    assert c2[0] > 0, "restart_continue=1 wiped the tracer at the seam"


@pytest.mark.parametrize("name,t_request", [("cold", 8 * T), ("half", 4 * T), ("cont", 8 * T),
                                            ("kill_resume", 8 * T), ("rec", 2 * T)])
def test_s12_run_stops_on_the_period_boundary(bundle, name, t_request):
    t = s.dump_time(rd(bundle, name) / "checkpoint.dump")
    want = s.stop_time(t_request)
    assert abs(t - want) <= 1e-9 * want, f"{name} stopped at t={t!r}, periods.py says {want!r}"


# S13 restart label (internal/docs/plans/2026-11-13-restart-label-plan.md): the dump's label carries the forcing that made
# the state; the restart ramps from it to params.json. The angle and its rate must not jump at the restart,
# for a continuation (same forcing), an early dump, a new level, and a new condition (amplitude AND frequency).
FORCED = ["cont", "kill_resume", "xl6", "caseB_warm"]


@pytest.mark.parametrize("name", FORCED)
def test_s13_restart_continues_the_angle_and_its_rate(bundle, name):
    """The parent's angle at the dump, from its label (Th, Th_d, Th_2d at t_forcing, one step before the dump,
    Taylor-extended), vs the restarted run's angle at its first acceleration step (restart_forcing.json).
    That step is one short step after the dump: restore() marks the solver events of the restored step as done
    (src/output.h), so that step moves t with frozen fields (seen with qcc -events, 2026-11-13; pitfall
    "restore-first-step"). Bar: at most 1e-3 period. Tolerances 1e-6 rad and 1e-3 of the peak rate: the
    Taylor rest is ~5e-7 rad / ~2e-4 at L6, while a ramp from rest jumps by ~0.1 rad and the full rate."""
    import json
    lab = json.loads((rd(bundle, bc.CASES[name].parent) / "checkpoint.phase").read_text())
    new = json.loads((rd(bundle, name) / "restart_forcing.json").read_text())
    # seconds from t_forcing to the restart's first step: parent units up to the dump, the child's after it
    # (a new condition has its own time unit T_bio)
    dt = (lab["t"] - lab["t_forcing"]) * lab["T_bio"] + (new["t"] - lab["t"]) * new["T_bio"]
    th = lab["Th"] + lab["Th_d"]*dt + 0.5*lab["Th_2d"]*dt*dt
    th_d = lab["Th_d"] + lab["Th_2d"]*dt
    peak = abs(lab["omega_b"]) * lab["theta_max"][0] * math.pi / 180
    lag = new["cycles"] - lab["cycles"]   # phase advance from the dump to the first step, in periods
    assert 0 <= lag <= 1e-3, f"first forcing step {lag:.2e} periods after the dump (bar 1e-3)"
    assert abs(new["Th"] - th) < 1e-6, f"angle jumps: {th!r} -> {new['Th']!r} rad"
    assert abs(new["Th_d"] - th_d) <= 1e-3 * peak, f"angular rate jumps: {th_d!r} -> {new['Th_d']!r} rad/s"


@pytest.mark.parametrize("name,msg", [("refuse_prev", "restart refused: params.json sets *_prev keys"),
                                      ("refuse_nolabel", "restart refused: no checkpoint.phase"),
                                      ("refuse_stale", "restart refused: checkpoint.phase is the label of t=")])
def test_s13_restart_refuses_an_ambiguous_source(bundle, name, msg):
    """*_prev next to a label (two sources), a restart with no label (no source), and a label written at
    another time than the dump (a stale label: wrong source) must stop the solver."""
    case(bundle, name)   # fails unless the solver stopped (status "refused")
    err = (rd(bundle, name) / "stderr.txt").read_text()
    assert "ERROR: " + msg in err, err[-500:]
