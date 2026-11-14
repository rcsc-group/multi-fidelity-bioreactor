"""tools/bundle/cases.py: the L6 case graph of the feature bundle (no solver, no SLURM)."""
import math
from pathlib import Path

import pytest

from tools.bundle import cases as bc

RUNS = Path("/tmp/bundle_runs")


def test_case_names_are_the_plan_list():
    assert set(bc.CASES) == {"cold", "cold2", "rec", "half", "cont", "kill", "kill_resume", "replay", "loop",
                             "lean", "l5", "xl6", "caseB_cold", "caseB_warm", "mpi4", "refuse_prev", "refuse_nolabel",
                             "refuse_stale"}


def test_order_puts_parents_first():
    seq = bc.order(bc.CASES)
    assert sorted(seq) == sorted(bc.CASES)
    for name, c in bc.CASES.items():
        if c.parent:
            assert seq.index(c.parent) < seq.index(name)


def test_cycle_raises():
    a = bc.Case("a", "default", "b", "checkpoint.dump", 6, 1)
    b = bc.Case("b", "default", "a", "checkpoint.dump", 6, 1)
    with pytest.raises(ValueError, match="cycle"):
        bc.order({"a": a, "b": b})


def test_unknown_parent_raises():
    a = bc.Case("a", "default", "ghost", "checkpoint.dump", 6, 1)
    with pytest.raises(ValueError, match="ghost"):
        bc.order({"a": a})


def test_every_case_has_a_known_build():
    for c in bc.CASES.values():
        assert c.build in bc.BUILDS, c.name


def test_restart_cases_set_t_checkpoint():
    for name, c in bc.CASES.items():
        p = bc.params(name, RUNS)
        if c.restart:
            assert p.get("t_checkpoint", 0) > 0, name   # the solver restarts only when t_checkpoint > 0
        else:
            assert "t_checkpoint" not in p, name


def test_caseB_warm_rescales_from_32p5():
    p = bc.params("caseB_warm", RUNS)
    assert p["omega_b"] == pytest.approx(25.0 * 2 * math.pi / 60)
    assert bc.CASES["caseB_warm"].parent == "half"   # the source forcing (32.5 rpm) comes from half's label
    assert p["n_transition_cycles"] == 1


def test_replay_reads_the_cold_record():
    p = bc.params("replay", RUNS)
    assert p["rec_dir"] == str(RUNS / "rec")
    assert bc.CASES["replay"].parent == "rec"
    assert bc.CASES["replay"].restart == "replay_start.dump"


def test_kill_sets_a_wall_clock_limit():
    assert "BIOREACTOR_MAXRUNTIME_S" in bc.CASES["kill"].env


# v2 (first bundle, 2026-11-13): the 3-cycle ramp, *_prev, and the record build's extra step landing

def test_references_are_default_builds():
    """The record build lands a step on the record end (replay.h:228), so after it its step sequence differs
    from the default and lean builds; the reference runs must be default builds."""
    for name in ("cold", "cold2", "half", "kill"):
        assert bc.CASES[name].build == "default", name
    assert bc.CASES["rec"].build == "record"


def test_restarts_start_after_the_ramp():
    T = bc.period(bc.RPM_A)
    for name in ("cont", "xl6", "caseB_warm"):
        assert bc.params(name, RUNS)["t_checkpoint"] >= bc.N_RAMP * T, name


def test_no_case_sets_prev_keys_but_the_refusal_case():
    """The source forcing of a restart comes from the dump's label (checkpoint.phase), never from params.json
    (plan docs/plans/2026-11-13-restart-label-plan.md). Only refuse_prev sets *_prev, to test the refusal."""
    for name in bc.CASES:
        keys = [k for k in bc.params(name, RUNS) if k.endswith("_prev")]
        assert bool(keys) == (name == "refuse_prev"), (name, keys)


def test_refusal_cases_expect_the_solver_to_stop():
    for name in ("refuse_prev", "refuse_nolabel", "refuse_stale"):
        c = bc.CASES[name]
        assert c.expect_refusal and c.parent == "half" and c.restart == "checkpoint.dump", name
    assert bc.CASES["refuse_prev"].stage_label and not bc.CASES["refuse_nolabel"].stage_label
    assert bc.CASES["refuse_stale"].label_from == "cold"   # cold's label (8T) next to half's dump (4T)
    assert not any(c.expect_refusal for n, c in bc.CASES.items() if not n.startswith("refuse_"))


def test_continuations_keep_the_tracer_and_the_oxygen_event():
    for name in ("cont", "kill_resume"):
        p = bc.params(name, RUNS)
        assert p["restart_continue"] == 1 and p["n_mix_cycles"] == 0, name


def test_caseB_runs_end_after_the_warm_start_transient():
    """S4 extension (lean diary 2026-11-14, pre-registered, TRANSIENT): warm vs cold settles within the bar
    after ~10 periods at B, so both caseB runs end at CASEB_END = 14 periods of B (8 before)."""
    TB = bc.period(bc.RPM_B)
    assert bc.CASEB_END == 14
    assert abs(bc.params("caseB_cold", RUNS)["t_end"] - 14 * TB) < 1e-12
    w = bc.params("caseB_warm", RUNS)
    assert abs(w["t_end"] - 14 * TB) < 1e-12
