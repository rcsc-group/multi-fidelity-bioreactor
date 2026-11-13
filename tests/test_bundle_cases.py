"""tools/bundle/cases.py: the L6 case graph of the feature bundle (no solver, no SLURM)."""
import math
from pathlib import Path

import pytest

from tools.bundle import cases as bc

RUNS = Path("/tmp/bundle_runs")


def test_case_names_are_the_plan_list():
    assert set(bc.CASES) == {"cold", "cold2", "rec", "half", "cont", "kill", "kill_resume", "replay", "loop",
                             "lean", "l5", "xl6", "caseB_cold", "caseB_warm", "mpi4"}


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
            assert p.get("restart_continue") == 1 or "omega_b_prev" in p, name
        else:
            assert "t_checkpoint" not in p, name


def test_caseB_warm_rescales_from_32p5():
    p = bc.params("caseB_warm", RUNS)
    assert p["omega_b"] == pytest.approx(25.0 * 2 * math.pi / 60)
    assert p["omega_b_prev"] == pytest.approx(32.5 * 2 * math.pi / 60)
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


def test_restarts_carry_the_current_forcing():
    """The solver ramps *_prev -> current over 3 cycles after every restart; without *_prev it ramps from 0."""
    for name in ("cont", "kill_resume", "xl6"):
        p = bc.params(name, RUNS)
        assert p["omega_b_prev"] == p["omega_b"] and p["theta_max_prev"] == p["theta_max"], name
        assert p["amplitude_h_prev"] == p["amplitude_h"] and p["omega_h_prev"] == p["omega_h"], name


def test_continuations_keep_the_tracer_and_the_oxygen_event():
    for name in ("cont", "kill_resume"):
        p = bc.params(name, RUNS)
        assert p["restart_continue"] == 1 and p["n_mix_cycles"] == 0, name
