"""tools/bundle/cases.py: the L6 case graph of the feature bundle (no solver, no SLURM)."""
import math
from pathlib import Path

import pytest

from tools.bundle import cases as bc

RUNS = Path("/tmp/bundle_runs")


def test_case_names_are_the_plan_list():
    assert set(bc.CASES) == {"cold", "cold2", "half", "cont", "kill", "kill_resume", "replay", "loop", "lean",
                             "l5", "xl6", "caseB_cold", "caseB_warm", "mpi4"}


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
    assert p["rec_dir"] == str(RUNS / "cold")
    assert bc.CASES["replay"].restart == "replay_start.dump"


def test_kill_sets_a_wall_clock_limit():
    assert "BIOREACTOR_MAXRUNTIME_S" in bc.CASES["kill"].env
