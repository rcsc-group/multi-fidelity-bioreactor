"""Parameter construction of the low-fidelity looped-replay pipeline (scripts/replay_lf.py); no solver, no Slurm."""
import json
from pathlib import Path

import pytest

from scripts import postprocess as pp
from scripts import replay_lf as lf

SRC = lf.SRC_RUNS / "kmix_l8_rpm32.5"
pytestmark = pytest.mark.skipif(not (SRC / "params.json").exists(), reason="source run not on this machine")


def test_source_ids():
    assert lf.source_id(32.5, 7.0) == "kmix_l8_rpm32.5"
    assert lf.source_id(15.0, 5.5) == "t6_l8_rpm15_th5.5"
    assert len(lf.conditions()) == 40


def test_record_and_replay_params():
    sp = json.loads((SRC / "params.json").read_text())
    _, T = pp._t_scales(sp)
    rec = lf.record_params(sp, "r", 123.0, Path("/x"))
    assert rec["n_mix_cycles"] == 1 and "restart_continue" not in rec and "t_rec0" not in rec
    assert rec["n_rec_periods"] == 1 and rec["t_end"] == pytest.approx(2 * T, rel=1e-5)
    assert rec["omega_b_prev"] == sp["omega_b"] and rec["theta_max_prev"] == sp["theta_max"]
    rep = lf.replay_params(sp, "p", 123.0 + T, Path("/x"))
    assert rep["restart_continue"] == 1 and rep["n_mix_cycles"] == 0 and rep["replay_loop"] == 1
    assert rep["t_end"] == pytest.approx(22 * T) and rep["replay_scalar_level"] == 8


def test_checkpoint_time_matches_cycle_230():
    _, T = pp._t_scales(json.loads((SRC / "params.json").read_text()))
    assert lf.checkpoint_time(SRC / "checkpoint.dump") / T == pytest.approx(230.0, abs=0.01)
