"""Tests for scripts/chain.py — chained simulation via checkpoint restart."""
import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1]))
from scripts.chain import build_chain

# ── minimal valid config ──────────────────────────────────────────────────────
_CFG = {
    "fidelity": 5,
    "geometry": {"a": 0.25, "b": 0.071, "n": 2.0},
    "fill_level": 0.7,
    "n_mix_cycles": 80,
    "n_transition_cycles": 10,
    "t_buffer": 150.0,
    "sweep": {
        "parameter": "omega_b",
        "values": [1.0, 1.1, 1.2],
    },
    "motion": {
        "n_harmonics": 1,
        "theta_max": [5.0, 0.0, 0.0],
        "phi_angular": [0.0, 0.0, 0.0],
        "omega_h": 0.0,
        "amplitude_h": [0.0, 0.0, 0.0],
        "phi_horizontal": [0.0, 0.0, 0.0],
    },
    "submit": False,
    "walltime": "01:00:00",
}


def test_build_chain_length():
    """3 sweep values → 3 param dicts."""
    chain = build_chain(_CFG)
    assert len(chain) == 3


def test_segment0_has_full_mix_cycles():
    """First segment uses the full n_mix_cycles (80), not the transition value."""
    chain = build_chain(_CFG)
    assert chain[0]["n_mix_cycles"] == _CFG["n_mix_cycles"]


def test_restart_segments_have_transition_cycles():
    """Segments 1+ use n_transition_cycles as n_mix_cycles for faster settling."""
    chain = build_chain(_CFG)
    for k in range(1, len(chain)):
        assert chain[k]["n_mix_cycles"] == _CFG["n_transition_cycles"], (
            f"segment {k} n_mix_cycles={chain[k]['n_mix_cycles']}, "
            f"expected {_CFG['n_transition_cycles']}"
        )


def test_restart_t_end_is_relative_and_shorter():
    """Restart segments have a shorter t_end than segment 0 (relative, not absolute)."""
    chain = build_chain(_CFG)
    t_end_0 = chain[0]["t_end"]
    for k in range(1, len(chain)):
        assert chain[k]["t_end"] < t_end_0, (
            f"segment {k} t_end={chain[k]['t_end']:.2f} should be < "
            f"segment 0 t_end={t_end_0:.2f}"
        )


def test_sweep_parameter_varies():
    """Each segment has the sweep parameter value from cfg['sweep']['values']."""
    chain = build_chain(_CFG)
    param = _CFG["sweep"]["parameter"]
    values = _CFG["sweep"]["values"]
    for k, expected in enumerate(values):
        assert math.isclose(chain[k][param], expected), (
            f"segment {k} {param}={chain[k][param]}, expected {expected}"
        )


def test_non_swept_motion_params_fixed():
    """Motion params not being swept are the same across all segments."""
    chain = build_chain(_CFG)
    for k in range(1, len(chain)):
        assert chain[k]["theta_max"] == _CFG["motion"]["theta_max"]
        assert chain[k]["n_harmonics"] == _CFG["motion"]["n_harmonics"]
        assert chain[k]["omega_h"] == _CFG["motion"]["omega_h"]


def test_fixed_params_propagate():
    """Fidelity, geometry, fill_level are present and correct in every segment."""
    chain = build_chain(_CFG)
    for k, p in enumerate(chain):
        assert p["fidelity"] == _CFG["fidelity"], f"segment {k} fidelity"
        assert p["geometry"] == _CFG["geometry"], f"segment {k} geometry"
        assert math.isclose(p["fill_level"], _CFG["fill_level"]), f"segment {k} fill_level"


def test_each_segment_has_unique_run_id():
    """Every param dict gets a distinct run_id."""
    chain = build_chain(_CFG)
    ids = [p["run_id"] for p in chain]
    assert len(set(ids)) == len(ids)


def test_single_value_sweep_produces_one_segment():
    """A sweep with one value is still valid."""
    cfg = {**_CFG, "sweep": {"parameter": "omega_b", "values": [1.0]}}
    chain = build_chain(cfg)
    assert len(chain) == 1
    assert chain[0]["n_mix_cycles"] == cfg["n_mix_cycles"]


def test_theta_max_sweep_parameter():
    """Sweeping theta_max_0 sets theta_max[0] per segment."""
    cfg = {
        **_CFG,
        "sweep": {"parameter": "theta_max_0", "values": [3.0, 5.0, 7.0]},
    }
    chain = build_chain(cfg)
    assert len(chain) == 3
    for k, val in enumerate([3.0, 5.0, 7.0]):
        assert math.isclose(chain[k]["theta_max"][0], val), (
            f"segment {k} theta_max[0]={chain[k]['theta_max'][0]}, expected {val}"
        )


# ── _prev motion param propagation ───────────────────────────────────────────

_MOTION_PREV_FIELDS = [
    "theta_max_prev",
    "phi_angular_prev",
    "amplitude_h_prev",
    "phi_horizontal_prev",
    "omega_h_prev",
]


def test_fresh_segment_has_no_prev_motion_params():
    """Segment 0 must not carry _prev motion fields (it has no predecessor)."""
    chain = build_chain(_CFG)
    for field in _MOTION_PREV_FIELDS:
        assert field not in chain[0], f"seg0 should not have {field}"


def test_restart_segments_have_all_prev_motion_params():
    """Every restart segment carries all five _prev motion fields."""
    chain = build_chain(_CFG)
    for k in range(1, len(chain)):
        for field in _MOTION_PREV_FIELDS:
            assert field in chain[k], f"segment {k} missing {field}"


def test_prev_motion_params_match_previous_segment():
    """_prev fields in segment k equal the live motion params of segment k-1."""
    chain = build_chain(_CFG)
    for k in range(1, len(chain)):
        prev = chain[k - 1]
        cur  = chain[k]
        assert cur["theta_max_prev"]      == prev["theta_max"]
        assert cur["phi_angular_prev"]    == prev["phi_angular"]
        assert cur["amplitude_h_prev"]    == prev["amplitude_h"]
        assert cur["phi_horizontal_prev"] == prev["phi_horizontal"]
        assert math.isclose(cur["omega_h_prev"], prev["omega_h"])


def test_prev_motion_params_reflect_swept_theta_max():
    """When theta_max_0 is swept, theta_max_prev[0] in seg k mirrors seg k-1's value."""
    cfg = {
        **_CFG,
        "sweep": {"parameter": "theta_max_0", "values": [3.0, 5.0, 8.0]},
    }
    chain = build_chain(cfg)
    assert math.isclose(chain[1]["theta_max_prev"][0], 3.0)
    assert math.isclose(chain[2]["theta_max_prev"][0], 5.0)


def test_prev_omega_b_still_set_on_restart():
    """omega_b_prev is still populated on restart segments (existing behaviour)."""
    chain = build_chain(_CFG)
    for k in range(1, len(chain)):
        assert "omega_b_prev" in chain[k]


# ── initial_checkpoint support ────────────────────────────────────────────────

_INITIAL_CK = {
    "run_id":   "deadbeef",
    "t_dump":    50.0,
    "omega_b":   1.0,
    "theta_max": [5.0, 0.0, 0.0],
    "phi_angular": [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h": [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
}

_CFG_WITH_CK = {
    **_CFG,
    "n_transition_cycles": 5,
    "initial_checkpoint": _INITIAL_CK,
}


def test_initial_checkpoint_all_segments_are_restarts():
    """With initial_checkpoint every segment (including k=0) must have t_checkpoint."""
    chain = build_chain(_CFG_WITH_CK)
    for k, p in enumerate(chain):
        assert "t_checkpoint" in p, f"segment {k} missing t_checkpoint"


def test_initial_checkpoint_seg0_uses_t_dump():
    """Segment 0 t_checkpoint equals initial_checkpoint.t_dump."""
    chain = build_chain(_CFG_WITH_CK)
    assert math.isclose(chain[0]["t_checkpoint"], _INITIAL_CK["t_dump"])


def test_initial_checkpoint_seg0_prev_fields_from_initial_ck():
    """Segment 0 _prev motion fields come from initial_checkpoint params."""
    chain = build_chain(_CFG_WITH_CK)
    assert chain[0]["omega_b_prev"]    == _INITIAL_CK["omega_b"]
    assert chain[0]["theta_max_prev"]  == _INITIAL_CK["theta_max"]
    assert chain[0]["omega_h_prev"]    == _INITIAL_CK["omega_h"]


def test_initial_checkpoint_all_use_transition_cycles():
    """With initial_checkpoint every segment uses n_transition_cycles, not n_mix_cycles."""
    chain = build_chain(_CFG_WITH_CK)
    for k, p in enumerate(chain):
        assert p["n_mix_cycles"] == _CFG_WITH_CK["n_transition_cycles"], (
            f"segment {k} should use n_transition_cycles"
        )


def test_initial_checkpoint_checkpoint_times_advance():
    """Each segment's t_checkpoint is strictly greater than the previous."""
    chain = build_chain(_CFG_WITH_CK)
    for k in range(1, len(chain)):
        assert chain[k]["t_checkpoint"] > chain[k - 1]["t_checkpoint"], (
            f"segment {k} t_checkpoint should exceed segment {k-1}"
        )
