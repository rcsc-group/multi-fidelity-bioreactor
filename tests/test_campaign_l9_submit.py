"""Campaign 2D (lead decision 2026-11-05): L9 warm submission. MPI jobs run on partition batch only (no batch2:
UCX wireup failure on gpu4001). L9 params = the completed c2d L8 source at the same (rpm, theta), fidelity 9,
same binary, 12 warm-up cycles, release, 40 post cycles (plan 1.3). Synthetic data only.
"""
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))
sys.path.insert(0, str(ROOT))

import campaign as c  # noqa: E402


def test_mpi_partition_batch_is_allowed():
    c.check_mpi_partition("batch")


def test_mpi_partition_batch2_is_refused():
    with pytest.raises(ValueError, match="batch2"):
        c.check_mpi_partition("batch2")


def test_l9_warm_params_are_warm_fidelity9_same_binary():
    src = {"fidelity": 8, "omega_b": 2 * math.pi * 20 / 60, "theta_max": [4.0, 0.0, 0.0], "_binary": "/b/2ef5b07",
           "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5, "n_harmonics": 1}
    p = c.l9_warm_params(src, "c2d_l8_rpm20.00_th4.00", t_ck=140.0, rpm=20.0, theta=4.0, binary="/b/2ef5b07")
    assert p["fidelity"] == 9
    assert p["t_checkpoint"] == 140.0
    assert p["n_mix_cycles"] == 12
    assert p["_binary"] == "/b/2ef5b07"
    assert p["_source_run"] == "c2d_l8_rpm20.00_th4.00"
    assert p["run_id"].startswith("c2d_l9_rpm20.00_th4.00")
    assert p["t_end"] > 0


def test_l9_warm_params_continue_the_source_forcing():
    """Pitfall restart-ramp (diary 2026-11-13): without *_prev the solver ramps the rocking again from zero after
    the restart. The W9 validation runs had *_prev; the first 8 campaign L9 runs did not."""
    src = {"fidelity": 8, "omega_b": 2 * math.pi * 20 / 60, "theta_max": [4.0, 0.0, 0.0], "_binary": "/b/2ef5b07",
           "geometry": {"a": 0.25, "b": 0.03575, "n": 8.0}, "fill_level": 0.5, "n_harmonics": 1}
    p = c.l9_warm_params(src, "c2d_l8_rpm20.00_th4.00", t_ck=140.0, rpm=20.0, theta=4.0, binary="/b/2ef5b07")
    assert p["omega_b_prev"] == src["omega_b"]
    assert p["theta_max_prev"] == src["theta_max"]
    for k in ("phi_angular_prev", "amplitude_h_prev", "phi_horizontal_prev"):
        assert p[k] == [0.0, 0.0, 0.0], k
    assert p["omega_h_prev"] == 0.0
