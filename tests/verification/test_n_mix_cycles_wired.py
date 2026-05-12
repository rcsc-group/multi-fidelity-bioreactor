"""Verify that n_mix_cycles is wired to t_mix (the oxygen/tracer start time).

n_mix_cycles is hardcoded to 80 in the upstream solver. If not wired from
params.json, setting it to a small value has no effect: the oxygen event
still fires at t = 80 * T_per_st ≈ 48.65 non-dim, which is past any
reasonable short test run. tr_oxy.dat would be all zeros.

After the fix: the oxygen event fires at t = n_mix_cycles * T_per_st,
so a run with n_mix_cycles=3 and t_end=5.0 will produce non-zero oxygen data.
"""
import math
import sys
import pathlib
import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).parents[2]))
from tests.conftest import CANONICAL_PARAMS, run_bioreactor


def _t_per_st(params: dict) -> float:
    """Non-dim rocking period for given params."""
    omega_b = params["omega_b"]
    L = params["geometry"]["a"]
    H = params["geometry"]["b"]
    th = math.radians(params["theta_max"][0])
    T_per = 2 * math.pi / omega_b
    V = L / 4 * (H + 0.5 * L * math.tan(th))
    U = V / (H * 0.5) / T_per
    T_bio = L / U
    return T_per / T_bio


def _first_nonzero_oxy_time(run_dir: pathlib.Path) -> float:
    """Return non-dim time of first non-zero oxy_liq_sum row, or inf."""
    path = run_dir / "tr_oxy.dat"
    for line in path.read_text().splitlines():
        parts = line.split()
        if not parts or parts[0] == "i":
            continue
        t, oxy = float(parts[1]), float(parts[2])
        if oxy > 0:
            return t
    return float("inf")


@pytest.mark.medium
def test_n_mix_cycles_controls_oxygen_start(tmp_path):
    """Oxygen transfer must begin at n_mix_cycles * T_per_st, not at 80 * T_per_st.

    Before the fix: n_mix_cycles=3 is ignored; t_mix = 80*T_per_st ≈ 48.65.
    With t_end=5.0, no oxygen data appears (first_oxy_time = inf).
    After the fix: t_mix = 3*T_per_st ≈ 1.82; oxygen data appears before t=5.0.
    """
    n_mix = 3
    params = {
        **CANONICAL_PARAMS,
        "run_id": "nmix_wire",
        "n_mix_cycles": n_mix,
        "t_end": 5.0,
    }
    T_per_st = _t_per_st(params)
    expected_t_mix = n_mix * T_per_st   # ≈ 1.82 non-dim

    run_dir = run_bioreactor(params, tmp_path, timeout=60)
    first_oxy = _first_nonzero_oxy_time(run_dir)

    assert first_oxy < float("inf"), (
        "No oxygen data in tr_oxy.dat — n_mix_cycles is likely not wired; "
        f"t_mix is still ~{80 * T_per_st:.1f} (80 cycles) instead of "
        f"{expected_t_mix:.2f} ({n_mix} cycles)"
    )
    assert first_oxy < params["t_end"], (
        f"First oxygen at t={first_oxy:.2f} is past t_end={params['t_end']}"
    )
