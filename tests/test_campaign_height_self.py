"""Height check compares a warm run with ITSELF (2026-11-07). The old check compared the L9 wave with the L8
source wave; near the 35.2 rpm mode the two grids give different waves, so steady runs failed (33.70/5.35: drift
0.1 %, height vs L8 11 %)."""
import numpy as np
from scripts.campaign2d.campaign import last_cycle_height


def _vol(h_of_t, period=1.0, tend=50.0, dt=0.01):
    t = np.arange(0.0, tend, dt)
    v = np.zeros((len(t), 6))
    v[:, 1] = t
    v[:, 4] = h_of_t(t)          # posY_max
    v[:, 5] = 0.0                # posY_min
    return v


def test_last_cycle_height_is_mean_over_last_complete_cycle():
    v = _vol(lambda t: np.where(t >= 48.5, 2.0, 1.0))
    assert abs(last_cycle_height(v, period=1.0) - 2.0) < 1e-9


def test_steady_run_has_zero_self_change():
    v = _vol(lambda t: 0.3 + 0.0 * t)
    assert last_cycle_height(v, period=1.0) == 0.3
