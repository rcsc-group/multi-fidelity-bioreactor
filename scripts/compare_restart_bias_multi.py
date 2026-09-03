"""Does the same-condition restart-chain bias result (diary.md 2026-09-03
(2): tau_mean_max +0.1%, tau_100_max -4.0% at 32.5rpm) hold across other
conditions? Compares fresh vs. chained arms at 17.5/25/32.5/37.5rpm,
each restricted to its own genuine QSS window (t > ramp end -- see
_ramp_end_nd/postprocess.py), same methodology as
compare_restart_bias_l6.py.

Usage:
    uv run python scripts/compare_restart_bias_multi.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.postprocess import _ramp_end_nd, _load_params

RUNS = Path(__file__).parent.parent / "runs"

CONDITIONS = {
    "17.5": {"fresh": "restart_bias_fresh_l6_rpm17.5",
             "chain": ["acbb1185", "a0ef9cfc", "14224e98", "941739f8"]},
    "25": {"fresh": "restart_bias_fresh_l6_rpm25",
           "chain": ["3275f001", "5aa1be08", "cb3e58a5", "71f23c2f"]},
    "32.5": {"fresh": "restart_bias_fresh_l6",
             "chain": ["459e1b76", "d74d35b7", "a3070b72", "7beefe12"]},
    "37.5": {"fresh": "restart_bias_fresh_l6_rpm37.5",
             "chain": ["74cdf0f8", "95869107", "d25d35ed", "3097f907"]},
}


def load(run_id):
    return np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)


print(f"{'rpm':>6} {'t_ramp':>8} {'t_common':>9} {'n_qss':>7}  "
      f"{'tau_mean fresh':>15} {'chain':>10} {'reldiff':>9}  "
      f"{'tau_100 fresh':>14} {'chain':>10} {'reldiff':>9}  {'peak_t f/c':>12}")

for rpm, ids in CONDITIONS.items():
    fresh = load(ids["fresh"])
    chain = np.vstack([load(seg) for seg in ids["chain"]])

    t_ramp = _ramp_end_nd(_load_params(RUNS / ids["fresh"]))
    t_common = chain[:, 1].max()

    fresh_qss = fresh[(fresh[:, 1] > t_ramp) & (fresh[:, 1] <= t_common)]
    chain_qss = chain[chain[:, 1] > t_ramp]

    row = [f"{rpm:>6}", f"{t_ramp:8.3f}", f"{t_common:9.3f}", f"{len(chain_qss):7d}"]
    peaks = []
    for col in (5, 4):  # tau_mean, tau_100
        f_max = fresh_qss[:, col].max()
        c_max = chain_qss[:, col].max()
        rel = 100 * (c_max - f_max) / f_max
        row += [f"{f_max:15.6f}" if col == 5 else f"{f_max:14.6f}",
                f"{c_max:10.6f}", f"{rel:+9.2f}%"]
        if col == 4:
            i_f = fresh_qss[:, 4].argmax()
            i_c = chain_qss[:, 4].argmax()
            peaks = [fresh_qss[i_f, 1], chain_qss[i_c, 1]]
    row.append(f"{peaks[0]:.2f}/{peaks[1]:.2f}")
    print(" ".join(row))
