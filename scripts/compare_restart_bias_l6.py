"""Fresh cold start vs. 4-segment same-condition checkpoint chain, same
condition (32.5rpm/theta=7deg, L6), compared over their common real
elapsed time window rather than assumed-matching cycle counts (diary.md
2026-09-03) -- AND restricted to genuine QSS (t > ramp end), not the
postprocess.py-hardcoded 3-cycle cutoff, which for the ramp-matched
binary used here is wrong by 8+ non-dim time (see
_ramp_end_nd/postprocess.py, diary.md 2026-09-03). Both fixes matter:
without the QSS restriction, the comparison window includes each arm's
own ramp-up transient, which is exactly the effect under test.

Usage:
    uv run python scripts/compare_restart_bias_l6.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from scripts.postprocess import _ramp_end_nd, _load_params

RUNS = Path(__file__).parent.parent / "runs"
CHAIN_SEGMENTS = ["459e1b76", "d74d35b7", "a3070b72", "7beefe12"]


def load(run_id):
    return np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)


fresh = load("restart_bias_fresh_l6")
chain = np.vstack([load(seg) for seg in CHAIN_SEGMENTS])

t_ramp = _ramp_end_nd(_load_params(RUNS / "restart_bias_fresh_l6"))
t_common = chain[:, 1].max()
print(f"ramp ends at t={t_ramp:.3f}; chain covers t=0..{t_common:.3f} "
      f"({len(CHAIN_SEGMENTS)} segments); QSS comparison window: t={t_ramp:.3f}..{t_common:.3f}")

fresh_qss = fresh[(fresh[:, 1] > t_ramp) & (fresh[:, 1] <= t_common)]
chain_qss = chain[chain[:, 1] > t_ramp]
print(f"QSS samples: fresh={len(fresh_qss)}  chain={len(chain_qss)}")

# columns: i t tau_95 tau_98 tau_100 tau_mean tau_100_strict tau_mean_strict tau_100_signed ediss_mean tau_mean_signed
for name, col in [("tau_100", 4), ("tau_mean", 5)]:
    f_max = fresh_qss[:, col].max()
    c_max = chain_qss[:, col].max()
    rel = 100 * (c_max - f_max) / f_max
    print(f"{name}_max: fresh={f_max:.6f}  chain={c_max:.6f}  rel_diff={rel:+.2f}%")

i_f = fresh_qss[:, 4].argmax()
i_c = chain_qss[:, 4].argmax()
print(f"tau_100 peak time: fresh t={fresh_qss[i_f, 1]:.3f}  chain t={chain_qss[i_c, 1]:.3f}")
