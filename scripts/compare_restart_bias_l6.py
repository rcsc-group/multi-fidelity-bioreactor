"""Fresh cold start vs. 4-segment same-condition checkpoint chain, same
condition (32.5rpm/theta=7deg, L6), compared over their common real
elapsed time window rather than assumed-matching cycle counts (diary.md
2026-09-03).

Usage:
    uv run python scripts/compare_restart_bias_l6.py
"""
from pathlib import Path

import numpy as np

RUNS = Path(__file__).parent.parent / "runs"
CHAIN_SEGMENTS = ["459e1b76", "d74d35b7", "a3070b72", "7beefe12"]


def load(run_id):
    return np.loadtxt(RUNS / run_id / "shear_stress.dat", skiprows=1)


fresh = load("restart_bias_fresh_l6")
chain = np.vstack([load(seg) for seg in CHAIN_SEGMENTS])

t_common = chain[:, 1].max()
print(f"chain covers t=0..{t_common:.3f} ({len(CHAIN_SEGMENTS)} segments); "
      f"fresh truncated to the same window for a fair comparison")

fresh_trunc = fresh[fresh[:, 1] <= t_common]

# columns: i t tau_95 tau_98 tau_100 tau_mean tau_100_strict tau_mean_strict tau_100_signed ediss_mean tau_mean_signed
for name, col in [("tau_100", 4), ("tau_mean", 5)]:
    f_max = fresh_trunc[:, col].max()
    c_max = chain[:, col].max()
    rel = 100 * (c_max - f_max) / f_max
    print(f"{name}_max: fresh={f_max:.6f}  chain={c_max:.6f}  rel_diff={rel:+.2f}%")

i_f = fresh_trunc[:, 4].argmax()
i_c = chain[:, 4].argmax()
print(f"tau_100 peak time: fresh t={fresh_trunc[i_f, 1]:.3f}  chain t={chain[i_c, 1]:.3f}")
