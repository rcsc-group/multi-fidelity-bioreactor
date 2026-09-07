"""How many cycles did condition 2 (20rpm, checkpointed from 17.5rpm)
actually take to stabilize after the t=24.9 handoff? Measured directly
from the data, not estimated by analogy (diary.md 2026-09-06).

Uses tau_mean (bulk, spatially-averaged), not tau_100 (pointwise), as
the settling indicator -- 2026-09-06's rank-invariance check found
tau_100 carries ~5% MPI-rank-count noise even deep in QSS (up to 16% at
individual timesteps), while tau_mean stays under ~0.2%. Using the
noisier statistic risked reading numerical noise as "still settling."

Looks at cycle-to-cycle peak tau_mean (per-cycle max) and asks when it
stops drifting monotonically and settles into a stable oscillation
(std/mean of the last few cycles' peaks below a tolerance).

Usage:
    uv run python scripts/measure_checkpoint_settling.py
"""
import math
from pathlib import Path

import numpy as np

RUN_DIR = Path(__file__).parent.parent / "runs" / "dcd481ed"
T_CHECKPOINT = 24.9

omega_b = 20.0 * 2 * math.pi / 60.0
L_bio, H_bio, th = 0.25, 2 * 0.03575, math.radians(7.0)
T_per = 2 * math.pi / omega_b
V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(th))
U_bio = V_bio / (H_bio * 0.5) / T_per
T_bio = L_bio / U_bio
T_per_nd = T_per / T_bio

d = np.loadtxt(RUN_DIR / "shear_stress.dat", skiprows=1)
t = d[:, 1]
tau_mean = d[:, 5]

cycles_since_ckpt = (t - T_CHECKPOINT) / T_per_nd
n_cycles = int(cycles_since_ckpt.max())
print(f"T_per_nd={T_per_nd:.4f}, {n_cycles} full cycles since checkpoint")

peak_per_cycle = []
for c in range(n_cycles):
    mask = (cycles_since_ckpt >= c) & (cycles_since_ckpt < c + 1)
    if mask.sum() > 0:
        peak_per_cycle.append(tau_mean[mask].max())

peak_per_cycle = np.array(peak_per_cycle)
print("cycle : peak tau_mean")
for c, p in enumerate(peak_per_cycle):
    print(f"{c:5d} : {p:.6f}")

# Settled = last-5-cycle peaks within 5% of their own mean, find first cycle
# where this holds for the remainder of the run.
tol = 0.05
for start in range(len(peak_per_cycle) - 5):
    tail = peak_per_cycle[start:]
    if tail.std() / tail.mean() < tol:
        print(f"\nSettles (tail std/mean < {tol:.0%}) starting at cycle {start} since checkpoint")
        break
else:
    print("\nNever settles within tolerance over the recorded run")
