"""Offset-from-reference vs Delta_theta_max, both resolutions, all from
identical-checkpoint (Test A style) comparisons (diary.md 2026-09-07).
Shows the corrected finding: continuous convergence to zero at both L6
and L7, with L7 steeper -- not a persistent floor.

Usage:
    uv run python scripts/plot_offset_vs_dtheta_convergence.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (d_theta, offset_pct), all identical-checkpoint comparisons (fixed
# source, varying only the requested target) -- see diary.md 2026-09-07
# Test A (L6, source dtheta_src_L6_th4.1) and Test A at L7 (source
# dtheta_src_L7_th6.9; the d_theta=0.1 point reuses the L6-vs-L7 check's
# theta_source=6.9->target=7 run, same source checkpoint).
L6 = [(0.0, 1.1), (0.05, 4.9), (0.1, 16.6)]
L7 = [(0.0, 0.1), (0.02, 4.9), (0.05, 6.7), (0.1, 12.1)]

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(*zip(*L6), marker="o", color="#2166ac", label="L6", linewidth=1.5, markersize=8)
ax.plot(*zip(*L7), marker="D", color="#b2182b", label="L7", linewidth=1.5, markersize=8)

ax.set_xlabel(r"$\Delta\theta_\mathrm{max}$ (deg)")
ax.set_ylabel("Offset from independent reference (%)")
ax.legend(frameon=False)
fig.tight_layout()

out_path = Path(__file__).parent.parent / "experiments" / "dtheta_monotonicity" / "offset_vs_dtheta_convergence.png"
fig.savefig(out_path, dpi=150)
print(f"Saved {out_path}")
