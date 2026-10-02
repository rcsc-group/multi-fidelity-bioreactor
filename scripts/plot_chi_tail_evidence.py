"""Evidence check: does the chi(t) curve up to 0.50 contain the late decay rate?

Plots ln(1 - chi) against time for each L9 rpm (one panel per rpm).
A single exponential tail is a straight line. Red: the line fitted on
chi in [0.25, 0.50] (what a truncated run has), extended. The horizontal
lines mark chi = 0.75 and 0.95. If the claim is true, the data after 0.50
stay on the red line.
"""
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs                     # noqa: E402
from scripts.pilot_truncated_chi import chi_curve, RPMS  # noqa: E402

plt.rcParams.update(fs.rcparams())
fig, axes = plt.subplots(2, 5, figsize=(14, 5.4), sharey=True)
for ax, r in zip(axes.ravel(), RPMS):
    t, chi = chi_curve(f"fig9_l9_rpm{r:g}")
    y = np.log(np.clip(1 - chi, 1e-6, 1))
    ax.plot(t, y, color="0.35", lw=0.9)
    k50 = int(np.argmax(chi >= 0.50))
    m = (chi >= 0.25) & (np.arange(len(t)) <= k50)
    b, a = np.polyfit(t[m], y[m], 1)
    tg = np.linspace(0, t.max(), 50)
    ax.plot(tg, a + b * tg, color="#d62728", lw=1.2, ls="--")
    ax.plot(t[m], y[m], color="#d62728", lw=2.2)
    for c in (0.75, 0.95):
        ax.axhline(np.log(1 - c), color="0.6", lw=0.7, ls=":")
    ax.set_title(f"{r:g} rpm", fontsize=10, loc="left")
    ax.set_ylim(-4.0, 0.1)
    ax.grid(**fs.GRID_KW)
for ax in axes[-1]:
    ax.set_xlabel("time after release (s)")
for ax in axes[:, 0]:
    ax.set_ylabel(r"$\ln(1-\chi)$")
fig.tight_layout()
out = ROOT / "experiments/multifidelity/chi_tail_evidence_l9.png"
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved", out)
