"""Grid convergence of tau_mean_max and tau_95_qss on L6-L10 at theta = 7 deg (data: scripts/hydro_dataset.py).

Left column: QoI against rpm, one curve per level. Right column: relative change between successive levels,
|f(L) - f(L-1)| / |f(L)|, per rpm (median over rpm as a thick line); convergence shows as a falling curve.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

QOIS = [("tau_mean_max", r"$\tau_{mean}$ (Pa)"), ("tau_95_qss", r"$\tau_{95}$ (Pa)")]
COLS = {6: "#c6dbef", 7: "#9ecae1", 8: "#4292c6", 9: "#2171b5", 10: "#08306b"}


def main():
    data = json.load(open(ROOT / "experiments/multifidelity/hydro_dataset.json"))
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    for row, (key, lab) in enumerate(QOIS):
        ax = axes[row, 0]
        table = {}
        for lev in sorted(COLS):
            pts = sorted((d["rpm"], d[key]) for d in data if d["level"] == lev)
            if not pts:
                continue
            rpm, val = np.array(pts).T
            table[lev] = dict(zip(rpm, val))
            ax.plot(rpm, val, "o-", color=COLS[lev], lw=1.5, ms=4, label=f"L{lev}")
        ax.set_xlabel("rocking speed (rpm)")
        ax.set_ylabel(lab)
        ax.grid(**fs.GRID_KW)
        ax2 = axes[row, 1]
        meds = []
        for lev in range(7, 11):
            common = sorted(set(table.get(lev, {})) & set(table.get(lev - 1, {})))
            rel = [abs(table[lev][r] - table[lev - 1][r]) / abs(table[lev][r]) for r in common]
            ax2.semilogy([lev] * len(rel), rel, "o", color="0.6", ms=3)
            meds.append(np.median(rel))
            print(f"{key}: L{lev - 1}->L{lev} median rel change {np.median(rel):.3f} (n={len(rel)})")
        ax2.semilogy(range(7, 11), meds, "k-", lw=2, label="median over rpm")
        ax2.set_xticks(range(7, 11), [f"L{l - 1}→L{l}" for l in range(7, 11)])
        ax2.set_ylabel("relative change between levels")
        ax2.grid(**fs.GRID_KW)
    axes[0, 0].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    axes[0, 1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    out = ROOT / "experiments/multifidelity/hydro_convergence.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
