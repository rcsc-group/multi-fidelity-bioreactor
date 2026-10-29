"""Test of the diagnostic hypothesis (diary 2026-10-29, pre-registered): tracer mass and liquid oxygen of each
DNS, measured with the f-weighted sums (S1, O = sum c f dv) and with the liquid part of the Henry split
(C1, OL = sum phi1 dv). Prediction: the half-cycle spikes and Kim's 0.9 % drift are in S1/O only.
Usage: uv run python scripts/fig_diag_test.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

T = 0.607329
RUNS = {"ours": (Path("/oscar/scratch/eaguerov/replay_gate/runs/pr28_A"), fs.level_colour(8)),
        "Kim": (Path("/oscar/scratch/eaguerov/kim_pipeline/prereg_out/pre_A_N256/off_m1"), fs.KIM)}
OUT = ROOT / "docs/meetings/prereg"


def main() -> None:
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(2, 2, figsize=(10, 5.5), sharex=True)
    for col, (name, (d, c)) in enumerate(RUNS.items()):
        a = np.genfromtxt(d / "prereg_diag.dat", names=True)
        a = a[int(np.argmax(a["C1"] > 0)):]
        a = a[a["t"] <= a["t"][0] + 3 * T + 1e-9]
        x = (a["t"] - a["t"][0]) / T
        ax[0, col].plot(x, 100 * (a["S1"] / a["S1"][0] - 1), color=c, lw=0.8, ls=":")
        ax[0, col].plot(x, 100 * (a["C1"] / a["C1"][0] - 1), color=c, lw=1.6)
        ax[1, col].plot(x, a["O"] / a["S0"], color=c, lw=0.8, ls=":")
        ax[1, col].plot(x, a["OL"] / a["S0"], color=c, lw=1.6)
        ax[0, col].set_title(name)
        ax[1, col].set_xlabel("time from release (cycles)")
    fig.legend(handles=[plt.Line2D([], [], color=fs.NEUTRAL, lw=0.8, ls=":", label=r"$\sum c f$"),
                        plt.Line2D([], [], color=fs.NEUTRAL, lw=1.6, label=r"$\sum \phi_1$")],
               loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    ax[0, 0].set_ylabel("tracer mass change (%)")
    ax[1, 0].set_ylabel("mean liquid oxygen (gas = 1)")
    for a in ax.flat:
        a.grid(**fs.GRID_KW)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"diag_test.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "diag_test.png")


if __name__ == "__main__":
    main()
