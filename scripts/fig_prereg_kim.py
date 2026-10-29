"""Kim's solver, N=256, 20 rpm, theta 7: full DNS against the hydrodynamics-only replays (diary 2026-10-29).
Reads prereg_out/<tag>/series.csv (compare_kim.py). Usage: uv run python scripts/fig_prereg_kim.py [tag]
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

TAG = sys.argv[1] if len(sys.argv) > 1 else "pre_A_N256"
SRC = Path("/oscar/scratch/eaguerov/kim_pipeline/prereg_out") / TAG / "series.csv"
OUT = ROOT / "docs/meetings/prereg"

# (column prefix, label, colour, style): the label says which face flux the replay uses and how often it changes
RUNS = [("A", "DNS", "black", "-"),
        ("floor_B", "DNS, 8 ranks", "0.6", "-"),
        ("off_m1", "stored flux, every step", "#0072B2", "--"),
        ("off_m3", "stored flux, held 1/100 cycle", "#56B4E9", "--"),
        ("his_S", "flux from u, every step", "#D55E00", "-."),
        ("his_H", "flux from u, held 1/100 cycle", "#E69F00", "-."),
        ("his_ship", "flux from u, held, cycle 1 looped", "#CC79A7", ":")]


def main() -> None:
    plt.rcParams.update(fs.rcparams())
    d = np.genfromtxt(SRC, delimiter=",", names=True)
    t = d["t_cyc"]
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 7), sharex=True)
    for key, lab, c, ls in RUNS:
        if f"{key}_m" not in d.dtype.names:
            continue
        ax[0].plot(t, 100 * (d[f"{key}_m"] - 1), color=c, ls=ls, lw=1.4, label=lab)
        ax[1].plot(t, d[f"{key}_chi"], color=c, ls=ls, lw=1.4)
        ax[2].plot(t, d[f"{key}_oxy"], color=c, ls=ls, lw=1.4)
    ax[0].set_ylabel("tracer mass change (%)")
    ax[1].set_ylabel(r"mixing $\chi$")
    ax[2].set_ylabel("liquid oxygen (code units)")
    ax[2].set_xlabel("time from release (cycles)")
    for a in ax:
        a.grid(**fs.GRID_KW)
    ax[2].set_xlim(0, 3)
    fig.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"kim_{TAG}_series.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / f"kim_{TAG}_series.png")


if __name__ == "__main__":
    main()
