"""Looped hydrodynamics-only workflows against the full DNS, ours vs Kim (diary 2026-10-29).
20 rpm, theta 7, grid 256^2 in both solvers, tracer released at cycle 20, recorded cycle 1 looped 3 times.
Curves per solver: DNS; stored projected flux (every step, looped); flux rebuilt from u (every step, looped).
Every curve is read from prereg_diag.dat, written every step by the SAME code in both solvers (pre_diag),
with the liquid parts of the Henry split (C1, C2, OL), not the f-weighted sums.
Kim's DNS curve is his exact replay (off_m1), equal to his producer to 1e-14 (diary, Kim N=256 result).
Usage: uv run python scripts/fig_prereg_loop.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

T = 0.607329                                  # T_per_st, same formula in both codes
OURS = Path("/oscar/scratch/eaguerov/replay_gate/runs")
KIM = Path("/oscar/scratch/eaguerov/kim_pipeline/prereg_out/pre_A_N256")
OUT = ROOT / "docs/meetings/prereg"

# workflow -> (line width, style); drawn in this order so the thick lines sit under the thin ones
WORK = [("DNS", 4.0, "-"), ("stored flux", 2.2, "--"), ("flux from u", 1.1, "-")]
DIRS = {"ours": {"DNS": OURS / "pr28_A", "stored flux": OURS / "pr28_repW", "flux from u": OURS / "pr28_repWS"},
        "Kim": {"DNS": KIM / "off_m1", "stored flux": KIM / "off_m4", "flux from u": KIM / "his_loop"}}
COLOUR = {"ours": fs.level_colour(8), "Kim": fs.KIM}


def load(d: Path) -> dict[str, np.ndarray]:
    a = np.genfromtxt(d / "prereg_diag.dat", names=True)
    a = a[int(np.argmax(a["C1"] > 0)):]      # from the release step on
    mean = a["B1"] / a["B0"]                  # chi over bulk-liquid cells (f >= 0.5): phi1/f is bounded there
    var = a["B2"] / a["B0"] - mean ** 2
    return {"t": a["t"], "m": a["C1"] / a["C1"][0], "chi": 1 - var / var[0], "oxy": a["OL"] / a["S0"]}


def main() -> None:
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 7), sharex=True)
    t0 = {s: load(DIRS[s]["DNS"])["t"][0] for s in DIRS}
    for wf, lw, ls in WORK:
        for s in DIRS:
            r = load(DIRS[s][wf])
            x = (r["t"] - t0[s]) / T
            kw = dict(color=COLOUR[s], lw=lw, ls=ls, alpha=0.45 if wf == "DNS" else 1.0)
            ax[0].plot(x, 100 * (r["m"] - 1), label=f"{s}, {wf}", **kw)
            ax[1].plot(x, r["chi"], **kw)
            ax[2].plot(x, r["oxy"], **kw)
    ax[0].set_ylabel("tracer mass change (%)")
    ax[1].set_ylabel(r"mixing $\chi$")
    ax[2].set_ylabel("mean liquid oxygen (gas = 1)")
    ax[2].set_xlabel("time from release (cycles)")
    for a in ax:
        a.grid(**fs.GRID_KW)
    ax[2].set_xlim(0, 3)
    fig.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"loop_ours_vs_kim.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "loop_ours_vs_kim.png")


if __name__ == "__main__":
    main()
