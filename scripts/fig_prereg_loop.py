"""Looped hydrodynamics-only workflows against the full DNS, ours vs Kim (diary 2026-10-29).
20 rpm, theta 7, grid 256^2 in both solvers, tracer released at cycle 20, recorded cycle 1 looped 3 times.
Curves per solver: DNS; stored projected flux (every step, looped); flux rebuilt from u (every step, looped).
Reads docs/meetings/prereg/ours_L8_series.csv (prereg_offline.py report 8) and
prereg_out/pre_A_N256/series.csv (compare_kim.py). Usage: uv run python scripts/fig_prereg_loop.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

OURS = ROOT / "docs/meetings/prereg/ours_L8_series.csv"
KIM = Path("/oscar/scratch/eaguerov/kim_pipeline/prereg_out/pre_A_N256/series.csv")
OUT = ROOT / "docs/meetings/prereg"

# workflow -> (line width, style); drawn in this order so the thick lines sit under the thin ones
WORK = [("DNS", 4.0, "-"), ("stored flux", 2.2, "--"), ("flux from u", 1.1, "-")]
COLS = {"ours": {"t": "t", "DNS": "A", "stored flux": "repW", "flux from u": "repWS"},
        "Kim": {"t": "t_cyc", "DNS": "A", "stored flux": "off_m4", "flux from u": "his_loop"}}
COLOUR = {"ours": fs.level_colour(8), "Kim": fs.KIM}


def main() -> None:
    plt.rcParams.update(fs.rcparams())
    data = {"ours": np.genfromtxt(OURS, delimiter=",", names=True),
            "Kim": np.genfromtxt(KIM, delimiter=",", names=True)}
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 7), sharex=True)
    for wf, lw, ls in WORK:
        for solver, d in data.items():
            c = COLS[solver]
            if f"{c[wf]}_m" not in d.dtype.names:
                print(f"missing: {solver} {wf}")
                continue
            t = d[c["t"]]
            ref = d[f"{c['DNS']}_oxy"][(t >= 2) & (t <= 3)]
            ref = np.nanmean(ref)                              # own DNS, cycle 3: the two solvers use other oxygen units
            kw = dict(color=COLOUR[solver], lw=lw, ls=ls, alpha=0.45 if wf == "DNS" else 1.0)
            ax[0].plot(t, 100 * (d[f"{c[wf]}_m"] - 1), label=f"{solver}, {wf}", **kw)
            ax[1].plot(t, d[f"{c[wf]}_chi"], **kw)
            ax[2].plot(t, d[f"{c[wf]}_oxy"] / ref, **kw)
    ax[0].set_ylabel("tracer mass change (%)")
    ax[1].set_ylabel(r"mixing $\chi$")
    ax[2].set_ylabel("liquid oxygen / DNS, cycle 3")
    ax[2].set_xlabel("time from release (cycles)")
    for a in ax:
        a.grid(**fs.GRID_KW)
    ax[2].set_xlim(0, 3)
    fig.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"loop_ours_vs_kim.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "loop_ours_vs_kim.png")


if __name__ == "__main__":
    main()
