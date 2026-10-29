"""Kim's workflow as shipped (flux rebuilt from u, flow held 1/100 cycle, cycle 1 looped) against his DNS, and
the two steps that separate them, with the corrected sums (diary 2026-10-29). Same three panels as
fig_prereg_loop.py. Usage: uv run python scripts/fig_kim_shipped.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts.fig_prereg_loop import KIM, T, load  # noqa: E402

OUT = ROOT / "docs/meetings/prereg"
# thick first
RUNS = [(KIM / "off_m1", "DNS", "0.55", 4.0, "-"),
        (KIM / "his_loop", "flux from u, looped", fs.KIM, 2.0, "--"),
        (KIM / "his_ship", "flux from u, held, looped", fs.KIM, 1.1, "-")]


def main() -> None:
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(3, 1, figsize=(7.5, 7), sharex=True)
    t0 = load(KIM / "off_m1")["t"][0]
    for d, lab, c, lw, ls in RUNS:
        r = load(d)
        x = (r["t"] - t0) / T
        ax[0].plot(x, 100 * (r["m"] - 1), color=c, lw=lw, ls=ls, label=lab)
        ax[1].plot(x, r["chi"], color=c, lw=lw, ls=ls)
        ax[2].plot(x, r["oxy"], color=c, lw=lw, ls=ls)
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
        fig.savefig(OUT / f"kim_shipped.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "kim_shipped.png")


if __name__ == "__main__":
    main()
