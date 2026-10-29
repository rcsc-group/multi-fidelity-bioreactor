"""Do the two DNS solve the same flow? Ours L8 vs Kim N=256, 20 rpm, theta 7, the 3 cycles after the release
(diary 2026-10-29). Free-surface extremes (vol_frac_interf.dat) and liquid rms velocities (normf.dat).
Usage: uv run python scripts/fig_dns_flow_compare.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

T = 0.607329          # T_per_st, same formula in both codes (20 rpm, theta 7, L_bio 0.25)
T_REL = 20 * T
RUNS = {"ours": (Path("/oscar/scratch/eaguerov/replay_gate/runs/pr8_A"), fs.level_colour(8), 3.5),
        "Kim": (Path("/oscar/scratch/eaguerov/kim_pipeline/prereg_out/pre_A_N256/producer"), fs.KIM, 1.2)}
PANELS = [("vol_frac_interf.dat", "posY_max", "highest free surface"),
          ("vol_frac_interf.dat", "posY_min", "lowest free surface"),
          ("normf.dat", "ux_liq_rms", r"liquid $u_x$ rms"),
          ("normf.dat", "uy_liq_rms", r"liquid $u_y$ rms")]
OUT = ROOT / "docs/meetings/prereg"


def col(path: Path, name: str) -> tuple[np.ndarray, np.ndarray]:
    head = path.read_text().splitlines()[0].split()
    a = np.loadtxt(path, skiprows=1, usecols=(head.index("t"), head.index(name)))
    m = (a[:, 0] >= T_REL - 3 * T) & (a[:, 0] <= T_REL + 3 * T)
    return (a[m, 0] - T_REL) / T, a[m, 1]


def main() -> None:
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(len(PANELS), 1, figsize=(7.5, 8), sharex=True)
    for name, (d, c, lw) in RUNS.items():
        for k, (fname, cname, lab) in enumerate(PANELS):
            t, v = col(d / fname, cname)
            ax[k].plot(t, v, color=c, lw=lw, alpha=0.5 if lw > 2 else 1.0, label=name if k == 0 else None)
            ax[k].set_ylabel(lab)
    for a in ax:
        a.grid(**fs.GRID_KW)
        a.axvline(0, color=fs.GUIDE, lw=0.8)
    ax[-1].set_xlabel("time from release (cycles)")
    fig.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"dns_flow_compare.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "dns_flow_compare.png")


if __name__ == "__main__":
    main()
