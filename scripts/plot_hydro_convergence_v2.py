"""Grid convergence of tau_mean_t (consistent window, hydro_dataset_v2.json) on L6-L9, theta = 7 deg.

Left:  tau_mean against rpm, one series per level (figstyle: colour = level, hollow = time mean, dotted = ours).
Right: relative change between successive levels, |f(L) - f(L-1)| / f(L), per rpm (grey), median over rpm (colour of
       the finer level, joined by a dotted line). A power law with order p predicts a fall by 2^-p per level.
Run: uv run --project /oscar/data/dharri15/eaguerov/Github/gcbml python scripts/plot_hydro_convergence_v2.py [LMAX]
"""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

KEY, LAB = "tau_mean_t", r"$\tau_{mean}$ (Pa)"
LMIN, LMAX = 6, int(sys.argv[1]) if len(sys.argv) > 1 else 9


def main():
    rows = [r for r in json.load(open(ROOT / "experiments/multifidelity/hydro_dataset_v2.json"))
            if LMIN <= r["level"] <= LMAX]
    plt.rcParams.update(fs.rcparams())
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    table = {}
    for lev in range(LMIN, LMAX + 1):
        pts = sorted((r["rpm"], r[KEY]) for r in rows if r["level"] == lev)
        a = np.array(pts)
        table[lev] = dict(pts)
        ax.plot(a[:, 0], a[:, 1], **fs.series_kw(fs.level_colour(lev), fs.MK["tau"], stat="mean", ours=True),
                label=f"L{lev}")
    ax.set_xlabel("rocking speed (rpm)")
    ax.set_ylabel(LAB)
    ax.grid(**fs.GRID_KW)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    levs, meds = [], []
    for lev in range(LMIN + 1, LMAX + 1):
        common = sorted(set(table[lev]) & set(table[lev - 1]))
        rel = np.array([abs(table[lev][r] - table[lev - 1][r]) / table[lev][r] for r in common])
        ax2.semilogy([lev] * len(rel), rel, "o", color=fs.GUIDE, ms=3)
        levs.append(lev)
        meds.append(np.median(rel))
        print(f"L{lev - 1}->L{lev}: median {np.median(rel):.1%}, min {rel.min():.1%}, max {rel.max():.1%} (n={len(rel)})")
    ax2.plot(levs, meds, ls=":", color=fs.NEUTRAL, lw=1.1)
    for lev, m in zip(levs, meds):
        ax2.plot([lev], [m], **fs.series_kw(fs.level_colour(lev), fs.MK["tau"], stat="mean", ours=True, ls="none"))
    ax2.set_xticks(levs, [f"L{l - 1}→L{l}" for l in levs])
    ax2.set_ylabel("relative change between levels")
    ax2.grid(**fs.GRID_KW)
    fig.tight_layout()
    out = ROOT / f"experiments/multifidelity/hydro_convergence_v2_L{LMIN}_{LMAX}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
