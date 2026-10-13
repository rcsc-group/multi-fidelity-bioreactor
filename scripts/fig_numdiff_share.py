"""Figure: the share of scalar transport that physical diffusion accounts for, per mesh level (L6-L10, 32.5 rpm).

Tracer (exact, from the variance budget d sigma^2/dt = -2 D <|grad c|^2>, tests V / V-L10): physical share = 1/r,
r = measured variance loss / loss the physical D can cause with the gradients on the grid (median, chi 0.1-0.4).
Oxygen (upper bound): measured kL = kLa_1T_25 / a; a resolved diffusive flux cannot exceed kL = 2 D/dx (whole
saturation jump across half a cell), so the physical share of the measured flux is at most (2 D/dx) / kL.
100% = transport fully explained by physical diffusion.
Data: test_v_l10.json, diag_kl_bound.json. Usage: uv run python scripts/fig_numdiff_share.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

MF = ROOT / "experiments/multifidelity"


def shares():
    v = json.loads((MF / "test_v_l10.json").read_text())["levels"]
    kb = json.loads((MF / "diag_kl_bound.json").read_text())
    tracer = {int(L): 100.0 / d["window"] for L, d in v.items()}
    oxygen = {int(r["level"]): 100.0 / r["ratio"] for r in kb}
    return tracer, oxygen


def main():
    tracer, oxygen = shares()
    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    for data, mk in ((tracer, fs.MK["deff"]), (oxygen, fs.MK["cstar_25"])):
        Ls = sorted(data)
        ax.plot(Ls, [data[L] for L in Ls], ls=":", color=fs.NEUTRAL, lw=fs.LW_OURS, zorder=1)
        for L in Ls:
            ax.plot([L], [data[L]], **fs.series_kw(fs.level_colour(L), mk, stat="mean", ours=True))
    ax.axhline(100.0, color=fs.NEUTRAL, ls="--", lw=1.0)
    ax.set_yscale("log")
    ax.set_ylim(3e-3, 300)
    ax.set_xticks(range(6, 11), [f"L{L}" for L in range(6, 11)])
    ax.set_xlabel(r"mesh level ($\Delta x$ = 3.9 mm $\to$ 0.24 mm)")
    ax.set_ylabel("physical share of transport (%)")
    ax.grid(**fs.GRID_KW)

    levels = sorted(set(tracer) | set(oxygen))
    datasets = [Line2D([], [], color=fs.level_colour(L), lw=fs.LW_OURS, ls=":", label=f"L{L}") for L in levels]
    marks = [Line2D([], [], color="0.35", marker=fs.MK["deff"], mfc="w", ls="none", ms=6,
                    label="tracer mixing"),
             Line2D([], [], color="0.35", marker=fs.MK["cstar_25"], mfc="w", ls="none", ms=6,
                    label=r"O$_2$ transfer (upper bound)")]
    guides = [Line2D([], [], color=fs.NEUTRAL, ls="--", lw=1.0, label="all physical")]
    blocks = []
    for handles, title, y in ((datasets, "dataset", 1.0), (marks, "marker", 0.50), (guides, "line", 0.22)):
        lg = ax.legend(handles=handles, fontsize=8, frameon=False, loc="upper left",
                       bbox_to_anchor=(1.02, y), title=title)
        lg.get_title().set_fontsize(8)
        blocks.append(lg)
    for lg in blocks[:-1]:
        ax.add_artist(lg)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = MF / f"fig_numdiff_share.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight", bbox_extra_artists=blocks)
        print("wrote", out)
    for L in levels:
        print(f"L{L}: tracer {tracer.get(L, float('nan')):.3g}%  O2 <= {oxygen.get(L, float('nan')):.3g}%")


if __name__ == "__main__":
    main()
