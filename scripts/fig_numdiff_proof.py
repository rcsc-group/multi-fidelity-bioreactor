"""Figure: mixing and oxygen transfer in our solver are carried by numerical, not physical, diffusion (L6-L10).

(a) Tracer: D_eff/D_phys from the variance budget (tests V and V-L10; median over chi 0.1-0.4, 32.5 rpm) vs dx.
    Physically controlled decay: D_eff/D_phys = 1.
(b) Oxygen: measured kL = kLa_1T_25 / a vs grid spacing, against the most a resolved diffusive layer can carry,
    2 D/dx (whole saturation jump across one cell). Data: diag_kl_bound.json (32.5 rpm, physical D).
Legend: the three decoder blocks of the figstyle grammar (dataset = colour, marker = quantity, line = guide), as in
plot_fig11_fig12.py. A level without data in a panel is absent from that panel.
Usage: uv run python scripts/fig_numdiff_proof.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

L_BIO, D_OXY = 0.25, 1.9e-9
V_JSON = "experiments/multifidelity/test_v_l10.json"  # all levels binned at equal chi
LS_A, LS_B = "--", "-."                               # the two guides mean different things, so differ in dash


def main():
    v = json.loads((ROOT / V_JSON).read_text())["levels"]
    kb = json.loads((ROOT / "experiments/multifidelity/diag_kl_bound.json").read_text())
    plt.rcParams.update(fs.rcparams())
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    for L in sorted(int(k) for k in v):
        a1.plot([L_BIO / 2 ** L * 1e3], [v[str(L)]["window"]],
                **fs.series_kw(fs.level_colour(L), fs.MK["deff"], stat="mean", ours=True))
    a1.axhline(1.0, color=fs.NEUTRAL, ls=LS_A, lw=1.0)
    a1.set_ylabel(r"$D_{\mathrm{eff}}/D_{\mathrm{phys}}$")
    a1.set_ylim(0.5, 3e4)
    dxs = np.array([r["dx"] for r in kb])
    for r in kb:
        a2.plot([r["dx"] * 1e3], [r["kL"]],
                **fs.series_kw(fs.level_colour(r["level"]), fs.MK["cstar_25"], stat="mean", ours=True))
    xs = np.geomspace(dxs.min() / 4, dxs.max() * 1.3, 50)
    a2.plot(xs * 1e3, 2 * D_OXY / xs, ls=LS_B, color=fs.NEUTRAL, lw=1.0)
    a2.set_ylabel(r"$k_L$ (m/s)")
    for ax in (a1, a2):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\Delta x$ (mm)")
        ax.grid(**fs.GRID_KW)
    a1.set_xlim(a2.get_xlim())

    levels = sorted({int(k) for k in v} | {int(r["level"]) for r in kb})
    datasets = [Line2D([], [], color=fs.level_colour(L), lw=fs.LW_OURS, ls=":", label=f"L{L}") for L in levels]
    marks = [Line2D([], [], color="0.35", marker=fs.MK["deff"], mfc="w", ls="none", ms=6,
                    label=r"$D_{\mathrm{eff}}/D_{\mathrm{phys}}$"),
             Line2D([], [], color="0.35", marker=fs.MK["cstar_25"], mfc="w", ls="none", ms=6,
                    label=r"$k_L$, $C^*=25\%$")]
    guides = [Line2D([], [], color=fs.NEUTRAL, ls=LS_A, lw=1.0, label=r"$D_{\mathrm{eff}}=D_{\mathrm{phys}}$"),
              Line2D([], [], color=fs.NEUTRAL, ls=LS_B, lw=1.0, label=r"$2D/\Delta x$")]
    blocks = []
    for handles, title, y in ((datasets, "dataset", 1.0), (marks, "marker", 0.50), (guides, "line", 0.22)):
        lg = a2.legend(handles=handles, fontsize=8, frameon=False, loc="upper left",
                       bbox_to_anchor=(1.02, y), title=title)
        lg.get_title().set_fontsize(8)
        blocks.append(lg)
    for lg in blocks[:-1]:
        a2.add_artist(lg)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = ROOT / f"experiments/multifidelity/fig_numdiff_proof.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print("wrote", out)


if __name__ == "__main__":
    main()
