"""Figure: mixing and oxygen transfer in our solver are carried by numerical, not physical, diffusion (L6-L10).

(a) Tracer: D_eff/D_phys from the variance budget (test V; median over chi 0.5-0.95, 32.5 rpm) vs grid spacing.
    Physically controlled decay: D_eff/D_phys = 1.
(b) Oxygen: measured kL = kLa_1T_25 / a vs grid spacing, against the most a resolved diffusive layer can carry,
    2 D/dx (whole saturation jump across one cell). Data: diag_kl_bound.json (32.5 rpm, physical D).
Usage: uv run python scripts/fig_numdiff_proof.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

L_BIO = 0.25
V_RUNS = {8: "scal_l8_rpm32.5_th7", 9: "scal_l9_rpm32.5_th7", 10: "v_l10_rpm32.5_th7"}


def main():
    v = json.loads((ROOT / "experiments/multifidelity/test_v_variance_budget.json").read_text())["runs"]
    kb = json.loads((ROOT / "experiments/multifidelity/diag_kl_bound.json").read_text())
    plt.rcParams.update(fs.rcparams())
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    for L, run in V_RUNS.items():
        if run in v:
            a1.plot([L_BIO / 2 ** L * 1e3], [v[run]["r_median"]], label=f"L{L}",
                    **fs.series_kw(fs.level_colour(L), fs.MK["chi_0.95"], stat="mean", ours=True))
    a1.axhline(1.0, color=fs.NEUTRAL, ls="--", lw=1.0)
    a1.set_ylabel(r"$D_{\mathrm{eff}}/D_{\mathrm{phys}}$")
    a1.set_ylim(0.5, 3e3)
    dxs = np.array([r["dx"] for r in kb])
    for r in kb:
        a2.plot([r["dx"] * 1e3], [r["kL"]], label=f"L{r['level']}",
                **fs.series_kw(fs.level_colour(r["level"]), fs.MK["cstar_25"], stat="mean", ours=True))
    xs = np.geomspace(dxs.min() / 4, dxs.max() * 1.3, 50)
    a2.plot(xs * 1e3, 2 * 1.9e-9 / xs, ls="--", color=fs.NEUTRAL, lw=1.0, label=r"$2D/\Delta x$")
    a2.set_ylabel(r"$k_L$ (m/s)")
    for ax in (a1, a2):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\Delta x$ (mm)")
        ax.grid(**fs.GRID_KW)
    a1.set_xlim(a2.get_xlim())
    h, lab = a2.get_legend_handles_labels()
    fig.legend(h, lab, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        out = ROOT / f"experiments/multifidelity/fig_numdiff_proof.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print("wrote", out)


if __name__ == "__main__":
    main()
