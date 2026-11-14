"""Kim et al. (2024) Fig. 9(a) (mixing time vs rpm) and Fig. 11 (kLa vs rpm), theta 7, with the looped-replay L8 points
(rc8_rpm*_rep, diary 2026-10-31) overlaid on the L8 DNS they replace (kmix_l8_rpm*). Same encoding as the replicas
(plot_fig9.py, plot_fig11_fig12.py): colour = dataset/level, marker = threshold, solid = Kim, dotted = ours.
Replay: L8 colour, no line, large hollow marker drawn on top, so a DNS point inside it means agreement.
Usage: uv run python scripts/figures/fig_kim_replay_overlay.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402
from scripts import replay_gate as rg  # noqa: E402
from scripts.autoextend import tip_results  # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
KIM = ROOT / "experiments/kimetal2024/csv_raw/mixing_kla_vs_frequency.csv"
OUT = ROOT / "internal/docs/meetings/prereg"
PANELS = [  # (our key, Kim column, marker)
    [("dtmix_0.95", "dtmix_strict_0.95", fs.threshold_marker("chi", 0.95)),
     ("dtmix_0.75", "dtmix_strict_0.75", fs.threshold_marker("chi", 0.75)),
     ("dtmix_0.50", "dtmix_strict_0.5", fs.threshold_marker("chi", 0.50))],
    [("kLa_1T_10", "kLa_exp5pts_10", fs.threshold_marker("cstar", 10)),
     ("kLa_1T_25", "kLa_exp5pts_25", fs.threshold_marker("cstar", 25)),
     ("kLa_1T_50", "kLa_exp5pts_50", fs.threshold_marker("cstar", 50))]]
SERIES = [(8, "kmix_l8_rpm{:g}"), (9, "fig9_l9_rpm{:g}"), (10, "fig11_l10_rpm{:g}")]


def ours(tmpl: str, root: Path = ROOT / "runs") -> pd.DataFrame:
    rows = []
    for r in RPMS:
        d = tip_results(root, tmpl.format(r)) if root == ROOT / "runs" else \
            (json.loads((root / tmpl.format(r) / "results.json").read_text())
             if (root / tmpl.format(r) / "results.json").exists() else None)
        if d:
            rows.append({"rpm": r, **d})
    return pd.DataFrame(rows)


def main() -> None:
    plt.rcParams.update(fs.rcparams())
    kim = pd.read_csv(KIM).sort_values("RPM")
    rep = ours("rc8_rpm{:g}_rep", rg.RUNS)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for a, panel in zip(ax, PANELS):
        for key, kcol, mk in panel:
            a.plot(kim["RPM"], kim[kcol], **fs.series_kw(fs.KIM, mk))
            for lev, tmpl in SERIES:
                d = ours(tmpl).dropna(subset=[key]) if True else None
                if key in d:
                    a.plot(d["rpm"], d[key], **fs.series_kw(fs.level_colour(lev), mk, ours=True))
            d = rep.dropna(subset=[key])
            a.plot(d["rpm"], d[key], ls="none", marker=mk, ms=11, mfc="none", mec=fs.level_colour(8), mew=1.4, zorder=5)
        a.set_xlabel("rocking speed (rpm)")
        a.set_yscale("log")
        a.grid(**fs.GRID_KW)
    ax[0].set_ylabel("mixing time (s)")
    ax[1].set_ylabel(r"$k_La$ (1/h)")
    lc = fs.level_colour
    handles = ([Line2D([], [], color=fs.KIM, lw=fs.LW_KIM, label="Kim et al.")]
               + [Line2D([], [], color=lc(lv), ls=":", lw=fs.LW_OURS, label=f"L{lv}") for lv, _ in SERIES]
               + [Line2D([], [], ls="none", marker="o", ms=11, mfc="none", mec=lc(8), label="L8 replay")]
               + [Line2D([], [], ls="none", marker=mk, color="0.3", label=lab) for (_, _, mk), lab in
                  zip(PANELS[0], (r"$\chi=0.95$", r"$\chi=0.75$", r"$\chi=0.50$"))]
               + [Line2D([], [], ls="none", marker=mk, color="0.3", label=lab) for (_, _, mk), lab in
                  zip(PANELS[1], (r"$C^*=10\%$", r"$C^*=25\%$", r"$C^*=50\%$"))])
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"kim_replay_overlay.{ext}", bbox_inches="tight", dpi=150)
    print(OUT / "kim_replay_overlay.png")


if __name__ == "__main__":
    main()
