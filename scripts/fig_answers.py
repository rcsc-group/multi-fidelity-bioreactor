"""Answer figures for the 2026-10-21 questions: (1) is kLa / mixing time well defined, (2) does warm starting work, (3) the two ideas.
Data: diary 2026-10-15/17/20 numbers (definition checks), experiments/settling_model_fit.json, test_warmstart_x8.json, h3 settled check
(scripts/check_settled_t10.py), kla_threshold_extrapolation.json, replay loop runs (when present).
Usage: uv run python scripts/fig_answers.py
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
from scripts.check_settled_t10 import POINTS, late  # noqa: E402
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.test_h1_mfbml import hf_run  # noqa: E402

MF = ROOT / "experiments/multifidelity"


def save(fig, name, extra=()):
    for ext in ("png", "pdf"):
        fig.savefig(MF / f"{name}.{ext}", dpi=200, bbox_inches="tight", bbox_extra_artists=list(extra))
    plt.close(fig)


def definition():
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.8), sharey=True)
    b = fs.level_colour(8)
    g = fs.level_colour(9)
    # 1: the same run again
    ax = axes[0]
    rows = [("L9, 17.5 rpm", 9, [-0.0077, 0.0, 0.06]), ("L9, 25 rpm", 9, [-0.0001, 0.0, 0.0]), ("L9, 32.5 rpm", 9, [0.0027, 0.0, 0.0])]
    for i, (lab, L, v) in enumerate(rows):
        for x, y, mk in zip((i - 0.2, i, i + 0.2), v, (fs.MK["cstar_10"], fs.MK["chi_0.50"], fs.MK["chi_0.95"])):
            ax.plot(x, y, ls="none", marker=mk, mfc="w", mec=fs.level_colour(L), ms=7, mew=1.2)
    for i, (r, v) in enumerate(((22.4, 0.08), (29.9, 0.5), (22.6, 1.0))):
        ax.plot(3 + i * 0.4 - 0.4, v, ls="none", marker=fs.MK["chi_0.95"], mfc="w", mec=b, ms=7, mew=1.2)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(["17.5", "25", "32.5", "22.4 / 29.9 / 22.6"], fontsize=8)
    ax.set_xlabel("rocking speed (rpm)")
    ax.set_title("same run, run again", fontsize=10)
    ax.set_ylabel("difference (%)")
    # 2: other code on the same fields
    ax = axes[1]
    pts = [(0, 8, fs.MK["chi_0.95"], 8.3), (0.25, 8, fs.MK["chi_0.95"], -1.1), (1, 9, fs.MK["chi_0.95"], 3.3), (1.25, 9, fs.MK["chi_0.95"], -1.8),
           (2, 8, fs.MK["cstar_25"], -1.0), (3, 8, fs.MK["cstar_50"], 1.6), (4, 8, fs.MK["cstar_10"], 13.6)]
    for x, L, mk, y in pts:
        ax.plot(x, y, ls="none", marker=mk, mfc="w", mec=fs.level_colour(L), ms=7, mew=1.2)
    ax.set_xticks([0.125, 1.125, 2, 3, 4])
    ax.set_xticklabels(["$t_{mix}$ L8", "$t_{mix}$ L9", "kLa 0.25", "kLa 0.50", "kLa 0.10"], fontsize=8)
    ax.set_title("Kim's code on our fields", fontsize=10)
    # 3: other estimator
    ax = axes[2]
    for x, y, mk in ((0, 3, fs.MK["chi_0.95"]), (1, 14, fs.MK["cstar_25"])):
        ax.plot(x, y, ls="none", marker=mk, mfc="w", mec=b, ms=7, mew=1.2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["$t_{mix}$ 0.95", "kLa 0.25"], fontsize=8)
    ax.set_xlim(-0.6, 1.6)
    ax.set_title("other estimator: scatter, new / old", fontsize=10)
    for ax in axes:
        ax.axhline(0, color=fs.GUIDE, lw=1.0)
        ax.grid(**fs.GRID_KW, axis="y")
    axes[0].set_ylim(-5, 16)
    fig.tight_layout()
    save(fig, "fig_answer_definition")


def warmstart():
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.2))
    fit = json.loads((ROOT / "experiments/settling_model_fit.json").read_text())
    ax = axes[0, 0]
    rng = np.random.default_rng(3)
    for e in fit:
        if e.get("confidence") == "flagged_anomaly":
            continue
        kind = 0 if e["d_rpm"] == 0 else 1
        size = abs(e["d_theta"]) if kind == 0 else abs(e["d_rpm"])
        x = kind + rng.uniform(-0.18, 0.18)
        col = fs.level_colour(e["level"])
        if e["settle_cycle"] is None:
            ax.errorbar(x, e["n_cycles_available"], yerr=18, lolims=True, ls="none", color=col, marker="o", mfc="w", ms=6, capsize=3)
        else:
            ax.plot(x, e["settle_cycle"], ls="none", marker="o", mfc="w", mec=col, ms=6, mew=1.2)
    ax.axhline(80, color="0.35", ls=":", lw=1.0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["angle change", "rocking speed change"])
    ax.set_ylabel("cycles to settle")
    ax.set_title("A  change of condition, same level", fontsize=10)
    ax.grid(**fs.GRID_KW, axis="y")
    # B: accuracy against a cold L10 run at 32.5 rpm
    ax = axes[0, 1]
    for k, (lab, v) in enumerate((("from L9", [0.2, 1.0, 0.6]), ("from L8", [2.6, 3.0, 0.1]))):
        for x, y, mk, fill in zip((k - 0.2, k, k + 0.2), v, (fs.MK["tau"], fs.MK["tau"], fs.MK["ediss"]), ("w", fs.level_colour(10), "w")):
            ax.plot(x, y, ls="none", marker=mk, mfc=fill if fill != "w" else "w", mec=fs.level_colour(10), ms=7, mew=1.2)
    ax.axhline(5, color="0.35", ls=":", lw=1.0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["start from L9", "start from L8"])
    ax.set_ylabel("difference to cold L10 (%)")
    ax.set_title("B  change of level, 32.5 rpm", fontsize=10)
    ax.grid(**fs.GRID_KW, axis="y")
    # C: settled check
    ax = axes[1, 0]
    names, vals = [], []
    for r, t in POINTS:
        name = "v_l10_rpm32.5_th7" if (r, t) == (32.5, 7.0) else f"t10_l10_rpm{r:g}_th{t:g}"
        run_tr = hf_run(r, t) if t != 5.5 else f"x8_l10_rpm{r:g}_th5.5"
        if not (ROOT / "runs" / name / "shear_stress.dat").exists() or not (ROOT / "runs" / run_tr / "shear_stress.dat").exists():
            continue
        try:
            a, bq = late(name), qois(run_tr)
        except Exception:
            continue
        vals.append([100 * (a[q] / bq[q] - 1) for q in ("tau_95_t", "tau_mean_t", "edr_mean_t")])
        names.append(f"{r:g}/{t:g}")
    vals = np.array(vals)
    for j, (mk, fill) in enumerate(((fs.MK["tau"], fs.level_colour(10)), (fs.MK["tau"], "w"), (fs.MK["ediss"], "w"))):
        ax.plot(np.arange(len(vals)) + (j - 1) * 0.2, vals[:, j], ls="none", marker=mk, mfc=fill, mec=fs.level_colour(10), ms=6, mew=1.1)
    ax.axhline(0, color=fs.GUIDE, lw=1.0)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=60, fontsize=7)
    ax.set_xlabel("rocking speed (rpm) / angle (deg)")
    ax.set_ylabel("late minus training window (%)")
    ax.set_title("C  is the training flow settled?", fontsize=10)
    ax.grid(**fs.GRID_KW, axis="y")
    # D: cycles needed
    ax = axes[1, 1]
    bars = [("cold", 80), ("from L8", 12), ("from L9", 6)]
    ax.bar(range(3), [b[1] for b in bars], color=fs.level_colour(10), alpha=0.8, width=0.55)
    ax.set_xticks(range(3))
    ax.set_xticklabels([b[0] for b in bars])
    ax.set_ylabel("L10 cycles before release")
    ax.set_title("D  L10 cycles needed", fontsize=10)
    ax.grid(**fs.GRID_KW, axis="y")
    h = [Line2D([], [], ls="none", marker="o", mfc="w", mec=fs.level_colour(L), label=f"L{L}") for L in (6, 7, 8)]
    lg = axes[0, 0].legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    c10 = fs.level_colour(10)
    hq = [Line2D([], [], ls="none", marker="o", mfc=c10, mec=c10, label=r"$\tau_{95}$"), Line2D([], [], ls="none", marker="o", mfc="w", mec=c10, label=r"$\overline{\tau}$"),
          Line2D([], [], ls="none", marker="s", mfc="w", mec=c10, label="EDR")]
    lg2 = axes[1, 1].legend(handles=hq, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    save(fig, "fig_answer_warmstart", [lg, lg2])


def ideas():
    d = json.loads((MF / "kla_threshold_extrapolation.json").read_text())
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    b = fs.level_colour(8)
    for i, q in enumerate(("kLa_1T_25", "kLa_1T_50")):
        r = d[q]["rms"]
        ax.plot([i - 0.12], [r["M1"]], ls="none", marker="o", mfc="w", mec="0.25", ms=8, mew=1.3)
        ax.plot([i + 0.12], [r["M2"]], ls="none", marker=fs.MK["cstar_10"], mfc="w", mec=b, ms=8, mew=1.3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([r"$k_La$, $C^*\!=\!0.25$", r"$k_La$, $C^*\!=\!0.50$"])
    ax.set_ylabel("error on unseen runs (RMS of ln ratio)")
    ax.set_ylim(0, 0.45)
    ax.grid(**fs.GRID_KW, axis="y")
    h = [Line2D([], [], ls="none", marker="o", mfc="w", mec="0.25", label="regression on rpm and angle"),
         Line2D([], [], ls="none", marker=fs.MK["cstar_10"], mfc="w", mec=b, label=r"line fit on $k_La$ at $C^*\!=\!0.10$")]
    lg = ax.legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    save(fig, "fig_answer_threshold", [lg])


if __name__ == "__main__":
    plt.rcParams.update(fs.rcparams())
    definition()
    warmstart()
    ideas()


def periodic():
    d = json.loads((MF / "f3b_loop.json").read_text())["deviation_percent"]
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    b = fs.level_colour(8)
    for j, (q, mk) in enumerate((("kLa_1T_10", fs.MK["cstar_10"]), ("dtmix_0.50", fs.MK["chi_0.50"]), ("kLa_1T_25", fs.MK["cstar_25"]))):
        xs = [i + (j - 1) * 0.18 for i, k in enumerate(("1", "2", "4"))]
        ax.plot(xs, [d[k][q] for k in ("1", "2", "4")], ls="none", marker=mk, mfc="w", mec=b, ms=8, mew=1.3)
    ax.axhspan(-5, 5, color=fs.GUIDE, alpha=0.25, lw=0)
    ax.axhline(0, color=fs.GUIDE, lw=1.0)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["1", "2", "4"])
    ax.set_xlabel("recorded periods repeated")
    ax.set_ylabel("repeated flow minus real flow (%)")
    ax.grid(**fs.GRID_KW, axis="y")
    h = [Line2D([], [], ls="none", marker=fs.MK["cstar_10"], mfc="w", mec=b, label=r"$k_La$, $C^*\!=\!0.10$"),
         Line2D([], [], ls="none", marker=fs.MK["chi_0.50"], mfc="w", mec=b, label=r"$t_{\mathrm{mix}}$, $\chi\!=\!0.50$"),
         Line2D([], [], ls="none", marker=fs.MK["cstar_25"], mfc="w", mec=b, label=r"$k_La$, $C^*\!=\!0.25$")]
    lg = ax.legend(handles=h, frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.03, 1.0))
    fig.tight_layout()
    save(fig, "fig_answer_periodic", [lg])


if __name__ == "__main__":
    periodic()
