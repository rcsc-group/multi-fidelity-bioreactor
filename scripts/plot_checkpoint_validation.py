"""Generate checkpoint-restart validation figure.

Two-panel method-comparison plot:
  Left  — identity scatter: baseline kLa (x) vs checkpointed kLa (y)
  Right — residual scatter: baseline kLa (x) vs |Δ|/baseline (%)

Outputs: experiments/figures/checkpoint_validation.pdf
         experiments/figures/checkpoint_validation.png
"""
from __future__ import annotations
import json, math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import numpy as np

_PROJECT_ROOT = Path(__file__).parents[1]
_RUNS_ROOT    = _PROJECT_ROOT / "runs"
_FIG_DIR      = _PROJECT_ROOT / "experiments" / "figures"

# ── colours ────────────────────────────────────────────────────────────────
C_BG       = "#EBF0F5"
C_INK      = "#14243A"
C_ACCENT   = "#1C62B7"   # instrument blue — fill sweep
C_AMBER    = "#C47F17"   # amber — theta sweep
C_MUTED    = "#8BA4BC"   # gridlines / reference lines
C_BAND     = "#D4E3F5"   # ±0.5 % tolerance band fill


def _load_sweep(exp_substr, key_fn, label_fn):
    records = {}
    for p in sorted(_RUNS_ROOT.glob("*/params.json"),
                    key=lambda f: f.stat().st_mtime):
        try:
            d = json.loads(p.read_text())
            if exp_substr not in d.get("_experiment_dir", ""):
                continue
            res = p.parent / "results.json"
            if not res.exists():
                continue
            r = json.loads(res.read_text())
            kla = r.get("kLa_25")
            if kla is None or not math.isfinite(float(kla)):
                continue
            key = key_fn(d)
            records[key] = {"kla": float(kla), "label": label_fn(d)}
        except Exception:
            pass
    return records


def _collect():
    om = lambda d: round(float(d.get("omega_b", 0)), 4)
    fl = lambda d: round(float(d.get("fill_level", 0)), 2)
    th = lambda d: round(float(
        (d.get("theta_max") or [0])[0]
        if isinstance(d.get("theta_max"), list)
        else d.get("theta_max", 0)), 1)

    old_fill  = _load_sweep("sweep_fb_fill_l7_v2",
                            lambda d: (om(d), fl(d)),
                            lambda d: f"ω={om(d):.2f} f={fl(d)}")
    new_fill  = _load_sweep("fill_l7_mpi_ckpt",
                            lambda d: (om(d), fl(d)),
                            lambda d: f"ω={om(d):.2f} f={fl(d)}")
    old_theta = _load_sweep("sweep_fb_theta_l7",
                            lambda d: (om(d), th(d)),
                            lambda d: f"ω={om(d):.2f} θ={th(d):.0f}°")
    new_theta = _load_sweep("theta_l7_mpi_ckpt",
                            lambda d: (om(d), th(d)),
                            lambda d: f"ω={om(d):.2f} θ={th(d):.0f}°")

    pts = {"fill": [], "theta": []}
    for k in sorted(set(old_fill) & set(new_fill)):
        o, n = old_fill[k]["kla"], new_fill[k]["kla"]
        pts["fill"].append((o, n, abs(n - o) / o * 100, old_fill[k]["label"]))
    for k in sorted(set(old_theta) & set(new_theta)):
        o, n = old_theta[k]["kla"], new_theta[k]["kla"]
        pts["theta"].append((o, n, abs(n - o) / o * 100, old_theta[k]["label"]))
    return pts


def main():
    pts = _collect()
    all_pts = pts["fill"] + pts["theta"]
    n_total  = len(all_pts)
    max_dev  = max(p[2] for p in all_pts)
    mean_dev = sum(p[2] for p in all_pts) / n_total

    # ── figure layout ───────────────────────────────────────────────────────
    fig = plt.figure(figsize=(11, 5.2), facecolor=C_BG)
    fig.patch.set_facecolor(C_BG)

    gs = fig.add_gridspec(1, 2, left=0.07, right=0.97,
                          bottom=0.13, top=0.82, wspace=0.38)
    ax_scatter  = fig.add_subplot(gs[0])
    ax_residual = fig.add_subplot(gs[1])

    for ax in (ax_scatter, ax_residual):
        ax.set_facecolor(C_BG)
        ax.tick_params(colors=C_INK, labelsize=8)
        for sp in ax.spines.values():
            sp.set_color(C_MUTED)
            sp.set_linewidth(0.6)

    # ── helper: draw fine gridlines ────────────────────────────────────────
    def _grid(ax):
        ax.grid(True, color=C_MUTED, linewidth=0.35, linestyle="--", alpha=0.55)
        ax.set_axisbelow(True)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Panel A — identity scatter
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    ax = ax_scatter
    _grid(ax)

    kmax = 0.155
    # ±0.5 % tolerance band
    xs = np.linspace(0, kmax, 200)
    ax.fill_between(xs, xs * 0.995, xs * 1.005,
                    color=C_BAND, alpha=0.7, zorder=1,
                    label="±0.5 % band")
    # identity line
    ax.plot([0, kmax], [0, kmax],
            color=C_MUTED, linewidth=0.9, linestyle="-", zorder=2)

    for sweep, color in (("fill", C_ACCENT), ("theta", C_AMBER)):
        data = pts[sweep]
        ox = [p[0] for p in data]
        ny = [p[1] for p in data]
        ax.scatter(ox, ny, s=28, color=color, edgecolors="white",
                   linewidths=0.4, zorder=4,
                   label=f"Fill sweep (n={len(data)})" if sweep == "fill"
                         else f"θ sweep (n={len(data)}, segs 0–1)")

    ax.set_xlim(0, kmax)
    ax.set_ylim(0, kmax)
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(mticker.MultipleLocator(0.02))
    ax.yaxis.set_major_locator(mticker.MultipleLocator(0.02))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.set_xlabel(r"$k_La$ — baseline (no checkpoint)  $[\mathrm{nd}^{-1}]$",
                  fontsize=8.5, color=C_INK, labelpad=6)
    ax.set_ylabel(r"$k_La$ — MPI checkpoint restart  $[\mathrm{nd}^{-1}]$",
                  fontsize=8.5, color=C_INK, labelpad=6)
    ax.set_title("A   Identity plot", loc="left",
                 fontsize=9, fontweight="bold", color=C_INK, pad=6)

    legend = ax.legend(fontsize=7.5, framealpha=0.85, edgecolor=C_MUTED,
                       loc="upper left", handlelength=1.2)
    legend.get_frame().set_facecolor(C_BG)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Panel B — residual scatter
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    ax = ax_residual
    _grid(ax)

    # zero reference line + 0.5 % ceiling
    ax.axhline(0,   color=C_MUTED, linewidth=0.9, linestyle="-", zorder=2)
    ax.axhline(0.5, color=C_MUTED, linewidth=0.6, linestyle=":",
               zorder=2, alpha=0.7, label="0.5 % ref")
    ax.fill_between([0, kmax], 0, 0.5, color=C_BAND, alpha=0.45, zorder=1)

    for sweep, color in (("fill", C_ACCENT), ("theta", C_AMBER)):
        data = pts[sweep]
        ox   = [p[0] for p in data]
        dev  = [p[2] for p in data]
        ax.scatter(ox, dev, s=28, color=color, edgecolors="white",
                   linewidths=0.4, zorder=4)

    ax.set_xlim(0, kmax)
    ax.set_ylim(-0.005, 0.14)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(0.02))
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.set_xlabel(r"$k_La$ — baseline  $[\mathrm{nd}^{-1}]$",
                  fontsize=8.5, color=C_INK, labelpad=6)
    ax.set_ylabel(r"$|k_La^{\,\mathrm{ckpt}} - k_La^{\,\mathrm{base}}| \;/\; k_La^{\,\mathrm{base}}$  [%]",
                  fontsize=8.5, color=C_INK, labelpad=6)
    ax.set_title("B   Relative deviation", loc="left",
                 fontsize=9, fontweight="bold", color=C_INK, pad=6)

    # annotate max point
    worst = max(all_pts, key=lambda p: p[2])
    ax.annotate(f"max {worst[2]:.3f}%",
                xy=(worst[0], worst[2]),
                xytext=(worst[0] + 0.012, worst[2] + 0.008),
                fontsize=7, color=C_INK,
                arrowprops=dict(arrowstyle="-", color=C_MUTED, lw=0.7))

    # ── summary header ──────────────────────────────────────────────────────
    fig.text(0.52, 0.935,
             f"MAX DEVIATION  {max_dev:.3f} %     MEAN  {mean_dev:.3f} %     "
             f"N = {n_total} conditions     FIDELITY L7 · 16 MPI RANKS",
             ha="center", va="center",
             fontsize=8.5, color="white", fontfamily="monospace",
             fontweight="bold",
             bbox=dict(boxstyle="round,pad=0.45", facecolor=C_ACCENT,
                       edgecolor="none"))

    fig.suptitle(
        "MPI Checkpoint Restart — $k_La$ Agreement with Baseline",
        fontsize=12, fontweight="bold", color=C_INK, y=0.975,
        fontfamily="monospace"
    )
    fig.text(0.52, 0.03,
             "Baseline: sweep_fb_fill_l7_v2 / sweep_fb_theta_l7  (no checkpoint, single SLURM job)  ·  "
             "Checkpointed: BioReactor-mpi-stripped, chain restart via dump/restore",
             ha="center", fontsize=6.5, color=C_MUTED)

    _FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(_FIG_DIR / "checkpoint_validation.pdf", bbox_inches="tight", dpi=150)
    fig.savefig(_FIG_DIR / "checkpoint_validation.png", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved: {_FIG_DIR / 'checkpoint_validation.pdf'}")
    print(f"Saved: {_FIG_DIR / 'checkpoint_validation.png'}")
    print(f"\nSummary: n={n_total}  max|Δ|={max_dev:.4f}%  mean|Δ|={mean_dev:.4f}%")


if __name__ == "__main__":
    main()
