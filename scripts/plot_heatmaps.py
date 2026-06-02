"""Generate heatmap figures for both sweep configurations.

Reads all completed results from runs/ and produces:
  experiments/figures/heatmap_theta_sweep.pdf   — theta x omega_b
  experiments/figures/heatmap_fill_sweep.pdf    — fill_level x omega_b

Usage
-----
    python scripts/plot_heatmaps.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

_PROJECT_ROOT = Path(__file__).parents[1]
_RUNS_ROOT    = _PROJECT_ROOT / "runs"
_FIG_DIR      = _PROJECT_ROOT / "experiments" / "figures"

_KPIS = [
    ("kLa_10",      r"$k_La$ at $C^*=0.10$ (5pt fit)",      "YlOrRd"),
    ("kLa_25",      r"$k_La$ at $C^*=0.25$ (5pt fit)",      "YlOrRd"),
    ("kLa_50",      r"$k_La$ at $C^*=0.50$ (5pt fit)",      "YlOrRd"),
    ("kLa_inst_10", r"$k_La$ at $C^*=0.10$ (inst)",         "YlOrBr"),
    ("kLa_inst_25", r"$k_La$ at $C^*=0.25$ (inst)",         "YlOrBr"),
    ("kLa_inst_50", r"$k_La$ at $C^*=0.50$ (inst)",         "YlOrBr"),
    ("dtmix_0.50",  r"$\Delta t_{mix}$ at $\chi=0.50$ (s)", "Blues_r"),
    ("dtmix_0.75",  r"$\Delta t_{mix}$ at $\chi=0.75$ (s)", "Blues_r"),
    ("dtmix_0.95",  r"$\Delta t_{mix}$ at $\chi=0.95$ (s)", "Blues_r"),
    ("vor_mean",    r"$\langle|\xi|\rangle$ (1/s)",          "Purples"),
]


def _load_results() -> list[dict]:
    """Load all fidelity-5 results with vor_mean (new postprocess format)."""
    records = {}
    for f in sorted(_RUNS_ROOT.glob("*/results.json"), key=lambda p: p.stat().st_mtime):
        try:
            r = json.loads(f.read_text())
            p = json.loads((f.parent / "params.json").read_text())
            if "vor_mean" not in r or p.get("fidelity") != 5:
                continue
            th = round(float(p.get("theta_max", [0])[0]), 1)
            fl = round(float(p.get("fill_level", 0.5)), 2)
            ob = round(float(p.get("omega_b", 0)) * 60 / (2 * math.pi), 2)
            key = (th, fl, ob)
            records[key] = {**r, "theta": th, "fill": fl, "rpm": ob}
        except Exception:
            continue
    return list(records.values())


def _pivot(records: list[dict], row_key: str, col_key: str,
           row_vals: list, col_vals: list, kpi: str) -> np.ndarray:
    """Build a 2-D array [rows x cols] for the given KPI."""
    lookup = {(r[row_key], r[col_key]): r.get(kpi, float("nan"))
              for r in records}
    arr = np.full((len(row_vals), len(col_vals)), float("nan"))
    for i, rv in enumerate(row_vals):
        for j, cv in enumerate(col_vals):
            arr[i, j] = lookup.get((rv, cv), float("nan"))
    return arr


def _make_figure(records: list[dict], row_key: str, col_key: str,
                 row_label: str, col_label: str, title: str,
                 out_path: Path) -> None:
    row_vals = sorted({r[row_key] for r in records})
    col_vals = sorted({r[col_key] for r in records})

    ncols = 5
    nrows = math.ceil(len(_KPIS) / ncols)
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 3.2, nrows * 2.6),
                             constrained_layout=True)
    axes = np.array(axes).flatten()

    for ax, (kpi, label, cmap) in zip(axes, _KPIS):
        data = _pivot(records, row_key, col_key, row_vals, col_vals, kpi)
        # mask NaN
        masked = np.ma.masked_invalid(data)
        im = ax.pcolormesh(
            range(len(col_vals) + 1),
            range(len(row_vals) + 1),
            masked,
            cmap=cmap,
            shading="flat",
        )
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_xticks(np.arange(len(col_vals)) + 0.5)
        ax.set_xticklabels([f"{v:.1f}" for v in col_vals],
                           rotation=45, ha="right", fontsize=7)
        ax.set_yticks(np.arange(len(row_vals)) + 0.5)
        ax.set_yticklabels([str(v) for v in row_vals], fontsize=7)
        ax.set_xlabel(col_label, fontsize=8)
        ax.set_ylabel(row_label, fontsize=8)
        ax.set_title(label, fontsize=8)

    # hide unused axes
    for ax in axes[len(_KPIS):]:
        ax.set_visible(False)

    fig.suptitle(title, fontsize=11, fontweight="bold")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


def main() -> None:
    records = _load_results()
    print(f"Loaded {len(records)} unique (theta, fill, omega_b) data points")

    # Theta sweep: theta_max vs rpm, fill=0.5
    theta_recs = [r for r in records if r["fill"] == 0.5]
    _make_figure(
        theta_recs,
        row_key="theta",
        col_key="rpm",
        row_label=r"$\theta_{max}$ (deg)",
        col_label="$f_b$ (rpm)",
        title="Theta sweep — KPI heatmaps (fill level = 0.5)",
        out_path=_FIG_DIR / "heatmap_theta_sweep.pdf",
    )

    # Fill sweep: fill_level vs rpm, theta=7
    fill_recs = [r for r in records if r["theta"] == 7.0]
    _make_figure(
        fill_recs,
        row_key="fill",
        col_key="rpm",
        row_label="Fill level",
        col_label="$f_b$ (rpm)",
        title="Fill sweep — KPI heatmaps ($\\theta_{max}$ = 7°)",
        out_path=_FIG_DIR / "heatmap_fill_sweep.pdf",
    )


if __name__ == "__main__":
    main()
