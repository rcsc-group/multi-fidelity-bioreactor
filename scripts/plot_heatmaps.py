"""Generate heatmap figures for both sweep configurations.

Reads all completed results from runs/ and produces:
  experiments/<exp_dir>/figures/heatmap_theta_sweep_l<N>.pdf
  experiments/<exp_dir>/figures/heatmap_fill_sweep_l<N>.pdf

The output directory is derived from --exp-suffix: the script searches
experiments/ for a unique subdirectory whose name contains the suffix and
writes figures there.  Falls back to experiments/figures/ if no unique match.

Usage
-----
    python scripts/plot_heatmaps.py [--fidelity 7] [--exp-suffix theta_l7_mpi_ckpt]
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

_PROJECT_ROOT  = Path(__file__).parents[1]
_RUNS_ROOT     = _PROJECT_ROOT / "runs"
_EXPERIMENTS   = _PROJECT_ROOT / "experiments"
_FIG_DIR       = _EXPERIMENTS / "figures"   # fallback for cross-experiment figures


def _fig_dir_for_suffix(exp_suffix: str | None) -> Path:
    """Return the figures/ dir inside the matching experiment folder.

    Searches experiments/ for a unique subdirectory whose name contains
    exp_suffix.  Falls back to the shared experiments/figures/ if zero or
    multiple matches are found.
    """
    if not exp_suffix:
        return _FIG_DIR
    matches = [p for p in _EXPERIMENTS.iterdir()
               if p.is_dir() and exp_suffix in p.name]
    if len(matches) == 1:
        return matches[0] / "figures"
    return _FIG_DIR

# 15 KPIs arranged in 5 rows x 3 columns — each row is one semantic group.
# dtmix_0.95 omitted: it is ~95% NaN at fidelity 8 (80-cycle budget does not
# reach 95% uniformity for most conditions) and produces blank columns there.
_KPIS = [
    # row 1 — 5-pt fit kLa
    ("kLa_10",      r"$k_La$ at $C^*=0.10$ (h$^{-1}$, 5pt)",      "YlOrRd"),
    ("kLa_25",      r"$k_La$ at $C^*=0.25$ (h$^{-1}$, 5pt)",      "YlOrRd"),
    ("kLa_50",      r"$k_La$ at $C^*=0.50$ (h$^{-1}$, 5pt)",      "YlOrRd"),
    # row 2 — instantaneous kLa
    ("kLa_inst_10", r"$k_La$ at $C^*=0.10$ (h$^{-1}$, inst)",         "YlOrBr"),
    ("kLa_inst_25", r"$k_La$ at $C^*=0.25$ (h$^{-1}$, inst)",         "YlOrBr"),
    ("kLa_inst_50", r"$k_La$ at $C^*=0.50$ (h$^{-1}$, inst)",         "YlOrBr"),
    # row 3 — mixing timescale + bulk flow
    ("dtmix_0.50",  r"$\Delta t_{mix}$ at $\chi=0.50$ (s)", "Blues_r"),
    ("dtmix_0.75",  r"$\Delta t_{mix}$ at $\chi=0.75$ (s)", "Blues_r"),
    ("vor_mean",    r"$\langle|\xi|\rangle$ (1/s)",          "Purples"),
    # row 4 — QSS shear stress percentiles
    ("tau_95_qss",  r"$\tau_{95}$ QSS (Pa)",                 "RdPu"),
    ("tau_98_qss",  r"$\tau_{98}$ QSS (Pa)",                 "RdPu"),
    ("tau_100_qss", r"$\tau_{100}$ QSS (Pa)",                "RdPu"),
    # row 5 — peak shear stress percentiles
    ("tau_95_max",  r"$\tau_{95}$ max (Pa)",                 "Reds"),
    ("tau_98_max",  r"$\tau_{98}$ max (Pa)",                 "Reds"),
    ("tau_100_max", r"$\tau_{100}$ max (Pa)",                "Reds"),
]


def _load_results(fidelity: int = 7, exp_suffix: str | None = None) -> list[dict]:
    """Load results matching fidelity and optional experiment-dir suffix."""
    records = {}
    for f in sorted(_RUNS_ROOT.glob("*/results.json"), key=lambda p: p.stat().st_mtime):
        try:
            r = json.loads(f.read_text())
            p = json.loads((f.parent / "params.json").read_text())
            if "vor_mean" not in r:
                continue
            if p.get("fidelity") != fidelity:
                continue
            if exp_suffix:
                exp_dir = p.get("_experiment_dir", "")
                if exp_suffix not in exp_dir:
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

    ncols = 3
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--fidelity", type=int, default=7)
    ap.add_argument("--exp-suffix", default="mpi_ckpt",
                    help="Only include runs whose _experiment_dir contains this string")
    ap.add_argument("--theta-fill", type=float, default=0.5)
    ap.add_argument("--fill-theta", type=float, default=7.0)
    ap.add_argument("--sweep-type", choices=["theta", "fill", "both"], default="both",
                    help="Which heatmap to generate: theta-vs-rpm, fill-vs-rpm, or both")
    args = ap.parse_args()

    records = _load_results(fidelity=args.fidelity, exp_suffix=args.exp_suffix)
    print(f"Loaded {len(records)} unique (theta, fill, omega_b) data points "
          f"[fidelity={args.fidelity}, exp_suffix={args.exp_suffix!r}]")

    tag     = f"l{args.fidelity}"
    fig_dir = _fig_dir_for_suffix(args.exp_suffix)

    # Theta sweep: theta_max vs rpm, at the specified fill level
    if args.sweep_type in ("theta", "both"):
        theta_recs = [r for r in records if r["fill"] == args.theta_fill]
        if theta_recs:
            _make_figure(
                theta_recs,
                row_key="theta",
                col_key="rpm",
                row_label=r"$\theta_{max}$ (deg)",
                col_label="$f_b$ (rpm)",
                title=f"Theta sweep — KPI heatmaps ({tag}, fill={args.theta_fill})",
                out_path=fig_dir / f"heatmap_theta_sweep_{tag}.pdf",
            )
        else:
            print(f"No theta-sweep records at fill={args.theta_fill}")

    # Fill sweep: fill_level vs rpm, at the specified theta
    if args.sweep_type in ("fill", "both"):
        fill_recs = [r for r in records if r["theta"] == args.fill_theta]
        if fill_recs:
            _make_figure(
                fill_recs,
                row_key="fill",
                col_key="rpm",
                row_label="Fill level",
                col_label="$f_b$ (rpm)",
                title=f"Fill sweep — KPI heatmaps ({tag}, $\\theta_{{max}}$={args.fill_theta}°)",
                out_path=fig_dir / f"heatmap_fill_sweep_{tag}.pdf",
            )
        else:
            print(f"No fill-sweep records at theta={args.fill_theta}")


if __name__ == "__main__":
    main()
