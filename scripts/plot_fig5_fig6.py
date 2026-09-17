"""Replicas of Kim et al. (2024) Figs. 5 and 6 -- initial tracer configuration.

Fig. 5: the degree of mixing chi(t) for four initial tracer configurations --
top half, left half, a circle at the centre, and a line across the middle --
with chi = 1.00, 0.95 and 0.50 marked.

Fig. 6: the tracer fields at chi = 0 and at chi = 0.5 for four configurations,
Kim's fourth being the BOTTOM half rather than the left.

One run gives all of it. The four configurations are four independent soluble
tracers advected by the same flow (BioReactor.c, EXTRA_TRACERS), and tr_oxy.dat
logs the first and second moments of each, so chi(t) for all four comes from
one time series. The bottom-half panel needs no tracer of its own: sigma^2 is
invariant under c -> 1 - c, so its chi(t) is identical to the top half's and
its field is 1 - c2 wherever there is water.

chi is normalised per tracer, by that tracer's OWN sigma^2 at injection.
sigma^2_max is 0.25 only for a configuration that fills exactly half the
water; a circle covering a fraction phi gives phi(1 - phi), and normalising
all four against 0.25 would scale the circle's and the line's curves by an
arbitrary factor while leaving the half-fills untouched.

Usage:  uv run python scripts/plot_fig5_fig6.py <run_id>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
OUT_DIR = ROOT / "experiments/kimetal2024/figure_replicas"
sys.path.insert(0, str(ROOT))
from scripts.fields import (  # noqa: E402
    cell_centres, liquid_mask, load_snapshot, nearest_snapshot,
)
from scripts.postprocess import _COL_F_LIQ, _joined_cols, _t_scales  # noqa: E402

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

# tracer key -> (sum column, sum^2 column, Kim's label, colour, linestyle)
TRACERS = {
    "c2": (8,  9,  "top half",  "#0072B2", "--"),
    "c1": (6,  7,  "left half", "#D55E00", "-."),
    "c3": (10, 11, "circle",    "#009E73", (0, (3, 1, 1, 1, 1, 1))),
    "c":  (4,  5,  "line",      "#E69F00", (0, (5, 2))),
}
MARKS = [(1.00, "-", "perfectly mixed"), (0.95, "--", "mixing time"),
         (0.50, ":", "partially mixed")]
BAG_HALF_HEIGHT = 0.143
# Fig. 6's panels, in Kim's order. "1-c2" is the bottom half by complement.
SNAP_PANELS = [("c2", "top half"), ("1-c2", "bottom half"),
               ("c3", "circle"), ("c", "line")]


def chi_series(run_dir: Path, col_sum: int, col_sum2: int):
    """(t non-dim, chi) for one tracer, normalised by its own sigma^2 at t_inj."""
    t, s1, s2 = _joined_cols(run_dir, "tr_oxy.dat", [col_sum, col_sum2])
    _, f_liq = _joined_cols(run_dir, "vol_frac_interf.dat", [_COL_F_LIQ])
    f_mean = float(f_liq.mean())
    if f_mean <= 0:
        raise ValueError("no liquid")
    mean, mean2 = s1 / f_mean, s2 / f_mean
    sigma2 = mean2 - mean ** 2
    nz = np.where(s1 > 1e-10 * f_mean)[0]
    if nz.size == 0:
        return None, None, np.nan
    i0 = int(nz[0])
    s0 = float(sigma2[i0])
    if s0 <= 0:
        return None, None, s0
    chi = np.clip(1.0 - sigma2 / s0, 0.0, 1.0)
    return t[i0:] - t[i0], chi[i0:], s0


def fig5(run_dir: Path, T_per_nd: float) -> dict:
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    t_half = {}
    for key, (cs_, cs2, label, colour, ls) in TRACERS.items():
        t, chi, s0 = chi_series(run_dir, cs_, cs2)
        if t is None:
            print(f"  {key} ({label}): never injected -- skipped")
            continue
        ax.plot(t / T_per_nd, chi, color=colour, ls=ls, lw=1.4, label=label)
        j = int(np.argmax(chi >= 0.5))
        t_half[key] = float(t[j]) if chi[j] >= 0.5 else np.nan
        print(f"  {key:>2} {label:<10} sigma^2_max={s0:.4f}  "
              f"chi_end={chi[-1]:.4f}  t(chi=0.5)="
              f"{t_half[key] / T_per_nd if np.isfinite(t_half[key]) else float('nan'):.2f} T_p")
    for level, ls, label in MARKS:
        ax.axhline(level, color="0.6", ls=ls, lw=0.9, zorder=0)
        ax.text(1.005, level, rf"$\chi={level:.2f}$", transform=
                ax.get_yaxis_transform(), fontsize=7.5, va="center",
                color="0.4")
    ax.set_xlabel(r"$(t-t_{inj})/T_p$", fontsize=12)
    ax.set_ylabel(r"$\chi$", fontsize=12)
    ax.set_ylim(0, 1.05)
    ax.tick_params(which="both", direction="in")
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")
    fig.tight_layout()
    out = OUT_DIR / "replicated_Fig5.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved {out}")
    return t_half


def _panel(ax, X, Y, field, mask):
    ax.pcolormesh(X, Y, np.where(mask, field, np.nan), cmap="viridis",
                  vmin=0, vmax=1, shading="auto", rasterized=True)
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-BAG_HALF_HEIGHT, BAG_HALF_HEIGHT)
    ax.set_aspect("equal")
    ax.tick_params(which="both", direction="in", labelsize=7)


def fig6(run_dir: Path, t_half: dict, t_inj_nd: float) -> None:
    first = nearest_snapshot(run_dir, t_inj_nd)
    if first is None:
        print("  no snapshots recorded -- Fig 6 needs a run with a snapshot "
              "window covering the tracer release")
        return
    fig, axes = plt.subplots(len(SNAP_PANELS), 2, figsize=(8.0, 6.4))
    for row, (key, label) in enumerate(SNAP_PANELS):
        src = "c2" if key == "1-c2" else key
        # chi = 0 is the release instant; chi = 0.5 is read from THIS tracer's
        # own crossing, which differs between configurations by design.
        t_mid = t_inj_nd + t_half.get(src, np.nan)
        for col, t_want in enumerate((t_inj_nd, t_mid)):
            ax = axes[row, col]
            if not np.isfinite(t_want):
                ax.text(0.5, 0.5, r"$\chi=0.5$ not reached",
                        transform=ax.transAxes, ha="center", fontsize=8)
                ax.set_axis_off()
                continue
            s = load_snapshot(nearest_snapshot(run_dir, t_want))
            X, Y = cell_centres(s["n"])
            m = liquid_mask(s)
            field = 1.0 - s[src] if key == "1-c2" else s[src]
            _panel(ax, X, Y, field, m)
            if col == 0:
                ax.set_ylabel(label, fontsize=9)
    axes[0, 0].set_title(r"$\chi=0$", fontsize=10)
    axes[0, 1].set_title(r"$\chi=0.5$", fontsize=10)
    fig.tight_layout()
    out = OUT_DIR / "replicated_Fig6.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"saved {out}")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    run_dir = ROOT / "runs" / sys.argv[1]
    params = json.loads((run_dir / "params.json").read_text())
    _, T_per_nd = _t_scales(params)
    t_inj_nd = params["n_mix_cycles"] * T_per_nd + params.get("t_checkpoint", 0.0)
    t_half = fig5(run_dir, T_per_nd)
    fig6(run_dir, t_half, t_inj_nd)


if __name__ == "__main__":
    main()
