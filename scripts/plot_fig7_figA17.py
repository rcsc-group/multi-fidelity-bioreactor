"""Replicas of Kim et al. (2024) Fig. 7 and Fig. A.17 -- oxygen transfer.

Fig. 7: the normalised oxygen level in the water C*(t), the mass transfer
coefficient kLa(t), and the oxygen field at C* = 0.10, 0.25 and 0.50.

Fig. A.17: kLa computed three ways -- a global fit from t = 0, and local fits
over 5 and over 11 consecutive points -- with beta_0 and the RMSE of each.

The three estimators are the appendix's own comparison, and they are not
interchangeable. From dC*/dt = kLa (1 - C*),

    global:   kLa(t0) = -(1/t0) ln[1 - C*(t0)]          (beta_0 = 1 by
                                                          construction)
    local:    ln(1 - C*) = -kLa t + ln(beta_0), fitted over a window

The global fit assumes the first-order model has held since t = 0 with a
single rate; the local fits do not, and beta_0 is exactly the measure of how
far that assumption has drifted. Reporting one number as "the" kLa hides
which of the three it was -- which is why Kim publishes both `kLa_inst_*` and
`kLa_exp5pts_*` columns, and why this project stores both.

The field panels are picked from the recorded snapshots by finding the
C* crossings in the time series; the solver is never told a threshold.

Usage:  uv run python scripts/plot_fig7_figA17.py <run_id>
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
from scripts.postprocess import _compute_c_star, _t_scales  # noqa: E402

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

THRESHOLDS = [0.10, 0.25, 0.50]
BAG_HALF_HEIGHT = 0.143
SEC_PER_HOUR = 3600.0


def local_fit(t_s, c_star, npts: int):
    """Rolling log-linear fit: (t centres, kLa [1/h], beta_0, RMSE).

    ln(1 - C*) = -kLa t + ln(beta_0) over `npts` consecutive samples.
    """
    y = np.log(np.clip(1.0 - c_star, 1e-12, None))
    half = npts // 2
    ts, kla, beta, rmse = [], [], [], []
    for i in range(half, len(t_s) - half):
        sl = slice(i - half, i + half + 1)
        slope, intercept = np.polyfit(t_s[sl], y[sl], 1)
        resid = y[sl] - (slope * t_s[sl] + intercept)
        ts.append(t_s[i])
        kla.append(-slope * SEC_PER_HOUR)
        beta.append(float(np.exp(intercept)))
        rmse.append(float(np.sqrt((resid ** 2).mean())))
    return (np.array(ts), np.array(kla), np.array(beta), np.array(rmse))


def global_fit(t_s, c_star):
    """kLa(t0) = -(1/t0) ln[1 - C*(t0)], in 1/h. beta_0 is 1 by construction."""
    with np.errstate(divide="ignore", invalid="ignore"):
        k = -np.log(np.clip(1.0 - c_star, 1e-12, None)) / t_s * SEC_PER_HOUR
    return k


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    run_dir = ROOT / "runs" / sys.argv[1]
    params = json.loads((run_dir / "params.json").read_text())
    T_bio, T_per_nd = _t_scales(params)
    t_nd, c_star = _compute_c_star(run_dir)
    if t_nd.size < 11:
        raise SystemExit(f"{run_dir.name}: only {t_nd.size} oxygen samples")
    # The series starts at the run's own t = 0, and oxygen is not released
    # until n_mix_cycles later -- ~80 periods of C* = 0 sit in front of it.
    # The global fit is kLa(t0) = -ln(1 - C*)/t0 with t0 measured FROM THE
    # RELEASE; anchoring it at the run start divides by the wrong elapsed time
    # and understates kLa by the ratio of the two (a factor of ~5 here). It
    # also leaves the local fits fitting a flat zero for most of their range.
    released = np.where(c_star > 1e-9)[0]
    if released.size == 0:
        raise SystemExit(f"{run_dir.name}: oxygen never released")
    i0 = int(released[0])
    t_nd, c_star = t_nd[i0:], c_star[i0:]
    t_rel_s = (t_nd - t_nd[0]) * T_bio
    # t0 = 0 makes the global fit singular at its first sample; step to the
    # next one rather than dividing by zero.
    t_rel_s[0] = t_rel_s[1] if t_rel_s[0] == 0 else t_rel_s[0]
    print(f"  oxygen released at t={t_nd[0]:.4f} nd; "
          f"{t_nd.size} samples over {t_rel_s[-1]:.1f} s")

    k_glob = global_fit(t_rel_s, c_star)
    t5, k5, b5, r5 = local_fit(t_rel_s, c_star, 5)
    t11, k11, b11, r11 = local_fit(t_rel_s, c_star, 11)

    crossings = {}
    for thr in THRESHOLDS:
        idx = int(np.argmax(c_star >= thr))
        crossings[thr] = t_nd[idx] if c_star[idx] >= thr else np.nan
        print(f"  C*={thr:.2f} crossed at t={crossings[thr]:.4f} nd "
              f"({(crossings[thr] - t_nd[0]) * T_bio:.1f} s after release)"
              if np.isfinite(crossings[thr]) else
              f"  C*={thr:.2f} never reached (max {c_star.max():.3f})")

    # ── Fig 7 ────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(9.4, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1.0])
    ax_c = fig.add_subplot(gs[0, :2])
    ax_c.plot((t_nd - t_nd[0]) / T_per_nd, c_star, color="#0072B2", lw=1.4)
    for thr, ls in zip(THRESHOLDS, ("-", "--", ":")):
        ax_c.axhline(thr, color="0.6", ls=ls, lw=0.9, zorder=0)
    ax_c.set_xlabel(r"$(t-t_{inj})/T_p$", fontsize=11)
    ax_c.set_ylabel(r"$C^*_{w,\mathrm{oxy}}$", fontsize=11, color="#0072B2")
    ax_c.tick_params(which="both", direction="in")
    ax_k = ax_c.twinx()
    ax_k.plot(t5 / (T_per_nd * T_bio), k5, color="#D55E00", lw=1.2)
    ax_k.set_ylabel(r"$k_La$ (h$^{-1}$)", fontsize=11, color="#D55E00")
    ax_k.tick_params(axis="y", colors="#D55E00")

    for col, thr in enumerate(THRESHOLDS):
        ax = fig.add_subplot(gs[1, col])
        p = (nearest_snapshot(run_dir, crossings[thr])
             if np.isfinite(crossings[thr]) else None)
        if p is None:
            ax.text(0.5, 0.5, "no snapshot", transform=ax.transAxes,
                    ha="center", fontsize=8)
            ax.set_axis_off()
            continue
        s = load_snapshot(p)
        X, Y = cell_centres(s["n"])
        m = liquid_mask(s)
        ax.pcolormesh(X, Y, np.where(m, s["oxy"], np.nan), cmap="viridis",
                      shading="auto", rasterized=True)
        ax.set_xlim(-0.5, 0.5)
        ax.set_ylim(-BAG_HALF_HEIGHT, BAG_HALF_HEIGHT)
        ax.set_aspect("equal")
        ax.set_title(rf"$C^*={thr:.2f}$", fontsize=9)
        ax.tick_params(which="both", direction="in", labelsize=7)
    fig.tight_layout()
    out7 = OUT_DIR / "replicated_Fig7.png"
    fig.savefig(out7, dpi=150, bbox_inches="tight")
    print(f"saved {out7}")

    # ── Fig A.17 ─────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))
    axes[0].plot(t_rel_s, k_glob, color="#0072B2", lw=1.3, label="global fit")
    axes[0].plot(t5, k5, color="#D55E00", lw=1.2, ls="--", label="local, 5 pts")
    axes[0].plot(t11, k11, color="#009E73", lw=1.2, ls="-.",
                 label="local, 11 pts")
    for thr in THRESHOLDS:
        idx = int(np.argmax(c_star >= thr))
        if c_star[idx] >= thr:
            axes[0].axvline(t_rel_s[idx], color="0.7", lw=0.8, zorder=0)
    axes[0].set_xlabel(r"$t-t_{inj}$ (s)", fontsize=11)
    axes[0].set_ylabel(r"$k_La$ (h$^{-1}$)", fontsize=11)
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].tick_params(which="both", direction="in")

    axes[1].plot(t5, b5, color="#D55E00", lw=1.2, label=r"$\beta_0$, 5 pts")
    axes[1].axhline(1.0, color="#0072B2", lw=1.2,
                    label=r"$\beta_0$, global $\equiv 1$")
    axes[1].set_xlabel(r"$t-t_{inj}$ (s)", fontsize=11)
    axes[1].set_ylabel(r"$\beta_0$", fontsize=11)
    ax_r = axes[1].twinx()
    ax_r.plot(t5, r5, color="0.45", lw=1.0, ls=":", label="RMSE, 5 pts")
    ax_r.set_ylabel("RMSE", fontsize=11, color="0.45")
    ax_r.tick_params(axis="y", colors="0.45")
    h1, l1 = axes[1].get_legend_handles_labels()
    h2, l2 = ax_r.get_legend_handles_labels()
    # beta_0 spikes by orders of magnitude in the first window after release,
    # where 1 - C* is still within rounding of 1 and the fitted intercept is
    # unconstrained. Clip to the settled range so the panel shows the drift it
    # exists to show rather than one transient.
    settled = b5[len(b5) // 10:]
    if settled.size:
        axes[1].set_ylim(0.0, float(np.percentile(settled, 99)) * 1.2)
    axes[1].legend(h1 + h2, l1 + l2, fontsize=8, frameon=False,
                   loc="center left", bbox_to_anchor=(1.22, 0.5))
    for ax, letter in zip(axes, "ab"):
        ax.text(-0.14, 1.03, rf"$({letter})$", transform=ax.transAxes,
                fontsize=13, style="italic")
    fig.tight_layout()
    outa = OUT_DIR / "replicated_FigA17.png"
    fig.savefig(outa, dpi=150, bbox_inches="tight")
    print(f"saved {outa}")


if __name__ == "__main__":
    main()
