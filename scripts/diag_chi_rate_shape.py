"""Test of the claim "ln(1 - chi) has a kink near chi = 0.75".

The local decay rate lam(chi) = -d ln(1 - chi)/dt sets the whole curve:
tau(chi) = integral of dchi / ((1 - chi) lam). A kink in ln(1 - chi) is a step
in lam. We average lam over one rocking period (removes the phase wobble) and
plot lam(chi) / lam(0.5) for L6-L9 at each rpm.

Falsifiers of "a fixed kink at 0.75":
  F1 the step location moves with rpm        -> the kink depends on x;
  F2 the step location moves with grid level -> the kink is numerical;
  F3 no step at all at some level/rpm.
Artifact check: tracer mass and liquid volume drift over the run (chi uses the
time-mean liquid volume, so a drift in either one bends chi).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs                                   # noqa: E402
from scripts.pilot_truncated_chi import chi_curve, RPMS              # noqa: E402
from scripts.postprocess import (_joined_cols, _COL_C2_LIQ_SUM,      # noqa: E402
                                 _COL_F_LIQ)

FAMILIES = {6: "kmix_l6_rpm{:g}", 7: "kmix_l7_rpm{:g}", 8: "kmix_l8_rpm{:g}",
            9: "fig9_l9_rpm{:g}"}
CHI_GRID = np.linspace(0.30, 0.95, 66)


def period_rate(run):
    """Period-averaged decay rate lam (1/s) against chi, on CHI_GRID."""
    t, chi = chi_curve(run)
    T_s = 60.0 / rpm_of(run)                                      # period in s
    y = np.log(np.clip(1 - chi, 1e-9, 1))
    tp, tm = t + T_s / 2, t - T_s / 2
    ok = (tm >= t[0]) & (tp <= t[-1])
    lam = -(np.interp(tp[ok], t, y) - np.interp(tm[ok], t, y)) / T_s
    c = chi[ok]
    # chi is not monotone inside a period; bin by chi of the window centre
    out = np.full(CHI_GRID.size, np.nan)
    for k, g in enumerate(CHI_GRID):
        m = np.abs(c - g) < 0.01
        if m.sum() >= 3:
            out[k] = np.median(lam[m])
    return out, chi.max()


def rpm_of(run):
    return float(run.split("rpm")[1].split("_")[0])


def drift(run):
    """Max relative drift of tracer mass and liquid volume after release."""
    d = ROOT / "runs" / run
    t, s = _joined_cols(d, "tr_oxy.dat", [_COL_C2_LIQ_SUM])
    tf, f = _joined_cols(d, "vol_frac_interf.dat", [_COL_F_LIQ])
    i0 = int(np.argmax(s > 1e-10 * f.mean()))
    s_rel = s[i0:] / s[i0] - 1
    f_rel = f[f > 0] / f.mean() - 1
    j = int(np.argmax(np.abs(s_rel)))
    return s_rel[j], np.abs(f_rel).max()


def main():
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(2, 5, figsize=(14, 5.6), sharex=True, sharey=True)
    cols = {6: "#9ecae1", 7: "#4292c6", 8: "#08519c", 9: "#08306b"}
    print(f"{'run':22s} {'chi_max':>7s} {'tracer drift':>12s} {'liq. drift':>10s} "
          f"{'lam(.85)/lam(.5)':>16s} {'lam(.6)/lam(.5)':>15s}")
    for ax, r in zip(axes.ravel(), RPMS):
        for lv, tmpl in FAMILIES.items():
            run = tmpl.format(r)
            if not (ROOT / "runs" / run).exists():
                continue
            lam, cmax = period_rate(run)
            k5 = np.argmin(np.abs(CHI_GRID - 0.5))
            ratio = lam / lam[k5]
            ax.plot(CHI_GRID, ratio, color=cols[lv], lw=1.2, label=f"L{lv}")
            ds, df = drift(run)
            k85 = np.argmin(np.abs(CHI_GRID - 0.85))
            k6 = np.argmin(np.abs(CHI_GRID - 0.60))
            print(f"{run:22s} {cmax:7.3f} {ds:+12.2%} {df:10.2%} "
                  f"{ratio[k85]:16.2f} {ratio[k6]:15.2f}")
        ax.axvline(0.75, color="0.6", lw=0.7, ls=":")
        ax.axhline(1.0, color="0.8", lw=0.6)
        ax.set_yscale("log")
        ax.set_yticks([0.25, 0.5, 1, 2, 4], ["0.25", "0.5", "1", "2", "4"])
        ax.minorticks_off()
        ax.set_title(f"{r:g} rpm", fontsize=10, loc="left")
        ax.grid(**fs.GRID_KW)
    for ax in axes[-1]:
        ax.set_xlabel(r"$\chi$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\lambda(\chi)/\lambda(0.5)$")
    axes[0, -1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    out = ROOT / "experiments/multifidelity/chi_decay_rate_shape.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main()
