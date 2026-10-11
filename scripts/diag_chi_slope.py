"""Why does dtmix_0.95 jump with tiny rpm offsets (test P) while dtmix_0.50 does not?

Hypothesis (diary 2026-10-11): chi(t) approaches 0.95 slowly and with per-period oscillations, so a small change in
the chi level moves the first crossing by a large time. Measured here, per run, from the solver's own chi(t)
(scripts/postprocess.py pipeline, columns of tr_oxy.dat / vol_frac_interf.dat, no re-implementation):
  - t_x: first crossing time (s after release) of each threshold (= results.json dtmix);
  - slope: period-averaged d chi/dt at the crossing (chi averaged over a sliding window of one period, slope by a
    linear fit over +-1 period around t_x), in 1/s;
  - dwell: total time (s) during which chi is within +-0.01 of the threshold;
  - sens = 1/slope: seconds of crossing shift per unit change of the chi level (dt_x = sens * d chi).
Writes experiments/multifidelity/diag_chi_slope.png (chi(t) of the test P families) and prints the table.

Usage: uv run python scripts/diag_chi_slope.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402
from scripts import figstyle as fs  # noqa: E402
from scripts.chi_scale import locate  # noqa: E402

FAMILIES = {
    "22.5": ["kmix_l8_rpm22.5", "p8_rpm22.5_d1e-6_th7", "p8_rpm22.5_d1e-3_th7", "p8_rpm22.5_d1e-2_th7",
             "t5_l8_rpm22.4_th7", "t5_l8_rpm22.6_th7"],
    "30": ["kmix_l8_rpm30", "p8_rpm30_d1e-6_th7", "p8_rpm30_d1e-3_th7", "p8_rpm30_d1e-2_th7",
           "t5_l8_rpm29.9_th7", "t5_l8_rpm30.1_th7"],
}
THR = (0.50, 0.75, 0.95)


def chi_of(run_id):
    loc = locate(run_id)
    d = loc["data_dir"]
    import json
    params = json.loads((d / "params.json").read_text())
    t, c_sum, c_sum2 = pp._joined_cols(d, "tr_oxy.dat", [pp._COL_C2_LIQ_SUM, pp._COL_C2_LIQ_SUM2])
    _, f_liq = pp._joined_cols(d, "vol_frac_interf.dat", [pp._COL_F_LIQ])
    fm = float(f_liq.mean())
    s2 = c_sum2 / fm - (c_sum / fm) ** 2
    i0 = int(np.where(c_sum > 1e-10 * fm)[0][0])
    T_bio, T_per_nd = pp._t_scales(params)
    ts = (t[i0:] - t[i0]) * T_bio
    chi = np.clip(1.0 - s2[i0:] / s2[i0], 0.0, 1.0)
    return ts, chi, T_per_nd * T_bio


def measures(ts, chi, period):
    dt = np.median(np.diff(ts))
    w = max(1, int(round(period / dt)))
    chi_avg = np.convolve(chi, np.ones(w) / w, mode="same")
    out = {}
    for thr in THR:
        idx = np.argmax(chi >= thr)
        if chi[idx] < thr:
            out[thr] = None
            continue
        tx = ts[idx]
        m = (ts > tx - period) & (ts < tx + period) & (ts > period / 2) & (ts < ts[-1] - period / 2)
        slope = np.polyfit(ts[m], chi_avg[m], 1)[0] if m.sum() > 3 else np.nan
        dwell = float(np.sum(np.abs(chi - thr) < 0.01) * dt)
        out[thr] = dict(t_x=float(tx), slope=float(slope), sens=float(1 / slope) if slope > 0 else np.inf,
                        dwell=dwell)
    return out


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    print(f"{'run':26s} " + " ".join(f"{'t_x' + str(th):>9s} {'slope':>8s} {'sens':>7s} {'dwell':>6s}" for th in THR))
    for ax, (base, runs) in zip(axes, FAMILIES.items()):
        for k, r in enumerate(runs):
            try:
                ts, chi, per = chi_of(r)
            except Exception as e:  # noqa: BLE001
                print(f"{r:26s} unavailable: {e}")
                continue
            ms = measures(ts, chi, per)
            row = []
            for th in THR:
                v = ms[th]
                row.append("      n/a" + " " * 24 if v is None else
                           f"{v['t_x']:9.2f} {v['slope']:8.5f} {v['sens']:7.1f} {v['dwell']:6.1f}")
            print(f"{r:26s} " + " ".join(row))
            ax.plot(ts, chi, lw=0.8, ls=":", color=plt.cm.viridis(k / max(1, len(runs) - 1)), label=r)
        ax.axhline(0.95, color="0.5", lw=0.6)
        ax.set_xlabel("t after release (s)")
        ax.set_title(f"{base} rpm")
        ax.set_ylim(0.85, 1.0)
    axes[0].set_ylabel(r"$\chi$")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7)
    fig.tight_layout()
    out = ROOT / "experiments/multifidelity/diag_chi_slope.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
