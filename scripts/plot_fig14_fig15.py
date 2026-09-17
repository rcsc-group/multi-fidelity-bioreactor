"""Replicas of Kim et al. (2024) Figs. 14 and 15 -- the free surface.

Fig. 14: the maximum surface elevation AT THE LEFT END, normalised by the
elevation of the equivalent static case, against rocking frequency (a) and
rocking angle (b).

Fig. 15: frequency spectra of the surface elevation, with frequency
normalised by the rocking frequency, at three rpm and three angles.

Both need the elevation at a FIXED END of the bag, not the global maximum.
`posY_max` is wherever the crest happens to be, and the crest travels from one
end to the other every half period, so its spectrum is dominated by that
traverse and its maximum is the same number at every condition. The solver
therefore probes the outermost 2% of the bag at each end
(`posY_left`, `posY_right` in vol_frac_interf.dat); runs recorded before
2026-09-16 have neither column and cannot produce either figure.

The static reference is the flat surface at the fill level, which for a
half-filled bag is the mid-plane: elevation measured from the bag floor is
then half the bag height, exactly, at every condition. It is computed rather
than simulated because a static case has no dynamics to get wrong.

Spectra are taken over whole settled cycles only. The transient after the
soft-start ramp is not periodic, and windowing it into the transform puts
power at every frequency.

Usage:
    uv run python scripts/plot_fig14_fig15.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
RUNS = ROOT / "runs"
OUT_DIR = ROOT / "experiments/kimetal2024/figure_replicas"
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _t_scales  # noqa: E402

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

_COL_T, _COL_POSY_LEFT = 1, 6
SETTLE_CYCLES = 5.0        # discard the soft start before measuring anything

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
ANGLES = [2, 3, 4, 5, 6, 7]
RPM_TEMPLATE = "fig14_l9_rpm{v:g}"
ANGLE_TEMPLATE = "fig14_l9_th{v:g}"
SPECTRA_RPM = [37.5, 32.5, 22.5]
SPECTRA_ANGLE = [5, 3, 2]


def elevation(run_id: str):
    """(t/T_p, left-end elevation / static elevation) over settled cycles."""
    run_dir = RUNS / run_id
    path = run_dir / "vol_frac_interf.dat"
    if not path.exists():
        return None
    data = np.loadtxt(path, skiprows=1)
    if data.ndim != 2 or data.shape[1] <= _COL_POSY_LEFT:
        return None                     # pre-2026-09-16 run: no end probes
    params = json.loads((run_dir / "params.json").read_text())
    _, T_per_nd = _t_scales(params)
    y_static = params["geometry"]["b"] / params["geometry"]["a"]   # half height
    t_tp = data[:, _COL_T] / T_per_nd
    y = data[:, _COL_POSY_LEFT]
    good = np.isfinite(y) & (t_tp >= SETTLE_CYCLES)
    if good.sum() < 32:
        return None
    return t_tp[good], (y[good] + y_static) / y_static


def sweep_panel(ax, values, template, colour, marker, xlabel):
    xs, ys = [], []
    for v in values:
        e = elevation(template.format(v=v))
        if e is None:
            continue
        xs.append(v)
        ys.append(float(np.max(e[1])))
    if xs:
        ax.plot(xs, ys, color=colour, marker=marker, ms=6, lw=1.3)
    ax.axhline(1.0, color="0.6", lw=1.0, zorder=0)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(r"$y_{left,max}/y_{static}$", fontsize=11)
    ax.set_xticks(values)
    ax.tick_params(which="both", direction="in")
    return len(xs)


def spectrum(t_tp, y):
    """(f/f_b, amplitude) of the elevation over a whole number of cycles."""
    n_whole = int(t_tp[-1] - t_tp[0])
    if n_whole < 2:
        return None
    m = t_tp <= t_tp[0] + n_whole
    t, s = t_tp[m], y[m] - y[m].mean()
    # Resample onto a uniform grid: the solver's output interval is not exactly
    # constant, and an FFT of unevenly sampled data smears every peak.
    n = 1 << int(math.floor(math.log2(len(t))))
    tu = np.linspace(t[0], t[-1], n)
    su = np.interp(tu, t, s)
    amp = np.abs(np.fft.rfft(su * np.hanning(n))) / n
    freq = np.fft.rfftfreq(n, d=(tu[1] - tu[0]))      # already in units of f_b
    return freq, amp


def main() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))
    n_rpm = sweep_panel(axes[0], RPMS, RPM_TEMPLATE, "#0072B2", "o",
                        r"Rocking frequency $f_b$ (rpm)")
    n_ang = sweep_panel(axes[1], ANGLES, ANGLE_TEMPLATE, "#D55E00", "^",
                        r"Rocking angle $\theta_{b,max}$ (deg)")
    for ax, letter in zip(axes, "ab"):
        ax.text(-0.14, 1.03, rf"$({letter})$", transform=ax.transAxes,
                fontsize=13, style="italic")
    fig.tight_layout()
    out14 = OUT_DIR / "replicated_Fig14.png"
    fig.savefig(out14, dpi=150, bbox_inches="tight")
    print(f"saved {out14}  ({n_rpm}/{len(RPMS)} rpm, {n_ang}/{len(ANGLES)} angle)")
    if n_rpm == 0 and n_ang == 0:
        print("  no run carries posY_left yet -- submit scripts/submit_fig14.py")

    fig, axes = plt.subplots(2, 3, figsize=(10.0, 5.2))
    panels = ([(v, RPM_TEMPLATE, f"{v:g} rpm") for v in SPECTRA_RPM]
              + [(v, ANGLE_TEMPLATE, rf"${v:g}^\circ$") for v in SPECTRA_ANGLE])
    drawn = 0
    for ax, (v, tmpl, label) in zip(axes.ravel(), panels):
        e = elevation(tmpl.format(v=v))
        sp = spectrum(*e) if e is not None else None
        if sp is None:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center",
                    fontsize=9)
        else:
            freq, amp = sp
            ax.semilogy(freq, np.clip(amp, 1e-12, None), color="0.2", lw=0.9)
            ax.set_xlim(0, 6)
            drawn += 1
        ax.set_title(label, fontsize=9, loc="left")
        ax.set_xlabel(r"$f/f_b$", fontsize=10)
        ax.set_ylabel("amplitude", fontsize=10)
        ax.tick_params(which="both", direction="in", labelsize=8)
    fig.tight_layout()
    out15 = OUT_DIR / "replicated_Fig15.png"
    fig.savefig(out15, dpi=150, bbox_inches="tight")
    print(f"saved {out15}  ({drawn}/{len(panels)} panels)")


if __name__ == "__main__":
    main()
