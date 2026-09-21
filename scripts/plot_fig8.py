"""Replica of Kim et al. (2024) Fig. 8 (tau_Ediss_evol).

Panels, following Kim's own axis limits throughout:
  (a1) spatial mean SIGNED shear stress <tau_w'> over the rocking phase,
       tau in [-3e-3, 3e-3] Pa.
  (a2) spatial mean EDR <eps_w'> over the rocking phase, [0, 0.4] W/m^3.
  (b)  histogram of signed tau across the liquid at the instant <tau_w'>
       peaks, over Kim's [-15e-3, 15e-3] Pa window, linear axes.
  (c)  histogram of EDR at the instant <eps_w'> peaks, [0, 1] W/m^3.

Kim's (a) is one twin-axis panel; it is split here because overlaying three
resolutions on twin axes was unreadable. (b) and (c) are L10 only, as in Kim
-- they are single-instant distributions, not a convergence statement.

Every value is recomputed from the saved fields with the BAG mask
(`f > 0.5` AND inside the embedded geometry), never from the on-disk KPIs:
runs predating the mask fix (diary.md 2026-09-15 (2)) logged spatial means
over the whole lower half-domain, which is 3.51x the bag's liquid area and
diluted every mean accordingly. Maxima were never affected, means always
were.

Panels a1/a2 fold on the TRUE forcing phase (absolute simulation time mod
T_p), not on each run's own first frame: the three runs sit at different
absolute times, and zeroing each at its window start imposed an arbitrary
per-level phase offset. That artifact -- not physics -- is what made L8 look
"totally off" when the split was first produced.

Data: L8 runs/l8_coldstart_vid, L9 runs/a34fc4d4, L10 runs/57f68830 (the
late-time probe, t/Tp 37->47). Frame counts differ by level because these
predate the off-period frame cadence (params_read.h, T_p/(N+0.618)); runs
on the old cadence sample only N distinct phases however long they run.

runs/l10_kim_fig8_signed must NOT be used again: `_binary: None` in its
params.json is validate_run.py's first hard-fail condition, and it predates
the H_bio nondim and tau-histogram OpenMP-race fixes.

Usage:  uv run python scripts/plot_fig8.py
"""
import json
import math
import struct
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
from pathlib import Path as _Path
sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from scripts import figstyle as fs  # noqa: E402

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
OUT_DIR = ROOT / "experiments/kimetal2024/figure_replicas"

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

# Kim's reported peaks at 32.5 rpm, theta=7 deg (csv_raw/shear_ediss_vs_frequency.csv)
KIM_TAU_PEAK = 1.8026e-3      # Pa
KIM_EDISS_PEAK = 0.28655      # W/m^3

RUNS = [("L8", "l8_coldstart_vid", fs.level_colour(8)),
        ("L9", "a34fc4d4", fs.level_colour(9)),
        ("L10", ("fig8_hist_l10"
                 if (ROOT / "runs" / "fig8_hist_l10" / "frames_tau").is_dir()
                 else "57f68830"), fs.level_colour(10))]
# L10 follows the same switch as HIST_RUN below. 57f68830 samples exactly five
# phases however long it runs (frame interval T_p/5 on the nose), which in
# panel (a) draws as five points joined by straight lines against L8's smooth
# 150-frame curve -- the sampling artifact, not the physics.
# Panels (b)/(c) need MANY distinct phases to locate the peak instant, and
# 57f68830 samples exactly five however long it runs (its frame interval is
# T_p/5 on the nose). runs/fig8_hist_l10 re-records the same converged state
# on the off-period cadence -- ~130 phases -- and is used as soon as it exists.
HIST_RUN = ("fig8_hist_l10"
            if (ROOT / "runs" / "fig8_hist_l10" / "frames_tau").is_dir()
            else "57f68830")


def load_frame(path):
    with open(path, "rb") as fh:
        (n,) = struct.unpack("i", fh.read(4))
        (t,) = struct.unpack("d", fh.read(8))
        fh.read(16)
        f = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
        tau = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
        ediss = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
    return t, f, tau, ediss


def scales(run_id):
    """(bag half-height /L, tau scale, EDR scale, T_bio, omega_b)."""
    p = json.load(open(ROOT / "runs" / run_id / "params.json"))
    L = p["geometry"]["a"]
    b = p["geometry"]["b"] / L
    tm = p["theta_max"]
    th = math.radians(tm[0] if isinstance(tm, list) else tm)
    omega_b = p["omega_b"]
    T_per = 2 * math.pi / omega_b
    H = 2 * L * b
    V = L / 4 * (H + 0.5 * L * math.tan(th))
    U = V / (H * 0.5) / T_per
    return b, 1000.0 * U**2, 1000.0 * U**3 / L, L / U, omega_b


def bag_mask(f, b):
    """Liquid INSIDE the embedded bag. `f` alone fills the whole lower half-domain."""
    n = f.shape[0]
    y = (-0.5 + (np.arange(n) + 0.5) / n)[:, None] * np.ones((1, n))
    return (f > 0.5) & (np.abs(y) < b)


# ── panels (a1), (a2): phase-folded spatial means at three resolutions ─────
series = {}
for lvl, run, col in RUNS:
    b, tau_scale, ediss_scale, T_bio, omega_b = scales(run)
    files = sorted((ROOT / "runs" / run / "frames_tau").glob("frame_*.bin"))
    files = files[len(files) // 2:]          # settled half only
    times, tau_mean, ediss_mean = [], [], []
    for fp in files:
        t, f, tau, ediss = load_frame(fp)
        m = bag_mask(f, b)
        if not m.any():
            continue
        times.append(t)
        tau_mean.append(tau[m].mean() * tau_scale)
        ediss_mean.append(ediss[m].mean() * ediss_scale)
    T_per_nd = (2 * math.pi / omega_b) / T_bio
    phase = (np.array(times) / T_per_nd) % 1.0
    o = np.argsort(phase)
    series[lvl] = dict(phase=phase[o], tau=np.array(tau_mean)[o],
                       ediss=np.array(ediss_mean)[o], col=col, n=len(o))
    print(f"{lvl}: {len(o)} frames, {len(set(np.round(phase, 6)))} distinct phases")


def phase_panel(key, fname, ylabel, ylim, kim_peak, symmetric):
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    # Mark size and line weight follow the SAMPLE DENSITY, because these runs
    # predate the off-period frame cadence and so differ wildly in it:
    # L8's frame spacing is 0.2006 T_p (slightly OFF commensurate), so its
    # phase creeps ~0.0006/cycle and 150 frames fill in 150 distinct phases --
    # it genuinely traces the curve. L9 (0.2000) pins to 16 phases and L10 to
    # just 5, however long they run. Drawn at one size, L8's 150 markers
    # clumped into drifting arcs and read as a heavy dashed rule rather than
    # as the densest, most informative series on the panel.
    for lvl, _, col in RUNS:
        d = series[lvl]
        dense = d["n"] > 60
        ax.plot(d["phase"], d[key], color=col,
                lw=1.0 if dense else 0.9,
                marker="o", ms=1.8 if dense else 4.2,
                alpha=0.55 if dense else 0.95, label=lvl, zorder=2 if dense else 3)
    # NO Kim series on this panel, and no horizontal reference line either.
    #
    # Kim's Fig 8(a) is a TIME SERIES over the rocking phase, but we do not
    # have it as data: Figures/Fig_tau_Ediss.pdf embeds the panel as a RASTER
    # image (one 2.7 MB image stream, zero vector path operators), and his
    # published dataset carries only the single PEAK value, not the curve.
    # Plotting that lone number as a rule across the panel was a poor stand-in
    # for a series -- it read as a heavy dashed line with no x-dependence,
    # unlike every other mark here. Kim's peak belongs in the caption until
    # the curve is pixel-digitised; then it can be drawn as markers joined by
    # a line, exactly like ours.
    ax.set_xlabel(r"rocking phase  $t/T_p$  (mod 1)", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(*ylim)
    ax.legend(fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    fig.savefig(OUT_DIR / fname, dpi=150, bbox_inches="tight")
    print("saved", fname)


phase_panel("tau", "replicated_Fig8_a1.png", r"$\langle\tau_w'\rangle$ (Pa)",
            (-3e-3, 3e-3), KIM_TAU_PEAK, symmetric=True)
phase_panel("ediss", "replicated_Fig8_a2.png", r"$\langle\epsilon_w'\rangle$ (W/m$^3$)",
            (0.0, 0.4), KIM_EDISS_PEAK, symmetric=False)


# ── panels (b), (c): field histograms at the peak instants ────────────────
b_hi, tau_scale, ediss_scale, _, _ = scales(HIST_RUN)
frames = [load_frame(p) for p in
          sorted((ROOT / "runs" / HIST_RUN / "frames_tau").glob("frame_*.bin"))]
frames = frames[len(frames) // 2:]

tau_means = [tau[bag_mask(f, b_hi)].mean() for _, f, tau, _ in frames]
ediss_means = [ed[bag_mask(f, b_hi)].mean() for _, f, _, ed in frames]
i_tau = int(np.argmax(np.abs(tau_means)))
i_ediss = int(np.argmax(ediss_means))

_, f_t, tau_field, _ = frames[i_tau]
_, f_e, _, ediss_field = frames[i_ediss]
tau_liquid = tau_field[bag_mask(f_t, b_hi)] * tau_scale
ediss_liquid = ediss_field[bag_mask(f_e, b_hi)] * ediss_scale
print(f"tau_liquid at peak frame: [{tau_liquid.min():.3e}, {tau_liquid.max():.3e}] Pa")
print(f"ediss_liquid at peak frame: [{ediss_liquid.min():.3e}, {ediss_liquid.max():.3e}] W/m3")


def histogram_panel(values, lo, hi, fname, xlabel, color, n_bins=30):
    """Bins span Kim's OWN window. Counts normalize by the TOTAL sample count,
    not the in-window subset, so bars are not inflated by excluding the tail;
    the tail's share is printed instead."""
    counts, edges = np.histogram(values, bins=np.linspace(lo, hi, n_bins + 1))
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    ax.bar(edges[:-1], counts / len(values), width=np.diff(edges), align="edge",
           color=color, edgecolor="none", alpha=0.85)
    ax.set_xlim(lo, hi)
    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel("Normalized frequency", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT_DIR / fname, dpi=150)
    outside = float(((values < lo) | (values > hi)).mean())
    print(f"saved {fname}  ({outside*100:.2f}% of samples outside the window)")


# One dataset, one quantity per panel: colour has no dimension left to
# encode here, so both histograms are neutral and the axis label alone says
# which quantity it is. Tinting them would imply a contrast that is not in
# the figure.
histogram_panel(tau_liquid, -15e-3, 15e-3, "replicated_Fig8_b.png",
                r"$\tau_w'$ (Pa)", fs.NEUTRAL)
histogram_panel(ediss_liquid, 0.0, 1.0, "replicated_Fig8_c.png",
                r"$\epsilon_w'$ (W/m$^3$)", fs.NEUTRAL)
