"""Replica of Kim et al. (2024) Fig. 8 (tau_Ediss_evol).

Panels, following Kim's own axis limits throughout:
  (a)  time evolution of the spatial-mean SIGNED shear stress <tau_w'> (left
       axis, +-3e-3 Pa) and EDR <eps_w'> (right axis, 0-0.4 W/m^3) over three
       unfolded settled cycles, from shear_stress.dat. Arrows (b) and (c) mark
       the instants the histograms are taken at.
  (b)  histogram of signed tau across the liquid at the (b) instant, the
       positive crest of <tau_w'>, over Kim's [-15e-3, 15e-3] Pa window.
  (c)  histogram of EDR at the (c) instant, the peak of <eps_w'>, [0, 1] W/m^3.
All at L10, as in Kim. Colour is the level (project convention), and line
style separates the two quantities.

(b) and (c) are recomputed from the saved fields with the BAG mask
(`f > 0.5` AND inside the embedded geometry), never from the on-disk KPIs:
runs predating the mask fix (diary.md 2026-09-15 (2)) logged spatial means
over the whole lower half-domain, which is 3.51x the bag's liquid area and
diluted every mean accordingly.

A phase-folded (t/T_p mod 1) level comparison used to be drawn here too. It
is not part of Kim's figure and was deleted on 2026-10-01.

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
from matplotlib.lines import Line2D
import sys
from pathlib import Path as _Path
sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from scripts import figstyle as fs  # noqa: E402
from scripts.postprocess import _t_scales  # noqa: E402

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
OUT_DIR = ROOT / "experiments/kimetal2024/figure_replicas"

plt.rcParams.update({
    "mathtext.fontset": "cm", "font.family": "serif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
})

# Kim's reported peaks at 32.5 rpm, theta=7 deg (csv_raw/shear_ediss_vs_frequency.csv)
KIM_TAU_PEAK = 1.8026e-3      # Pa
KIM_EDISS_PEAK = 0.28655      # W/m^3

def _l10_source() -> str:
    """Newest L10 recording on the off-period cadence.

    fig8_hist_l10b re-records for 20 cycles because fig8_hist_l10 ended with
    its restart transient still present (2x/1x harmonic of signed <tau> 0.23;
    half-period antisymmetry forbids it in a settled flow). l10b decays to
    0.05 by cycles 17-20. 57f68830 samples exactly five phases (T_p/5 on the
    nose) and is the last resort.
    """
    for run in ("fig8_hist_l10b", "fig8_hist_l10"):
        if (ROOT / "runs" / run / "frames_tau").is_dir():
            return run
    return "57f68830"


RUNS = [("L8", "l8_coldstart_vid", fs.level_colour(8)),
        ("L9", "a34fc4d4", fs.level_colour(9)),
        ("L10", _l10_source(), fs.level_colour(10))]
HIST_RUN = _l10_source()
# Fraction of HIST_RUN's frames kept for (b)/(c), from the END. l10b's even
# harmonic is <= 0.09 only over its last ~6 of 21 cycles.
HIST_TAIL = 0.25 if HIST_RUN == "fig8_hist_l10b" else 0.5


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


def _shear_series(run):
    """(t/T_p, <tau>, <eps>) from shear_stress.dat -- the DENSE series.

    Panel (a) used to be built from frames_tau, which is the wrong source
    twice over. Kim's panel is a time series of SPATIAL MEANS, and that is
    exactly what this file logs; the frames exist for panels (b)/(c), which
    need the per-cell fields. And the frames are sparse in time by design --
    each is ~10 MB, so the cadence is a storage decision -- giving 13.6
    samples per cycle against this file's 30.4, which is why the curve read
    as a polygon. On L8 the gap is far wider: 5 frames per cycle against
    1823 rows here.

    Columns are read by NAME. tau_kim_mean / ediss_kim_mean use Kim's exact
    mask and exist only from 4b3a435; older runs fall back to the project's
    own signed mean, which on fig8_hist_l10 agrees with the Kim-masked
    column to three digits (2x/1x of 0.221 vs 0.222).
    """
    d = ROOT / "runs" / run
    f = d / "shear_stress.dat"
    if not f.exists():
        return None
    head = f.open().readline().split()
    a = np.loadtxt(f, skiprows=1)
    if a.ndim != 2 or len(a) < 10:
        return None

    def col(*names):
        for n in names:
            if n in head:
                return a[:, head.index(n)]
        return None

    t = col("t")
    tau = col("tau_kim_mean", "tau_mean_signed")
    eps = col("ediss_kim_mean", "ediss_mean")
    if t is None or tau is None or eps is None:
        return None
    params = json.loads((d / "params.json").read_text())
    _, T_per_nd = _t_scales(params)
    # shear_stress.dat is NON-DIMENSIONAL (postprocess.py:719). The frames
    # path applies these same factors, so omitting them here silently
    # rescaled the panel by rho*U^2 -- caught only because the tau peak moved
    # from 2.2e-3 to 3e-4 Pa when the data source changed.
    _, tau_scale, ediss_scale, _, _ = scales(run)
    return t / T_per_nd, tau * tau_scale, eps * ediss_scale


def panel_window(lvl="L10", n_cycles=3.0):
    """The last n_cycles of the reference level's shear_stress.dat, in t/T_p."""
    got = _shear_series(dict((n, r) for n, r, _ in RUNS)[lvl])
    c = got[0]
    t1 = float(c.max())
    return max(float(c.min()), t1 - n_cycles), t1


def time_panel(lvl, t_b, t_c, n_cycles=3.0):
    """Kim's Fig 8(a): the TIME EVOLUTION, both quantities on one axes.

    His panel runs t/T_p = 80 to 83 -- three unfolded cycles -- with
    <tau'_w> on the left axis in blue and <eps'_w> on the right in red, and
    the instants where each average peaks labelled (b) and (c), which is what
    panels (b)/(c) are histograms of.

    Our earlier version folded the x axis mod 1, split the two quantities
    across separate panels, and coloured by refinement level. The fold and
    the split were both departures from the paper, and the level comparison
    is what Fig 13 is for.

    The cycle NUMBERS differ and are not forced to match. fig8_hist_l10
    continues a converged state at t/T_p = 47, not 80; for a settled periodic
    flow any three consecutive settled cycles are the same three cycles, and
    relabelling them 80-83 would assert a correspondence with Kim's tracer
    protocol that this run does not have.

    COLOUR STAYS ON THE LEVEL. Kim draws one dataset and spends colour on the
    quantity (blue tau, red EDR), and copying that here would break the
    project convention that colour means refinement level -- an established
    convention worth more than matching one panel's palette, because a reader
    carries it across all eighteen figures. So every level present is drawn
    in its own hue and the QUANTITY is separated by the axis it belongs to
    and by line style: solid on the left axis is tau, dashed on the right is
    EDR. The axis labels are tinted neutral rather than per-quantity for the
    same reason.
    """
    dense = {}
    for name, run, col in RUNS:
        got = _shear_series(run)
        if got is not None:
            dense[name] = (got, col)
    if lvl not in dense:
        print(f"time_panel: no shear_stress.dat for {lvl}")
        return
    # Window set by the reference level, so every series shows the SAME three
    # cycles of its own record rather than three arbitrary different ones.
    t0, t1 = panel_window(lvl, n_cycles)

    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    ax2 = ax.twinx()
    shown = []
    for name, _, _ in RUNS:
        if name not in dense:
            continue
        (c, tau, eps), col = dense[name]
        m = (c >= t0) & (c <= t1)
        if m.sum() < 20:
            continue
        ax.plot(c[m], tau[m], color=col, lw=1.4, ls="-")
        ax2.plot(c[m], eps[m], color=col, lw=1.1, ls="--")
        shown.append((name, col))
        if name == lvl:
            # (b) and (c) mark the EXACT field frames the histograms are taken
            # from (t_b, t_c), not peaks of this curve, so the labels and the
            # distributions refer to the same instants. Arrows as in Kim's panel.
            for lab, tx, axis, y in (("(b)", t_b, ax, tau), ("(c)", t_c, ax2, eps)):
                yv = float(np.interp(tx, c, y))
                axis.annotate(lab, xy=(tx, yv), xytext=(0, 22),
                              textcoords="offset points", color=col,
                              fontsize=10, ha="center",
                              arrowprops=dict(arrowstyle="->", color=col, lw=0.9))
    if not shown:
        plt.close(fig)
        print("time_panel: no level had frames in the window")
        return

    (rc, rtau, reps), _ = dense[lvl]
    m = (rc >= t0) & (rc <= t1)
    tau, eps = rtau, reps
    ax.set_xlabel(r"$t/T_p$", fontsize=12)
    ax.set_ylabel(r"$\langle\tau'_w\rangle$ (Pa)   (solid)", fontsize=11)
    ax2.set_ylabel(r"$\langle\epsilon'_w\rangle$ (W/m$^3$)   (dashed)",
                   fontsize=11)
    # Only levels actually drawn in this window. A level whose recording
    # covers different cycles has no curve here, and listing it would promise
    # a series the reader then hunts for.
    lv = [Line2D([], [], color=c, lw=1.4, label=n) for n, c in shown]
    lg = ax.legend(handles=lv, fontsize=8, frameon=False, loc="upper left",
                   bbox_to_anchor=(1.10, 1.0), title="level")
    lg.get_title().set_fontsize(8)
    ax.set_xlim(t0, t1)
    # Kim's own limits, so the two panels can be laid side by side:
    # +-3e-3 Pa (scaled x10^-3, as he draws it) and 0-0.4 W/m^3.
    ax.set_ylim(-3e-3, 3e-3)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, -3), useMathText=True)
    ax2.set_ylim(0.0, 0.4)
    ax.tick_params(which="both", direction="in")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "replicated_Fig08_a.png", dpi=150,
                bbox_inches="tight")
    print(f"saved replicated_Fig08_a.png  ({lvl}, t/T_p {t0:.2f}-{t1:.2f}, "
          f"{m.sum()} samples, {m.sum()/n_cycles:.1f}/cycle)")


# Panel (a) as Kim draws it: one axes, three unfolded cycles, both
# quantities. It is drawn after the histogram frames are chosen (below), so
# its (b)/(c) arrows point at those exact instants.


# ── panels (b), (c): field histograms at the peak instants ────────────────
b_hi, tau_scale, ediss_scale, _, _ = scales(HIST_RUN)
frames = [load_frame(p) for p in
          sorted((ROOT / "runs" / HIST_RUN / "frames_tau").glob("frame_*.bin"))]
# Only frames inside panel (a)'s window, so (b) and (c) are instants the reader
# can see there. The window is the run's last 3 cycles, inside the settled tail
# that HIST_TAIL was chosen for.
_, _, _, T_bio_h, omega_h = scales(HIST_RUN)
T_per_h = (2 * math.pi / omega_h) / T_bio_h
w0, w1 = panel_window("L10")
frames = [fr for fr in frames[int(len(frames) * (1 - HIST_TAIL)):]
          if w0 <= fr[0] / T_per_h <= w1]

tau_means = [tau[bag_mask(f, b_hi)].mean() for _, f, tau, _ in frames]
ediss_means = [ed[bag_mask(f, b_hi)].mean() for _, f, _, ed in frames]
# Kim's (b) is the POSITIVE peak of the signed mean (his arrow sits on a crest).
i_tau = int(np.argmax(tau_means))
i_ediss = int(np.argmax(ediss_means))
time_panel("L10", frames[i_tau][0] / T_per_h, frames[i_ediss][0] / T_per_h)

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
histogram_panel(tau_liquid, -15e-3, 15e-3, "replicated_Fig08_b.png",
                r"$\tau_w'$ (Pa)", fs.NEUTRAL)
histogram_panel(ediss_liquid, 0.0, 1.0, "replicated_Fig08_c.png",
                r"$\epsilon_w'$ (W/m$^3$)", fs.NEUTRAL)
