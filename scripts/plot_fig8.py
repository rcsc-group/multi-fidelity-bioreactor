"""Replica of Kim et al. (2024) Fig. 8 (tau_Ediss_evol):
(a) time evolution of spatially-averaged shear stress (blue) and energy
    dissipation rate (red), dual y-axis, theta=7deg f_b=32.5rpm baseline.
(b) normalized histogram of shear stress across the liquid, at the instant
    tau_mean(t) peaks.
(c) normalized histogram of EDR across the liquid, at the instant
    ediss_mean(t) peaks.

Data: runs/l10_kim_fig8 (L10, warm-restarted from l10_kim_seg2's
checkpoint at t=20.65, +1.8 nondim time / ~3 rocking periods).
"""
import json
import math
import struct
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

mpl_rc = plt.rcParams
mpl_rc['mathtext.fontset'] = 'cm'
mpl_rc['font.family'] = 'serif'
mpl_rc['axes.linewidth'] = 1.4
mpl_rc['xtick.direction'] = 'in'
mpl_rc['ytick.direction'] = 'in'

RUN_DIR = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/l10_kim_fig8")
OUT_DIR = Path("/oscar/scratch/eaguerov/tmp/fig8")

params = json.load(open(RUN_DIR / "params.json"))
L = params["geometry"]["a"]; H = 2 * params["geometry"]["b"]
th = math.radians(params["theta_max"][0])
omega_b = params["omega_b"]
T_per = 2 * math.pi / omega_b
V = L / 4 * (H + 0.5 * L * math.tan(th))
U_bio = V / (H * 0.5) / T_per
T_bio = L / U_bio
T_per_nd = T_per / T_bio
rho_w = 1000.0
tau_scale = rho_w * U_bio**2
ediss_scale = rho_w * U_bio**3 / L

# ── panel (a): time series ────────────────────────────────────────────────
d = np.loadtxt(RUN_DIR / "shear_stress.dat", skiprows=1)
t_nd = d[:, 1]
tau_mean_nd = d[:, 5]
ediss_mean_nd = d[:, 9]
tau_mean_pa = tau_mean_nd * tau_scale
ediss_mean_wm3 = ediss_mean_nd * ediss_scale
t_over_Tp = t_nd / T_per_nd

i_tau_peak = int(np.argmax(tau_mean_nd))
i_ediss_peak = int(np.argmax(ediss_mean_nd))
t_tau_peak = t_nd[i_tau_peak]
t_ediss_peak = t_nd[i_ediss_peak]
print(f"tau_mean peak at t={t_tau_peak:.4f} (t/Tp={t_over_Tp[i_tau_peak]:.3f}), {tau_mean_pa[i_tau_peak]:.6f} Pa")
print(f"ediss_mean peak at t={t_ediss_peak:.4f} (t/Tp={t_over_Tp[i_ediss_peak]:.3f}), {ediss_mean_wm3[i_ediss_peak]:.4f} W/m3")

fig, ax1 = plt.subplots(figsize=(6.0, 4.2))
l1, = ax1.plot(t_over_Tp, tau_mean_pa, color="royalblue", lw=1.6, label=r"$\langle\tau_w'\rangle$")
ax1.axvline(t_over_Tp[i_tau_peak], color="royalblue", lw=0.8, ls=":")
ax1.set_xlabel(r"$t/T_p$", fontsize=13)
ax1.set_ylabel(r"$\langle\tau_w'\rangle$ (Pa)", fontsize=13, color="royalblue")
ax1.tick_params(axis="y", labelcolor="royalblue")

ax2 = ax1.twinx()
l2, = ax2.plot(t_over_Tp, ediss_mean_wm3, color="firebrick", lw=1.6, label=r"$\langle\epsilon_w'\rangle$")
ax2.axvline(t_over_Tp[i_ediss_peak], color="firebrick", lw=0.8, ls=":")
ax2.set_ylabel(r"$\langle\epsilon_w'\rangle$ (W/m$^3$)", fontsize=13, color="firebrick")
ax2.tick_params(axis="y", labelcolor="firebrick")

ax1.legend(handles=[l1, l2], fontsize=10, loc="upper right")
ax1.text(-0.14, 1.02, r'$(a)$', transform=ax1.transAxes, fontsize=16, style='italic')
fig.tight_layout()
fig.savefig(OUT_DIR / "fig8_a.png", dpi=150)
print("saved fig8_a.png")


# ── panels (b),(c): field histograms at the peak instants ──────────────────
def load_frame(path):
    with open(path, "rb") as fh:
        (n,) = struct.unpack("i", fh.read(4))
        (t,) = struct.unpack("d", fh.read(8))
        fh.read(16)
        f = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
        tau = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
        ediss = np.frombuffer(fh.read(n * n * 4), dtype=np.float32).reshape(n, n)
    return t, f, tau, ediss


frame_files = sorted((RUN_DIR / "frames_tau").glob("frame_*.bin"))
frame_times = []
frames = []
for p in frame_files:
    t, f, tau, ediss = load_frame(p)
    frame_times.append(t)
    frames.append((f, tau, ediss))
frame_times = np.array(frame_times)

idx_tau_frame = int(np.argmin(np.abs(frame_times - t_tau_peak)))
idx_ediss_frame = int(np.argmin(np.abs(frame_times - t_ediss_peak)))
print(f"nearest video frame to tau peak: t={frame_times[idx_tau_frame]:.4f} (target {t_tau_peak:.4f})")
print(f"nearest video frame to ediss peak: t={frame_times[idx_ediss_frame]:.4f} (target {t_ediss_peak:.4f})")

f_tau, tau_field, _ = frames[idx_tau_frame]
f_ediss, _, ediss_field = frames[idx_ediss_frame]

tau_liquid = tau_field[f_tau > 0.5] * tau_scale
ediss_liquid = ediss_field[f_ediss > 0.5] * ediss_scale

def _plot_loglog_hist(values, color, xlabel, panel_label, out_path):
    """Log-spaced bins + log-y counts: these tau/EDR distributions are
    heavily right-skewed (near-zero through the laminar bulk, a thin tail
    near the wall/interface -- same shape found earlier investigating the
    B.17 kLa fits and the tau_max non-convergence), so a linear histogram
    is >95% one bin and hides the entire tail. Floor at the smallest
    positive value (not 0) so log-spaced bins are well-defined."""
    positive = values[values > 0]
    floor = np.percentile(positive, 0.01)
    bins = np.logspace(np.log10(floor), np.log10(values.max()), 60)
    counts, edges = np.histogram(values, bins=bins)
    counts_norm = counts / counts.sum()
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    ax.bar(edges[:-1], counts_norm, width=np.diff(edges), align="edge",
           color=color, edgecolor="none", alpha=0.85)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(bottom=1e-6)
    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel("Normalized frequency", fontsize=13)
    ax.text(-0.16, 1.02, panel_label, transform=ax.transAxes, fontsize=16, style='italic')
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"saved {out_path.name}")


_plot_loglog_hist(tau_liquid, "royalblue", r"$\tau_w'$ (Pa)", r'$(b)$', OUT_DIR / "fig8_b.png")
_plot_loglog_hist(ediss_liquid, "firebrick", r"$\epsilon_w'$ (W/m$^3$)", r'$(c)$', OUT_DIR / "fig8_c.png")
