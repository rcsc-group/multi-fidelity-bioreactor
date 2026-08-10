"""Replica of Kim et al. (2024) Fig. 8 (tau_Ediss_evol), v2 -- corrected
after actually viewing docs/kimetal2024/Figures/Fig_tau_Ediss.pdf:

  (a) domain-mean SIGNED shear stress <tau_w'> (blue) and EDR <eps_w'>
      (red), dual y-axis, LINEAR scale, matching Kim's style. v1 plotted
      mean(|tau|) (always >=0) instead of Kim's actual plotted quantity,
      the signed mean (oscillates through zero) -- a real bug, not a
      styling choice, caught by comparing against Kim's real figure.
  (b) histogram of the SIGNED shear stress across the liquid at the
      instant <tau_w'> peaks. LINEAR axes (Kim's own figure uses linear,
      not log) -- v1 used log-log, which was a deviation introduced
      without first checking Kim's actual convention.
  (c) histogram of EDR (non-negative by construction, no sign issue)
      across the liquid at the instant <eps_w'> peaks. LINEAR axes.

Data: runs/l10_kim_fig8_signed (L10, warm-restarted from l10_kim_seg2's
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

RUN_DIR = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/l10_kim_fig8_signed")
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

# ── panel (a): time series, SIGNED tau_mean this time ──────────────────────
d = np.loadtxt(RUN_DIR / "shear_stress.dat", skiprows=1)
t_nd = d[:, 1]
tau_mean_signed_nd = d[:, 10]   # new column
ediss_mean_nd = d[:, 9]
tau_mean_pa = tau_mean_signed_nd * tau_scale
ediss_mean_wm3 = ediss_mean_nd * ediss_scale
t_over_Tp = t_nd / T_per_nd

# peak = max |signed mean| (matches "amplitude" reading off Kim's fig,
# not just the positive-going max) for picking the panel-(b) snapshot instant
i_tau_peak = int(np.argmax(np.abs(tau_mean_pa)))
i_ediss_peak = int(np.argmax(ediss_mean_wm3))
t_tau_peak = t_nd[i_tau_peak]
t_ediss_peak = t_nd[i_ediss_peak]
print(f"tau_mean_signed |peak| at t={t_tau_peak:.4f} (t/Tp={t_over_Tp[i_tau_peak]:.3f}), {tau_mean_pa[i_tau_peak]:.6f} Pa")
print(f"ediss_mean peak at t={t_ediss_peak:.4f} (t/Tp={t_over_Tp[i_ediss_peak]:.3f}), {ediss_mean_wm3[i_ediss_peak]:.4f} W/m3")
print(f"tau_mean_signed range: [{tau_mean_pa.min():.6f}, {tau_mean_pa.max():.6f}] Pa  (Kim: [-1.8e-3, +1.8e-3])")
print(f"ediss_mean range: [{ediss_mean_wm3.min():.4f}, {ediss_mean_wm3.max():.4f}] W/m3  (Kim: [0, ~0.35])")

fig, ax1 = plt.subplots(figsize=(6.0, 4.2))
l1, = ax1.plot(t_over_Tp, tau_mean_pa, color="royalblue", lw=1.6, label=r"$\langle\tau_w'\rangle$")
ax1.axhline(0, color="gray", lw=0.6, ls="-")
ax1.set_xlabel(r"$t/T_p$", fontsize=13)
ax1.set_ylabel(r"$\langle\tau_w'\rangle$ (Pa)", fontsize=13, color="royalblue")
ax1.set_ylim(-3e-3, 3e-3)
ax1.tick_params(axis="y", labelcolor="royalblue")

ax2 = ax1.twinx()
l2, = ax2.plot(t_over_Tp, ediss_mean_wm3, color="firebrick", lw=1.6, label=r"$\langle\epsilon_w'\rangle$")
ax2.set_ylim(0.0, 0.4)
ax2.set_ylabel(r"$\langle\epsilon_w'\rangle$ (W/m$^3$)", fontsize=13, color="firebrick")
ax2.tick_params(axis="y", labelcolor="firebrick")

ax1.legend(handles=[l1, l2], fontsize=10, loc="upper right")
ax1.text(-0.16, 1.02, r'$(a)$', transform=ax1.transAxes, fontsize=16, style='italic')
fig.tight_layout()
fig.savefig(OUT_DIR / "fig8_a_v2.png", dpi=150)
print("saved fig8_a_v2.png")


# ── panels (b),(c): field histograms at the peak instants, LINEAR axes ────
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
frame_times, frames = [], []
for p in frame_files:
    t, f, tau, ediss = load_frame(p)
    frame_times.append(t)
    frames.append((f, tau, ediss))
frame_times = np.array(frame_times)

idx_tau_frame = int(np.argmin(np.abs(frame_times - t_tau_peak)))
idx_ediss_frame = int(np.argmin(np.abs(frame_times - t_ediss_peak)))
print(f"nearest video frame to tau peak: t={frame_times[idx_tau_frame]:.4f} (target {t_tau_peak:.4f})")
print(f"nearest video frame to ediss peak: t={frame_times[idx_ediss_frame]:.4f} (target {t_ediss_peak:.4f})")

f_tau, tau_field, _ = frames[idx_tau_frame]      # tau_field is now SIGNED (2026-08-09 fix)
f_ediss, _, ediss_field = frames[idx_ediss_frame]

tau_liquid = tau_field[f_tau > 0.5] * tau_scale
ediss_liquid = ediss_field[f_ediss > 0.5] * ediss_scale
print(f"tau_liquid (signed) range at peak frame: [{tau_liquid.min():.6f}, {tau_liquid.max():.6f}] Pa")
print(f"ediss_liquid range at peak frame: [{ediss_liquid.min():.6f}, {ediss_liquid.max():.6f}] W/m3")

# Bins span Kim's OWN window, not our full data range -- binning over our
# much wider range (tau tail reaches +-0.25 Pa) and then zooming into his
# +-15e-3 window would leave only ~3 giant bins visible instead of a real
# distribution shape. Counts still normalize by the TOTAL sample count
# (not just the in-window subset), so the visible bars are not inflated by
# excluding the tail -- the tail's share is reported separately below.
def _hist_in_window(values, lo, hi, n_bins=30):
    bins = np.linspace(lo, hi, n_bins + 1)
    counts, edges = np.histogram(values, bins=bins)
    counts_norm = counts / len(values)  # normalize by ALL samples, not just in-window
    return counts_norm, edges

fig, ax = plt.subplots(figsize=(5.0, 4.0))
counts_norm, edges = _hist_in_window(tau_liquid, -15e-3, 15e-3)
ax.bar(edges[:-1], counts_norm, width=np.diff(edges), align="edge",
       color="royalblue", edgecolor="none", alpha=0.85)
ax.set_xlim(-15e-3, 15e-3)
ax.set_xlabel(r"$\tau_w'$ (Pa)", fontsize=13)
ax.set_ylabel("Normalized frequency", fontsize=13)
ax.text(-0.16, 1.02, r'$(b)$', transform=ax.transAxes, fontsize=16, style='italic')
fig.tight_layout()
fig.savefig(OUT_DIR / "fig8_b_v2.png", dpi=150)
print("saved fig8_b_v2.png")
frac_outside_b = float(((tau_liquid < -15e-3) | (tau_liquid > 15e-3)).mean())
print(f"fraction of tau_liquid outside Kim's [-15e-3,15e-3] window: {frac_outside_b*100:.2f}%")

fig, ax = plt.subplots(figsize=(5.0, 4.0))
counts_norm, edges = _hist_in_window(ediss_liquid, 0.0, 1.0)
ax.bar(edges[:-1], counts_norm, width=np.diff(edges), align="edge",
       color="firebrick", edgecolor="none", alpha=0.85)
ax.set_xlim(0.0, 1.0)
ax.set_xlabel(r"$\epsilon_w'$ (W/m$^3$)", fontsize=13)
ax.set_ylabel("Normalized frequency", fontsize=13)
ax.text(-0.16, 1.02, r'$(c)$', transform=ax.transAxes, fontsize=16, style='italic')
fig.tight_layout()
fig.savefig(OUT_DIR / "fig8_c_v2.png", dpi=150)
frac_outside_c = float((ediss_liquid > 1.0).mean())
print(f"fraction of ediss_liquid above Kim's 1.0 W/m3 window: {frac_outside_c*100:.2f}%")
print("saved fig8_c_v2.png")
