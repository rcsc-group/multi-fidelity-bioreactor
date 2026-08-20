"""Ours-vs-upstream velocity/shear-stress comparison stats across all 13
matching L10 snapshots, using the CORRECTED liquid mask (diary.md,
2026-08-19 "MAJOR CORRECTION" entry): the naive `f>0.5` mask includes
~71% dead solid-region cells (frozen VOF artifact from the half-space
fill IC, which ignores the embedded boundary). Corrected mask:
  - upstream: f>0.5 AND cs>0.5 (cs = upstream's real solid-indicator
    column, index 6 in Data_all_*.txt)
  - ours: f>0.5 AND fabs(y)<b_nd (analytic reconstruction -- our fork's
    periodic dump doesn't carry a cs column). b_nd=0.143 is the only
    binding embed constraint (a_nd=1 exceeds the domain half-width 0.5,
    so it never binds; verified against upstream's real cs at t=12.7596
    to within 0.02% cell-count agreement).

Zero new compute -- reuses fork_l10_periodic (job 5073228) and upstream
(job 4961226) data, same as the single-instant check this redoes.

Usage:
    uv run python scripts/compute_us_vs_upstream_stats_corrected.py
"""
import glob
import math

import numpy as np

N = 1024
B_ND = 0.143
TIMES = [0, 1.0633, 2.1266, 3.1899, 4.2532, 5.3165, 6.3798, 7.4431,
         8.5064, 9.5697, 10.633, 11.6963, 12.7596]
RAMP_DONE_T = 11.6132  # upstream's own ramp (t_change=30s literal) completes at t_change_st = 30/T_bio = 11.6132
                        # (recomputed directly from T_bio=L_bio/U_bio; corrects an earlier ~9.869 estimate).
                        # NOTE: our fork's N_RAMP_CYCLES ramp is INERT for this run (theta_max_prev==theta_max
                        # in params.json -> Ak interpolation is a no-op) -- our fork runs at FULL Th_max
                        # amplitude from t=0, zero ramp. Only upstream ramps in this comparison.

OURS_GLOB = "/oscar/scratch/eaguerov/tmp/fork_l10_periodic/rundir/DataOurs_1024_{t:g}_*.txt"
UPSTREAM_GLOB = "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/Data_all/Data_all_1024_{t:g}_*.txt"

rho_w, mu_w = 1.0e3, 1.0e-3
mu_a = 1.81e-5
L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
a_geom, b_geom = 0.25, 0.03575
Ly = b_geom / L_bio
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
H_bio = 2. * L_bio * Ly  # [FIX 2026-08-20] was missing the factor of 2 the production driver has had since 2026-08-03 (diary.md) -- Ly is a HALF-height ratio, H_bio must be the FULL bag height
V_bio = L_bio / 4 * (H_bio + 0.5 * L_bio * math.tan(Th_max))
U_bio = V_bio / (H_bio * 0.5) / T_per
Re_w = rho_w * U_bio * L_bio / mu_w
mur = mu_a / mu_w
mu1 = 1.0 / Re_w
mu2 = mur * mu1


def mu_of_f(f):
    fc = np.clip(f, 0, 1)
    return 1.0 / (fc * (1.0 / mu1 - 1.0 / mu2) + 1.0 / mu2)


def load(glob_pattern, ncols):
    files = sorted(glob.glob(glob_pattern))
    if not files:
        return None
    chunks = []
    for fp in files:
        try:
            arr = np.loadtxt(fp)
        except ValueError:
            arr = np.loadtxt(fp, skiprows=1)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        chunks.append(arr[:, :ncols])
    return np.vstack(chunks)


def to_grid(data, ncols, n=N):
    dx = 1.0 / n
    ix = np.clip(np.round((data[:, 0] + 0.5 - dx / 2) / dx).astype(int), 0, n - 1)
    iy = np.clip(np.round((data[:, 1] + 0.5 - dx / 2) / dx).astype(int), 0, n - 1)
    grids = [np.full((n, n), np.nan) for _ in range(ncols - 2)]
    for k, g in enumerate(grids):
        g[iy, ix] = data[:, 2 + k]
    return grids


def fields_and_mask(glob_pattern, ncols, is_upstream):
    data = load(glob_pattern, ncols)
    if data is None:
        return None
    grids = to_grid(data, ncols)
    ux, uy, f = grids[0], grids[1], grids[2]
    dx = 1.0 / N
    du_dy = np.full((N, N), np.nan)
    dv_dx = np.full((N, N), np.nan)
    du_dy[:, 1:-1] = (ux[:, 2:] - ux[:, :-2]) / (2 * dx)
    dv_dx[1:-1, :] = (uy[2:, :] - uy[:-2, :]) / (2 * dx)
    tau = mu_of_f(f) * (du_dy + dv_dx)
    speed = np.sqrt(ux ** 2 + uy ** 2)

    if is_upstream:
        # Data_all columns: x y ux uy vol_frac tracer solid oxygen ... -->
        # solid (cs) is column index 6 = grids[4] (grids[k] = data[:, 2+k]).
        cs = grids[4]
        in_bag = (f > 0.5) & (cs > 0.5)
    else:
        yy = (np.arange(N) + 0.5) / N - 0.5
        y2d = np.tile(yy.reshape(-1, 1), (1, N))
        in_bag = (f > 0.5) & (np.abs(y2d) < B_ND)

    interior = np.zeros((N, N), dtype=bool)
    interior[:, 1:-1] = True
    interior[1:-1, :] &= True
    mask = in_bag & interior & ~np.isnan(tau)
    return speed, tau, mask


results = []
for t in TIMES:
    ours = fields_and_mask(OURS_GLOB.format(t=t), 5, is_upstream=False)
    up = fields_and_mask(UPSTREAM_GLOB.format(t=t), 12, is_upstream=True)
    if ours is None or up is None:
        print(f"t={t:8.4f}  MISSING DATA")
        continue
    s1, tau1, m1 = ours
    s2, tau2, m2 = up
    valid = m1 & m2
    n_valid = valid.sum()
    if n_valid == 0:
        print(f"t={t:8.4f}  NO OVERLAPPING VALID CELLS")
        continue

    speed_scale = np.mean(s2[valid])
    speed_relerr = np.mean(np.abs(s1[valid] - s2[valid])) / speed_scale

    t1, t2 = tau1[valid], tau2[valid]
    corr = np.corrcoef(t1, t2)[0, 1]
    sign_agree = np.mean(np.sign(t1) == np.sign(t2))
    tau_scale = np.mean(np.abs(t2))
    tau_bulk_mean_relerr = (np.mean(np.abs(t1)) - np.mean(np.abs(t2))) / tau_scale

    flag = "" if t >= RAMP_DONE_T else "  [ramp-confounded]"
    results.append(dict(t=t, n_valid=int(n_valid), speed_relerr=speed_relerr,
                         corr=corr, sign_agree=sign_agree,
                         tau_bulk_mean_relerr=tau_bulk_mean_relerr))
    print(f"t={t:8.4f}  n_valid={n_valid:7d}  "
          f"speed_relerr={speed_relerr*100:6.2f}%  "
          f"tau_corr={corr:+.4f}  "
          f"tau_sign_agree={sign_agree*100:5.1f}%  "
          f"tau_bulk_mean_relerr={tau_bulk_mean_relerr*100:+8.1f}%{flag}")

settled = [r for r in results if r["t"] >= RAMP_DONE_T]
print()
print(f"Settled snapshots only (t >= {RAMP_DONE_T}, n={len(settled)}):")
print(f"  mean speed_relerr    = {np.mean([r['speed_relerr'] for r in settled])*100:.2f}%")
print(f"  mean tau_corr        = {np.mean([r['corr'] for r in settled]):+.4f}")
print(f"  mean tau_sign_agree  = {np.mean([r['sign_agree'] for r in settled])*100:.1f}%")
