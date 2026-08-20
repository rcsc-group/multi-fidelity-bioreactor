"""Direct check for a phase mismatch between the restart-recovered run
and the fresh baseline (user's concern, diary.md 2026-08-20).

The forcing is Th_max*sin(w_bio_st*t) as a function of ABSOLUTE
simulation time t -- not time since either run's own ramp start. Two
ramps that start at different absolute times (baseline: t=0; restart:
t=t_checkpoint~11.46) are offset by (t_checkpoint mod T_per_nd) in
forcing phase PURELY from when they started, before any flow dynamics
enter into it at all. So the right comparison folds each run's tail by
ABSOLUTE t mod T_per_nd (matching compare_restart_recovery.py's
convention) -- comparing by "time since ramp start" instead (as an
earlier version of this script did) bakes in that startup-time offset
and looks like a dynamical phase lag when it may just be arithmetic.

Bins BOTH runs' settled tails by absolute-t phase and compares the
PER-BIN MEAN (not a raw sample-to-sample cross-correlation, which is
dominated by tau_95's considerable sample-level noise -- see the
compare_restart_recovery.py scatter, which already showed a wide
vertical spread at every phase for both runs).

Usage:
    uv run python scripts/check_restart_recovery_phase_lag.py
"""
import numpy as np

RESTART_DIR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/c7e9eca7"
BASELINE_DIR = "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/runs/d054ff02"
T_PER_ND = 0.6073  # theta=7deg, 32.5rpm (corrected H_bio -- diary.md 2026-08-20 (6))
N_BINS = 16


def load_tail(run_dir, t_ramp_start, n_cycles_tail=8):
    t, tau95 = np.loadtxt(f"{run_dir}/shear_stress.dat", skiprows=1, usecols=(1, 2), unpack=True)
    t_end = t[-1]
    mask = t > t_end - n_cycles_tail * T_PER_ND
    return t[mask], tau95[mask]


t_r, tau_r = load_tail(RESTART_DIR, t_ramp_start=None)
t_b, tau_b = load_tail(BASELINE_DIR, t_ramp_start=None)

phase_r = np.mod(t_r, T_PER_ND) / T_PER_ND
phase_b = np.mod(t_b, T_PER_ND) / T_PER_ND

bins = np.linspace(0, 1, N_BINS + 1)
bin_idx_r = np.clip(np.digitize(phase_r, bins) - 1, 0, N_BINS - 1)
bin_idx_b = np.clip(np.digitize(phase_b, bins) - 1, 0, N_BINS - 1)

mean_r = np.array([tau_r[bin_idx_r == k].mean() if np.any(bin_idx_r == k) else np.nan for k in range(N_BINS)])
mean_b = np.array([tau_b[bin_idx_b == k].mean() if np.any(bin_idx_b == k) else np.nan for k in range(N_BINS)])
n_r = np.array([np.sum(bin_idx_r == k) for k in range(N_BINS)])
n_b = np.array([np.sum(bin_idx_b == k) for k in range(N_BINS)])

valid = ~np.isnan(mean_r) & ~np.isnan(mean_b)
corr = np.corrcoef(mean_r[valid], mean_b[valid])[0, 1]
rel_diff = np.abs(mean_r[valid] - mean_b[valid]) / mean_b[valid].mean()

print(f"phase-binned (N={N_BINS}) tau_95 mean profiles:")
for k in range(N_BINS):
    print(f"  bin {k:2d} [{bins[k]:.2f}-{bins[k+1]:.2f}]: restart={mean_r[k]:.5g} (n={n_r[k]:3d})  "
          f"baseline={mean_b[k]:.5g} (n={n_b[k]:3d})")
print(f"correlation between the two phase-binned mean profiles: {corr:.4f}")
print(f"mean relative diff between binned means: {rel_diff.mean()*100:.1f}%")

# Also check: does a small lag (in bin-index/phase units) improve the match?
best = (0, corr)
for shift in range(1, N_BINS // 2):
    shifted = np.roll(mean_r, shift)
    c = np.corrcoef(shifted[valid], mean_b[valid])[0, 1]
    if c > best[1]:
        best = (shift, c)
    shifted = np.roll(mean_r, -shift)
    c = np.corrcoef(shifted[valid], mean_b[valid])[0, 1]
    if c > best[1]:
        best = (-shift, c)
print(f"best phase-bin shift: {best[0]} bins ({best[0]/N_BINS:.3f} cycles), correlation there: {best[1]:.4f}")
