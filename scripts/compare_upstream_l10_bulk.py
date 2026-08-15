"""Compare the raw bulk field (u.x, u.y, f) between Kim et al.'s own upstream
driver and our fork, both run cold-start to the same nondimensional time
(~8 rocking cycles, t=4.85863) at matching L10 resolution.

Diagnostic script for diary.md 2026-08-11/12 -- answers "does the bulk look
identical to ours" directly, sidestepping the postprocessing-formula
questions raised by the abandoned metric-correction check.

Usage:
    uv run python scripts/compare_upstream_l10_bulk.py
"""
import glob
import sys

import numpy as np

UPSTREAM_GLOB = "/oscar/scratch/eaguerov/tmp/upstream_l10/rundir/DumpEarly_1024_12.14*_*.txt"
FORK_GLOB = "/oscar/scratch/eaguerov/tmp/fork_l10_coldstart/rundir/DumpEarlyFork_1024_12.14*_*.txt"


def load_dump(glob_pattern):
    """Load all per-rank dump files matching glob_pattern.

    Columns: x y ux uy f (6th column ignored -- semantics differ between
    upstream's tracer field and our fork's cs; not used in this comparison).
    """
    files = sorted(glob.glob(glob_pattern))
    if not files:
        raise SystemExit(f"No files matched: {glob_pattern}")
    chunks = []
    for fp in files:
        arr = np.loadtxt(fp)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        chunks.append(arr[:, :5])  # x y ux uy f
    data = np.vstack(chunks)
    return data, len(files)


def to_grid(data, n=1024):
    """Bin (x, y, ux, uy, f) point data onto an n x n regular grid, indexed
    by nearest-cell-center coordinates. Assumes L0=1, origin=(-0.5,-0.5)."""
    dx = 1.0 / n
    ix = np.clip(np.round((data[:, 0] + 0.5 - dx / 2) / dx).astype(int), 0, n - 1)
    iy = np.clip(np.round((data[:, 1] + 0.5 - dx / 2) / dx).astype(int), 0, n - 1)
    grid_ux = np.full((n, n), np.nan)
    grid_uy = np.full((n, n), np.nan)
    grid_f = np.full((n, n), np.nan)
    grid_ux[iy, ix] = data[:, 2]
    grid_uy[iy, ix] = data[:, 3]
    grid_f[iy, ix] = data[:, 4]
    return grid_ux, grid_uy, grid_f


def summarize(name, ux, uy, f):
    speed = np.sqrt(ux**2 + uy**2)
    liquid = f > 0.5
    print(f"\n=== {name} ===")
    print(f"  cells populated       : {np.sum(~np.isnan(ux))} / {ux.size}")
    print(f"  liquid cells (f>0.5)  : {np.sum(liquid[~np.isnan(f)])}")
    print(f"  mean |u| (all cells)  : {np.nanmean(speed):.6g}")
    print(f"  mean |u| (liquid only): {np.nanmean(speed[liquid]):.6g}")
    print(f"  max  |u| (all cells)  : {np.nanmax(speed):.6g}")
    print(f"  rms ux (liquid)       : {np.sqrt(np.nanmean(ux[liquid]**2)):.6g}")
    print(f"  rms uy (liquid)       : {np.sqrt(np.nanmean(uy[liquid]**2)):.6g}")
    print(f"  mean f                : {np.nanmean(f):.6g}")
    return speed


def main():
    print("Loading upstream dump...")
    up_data, up_nfiles = load_dump(UPSTREAM_GLOB)
    print(f"  {up_nfiles} rank files, {up_data.shape[0]} total points")

    print("Loading fork dump...")
    fk_data, fk_nfiles = load_dump(FORK_GLOB)
    print(f"  {fk_nfiles} rank files, {fk_data.shape[0]} total points")

    up_ux, up_uy, up_f = to_grid(up_data)
    fk_ux, fk_uy, fk_f = to_grid(fk_data)

    up_speed = summarize("UPSTREAM (Kim et al. own driver)", up_ux, up_uy, up_f)
    fk_speed = summarize("FORK (our driver)", fk_ux, fk_uy, fk_f)

    # Direct cell-by-cell diff (both on the identical N=1024 grid, same
    # origin/L0, so no interpolation needed).
    valid = ~np.isnan(up_speed) & ~np.isnan(fk_speed)
    diff_ux = fk_ux[valid] - up_ux[valid]
    diff_uy = fk_uy[valid] - up_uy[valid]
    diff_f = fk_f[valid] - up_f[valid]

    print("\n=== CELL-BY-CELL DIFF (fork - upstream, same grid) ===")
    print(f"  valid overlapping cells : {np.sum(valid)}")
    print(f"  mean |diff_u|           : {np.mean(np.sqrt(diff_ux**2 + diff_uy**2)):.6g}")
    print(f"  max  |diff_u|           : {np.max(np.sqrt(diff_ux**2 + diff_uy**2)):.6g}")
    print(f"  mean |diff_f|           : {np.mean(np.abs(diff_f)):.6g}")
    print(f"  max  |diff_f|           : {np.max(np.abs(diff_f)):.6g}")

    ratio = np.nanmean(fk_speed[valid]) / np.nanmean(up_speed[valid])
    print(f"\n  mean|u| ratio (fork/upstream): {ratio:.4g}")


if __name__ == "__main__":
    main()
