"""First DMD/POD reconstruction-error prototype on the rocking-bioreactor
velocity field.

Loads a settled-window time series of u.x/u.y snapshots (frames_dmd/*.bin,
written by a scratch-only event added to BioReactor.c for this experiment --
see diary.md 2026-08-14), fits both POD and exact DMD, and reports relative
reconstruction error vs. number of modes retained. Answers the concrete
question: how many modes does it take to reconstruct this flow "almost
losslessly"?

Usage:
    uv run python scripts/dmd_pod_experiment.py
"""
import glob
import struct

import numpy as np

FRAMES_DIR = "/oscar/scratch/eaguerov/tmp/dmd_experiment/rundir/frames_dmd"


def load_frame(fpath):
    with open(fpath, "rb") as fh:
        n = struct.unpack("i", fh.read(4))[0]
        t = struct.unpack("d", fh.read(8))[0]
        ux = np.frombuffer(fh.read(4 * n * n), dtype=np.float32).reshape(n, n)
        uy = np.frombuffer(fh.read(4 * n * n), dtype=np.float32).reshape(n, n)
    return t, ux, uy


def load_all_frames():
    files = sorted(glob.glob(f"{FRAMES_DIR}/frame_*.bin"))
    if not files:
        raise SystemExit(f"No frames found in {FRAMES_DIR}")
    ts, snaps = [], []
    for fp in files:
        t, ux, uy = load_frame(fp)
        ts.append(t)
        snaps.append(np.concatenate([ux.ravel(), uy.ravel()]))
    return np.array(ts), np.array(snaps).T  # X: (state_dim, n_snapshots)


def pod_reconstruction_error(X, max_modes):
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    errs = []
    xnorm = np.linalg.norm(X)
    for r in range(1, max_modes + 1):
        Xr = U[:, :r] @ np.diag(S[:r]) @ Vt[:r, :]
        errs.append(np.linalg.norm(X - Xr) / xnorm)
    return np.array(errs), S


def exact_dmd(X, dt, max_modes):
    """Standard exact DMD (Tu et al. 2014). X: (state_dim, n_snapshots)."""
    X1, X2 = X[:, :-1], X[:, 1:]
    U, S, Vt = np.linalg.svd(X1, full_matrices=False)
    V = Vt.conj().T
    n_snap = X1.shape[1]
    max_modes = min(max_modes, n_snap)

    errs = []
    xnorm = np.linalg.norm(X)
    for r in range(1, max_modes + 1):
        Ur, Sr, Vr = U[:, :r], S[:r], V[:, :r]
        Atilde = Ur.conj().T @ X2 @ Vr @ np.diag(1.0 / Sr)
        eigvals, W = np.linalg.eig(Atilde)
        Phi = X2 @ Vr @ np.diag(1.0 / Sr) @ W  # DMD modes, (state_dim, r)

        # Reconstruct the full snapshot sequence from the r-mode expansion.
        b = np.linalg.lstsq(Phi, X[:, 0], rcond=None)[0]
        n_snap_full = X.shape[1]
        time_dynamics = np.array(
            [eigvals**k for k in range(n_snap_full)]
        ).T  # (r, n_snap_full)
        Xr = (Phi @ (b[:, None] * time_dynamics)).real
        errs.append(np.linalg.norm(X - Xr) / xnorm)
    return np.array(errs)


def main():
    print(f"Loading frames from {FRAMES_DIR} ...")
    ts, X = load_all_frames()
    print(f"  {X.shape[1]} snapshots, state dim {X.shape[0]}, t in [{ts[0]:.4g}, {ts[-1]:.4g}]")
    dt = np.mean(np.diff(ts))
    print(f"  dt (mean) = {dt:.6g}")

    max_modes = min(30, X.shape[1] - 1)

    print("\nFitting POD...")
    pod_errs, singular_values = pod_reconstruction_error(X, max_modes)

    print("Fitting exact DMD...")
    dmd_errs = exact_dmd(X, dt, max_modes)

    print("\n=== Relative reconstruction error (Frobenius norm) vs. modes retained ===")
    print(f"{'modes':>6} {'POD rel err':>14} {'DMD rel err':>14}")
    for r in range(max_modes):
        print(f"{r+1:>6} {pod_errs[r]:>14.6g} {dmd_errs[r]:>14.6g}")

    print("\n=== Singular value spectrum (normalized, first 15) ===")
    sv_norm = singular_values / singular_values[0]
    for i, sv in enumerate(sv_norm[:15]):
        print(f"  mode {i+1}: {sv:.6g}")

    # Report the mode count needed to reach a few representative tolerances.
    for tol in (1e-1, 1e-2, 1e-3):
        pod_r = np.searchsorted(pod_errs < tol, True) + 1 if np.any(pod_errs < tol) else None
        dmd_r = np.searchsorted(dmd_errs < tol, True) + 1 if np.any(dmd_errs < tol) else None
        print(f"\nModes needed for rel err < {tol}: POD={pod_r}, DMD={dmd_r}")


if __name__ == "__main__":
    main()
