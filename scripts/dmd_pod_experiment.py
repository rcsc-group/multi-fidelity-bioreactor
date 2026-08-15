"""DMD/POD reconstruction-error prototype on the rocking-bioreactor velocity
field -- properly re-evaluated.

Loads a settled-window time series of u.x/u.y snapshots (frames_dmd/*.bin,
written by a scratch-only event added to BioReactor.c -- see diary.md
2026-08-14), and compares FOUR reconstruction methods at increasing mode
count:

  1. POD                    -- truncated SVD, baseline.
  2. exact DMD (t0-anchored) -- Tu et al. 2014, amplitudes fit to z(t=0) only
                                 and extrapolated forward (eq. 32-33/35 of
                                 Askham & Kutz 2018). This is what the first
                                 pass of this experiment used, and it
                                 performed WORSE than POD -- diagnosed as an
                                 artifact of anchoring only at t=0 combined
                                 with decaying eigenvalues (see diary.md
                                 2026-08-14).
  3. exact DMD (global b-fit) -- SAME eigenvalues/modes as (2), but the
                                 amplitudes b are fit by least squares
                                 against the WHOLE snapshot sequence at once
                                 (eq. 34 of Askham & Kutz 2018), not just
                                 z(t=0). Isolates whether the anchoring was
                                 the whole bug.
  4. optimized DMD           -- Askham & Kutz 2018 "Variable projection
                                 methods for an optimized dynamic mode
                                 decomposition" (arXiv:1704.02343), Algorithm
                                 3 (approximate/efficient optimized DMD).
                                 Jointly fits the continuous-time eigenvalues
                                 AND amplitudes via nonlinear least squares
                                 against the full sequence, in the rank-r
                                 POD-compressed coordinate system. This is
                                 the standard fix for exact DMD's bias/poor
                                 long-horizon reconstruction (Remark 5 of the
                                 same paper: exact DMD only fits pairwise
                                 one-step transitions, so its eigenvalues
                                 need not be optimal for reconstructing the
                                 whole sequence).

Usage:
    uv run python scripts/dmd_pod_experiment.py
"""
import glob
import struct

import numpy as np
from scipy.optimize import least_squares

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


def _exact_dmd_modes(X, r):
    """Algorithm 1 (Tu et al. 2014) -- discrete-time eigenvalues/modes at rank r."""
    X1, X2 = X[:, :-1], X[:, 1:]
    U, S, Vt = np.linalg.svd(X1, full_matrices=False)
    V = Vt.conj().T
    Ur, Sr, Vr = U[:, :r], S[:r], V[:, :r]
    Atilde = Ur.conj().T @ X2 @ Vr @ np.diag(1.0 / Sr)
    eigvals, W = np.linalg.eig(Atilde)
    Phi = X2 @ Vr @ np.diag(1.0 / Sr) @ W  # DMD modes (state_dim, r)
    return eigvals, Phi


def exact_dmd_t0_anchored(X, max_modes):
    """Method 2: amplitudes fit at t=0 only, extrapolated forward (the
    ORIGINAL, buggy evaluation from the first pass of this experiment)."""
    errs = []
    xnorm = np.linalg.norm(X)
    n_snap_full = X.shape[1]
    for r in range(1, max_modes + 1):
        eigvals, Phi = _exact_dmd_modes(X, r)
        b = np.linalg.lstsq(Phi, X[:, 0], rcond=None)[0]
        time_dynamics = np.array([eigvals**k for k in range(n_snap_full)]).T
        Xr = (Phi @ (b[:, None] * time_dynamics)).real
        errs.append(np.linalg.norm(X - Xr) / xnorm)
    return np.array(errs)


def exact_dmd_global_bfit(X, max_modes):
    """Method 3: SAME eigenvalues/modes as exact DMD, but amplitudes b fit
    by least squares against the WHOLE snapshot sequence (eq. 34, Askham &
    Kutz 2018) instead of anchored at t=0.

    X ~= Phi @ diag(b) @ V (V = Vandermonde in the discrete eigenvalues) is
    linear in b, but forming the (state_dim*n_snap) x r explicit system is
    wasteful. Solve the r x r normal equations directly instead:
      (M^H M)[i,k] = (Phi[:,i]^H Phi[:,k]) * (V[i,:] @ V[k,:].conj())
      (M^H x)[i]   = Phi[:,i]^H @ (X @ V[i,:].conj())
    """
    errs = []
    xnorm = np.linalg.norm(X)
    n_snap_full = X.shape[1]
    for r in range(1, max_modes + 1):
        eigvals, Phi = _exact_dmd_modes(X, r)
        V = np.array([eigvals**k for k in range(n_snap_full)]).T  # (r, n_snap)

        G_phi = Phi.conj().T @ Phi  # (r, r)
        G_v = V @ V.conj().T  # (r, r)
        MhM = G_phi * G_v  # Hadamard product, (r, r)
        Mhx = np.array([Phi[:, i].conj() @ (X @ V[i, :].conj()) for i in range(r)])

        b = np.linalg.solve(MhM, Mhx)
        Xr = (Phi @ (b[:, None] * V)).real
        errs.append(np.linalg.norm(X - Xr) / xnorm)
    return np.array(errs)


def optimized_dmd(X, ts, max_modes, verbose=False):
    """Method 4: Algorithm 3, "Approximate optimized DMD" (Askham & Kutz
    2018, arXiv:1704.02343). Jointly fits continuous-time eigenvalues alpha
    AND amplitudes B by nonlinear least squares against the full sequence,
    in the rank-r POD-compressed coordinate system (cheap: r x r problem,
    not state_dim x r). Initial guess for alpha comes from exact DMD's
    discrete eigenvalues via alpha = log(lambda)/dt_mean (the paper's own
    suggested initialization scheme)."""
    dt_mean = np.mean(np.diff(ts))
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    errs = []
    xnorm = np.linalg.norm(X)

    for r in range(1, max_modes + 1):
        Ur, Sr, Vr = U[:, :r], S[:r], Vt[:r, :].conj().T  # Vr: (n_snap, r)
        Y = Vr @ np.diag(Sr)  # (n_snap, r) -- target for the r x r fit

        # Initial guess: exact DMD's discrete eigenvalues -> continuous rates.
        eigvals0, _ = _exact_dmd_modes(X, r)
        alpha0 = np.log(eigvals0.astype(complex)) / dt_mean

        def build_Phi(alpha):
            return np.exp(np.outer(ts, alpha))  # (n_snap, r)

        def residual(params):
            alpha = params[:r] + 1j * params[r:]
            Phi = build_Phi(alpha)
            B = np.linalg.lstsq(Phi, Y, rcond=None)[0]  # variable projection
            res = Y - Phi @ B
            return np.concatenate([res.real.ravel(), res.imag.ravel()])

        # Bound the real part (growth/decay rate) to a physically sane range
        # -- the settled flow is confirmed periodic/non-exploding, so allow
        # generous decay but essentially no growth. 'lm' has no bound support
        # and an unconstrained LM step can send exp(alpha*t) to overflow
        # (hit this directly: NaNs from an unbounded real-part excursion
        # broke the lstsq inside residual() with 'SVD did not converge').
        # Imaginary part (oscillation frequency) is left unbounded.
        lower = np.concatenate([np.full(r, -50.0), np.full(r, -np.inf)])
        upper = np.concatenate([np.full(r, 0.5), np.full(r, np.inf)])
        x0 = np.concatenate([alpha0.real, alpha0.imag])
        x0_clipped = np.clip(x0, lower, upper)
        result = least_squares(
            residual, x0_clipped, method="trf", bounds=(lower, upper),
            max_nfev=2000 * (r + 1),
        )
        alpha_hat = result.x[:r] + 1j * result.x[r:]
        Phi_hat = build_Phi(alpha_hat)
        B_hat = np.linalg.lstsq(Phi_hat, Y, rcond=None)[0]  # (r, r)

        # Reconstruct in the FULL state space: z(t) = U_r @ B_hat.T-ish
        # combination -- per eq. 55, phi_i = normalize(U_r @ B_hat.T[:,i]),
        # b_i = norm. Equivalently, full reconstruction is:
        #   X_r = U_r @ (B_hat @ Phi_hat.T).T = U_r @ Phi_hat @ B_hat ... check dims
        # Y ~= Phi_hat @ B_hat, and X ~= U_r @ Vr.T @ diag(Sr) = U_r @ Y.T
        # so X_reconstructed ~= U_r @ (Phi_hat @ B_hat).T
        Xr = (Ur @ (Phi_hat @ B_hat).T).real
        err = np.linalg.norm(X - Xr) / xnorm
        errs.append(err)
        if verbose:
            print(f"  r={r}: cost={result.cost:.4g}, nfev={result.nfev}, rel_err={err:.6g}")

    return np.array(errs)


def main():
    print(f"Loading frames from {FRAMES_DIR} ...")
    ts, X = load_all_frames()
    print(f"  {X.shape[1]} snapshots, state dim {X.shape[0]}, t in [{ts[0]:.4g}, {ts[-1]:.4g}]")
    dt = np.mean(np.diff(ts))
    print(f"  dt (mean) = {dt:.6g}, dt std = {np.std(np.diff(ts)):.4g}")

    max_modes = min(20, X.shape[1] - 1)

    print("\nFitting POD...")
    pod_errs, singular_values = pod_reconstruction_error(X, max_modes)

    print("Fitting exact DMD (t0-anchored, ORIGINAL buggy method)...")
    dmd_t0_errs = exact_dmd_t0_anchored(X, max_modes)

    print("Fitting exact DMD (global b-fit, eq. 34)...")
    dmd_global_errs = exact_dmd_global_bfit(X, max_modes)

    print("Fitting optimized DMD (Algorithm 3, variable projection)...")
    opt_dmd_errs = optimized_dmd(X, ts, max_modes, verbose=True)

    print("\n=== Relative reconstruction error (Frobenius norm) vs. modes retained ===")
    header = f"{'modes':>6} {'POD':>12} {'DMD(t0)':>12} {'DMD(global b)':>14} {'Optimized DMD':>15}"
    print(header)
    for r in range(max_modes):
        print(
            f"{r+1:>6} {pod_errs[r]:>12.6g} {dmd_t0_errs[r]:>12.6g} "
            f"{dmd_global_errs[r]:>14.6g} {opt_dmd_errs[r]:>15.6g}"
        )

    print("\n=== Singular value spectrum (normalized, first 15) ===")
    sv_norm = singular_values / singular_values[0]
    for i, sv in enumerate(sv_norm[:15]):
        print(f"  mode {i+1}: {sv:.6g}")

    for tol in (1e-1, 1e-2, 1e-3):
        def modes_needed(errs):
            hit = np.where(errs < tol)[0]
            return int(hit[0] + 1) if len(hit) else None

        print(
            f"\nModes needed for rel err < {tol}: "
            f"POD={modes_needed(pod_errs)}, "
            f"DMD(t0)={modes_needed(dmd_t0_errs)}, "
            f"DMD(global b)={modes_needed(dmd_global_errs)}, "
            f"OptimizedDMD={modes_needed(opt_dmd_errs)}"
        )


if __name__ == "__main__":
    main()
