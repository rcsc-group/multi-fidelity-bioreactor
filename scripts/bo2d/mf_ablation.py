"""Model ablation on REAL data (2026-11-10, user: "use the one which performs better; start with 2-level L8-L9 not fully
Bayesian, then compare number of levels, fully Bayesian, basis").

PRE-REGISTERED (written before any result):
- Data: every real run in the frozen snapshots (scripts/bo2d/state/truth_snapshot_*.csv), ln units, mixing as the rate.
- Split by SITE (rounded rpm, theta): n9 L9 sites for training, all other L9 sites held out. L7 and L8 rows are all
  used (cheap levels are assumed available everywhere; stated).
- Variants (factor: levels x hyperparameters x basis), all with upstream mfbml KernelRidgeRegression as the lowest LF:
    U2-lin   upstream KRR-LR-GPR, L8 -> L9, ML hyperparameters (L-BFGS-B, 10 restarts), basis [1, f8]   <- the paper
    U2-quad  same, lf_poly_order="quadratic"
    U3-lin   upstream stacked: KRR-LR-GPR L7 -> L8, then its prediction as the LF of KRR-LR-GPR -> L9 (linear)
    U3-quad  same, quadratic at the L9 layer
    B2-lin   our FullBayes residual GP (isotropic, 2400-node hyper grid), LF = upstream KRR on L8, basis [1, f8]
    B2-quad  same, basis [1, f8, f8^2]
    B3-lin   KRR L7 -> FullBayes L8 -> FullBayes L9 (linear)
    B3-quad  same, quadratic at the L9 layer
  Baselines (L9 only): U1 = upstream GP (lf_poly_order="ordinary"), B1 = FullBayes constant basis.
- Score: RMSE (ln) on held-out L9 sites; 20 random site draws per n9 in {4, 8, 16}. Coverage of the 90 % interval is
  reported as a secondary number (upstream std as returned by predict(return_std=True)).
- Selection rule: per QoI, the variant with the lowest median RMSE at n9 = 8; ties (< 5 % apart) go to the simpler
  variant in the order U2-lin, U2-quad, U3-*, B2-*, B3-*. The selected variant is then frozen for the BO split study.
Usage: uv run python scripts/bo2d/mf_ablation.py QOI N9  -> CSV rows"""
import sys
import time
import warnings
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts")); sys.path.insert(0, str(ROOT / "scripts" / "campaign2d"))
warnings.filterwarnings("ignore")
from mfbml.methods.krr_lr_gpr import KernelRidgeLinearGaussianProcess as KLG  # noqa: E402  (upstream)
from mfbml.methods.kernel_ridge_regression import KernelRidgeRegression as KRR  # noqa: E402  (upstream)
from scripts.bo2d import bo  # noqa: E402
from scripts.bo2d import loop as LP  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402
import campaign as C  # noqa: E402

Q, N9 = sys.argv[1], int(sys.argv[2]); D = C.DOMAIN; Z90 = 1.645


class KLG(KLG):
    """Upstream KRR-LR-GPR, unchanged except ONE numpy-2 compatibility line: upstream writes
    `nll[i] = -logp.ravel()` (a 1-element array into a float slot), which numpy >= 2 rejects. Every other line of
    _logLikelihood is copied verbatim from mfbml 90d344a (methods/krr_lr_gpr.py L247-311), including its use of
    sigma2 (not log sigma2) in the noisy branch."""
    def _logLikelihood(self, params):
        from numpy.linalg import cholesky, solve
        params = np.atleast_2d(params); num_params = params.shape[1]; nll = np.zeros(params.shape[0])
        for i in range(params.shape[0]):
            if self.noise is None:
                param = params[i, 0: num_params - 1]; noise_sigma = params[i, -1]
            else:
                param = params[i, :]; noise_sigma = self.noise / self.yh_std
            K = self.kernel(self.sample_xh_scaled, self.sample_xh_scaled, param) + noise_sigma**2 * np.eye(self._num_xh)
            L = cholesky(K)
            alpha = solve(L.T, solve(L, self.sample_yh_scaled))
            KF = solve(L.T, solve(L, self.f))
            ld = cholesky(np.dot(self.f.T, KF))
            beta = solve(ld.T, solve(ld, np.dot(self.f.T, alpha)))
            gamma = solve(L.T, solve(L, (self.sample_yh_scaled - np.dot(self.f, beta))))
            sigma2 = np.dot((self.sample_yh_scaled - np.dot(self.f, beta)).T, gamma) / self._num_xh
            if self.noise == 0.0:
                logp = -0.5 * self._num_xh * np.log(sigma2) - np.sum(np.log(np.diag(L)))
            else:
                logp = -0.5 * self._num_xh * sigma2 - np.sum(np.log(np.diag(L)))
            nll[i] = float(-np.ravel(logp)[0])   # <- the only change
        return nll


class _Quad(FullBayesKRRLRGP):
    def _basis(self, X):
        F = super()._basis(X)
        return np.hstack([F, F[:, 1:2] ** 2])


class _Fixed:
    """A trained model used as the LF of the next upstream layer: train() is a no-op (already trained)."""
    def __init__(self, f):
        self.f = f

    def train(self, *a, **k):
        pass

    def predict(self, X):
        return np.asarray(self.f(X), float).reshape(-1, 1)


def krr(X, y):
    m = KRR(design_space=D, params_optimize=True, noise_data=True, optimizer_restart=10, seed=42)
    m.train(X, y.reshape(-1, 1), portion_test=0.2)
    return m


def upstream(Xh, yh, Xl, yl, order, lf=None):
    m = KLG(design_space=D, optimizer_restart=10, lf_model=lf, lf_poly_order=order, seed=42)
    m.train([Xh, Xl], [yh.reshape(-1, 1), yl.reshape(-1, 1)])
    return m


def up_pred(m, X):
    mu, sd = m.predict(X, return_std=True)
    return np.ravel(mu), np.ravel(sd)


def fb(Xh, yh, lf, quad, floor):
    K = _Quad if quad else FullBayesKRRLRGP
    m = K(D, lf_model=lf, basis="linear" if lf is not None else "ordinary", grid=C.hf_grid(floor, float(np.std(yh))))
    return m.train(Xh, yh)


def load(q):
    d = pd.read_csv(ROOT / "scripts/bo2d/state" / f"truth_snapshot_{q}.csv")
    d = d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]; d = LP.as_objective(d, q)
    d["ly"] = np.log(d.y.astype(float)); d["key"] = d.rpm.round(2).astype(str) + "_" + d.theta.round(2).astype(str)
    return d


d = load(Q); floor = C.noise_floor(Q)
xy = lambda s: (s[["rpm", "theta"]].to_numpy(float), s.ly.to_numpy(float))
X7, y7 = xy(d[d.level == 7]); X8, y8 = xy(d[d.level == 8]); d9 = d[d.level == 9]; sites = d9.key.unique()
print("qoi,n9,rep,variant,rmse,cover90,seconds", flush=True)
lf8 = krr(X8, y8); lf7 = krr(X7, y7)
up78 = {o: upstream(X8, y8, X7, y7, o, lf=_Fixed(lf7.predict)) for o in ("linear",)}
fb8 = fb(X8, y8, lf7, False, floor)
for rep in range(20):
    rng = np.random.default_rng(100 * rep + N9)
    tr = rng.choice(sites, N9, replace=False); A = d9[d9.key.isin(tr)]; B = d9[~d9.key.isin(tr)]
    Xa, ya = xy(A); Xb, yb = xy(B)
    V = {
        "U1": lambda: up_pred(upstream(Xa, ya, X8, y8, "ordinary", lf=_Fixed(lf8.predict)), Xb),
        "U2-lin": lambda: up_pred(upstream(Xa, ya, X8, y8, "linear", lf=_Fixed(lf8.predict)), Xb),
        "U2-quad": lambda: up_pred(upstream(Xa, ya, X8, y8, "quadratic", lf=_Fixed(lf8.predict)), Xb),
        "U3-lin": lambda: up_pred(upstream(Xa, ya, X8, y8, "linear", lf=_Fixed(lambda Z: up78["linear"].predict(Z))), Xb),
        "U3-quad": lambda: up_pred(upstream(Xa, ya, X8, y8, "quadratic", lf=_Fixed(lambda Z: up78["linear"].predict(Z))), Xb),
        "B1": lambda: bo.gauss(fb(Xa, ya, None, False, floor), Xb),
        "B2-lin": lambda: bo.gauss(fb(Xa, ya, lf8, False, floor), Xb),
        "B2-quad": lambda: bo.gauss(fb(Xa, ya, lf8, True, floor), Xb),
        "B3-lin": lambda: bo.gauss(fb(Xa, ya, bo._Mean(fb8), False, floor), Xb),
        "B3-quad": lambda: bo.gauss(fb(Xa, ya, bo._Mean(fb8), True, floor), Xb),
    }
    for name, f in V.items():
        t0 = time.time()
        try:
            mu, sd = f()
            rmse = float(np.sqrt(np.mean((mu - yb) ** 2))); cov = float(np.mean(np.abs(mu - yb) <= Z90 * sd))
        except Exception as e:  # report, never hide
            rmse, cov = np.nan, np.nan; print(f"# {name} rep {rep}: {type(e).__name__}: {e}", flush=True)
        print(f"{Q},{N9},{rep},{name},{rmse:.4f},{cov:.3f},{time.time() - t0:.1f}", flush=True)
