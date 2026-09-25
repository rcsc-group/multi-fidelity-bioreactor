"""Multi-fidelity demo on max wall shear stress vs rpm: L8 -> L10 (Yi et al. 2024).

Not a Kim replica. Tests KRR-LR-GPR (arXiv 2407.15110) on a quantity where
both fidelities exist at all nine rpms of Kim's Fig 13(a), so the model can be
trained on a few L10 points and SCORED on the rest.

  LF  = L8 at 9 rpm, fig13a_rm_mf_rpm*   (KRR, deterministic)
  HF  = L10 at 3 rpm, l10_fig13a_xlevel_rpm*   (linear transfer + GPR residual)
  test= the other 6 L10 points, never seen in training

Line style: solid = Kim, dashed = a fitted model; our data are markers.
Error bars and the model noise are both the standard error of a mean over
settled cycles, so they mean the same thing.

Response: signed max wall shear over one cycle, averaged over settled cycles
(Kim's tau_liq_max takes max() of the signed field). Per-point aleatoric
uncertainty is MEASURED as the cycle-to-cycle standard deviation of that
per-cycle max, independent of any model; it is compared against the single
homoscedastic noise sigma the GPR infers. Epistemic = GPR posterior std.

Baseline: the same model with lf_poly_order="ordinary", which drops the LF
basis and reduces to a single-fidelity GPR on the 3 L10 points. If MF does not
beat it on the held-out 6, the L8 data added nothing.

Usage:  uv run python scripts/plot_mf_tau_max.py [--hf 17.5 27.5 37.5]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import figstyle as fs                            # noqa: E402
from mfbml_local.krr_lr_gpr import KernelRidgeLinearGaussianProcess  # noqa: E402
from mfbml_local.kernel_ridge_regression import KernelRidgeRegression  # noqa: E402
from scipy.optimize import minimize                           # noqa: E402

RPMS = [17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
LF_PREFIX, HF_PREFIX = "fig13a_rm_mf_rpm", "l10_fig13a_xlevel_rpm"
KIM_CSV = ROOT / "experiments/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv"
OUT = ROOT / "experiments/multifidelity/mf_tau_max_rpm.png"
COL_TAU_SIGNED_MAX = 8          # shear_stress.dat column, as in plot_fig13.py


def per_cycle_max(run: str, settled_from: float) -> np.ndarray:
    """Signed max wall shear [Pa] in each whole settled cycle."""
    d = ROOT / "runs" / run
    p = json.load(open(d / "params.json"))
    L, H = p["geometry"]["a"], 2 * p["geometry"]["b"]
    th = math.radians(p["theta_max"][0])
    Tp = 2 * math.pi / p["omega_b"]
    U = L / 4 * (H + 0.5 * L * math.tan(th)) / (H * 0.5) / Tp
    T_nd = Tp / (L / U)
    a = np.loadtxt(d / "shear_stress.dat", skiprows=1)
    c = (a[:, 1] - a[0, 1]) / T_nd
    tau = a[:, COL_TAU_SIGNED_MAX] * 1000.0 * U ** 2
    k0, k1 = int(math.ceil(settled_from * c[-1])), int(math.floor(c[-1]))
    return np.array([tau[(c >= k) & (c < k + 1)].max() for k in range(k0, k1)])


def dataset(prefix: str, settled_from: float) -> pd.DataFrame:
    rows = []
    for r in RPMS:
        m = per_cycle_max(f"{prefix}{r:g}", settled_from)
        rows.append(dict(rpm=r, y=m.mean(), sd=m.std(ddof=1), n=len(m)))
    return pd.DataFrame(rows)


class LooKRR(KernelRidgeRegression):
    """Upstream KRR, hyperparameters by leave-one-out CV.

    Deviation from Yi et al.: upstream scores (theta, noise) on ONE random
    split holding out portion_test of the LF points -- 2 of 9 here -- which
    selected a flat fit (0.0918..0.0939 Pa across data spanning 0.05..0.25).
    Their benchmarks use 200*d LF points, where one split is enough.
    """

    def _optimize_kernel_params(self) -> None:
        X, Y, n = self.sample_x_scaled, self.sample_y_scaled, len(self.sample_x_scaled)

        def loo(params):
            self._set_kernel_params(params=params)
            err = 0.0
            for i in range(n):
                k = np.arange(n) != i
                K = self._training_kernel_matrix(scaled_x=X[k],
                                                 scaled_noise_std=self.noise_std / self.y_std)
                L = np.linalg.cholesky(K)
                W = np.linalg.solve(L.T, np.linalg.solve(L, Y[k]))
                pred = np.dot(W.T, self.kernel.get_kernel_matrix(X[k], X[i:i + 1])).ravel()
                err += float((Y[i] - pred[0]) ** 2)
            return err / n

        bounds = self._bound_definition_for_optimization()
        best, rng = None, np.random.default_rng(self.seed)
        for _ in range(self.optimizer_restart + 1):
            r = minimize(loo, rng.uniform(bounds[:, 0], bounds[:, 1]),
                         method="l-bfgs-b", bounds=bounds.tolist())
            if best is None or r.fun < best.fun:
                best = r
        self._set_kernel_params(params=best.x)


def fit(X_hf, y_hf, X_lf, y_lf, poly, noise=None):
    ds = np.array([[min(RPMS), max(RPMS)]])
    lf_model = LooKRR(design_space=ds, params_optimize=True, noise_data=True,
                      optimizer_restart=10, seed=42)
    mdl = KernelRidgeLinearGaussianProcess(
        design_space=ds, optimizer_restart=10, lf_poly_order=poly, seed=42,
        lf_model=lf_model, noise_prior=noise)
    mdl.train(X=[X_hf, X_lf], Y=[y_hf, y_lf])
    return mdl


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hf", type=float, nargs="+", default=[17.5, 27.5, 37.5])
    a = ap.parse_args()

    # L8 is a cold start (33 cycles): use the last 40%, as plot_fig13.py does.
    # L10 warm-starts settled (per-cycle amplitude flat to <6%, plot_fig13.py).
    lf = dataset(LF_PREFIX, 0.6)
    hf = dataset(HF_PREFIX, 0.0)
    train = hf[hf.rpm.isin(a.hf)]
    test = hf[~hf.rpm.isin(a.hf)]
    kim = pd.read_csv(KIM_CSV, skiprows=[1])
    kim["RPM"] = pd.to_numeric(kim["RPM"])
    kim = kim.sort_values("RPM")

    col = lambda s: s.to_numpy().reshape(-1, 1)
    X_lf, y_lf = col(lf.rpm), lf.y.to_numpy()
    X_hf, y_hf = col(train.rpm), train.y.to_numpy()
    # Three HF points cannot separate noise from signal (a free noise
    # hyperparameter collapses to 0 and 95% coverage to 0%). Fix it instead at
    # the MEASURED aleatoric level of the response: the standard error of a
    # mean over n settled cycles, pooled over the training points only.
    se = float(np.sqrt(np.mean(train.sd.to_numpy() ** 2 / train.n.to_numpy())))
    print(f"aleatoric noise fixed at measured L10 standard error {se:.4f} Pa")
    mf_free = fit(X_hf, y_hf, X_lf, y_lf, "linear")
    mf = fit(X_hf, y_hf, X_lf, y_lf, "linear", noise=se)
    sf = fit(X_hf, y_hf, X_lf, y_lf, "ordinary", noise=se)

    xg = np.linspace(min(RPMS), max(RPMS), 200).reshape(-1, 1)
    krr = mf.predict_lf(xg).ravel()
    mu, tot = (v.ravel() for v in mf.predict(xg, return_std=True))
    epi = mf.epistemic.ravel()

    def score(m, name):
        mt, st = (v.ravel() for v in m.predict(col(test.rpm), return_std=True))
        err = mt - test.y.to_numpy()
        cov = np.mean(np.abs(err) <= 1.96 * st)
        print(f"  {name:18s} held-out RMSE {np.sqrt(np.mean(err**2)):.4f} Pa  "
              f"max|err| {np.abs(err).max():.4f}  95% coverage {cov:.0%}")
    print(f"LF L8: {len(lf)} pts, HF L10 train {list(train.rpm)}, test {list(test.rpm)}")
    print(f"transfer rho = {np.ravel(mf.beta)}  (normalised units)")
    print(f"GPR noise sigma = {float(mf.noise):.4f} Pa")
    print(f"measured cycle-to-cycle sd: L8 {lf.sd.mean():.4f} Pa (median n={int(lf.n.median())}), "
          f"L10 {hf.sd.mean():.4f} Pa (n={int(hf.n.median())})")
    score(mf_free, "MF, noise free")
    score(mf, "MF KRR-LR-GPR")
    score(sf, "SF GPR (L10 only)")
    score_lf = np.interp(test.rpm, lf.rpm, lf.y) - test.y.to_numpy()
    print(f"  {'L8 as-is':18s} held-out RMSE {np.sqrt(np.mean(score_lf**2)):.4f} Pa")

    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    c8, c10 = fs.level_colour(8), fs.level_colour(10)
    ax.fill_between(xg.ravel(), mu - 1.96 * tot, mu + 1.96 * tot, color=c10,
                    alpha=0.12, lw=0, label="MF 95%, total")
    ax.fill_between(xg.ravel(), mu - 1.96 * epi, mu + 1.96 * epi, color=c10,
                    alpha=0.28, lw=0, label="MF 95%, epistemic")
    ax.plot(xg, krr, color=c8, lw=1.1, ls="--", label="KRR on L8")
    ax.plot(xg, mu, color=c10, lw=1.4, ls="--", label="MF KRR-LR-GPR")
    ax.errorbar(lf.rpm, lf.y, yerr=1.96 * lf.sd / np.sqrt(lf.n), fmt="none", ecolor=c8, lw=0.8, capsize=2)
    ax.plot(lf.rpm, lf.y, **fs.series_kw(c8, fs.MK["tau"], ours=True, ls="none"), label="L8")
    ax.errorbar(hf.rpm, hf.y, yerr=1.96 * hf.sd / np.sqrt(hf.n), fmt="none", ecolor=c10, lw=0.8, capsize=2)
    ax.plot(hf.rpm, hf.y, **fs.series_kw(c10, fs.MK["tau"], ours=True, ls="none"), label="L10")
    ax.plot(train.rpm, train.y, ls="none", marker="o", ms=2.3 * fs.MS_OURS, mfc="none",
            mec=c10, mew=1.0, label="L10 used for training")
    ax.plot(kim.RPM, kim.tau_liq_max, **fs.series_kw(fs.KIM, fs.MK["tau"]), label="Kim et al.")
    ax.set_xlabel("rocking frequency [rpm]")
    ax.set_ylabel(r"max wall shear stress $\tau_{\max}$ [Pa]")
    ax.grid(**fs.GRID_KW)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
