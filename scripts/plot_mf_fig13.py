"""Multi-fidelity demos on Kim's Fig 13 quantities (Yi et al. 2024, KRR-LR-GPR).

Not Kim replicas. Each case trains on many LF points and a few HF points and
is SCORED on HF points it never saw; those held-out points are reported in the
printout but not drawn.

  tau_rpm   Fig 13(a): max wall shear vs rpm at 7 deg.  LF L8, HF L10.
  edr_angle Fig 13(b): max EDR vs angle at 32.5 rpm.    LF L7, HF L9.
  edr_mean_angle  same, max of the spatially averaged EDR (Kim's hollow squares).

Response: per-cycle max of the quantity, averaged over settled cycles. Its
aleatoric uncertainty is MEASURED as the standard error over those cycles, with
no model involved, and the GPR noise is fixed at the pooled training-point
value: with three HF points a free noise hyperparameter collapses to zero
(measured: 95% coverage 0%). Epistemic = GPR posterior std.

KRR hyperparameters by leave-one-out CV, not upstream's single random split,
which held out 2 of 9 LF points and selected a flat fit (diary 2026-09-25).

Baseline: lf_poly_order="ordinary", i.e. a single-fidelity GPR on the HF
training points alone. If MF does not beat it, the LF data added nothing.

Line style: solid = Kim, dashed = a fitted model; our data are markers.

Usage:  uv run python scripts/plot_mf_fig13.py tau_rpm [--hf 17.5 27.5 37.5]
        uv run python scripts/plot_mf_fig13.py edr_angle [--hf 2 4 7]
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
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from scripts import figstyle as fs                                      # noqa: E402
from mfbml_local.krr_lr_gpr import KernelRidgeLinearGaussianProcess     # noqa: E402
from mfbml_local.kernel_ridge_regression import KernelRidgeRegression   # noqa: E402

CSV = ROOT / "experiments/kimetal2024/csv_raw"
CASES = {
    "tau_rpm": dict(
        xs=[17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5], lf_level=8, hf_level=10,
        lf_run="fig13a_rm_mf_rpm{:g}", hf_run="l10_fig13a_xlevel_rpm{:g}",
        # L8 cold-starts over 33 cycles: last 40%. L10 warm-starts settled.
        lf_settled=0.6, hf_settled=0.0, hf_train=[17.5, 27.5, 37.5],
        column="tau_100_signed", kind="tau", marker="tau",
        kim_csv=CSV / "shear_ediss_vs_frequency.csv", kim_skip=[1],
        kim_x="RPM", kim_y="tau_liq_max",
        xlabel="rocking frequency [rpm]",
        ylabel=r"max wall shear stress $\tau_{\max}$ [Pa]",
        out="mf_tau_max_rpm.png"),
    "edr_angle": dict(
        xs=[2.0, 3.0, 4.0, 5.0, 6.0, 7.0], lf_level=7, hf_level=9,
        lf_run="mf_l7_th{:g}", hf_run="fig10_l9_th{:g}",
        # both cold starts; last 40% of the L7 run, and of the L9 spin-up
        lf_settled=0.6, hf_settled=0.6, hf_train=[2.0, 4.0, 7.0],
        column="ediss_kim_max", kind="ediss", marker="ediss",
        kim_csv=CSV / "shear_ediss_vs_angle.csv", kim_skip=None,
        kim_x="theta_deg", kim_y="Ediss_liq_max",
        xlabel=r"rocking angle $\theta_{\max}$ [deg]",
        ylabel=r"max energy dissipation rate $\varepsilon_{\max}$ [W m$^{-3}$]",
        out="mf_edr_max_angle.png"),
}
# Same runs, Kim's hollow squares: the max over a cycle of the SPATIALLY
# AVERAGED EDR. L7/L9 correlate at r = 0.997 here against 0.46 for the pointwise
# max above, which is a grid-scale extreme L7 cannot resolve -- so the pair is
# the method working and its negative control.
CASES["edr_mean_angle"] = dict(
    CASES["edr_angle"], column="ediss_kim_mean", kim_y="Ediss_liq_mean", stat="mean",
    ylabel=r"max spatially averaged EDR $\langle\varepsilon\rangle_{\max}$ [W m$^{-3}$]",
    out="mf_edr_mean_angle.png")


def per_cycle_max(run: str, column: str, kind: str, settled_from: float) -> np.ndarray:
    """Per-cycle max of a shear_stress.dat column, dimensional, settled cycles.

    shear_stress.dat is non-dimensional: tau scales with rho U^2, EDR with
    rho U^3 / L. Only cycles before tracer/oxygen release are used; they are
    passive, so this is a choice of window, not a correction.
    """
    d = ROOT / "runs" / run
    p = json.load(open(d / "params.json"))
    L, H = p["geometry"]["a"], 2 * p["geometry"]["b"]
    th = math.radians(p["theta_max"][0])
    Tp = 2 * math.pi / p["omega_b"]
    U = L / 4 * (H + 0.5 * L * math.tan(th)) / (H * 0.5) / Tp
    T_nd = Tp / (L / U)
    df = pd.read_csv(d / "shear_stress.dat", sep=r"\s+")
    t = df["t"].to_numpy()
    t_rel = (p.get("t_checkpoint") or 0) + p["n_mix_cycles"] * T_nd
    keep = t < t_rel if t_rel < t[-1] else np.ones_like(t, bool)
    t = t[keep]
    scale = 1000.0 * U ** 2 if kind == "tau" else 1000.0 * U ** 3 / L
    y = df[column].to_numpy()[keep] * scale
    c = (t - t[0]) / T_nd
    k0, k1 = int(math.ceil(settled_from * c[-1])), int(math.floor(c[-1]))
    return np.array([y[(c >= k) & (c < k + 1)].max() for k in range(k0, k1)])


def dataset(case: dict, which: str) -> pd.DataFrame:
    rows = []
    for x in case["xs"]:
        m = per_cycle_max(case[f"{which}_run"].format(x), case["column"],
                          case["kind"], case[f"{which}_settled"])
        rows.append(dict(x=x, y=m.mean(), se=m.std(ddof=1) / np.sqrt(len(m)), n=len(m)))
    return pd.DataFrame(rows)


class LooKRR(KernelRidgeRegression):
    """Upstream KRR, hyperparameters by leave-one-out CV (see module docstring)."""

    def _optimize_kernel_params(self) -> None:
        X, Y, n = self.sample_x_scaled, self.sample_y_scaled, len(self.sample_x_scaled)

        def loo(params):
            self._set_kernel_params(params=params)
            err = 0.0
            for i in range(n):
                k = np.arange(n) != i
                K = self._training_kernel_matrix(
                    scaled_x=X[k], scaled_noise_std=self.noise_std / self.y_std)
                Lc = np.linalg.cholesky(K)
                W = np.linalg.solve(Lc.T, np.linalg.solve(Lc, Y[k]))
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


def fit(ds, X_hf, y_hf, X_lf, y_lf, poly, noise):
    lf_model = LooKRR(design_space=ds, params_optimize=True, noise_data=True,
                      optimizer_restart=10, seed=42)
    m = KernelRidgeLinearGaussianProcess(
        design_space=ds, optimizer_restart=10, lf_poly_order=poly, seed=42,
        lf_model=lf_model, noise_prior=noise)
    m.train(X=[X_hf, X_lf], Y=[y_hf, y_lf])
    return m


def lengthscale(kernel, span: float) -> float:
    """RBF here is exp(-theta d^2) on inputs scaled to [0, 1]."""
    return float(span / np.sqrt(2 * np.ravel(kernel.param)[0]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=sorted(CASES))
    ap.add_argument("--hf", type=float, nargs="+")
    a = ap.parse_args()
    cs = CASES[a.case]
    hf_train = a.hf or cs["hf_train"]

    lf, hf = dataset(cs, "lf"), dataset(cs, "hf")
    train, test = hf[hf.x.isin(hf_train)], hf[~hf.x.isin(hf_train)]
    kim = pd.read_csv(cs["kim_csv"], skiprows=cs["kim_skip"])
    kim[cs["kim_x"]] = pd.to_numeric(kim[cs["kim_x"]])
    kim = kim.sort_values(cs["kim_x"])

    col = lambda s: s.to_numpy().reshape(-1, 1)
    lo, hi = min(cs["xs"]), max(cs["xs"])
    ds = np.array([[lo, hi]])
    noise = float(np.sqrt(np.mean(train.se.to_numpy() ** 2)))
    args = (ds, col(train.x), train.y.to_numpy(), col(lf.x), lf.y.to_numpy())
    mf = fit(*args, "linear", noise)
    sf = fit(*args, "ordinary", noise)

    print(f"{a.case}: r(LF,HF) = {np.corrcoef(lf.y, hf.y)[0, 1]:+.3f}; LF L{cs['lf_level']} x{len(lf)}, HF L{cs['hf_level']} "
          f"train {list(train.x)}, held out {list(test.x)}")
    print(f"  aleatoric noise (pooled training SE) {noise:.4g}")
    print(f"  lengthscale KRR {lengthscale(mf.lf_model.kernel, hi - lo):.3g}, "
          f"GPR residual {lengthscale(mf.kernel, hi - lo):.3g}  (x units; range {hi - lo:g})")
    print(f"  transfer rho {np.round(np.ravel(mf.beta), 3)} (normalised)")
    for m, name in ((mf, "MF KRR-LR-GPR"), (sf, f"GPR on L{cs['hf_level']} only")):
        mt, st = (v.ravel() for v in m.predict(col(test.x), return_std=True))
        err = mt - test.y.to_numpy()
        print(f"  {name:18s} held-out RMSE {np.sqrt(np.mean(err ** 2)):.4g}  "
              f"max|err| {np.abs(err).max():.4g}  95% coverage "
              f"{np.mean(np.abs(err) <= 1.96 * st):.0%}")
    e_lf = np.interp(test.x, lf.x, lf.y) - test.y.to_numpy()
    print(f"  {'LF as-is':18s} held-out RMSE {np.sqrt(np.mean(e_lf ** 2)):.4g}")

    xg = np.linspace(lo, hi, 200).reshape(-1, 1)
    krr = mf.predict_lf(xg).ravel()
    mu, tot = (v.ravel() for v in mf.predict(xg, return_std=True))
    epi = mf.epistemic.ravel()

    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    cl, ch = fs.level_colour(cs["lf_level"]), fs.level_colour(cs["hf_level"])
    mk = fs.MK[cs["marker"]]
    stat = cs.get("stat", "max")
    ax.fill_between(xg.ravel(), mu - 1.96 * tot, mu + 1.96 * tot, color=ch,
                    alpha=0.12, lw=0, label="MF 95%, total")
    ax.fill_between(xg.ravel(), mu - 1.96 * epi, mu + 1.96 * epi, color=ch,
                    alpha=0.28, lw=0, label="MF 95%, epistemic")
    ax.plot(xg, krr, color=cl, lw=1.1, ls="--", label=f"KRR on L{cs['lf_level']}")
    ax.plot(xg, mu, color=ch, lw=1.4, ls="--", label="MF KRR-LR-GPR")
    for d, c, lab in ((lf, cl, f"L{cs['lf_level']}"), (train, ch, f"L{cs['hf_level']}")):
        ax.errorbar(d.x, d.y, yerr=1.96 * d.se, fmt="none", ecolor=c, lw=0.8, capsize=2)
        ax.plot(d.x, d.y, **fs.series_kw(c, mk, stat=stat, ours=True, ls="none"), label=lab)
    ax.plot(kim[cs["kim_x"]], kim[cs["kim_y"]], **fs.series_kw(fs.KIM, mk, stat=stat),
            label="Kim et al.")
    ax.set_xlabel(cs["xlabel"])
    ax.set_ylabel(cs["ylabel"])
    ax.grid(**fs.GRID_KW)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    out = ROOT / "experiments/multifidelity" / cs["out"]
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
