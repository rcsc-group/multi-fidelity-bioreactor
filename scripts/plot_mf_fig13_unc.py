"""Fig 13 multi-fidelity demos with numerical (grid) uncertainty on the HF data.

Not Kim replicas. Same structure as plot_mf_fig13.py (Yi et al. 2024: a KRR on
the LF data, then a linear transfer plus a GP residual on the HF data). What
changes is what the HF uncertainty means.

HF noise. Each HF point i gets its own variance
    sigma_i^2 = SE_i^2 + (U_i / 1.96)^2
where SE_i is the cycle-to-cycle standard error (sampling), and U_i is the
~95% numerical-uncertainty bound on |q_i - q_exact| from least-squares grid
fits (scripts/diag_num_uncertainty.py; simplified Eca & Hoekstra 2014). U is
10-100x SE for these statistics (diary 2026-09-29), so SE alone made the old
bands claim knowledge of L10, not of the physics.
  13(a): U/q from levels {8, 9, 10} at each rpm, pooled as an RMS over rpm
         (each single-rpm estimate has one degree of freedom).
  13(b): no L10 at the angles yet. U(L9)/q(L9) is measured at the one condition
         where it can be (32.5 rpm, 7 deg; levels 8, 9, 10) and applied at every
         angle. This is an ASSUMPTION, and it is replaced once L10 angle runs land.

Transfer step. Universal kriging with KNOWN heteroscedastic noise:
    y = F beta + r,   Cov = s^2 R_l + D,   F = [1, f_lf(x)],   D = diag(sigma_i^2).
Because D is not proportional to R, s^2 has no closed form. (s^2, l) are
estimated by REML (restricted likelihood; it has n - p = 1 degree of freedom
here, and ML would divide by n and bias s^2 low). beta is GLS. The upstream
code's scalar nugget and its likelihood term -0.5 n sigma2 are not used.
Because D contains the grid error, the GP mean estimates q_exact(x), and its
posterior std, including the uncertainty in beta, is the band plotted.

Treating the grid error as independent across x is an approximation: it is
smooth in x in reality, so the band is if anything wide at the HF points.

Score. The held-out HF points are covered when |mu - y| <= 1.96 sqrt(var + sigma_i^2).
They are printed, not drawn.

Usage:  uv run python scripts/plot_mf_fig13_unc.py tau_mean_rpm
        uv run python scripts/plot_mf_fig13_unc.py edr_mean_angle
"""
from __future__ import annotations

import argparse
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
from scripts import figstyle as fs                                   # noqa: E402
from scripts.plot_mf_fig13 import CSV, LooKRR, per_cycle_max        # noqa: E402
from scripts.diag_num_uncertainty import uncertainty                # noqa: E402

CONV = {8: ("fig13a_rm_mf_rpm{:g}", 0.6), 9: ("fig9_l9_rpm{:g}", 0.6),
        10: ("l10_fig13a_xlevel_rpm{:g}", 0.0)}
CASES = {
    "tau_mean_rpm": dict(
        xs=[17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5], lf_level=8, hf_level=10,
        lf_run="fig13a_rm_mf_rpm{:g}", hf_run="l10_fig13a_xlevel_rpm{:g}",
        lf_settled=0.6, hf_settled=0.0, hf_train=[17.5, 27.5, 37.5],
        column="tau_mean", kind="tau", conv_column="tau_mean", marker="tau",
        kim_csv=CSV / "shear_ediss_vs_frequency.csv", kim_skip=[1],
        kim_x="RPM", kim_y="tau_liq_mean",
        xlabel="rocking frequency [rpm]",
        ylabel=r"max spatially averaged wall shear $\langle\tau\rangle_{\max}$ [Pa]",
        out="mf_tau_mean_rpm_unc.png"),
    "edr_mean_angle": dict(
        xs=[2.0, 3.0, 4.0, 5.0, 6.0, 7.0], lf_level=7, hf_level=9,
        lf_run="mf_l7_th{:g}", hf_run="fig10_l9_th{:g}",
        lf_settled=0.6, hf_settled=0.6, hf_train=[2.0, 4.0, 7.0],
        column="ediss_kim_mean", kind="ediss", marker="ediss",
        # the L6-L10 rpm runs predate the ediss_kim_* columns; ediss_mean is
        # the same spatial mean over the liquid
        conv_column="ediss_mean", conv_rpm=32.5,
        kim_csv=CSV / "shear_ediss_vs_angle.csv", kim_skip=None,
        kim_x="theta_deg", kim_y="Ediss_liq_mean",
        xlabel=r"rocking angle $\theta_{\max}$ [deg]",
        ylabel=r"max spatially averaged EDR $\langle\varepsilon\rangle_{\max}$ [W m$^{-3}$]",
        out="mf_edr_mean_angle_unc.png"),
}


def stat(run, column, kind, settled):
    m = per_cycle_max(run, column, kind, settled)
    return m.mean(), m.std(ddof=1) / math.sqrt(len(m))


def rel_u(rpm, column, kind, target):
    """U / q at level `target`, from the L8/L9/L10 fit at one rpm (7 deg)."""
    q = [stat(t.format(rpm), column, kind, s)[0] for t, s in CONV.values()]
    _, _, U = uncertainty(list(CONV), q, target=target)
    return U / q[list(CONV).index(target)]


def dataset(cs, which):
    rows = []
    for x in cs["xs"]:
        y, se = stat(cs[f"{which}_run"].format(x), cs["column"], cs["kind"],
                     cs[f"{which}_settled"])
        rows.append(dict(x=x, y=y, se=se))
    d = pd.DataFrame(rows)
    if which == "hf":
        if "conv_rpm" in cs:
            r = rel_u(cs["conv_rpm"], cs["conv_column"], cs["kind"], cs["hf_level"])
            d["U"] = r * d.y
        else:
            # A per-point U comes from a fit with one degree of freedom and scatters
            # 2-44% between neighbouring rpm, so it is itself a noisy estimate.
            # Pool the relative U as an RMS over all rpm, fixed before scoring.
            r = np.array([rel_u(x, cs["conv_column"], cs["kind"], cs["hf_level"])
                          for x in d.x])
            d["U_point"] = r * d.y
            d["U"] = math.sqrt(np.mean(r ** 2)) * d.y
        d["sig"] = np.sqrt(d.se ** 2 + (d.U / 1.96) ** 2)
    return d


class HeteroUK:
    """Universal kriging, known diagonal noise D, (s^2, l) by REML."""

    def __init__(self, lo, hi, basis):
        self.lo, self.hi, self.basis = lo, hi, basis

    def _R(self, a, b, l):
        d = (a[:, None] - b[None, :]) / (self.hi - self.lo)
        return np.exp(-0.5 * (d / l) ** 2)

    def _solve(self, logs2, logl):
        s2, l = math.exp(logs2), math.exp(logl)
        C = s2 * self._R(self.x, self.x, l) + np.diag(self.D)
        Ci = np.linalg.inv(C)
        A = self.F.T @ Ci @ self.F
        beta = np.linalg.solve(A, self.F.T @ Ci @ self.y)
        return s2, l, C, Ci, A, beta

    def _nreml(self, p):
        _, _, C, Ci, A, beta = self._solve(*p)
        e = self.y - self.F @ beta
        return 0.5 * (np.linalg.slogdet(C)[1] + np.linalg.slogdet(A)[1] + e @ Ci @ e)

    def fit(self, x, y, sig):
        self.x, self.y, self.D = np.asarray(x, float), np.asarray(y, float), np.asarray(sig) ** 2
        self.F = self.basis(self.x)
        v = np.var(self.y) + 1e-30
        bounds = [(math.log(v * 1e-4), math.log(v * 1e2)), (math.log(0.1), math.log(2.0))]
        best = None
        for s0 in np.linspace(*bounds[0], 5):
            for l0 in np.linspace(*bounds[1], 4):
                r = minimize(self._nreml, [s0, l0], method="L-BFGS-B", bounds=bounds)
                if best is None or r.fun < best.fun:
                    best = r
        self.s2, self.l, self.C, self.Ci, self.A, self.beta = self._solve(*best.x)
        self.e = self.y - self.F @ self.beta
        return self

    def predict(self, xs):
        xs = np.asarray(xs, float)
        c = self.s2 * self._R(xs, self.x, self.l)
        f = self.basis(xs)
        mu = f @ self.beta + c @ self.Ci @ self.e
        u = f - c @ self.Ci @ self.F
        var = (self.s2 - np.einsum("ij,jk,ik->i", c, self.Ci, c)
               + np.einsum("ij,jk,ik->i", u, np.linalg.inv(self.A), u))
        return mu, np.sqrt(np.maximum(var, 0))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=sorted(CASES))
    cs = CASES[ap.parse_args().case]
    lf, hf = dataset(cs, "lf"), dataset(cs, "hf")
    train, test = hf[hf.x.isin(cs["hf_train"])], hf[~hf.x.isin(cs["hf_train"])]
    lo, hi = min(cs["xs"]), max(cs["xs"])

    krr = LooKRR(design_space=np.array([[lo, hi]]), params_optimize=True,
                 noise_data=True, optimizer_restart=10, seed=42)
    krr.train(lf.x.to_numpy().reshape(-1, 1), lf.y.to_numpy())
    f_lf = lambda x: np.ravel(krr.predict(np.asarray(x, float).reshape(-1, 1)))

    mf = HeteroUK(lo, hi, lambda x: np.column_stack([np.ones_like(x), f_lf(x)]))
    sf = HeteroUK(lo, hi, lambda x: np.ones((len(x), 1)))
    for m in (mf, sf):
        m.fit(train.x, train.y, train.sig)

    print(f"{cs['out']}: LF L{cs['lf_level']} x{len(lf)}, HF L{cs['hf_level']} train "
          f"{list(train.x)}")
    print(hf.assign(U_rel=hf.U / hf.y, se_rel=hf.se / hf.y).round(5).to_string(index=False))
    print(f"  MF beta {np.round(mf.beta, 4)}, s {math.sqrt(mf.s2):.3g}, "
          f"l {mf.l * (hi - lo):.3g} (x units)")
    print(f"  SF beta {np.round(sf.beta, 4)}, s {math.sqrt(sf.s2):.3g}, "
          f"l {sf.l * (hi - lo):.3g}")
    for m, name in ((mf, "MF"), (sf, "HF-only GP")):
        mu, sd = m.predict(test.x)
        err = mu - test.y.to_numpy()
        z = err / np.sqrt(sd ** 2 + test.sig.to_numpy() ** 2)
        print(f"  {name:10s} held-out RMSE {np.sqrt(np.mean(err ** 2)):.4g}  "
              f"95% coverage {np.mean(np.abs(z) <= 1.96):.0%}  |z| {np.round(np.abs(z), 2)}")
    e_lf = f_lf(test.x) - test.y.to_numpy()
    print(f"  {'KRR on LF':10s} held-out RMSE {np.sqrt(np.mean(e_lf ** 2)):.4g}")

    xg = np.linspace(lo, hi, 200)
    mu, sd = mf.predict(xg)
    kim = pd.read_csv(cs["kim_csv"], skiprows=cs["kim_skip"])
    kim[cs["kim_x"]] = pd.to_numeric(kim[cs["kim_x"]])
    kim = kim.sort_values(cs["kim_x"])

    plt.rcParams.update(fs.rcparams())
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    cl, ch, mk = fs.level_colour(cs["lf_level"]), fs.level_colour(cs["hf_level"]), fs.MK[cs["marker"]]
    ax.fill_between(xg, mu - 1.96 * sd, mu + 1.96 * sd, color=ch, alpha=0.18, lw=0,
                    label="MF 95%")
    ax.plot(xg, f_lf(xg), color=cl, lw=1.1, ls="--", label=f"KRR on L{cs['lf_level']}")
    ax.plot(xg, mu, color=ch, lw=1.4, ls="--", label="MF KRR-LR-GPR")
    ax.plot(lf.x, lf.y, **fs.series_kw(cl, mk, stat="mean", ours=True, ls="none"),
            label=f"L{cs['lf_level']}")
    ax.errorbar(train.x, train.y, yerr=1.96 * train.sig, fmt="none", ecolor=ch, lw=0.8, capsize=2)
    ax.plot(train.x, train.y, **fs.series_kw(ch, mk, stat="mean", ours=True, ls="none"),
            label=f"L{cs['hf_level']}")
    ax.plot(kim[cs["kim_x"]], kim[cs["kim_y"]], **fs.series_kw(fs.KIM, mk, stat="mean"),
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
