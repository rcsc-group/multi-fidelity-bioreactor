"""Kim's Fig 13 means, predicted multi-fidelity with grid uncertainty.

Not Kim replicas. 13(a) is vs rpm at 7 deg, 13(b) is vs angle at 32.5 rpm.
Each has two panels, Kim's hollow markers, defined exactly as in the replica
(plot_fig13.py):
  <tau> = half the peak-to-peak swing of the SIGNED spatial-mean shear
          (tau_mean_signed). Kim's <tau> oscillates through zero; Fig 13 plots its amplitude.
  <eps> = the peak of the spatial-mean EDR (ediss_mean).
Both columns are logged at every level L7-L10, so the grid fit below uses the
plotted statistic itself. The pointwise maxima are left out on purpose,
because they have no grid-converged limit here: tau_max grows 1.4-2.5x per
level (diag_richardson_tau.py), so a band around it would bound nothing.

LF L7 (mf_l7_rpm*, mf_l7_th*), HF L9 (fig9_l9_rpm*, fig10_l9_th*): same binary,
so they differ in level only. 3 HF points train, and the rest are scored,
printed but not drawn.
Response: the per-cycle amplitude (tau) or peak (eps), averaged over the settled cycles.

HF noise. sigma_i^2 = SE_i^2 + (U_i/1.96)^2. SE is the cycle-to-cycle standard error.
U is the ~95% bound on |q_L9 - q_exact| from least-squares grid fits
(diag_num_uncertainty.py, simplified Eca & Hoekstra 2014) over L8/L9/L10 at 7 deg.
It is evaluated at L9 per rpm, and U/q is pooled as an RMS over rpm, because
each per-rpm fit has one degree of freedom.
ASSUMPTION: 13(b) inherits the 7-deg U/q at every angle.

Transfer: universal kriging, y = F beta + r, Cov = s^2 R_l + D, F = [1, f_KRR],
D = diag(sigma_i^2) known. (s^2, l) are fitted by REML, beta by GLS. The band is the
posterior of q_exact, including the uncertainty in beta.
Baseline: the same GP with F = [1] on the HF points alone.

Usage:  uv run python scripts/plot_mf_fig13_unc.py 13a
        uv run python scripts/plot_mf_fig13_unc.py 13b
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
from scripts import figstyle as fs                                   # noqa: E402
from scripts.plot_mf_fig13 import CSV, LooKRR                       # noqa: E402
from scripts.diag_num_uncertainty import uncertainty                # noqa: E402

LF_LEVEL, HF_LEVEL = 7, 9
CONV = {8: ("fig13a_rm_mf_rpm{:g}", 0.6), 9: ("fig9_l9_rpm{:g}", 0.6),
        10: ("l10_fig13a_xlevel_rpm{:g}", 0.0)}
CONV_RPMS = [17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
QTY = {
    "tau": dict(column="tau_mean_signed", reduce="amp", kim_y="tau_liq_mean",
                ylabel=r"$\langle\tau\rangle$ [Pa]"),
    "ediss": dict(column="ediss_mean", reduce="max", kim_y="Ediss_liq_mean",
                  ylabel=r"$\langle\varepsilon\rangle$ [W m$^{-3}$]"),
}
FIGS = {
    "13a": dict(xs=CONV_RPMS, train=[17.5, 27.5, 37.5],
                lf_run="mf_l7_rpm{:g}", hf_run="fig9_l9_rpm{:g}",
                kim_csv=CSV / "shear_ediss_vs_frequency.csv", kim_skip=[1], kim_x="RPM",
                xlabel="rocking frequency [rpm]", title=r"$\theta_{\max}=7^\circ$",
                out="mf_fig13a_means.png"),
    "13b": dict(xs=[2.0, 3.0, 4.0, 5.0, 6.0, 7.0], train=[2.0, 4.0, 7.0],
                lf_run="mf_l7_th{:g}", hf_run="fig10_l9_th{:g}",
                kim_csv=CSV / "shear_ediss_vs_angle.csv", kim_skip=None, kim_x="theta_deg",
                xlabel=r"rocking angle $\theta_{\max}$ [deg]", title=r"$f_b=32.5$ rpm",
                out="mf_fig13b_means.png"),
}


def per_cycle(run, kind, settled):
    """Per-cycle statistic of QTY[kind], dimensional, settled cycles before release."""
    d = ROOT / "runs" / run
    p = json.load(open(d / "params.json"))
    L, H = p["geometry"]["a"], 2 * p["geometry"]["b"]
    Tp = 2 * math.pi / p["omega_b"]
    U = L / 4 * (H + 0.5 * L * math.tan(math.radians(p["theta_max"][0]))) / (H * 0.5) / Tp
    T_nd = Tp / (L / U)
    df = pd.read_csv(d / "shear_stress.dat", sep=r"\s+")
    t = df["t"].to_numpy()
    t_rel = (p.get("t_checkpoint") or 0) + p["n_mix_cycles"] * T_nd
    keep = t < t_rel if t_rel < t[-1] else np.ones_like(t, bool)
    t = t[keep]
    y = df[QTY[kind]["column"]].to_numpy()[keep]
    y = y * (1000 * U ** 2 if kind == "tau" else 1000 * U ** 3 / L)
    c = (t - t[0]) / T_nd
    k0, k1 = int(math.ceil(settled * c[-1])), int(math.floor(c[-1]))
    out = []
    for k in range(k0, k1):
        w = y[(c >= k) & (c < k + 1)]
        out.append(0.5 * (w.max() - w.min()) if QTY[kind]["reduce"] == "amp" else w.max())
    return np.array(out)


def stat(run, kind, settled=0.6):
    m = per_cycle(run, kind, settled)
    return m.mean(), m.std(ddof=1) / math.sqrt(len(m))


def pooled_rel_u(kind):
    """RMS over rpm of U(L9)/q(L9), from the L8/L9/L10 fit of the same statistic."""
    r = []
    for rpm in CONV_RPMS:
        q = [stat(t.format(rpm), kind, s)[0] for t, s in CONV.values()]
        _, _, U = uncertainty(list(CONV), q, target=HF_LEVEL)
        r.append(U / q[list(CONV).index(HF_LEVEL)])
    return math.sqrt(np.mean(np.square(r))), np.array(r)


def dataset(fig, kind, which, u_rel=None):
    rows = []
    for x in fig["xs"]:
        y, se = stat(fig[f"{which}_run"].format(x), kind)
        rows.append(dict(x=x, y=y, se=se))
    d = pd.DataFrame(rows)
    if u_rel is not None:
        d["U"] = u_rel * d.y
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


def panel(ax, fig, kind, u_rel, legend):
    lf, hf = dataset(fig, kind, "lf"), dataset(fig, kind, "hf", u_rel)
    train, test = hf[hf.x.isin(fig["train"])], hf[~hf.x.isin(fig["train"])]
    lo, hi = min(fig["xs"]), max(fig["xs"])

    krr = LooKRR(design_space=np.array([[lo, hi]]), params_optimize=True,
                 noise_data=True, optimizer_restart=10, seed=42)
    krr.train(lf.x.to_numpy().reshape(-1, 1), lf.y.to_numpy())
    f_lf = lambda x: np.ravel(krr.predict(np.asarray(x, float).reshape(-1, 1)))
    mf = HeteroUK(lo, hi, lambda x: np.column_stack([np.ones_like(x), f_lf(x)]))
    sf = HeteroUK(lo, hi, lambda x: np.ones((len(x), 1)))
    for m in (mf, sf):
        m.fit(train.x, train.y, train.sig)

    print(f"  [{kind}] r(L7,L9) {np.corrcoef(lf.y, hf.y)[0, 1]:+.3f}  "
          f"beta {np.round(mf.beta, 4)}  s {math.sqrt(mf.s2):.3g}  l {mf.l * (hi - lo):.3g}")
    for m, name in ((mf, "MF"), (sf, "L9-only GP")):
        mu, sd = m.predict(test.x)
        err = mu - test.y.to_numpy()
        z = err / np.sqrt(sd ** 2 + test.sig.to_numpy() ** 2)
        print(f"    {name:10s} held-out RMSE {np.sqrt(np.mean(err ** 2)):.4g} "
              f"({np.sqrt(np.mean((err / test.y) ** 2)):.1%})  "
              f"95% cov {np.mean(np.abs(z) <= 1.96):.0%}  |z| {np.round(np.abs(z), 2)}")
    e = f_lf(test.x) - test.y.to_numpy()
    print(f"    {'KRR on L7':10s} held-out RMSE {np.sqrt(np.mean(e ** 2)):.4g}")

    kim = pd.read_csv(fig["kim_csv"], skiprows=fig["kim_skip"])
    kim[fig["kim_x"]] = pd.to_numeric(kim[fig["kim_x"]])
    kim = kim.sort_values(fig["kim_x"])
    xg = np.linspace(lo, hi, 200)
    mu, sd = mf.predict(xg)
    cl, ch, mk = fs.level_colour(LF_LEVEL), fs.level_colour(HF_LEVEL), fs.MK[kind]
    ax.fill_between(xg, mu - 1.96 * sd, mu + 1.96 * sd, color=ch, alpha=0.18, lw=0,
                    label="MF 95%")
    ax.plot(xg, f_lf(xg), color=cl, lw=1.1, ls="--", label=f"KRR on L{LF_LEVEL}")
    ax.plot(xg, mu, color=ch, lw=1.4, ls="--", label="MF KRR-LR-GPR")
    ax.plot(lf.x, lf.y, **fs.series_kw(cl, mk, stat="mean", ours=True, ls="none"),
            label=f"L{LF_LEVEL}")
    ax.errorbar(train.x, train.y, yerr=1.96 * train.sig, fmt="none", ecolor=ch, lw=0.8,
                capsize=2)
    ax.plot(train.x, train.y, **fs.series_kw(ch, mk, stat="mean", ours=True, ls="none"),
            label=f"L{HF_LEVEL}")
    ax.plot(kim[fig["kim_x"]], kim[QTY[kind]["kim_y"]],
            **fs.series_kw(fs.KIM, mk, stat="mean"), label="Kim et al.")
    ax.set_ylabel(QTY[kind]["ylabel"])
    ax.grid(**fs.GRID_KW)
    if legend:
        ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("fig", choices=sorted(FIGS))
    name = ap.parse_args().fig
    fig = FIGS[name]
    plt.rcParams.update(fs.rcparams())
    f, axes = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True)
    print(f"Fig {name}: train {fig['train']}")
    for i, kind in enumerate(("tau", "ediss")):
        u, per = pooled_rel_u(kind)
        print(f"  [{kind}] U(L9)/q per rpm {np.round(per, 3)} -> pooled RMS {u:.3f}")
        panel(axes[i], fig, kind, u, legend=(i == 0))
    axes[0].set_title(fig["title"])
    axes[-1].set_xlabel(fig["xlabel"])
    out = ROOT / "experiments/multifidelity" / fig["out"]
    out.parent.mkdir(parents=True, exist_ok=True)
    f.savefig(out, dpi=200, bbox_inches="tight")
    print(f"  saved {out}")


if __name__ == "__main__":
    main()
