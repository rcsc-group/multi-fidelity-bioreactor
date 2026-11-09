"""Two figures on 'does L8 inform L9?' (user 2026-11-09: show figures, not prose).
 A  bo_level_raw.png: raw data. Row 1 kLa_1T_10, row 2 tau_95. Left: vs rpm at theta 7 (L7, L8, L9 lines).
    Right: L8 vs L9 at identical (rpm, theta), colour = theta, dashed y = x.
 B  bo_level_lcurve.png: learning curve. Train on n9 random L9 points (+ all L8 for MF), predict the rest;
    median held-out RMSE (ln) over 30 draws vs n9; one panel per QoI. Solid L9-only, dotted MF ([1, f8] + GP).
Usage: uv run python scripts/bo2d/fig_level_info.py"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import loop as L  # noqa: E402
from scripts.bo2d import bo  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import figstyle as fs  # noqa: E402

C = L.C
OUT = Path(__file__).resolve().parents[2] / "talks" / "dfd2026" / "public"
QN = {"kLa_1T_10": r"$k_La$ (1/h)", "tau_95": r"$\tau_{95}$ (-)", "edr_95": r"EDR$_{95}$ (-)",
      "dtmix_0.50": r"$t_{mix,50}$ (s)"}
MS = {7: 4, 8: 7, 9: 11}


class NAR(bo.FullBayesKRRLRGP):
    """Nonlinear autoregressive MF (NARGP form, Perdikaris et al. 2017): basis [1, f8] and a GP over (rpm, theta, f8),
    so the L8 -> L9 map may bend and change with (rpm, theta)."""
    f_lo, f_hi = 0.0, 1.0

    def normalize_input(self, X):
        X = np.atleast_2d(X)
        xn = (X - self.bounds[:, 0]) / (self.bounds[:, 1] - self.bounds[:, 0])
        f = np.asarray(self.lf_model.predict(X), float).reshape(-1, 1)
        return np.hstack([xn, (f - self.f_lo) / (self.f_hi - self.f_lo)])


def data(q):
    db = C.load_db(q)[["level", "rpm", "theta", "run", "y"]]
    extra = L.rows(q, ["bo1_", "bo4_", "bo_l7_"], {7, 8, 9})[["level", "rpm", "theta", "run", "y"]]
    d = pd.concat([db, extra]).drop_duplicates("run")
    return d[np.isfinite(d.y.astype(float)) & (d.y.astype(float) > 0)]


def fig_raw():
    plt.rcParams.update(fs.rcparams())
    f, ax = plt.subplots(2, 2, figsize=(8, 6.4), constrained_layout=True)
    for r, q in enumerate(("kLa_1T_10", "tau_95")):
        d = data(q)
        a = ax[r, 0]
        for lev in (7, 8, 9):
            s = d[(d.level == lev) & (abs(d.theta - 7) < 0.01)].groupby("rpm").y.mean()
            a.plot(s.index, s.values, "-o", color=fs.level_colour(lev), ms=MS[lev], label=f"L{lev}")
        a.set_yscale("log"); a.set_xlabel("rpm"); a.set_ylabel(QN[q])
        if r == 0:
            a.set_title(r"$\theta = 7^\circ$"); a.legend(frameon=False)
        a = ax[r, 1]
        d8, d9 = d[d.level == 8], d[d.level == 9]
        pr = []
        for _, x in d9.iterrows():
            m = d8[(abs(d8.rpm - x.rpm) < 0.05) & (abs(d8.theta - x.theta) < 0.02)]
            if len(m):
                pr.append((x.theta, float(m.y.mean()), x.y))
        p = np.array(pr)
        sc = a.scatter(p[:, 1], p[:, 2], c=p[:, 0], cmap="viridis", vmin=2, vmax=7, s=45, edgecolor="k", lw=0.4)
        lo, hi = p[:, 1:].min() * 0.8, p[:, 1:].max() * 1.25
        a.plot([lo, hi], [lo, hi], "--", color=fs.NEUTRAL, lw=1, label="L8 = L9")
        a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(lo, hi); a.set_ylim(lo, hi)
        a.set_xlabel(f"L8 {QN[q]}"); a.set_ylabel(f"L9 {QN[q]}")
        if r == 0:
            a.set_title("same (rpm, θ)"); a.legend(frameon=False, loc="upper left")
    f.colorbar(sc, ax=ax[:, 1], label=r"$\theta$ (deg)", shrink=0.6)
    f.savefig(OUT / "bo_level_raw.png", dpi=200); print("wrote", OUT / "bo_level_raw.png")


def lcurve(q, ns=(4, 6, 8, 12, 16), draws=30):
    d = data(q); floor = C.noise_floor(q)
    d8, d9 = d[d.level == 8], d[d.level == 9]
    X8, y8 = d8[["rpm", "theta"]].to_numpy(float), np.log(d8.y.to_numpy(float))
    X9, y9 = d9[["rpm", "theta"]].to_numpy(float), np.log(d9.y.to_numpy(float))
    m8 = bo.FullBayesKRRLRGP(C.DOMAIN, basis="ordinary", grid=C.hf_grid(floor, float(np.std(y8)))).train(X8, y8)
    lf = bo._Mean(m8)
    NAR.f_lo, NAR.f_hi = float(y8.min()), float(y8.max())
    rng = np.random.default_rng(0)
    out = {"L9-only": [], "MF": [], "NAR": []}
    for n in ns:
        e = {k: [] for k in out}
        for _ in range(draws):
            tr = rng.choice(len(y9), n, replace=False); te = np.setdiff1d(np.arange(len(y9)), tr)
            g = C.hf_grid(floor, max(float(np.std(y9[tr])), 1e-3))
            for k, m in (("L9-only", bo.FullBayesKRRLRGP(C.DOMAIN, basis="ordinary", grid=g)),
                         ("MF", bo.FullBayesKRRLRGP(C.DOMAIN, lf_model=lf, grid=g)),
                         ("NAR", NAR(C.DOMAIN, lf_model=lf, grid=g))):
                m.train(X9[tr], y9[tr])
                e[k].append(np.sqrt(np.mean((m.predict_mean(X9[te]) - y9[te]) ** 2)))
        for k in out:
            out[k].append(np.percentile(e[k], [25, 50, 75]))
        print(q, n, {k: round(float(np.median(e[k])), 3) for k in e})
    return np.array(ns), {k: np.array(v) for k, v in out.items()}, float(np.std(y9))


def fig_lcurve():
    plt.rcParams.update(fs.rcparams())
    qs = ["kLa_1T_10", "dtmix_0.50", "tau_95", "edr_95"]
    f, ax = plt.subplots(1, 4, figsize=(12, 3.2), constrained_layout=True)
    for a, q in zip(ax, qs):
        ns, out, sd = lcurve(q)
        lab = {"L9-only": "L9 only", "MF": "L8 + L9, linear map", "NAR": "L8 + L9, nonlinear map"}
        for k, ls in (("L9-only", "-"), ("MF", ":"), ("NAR", "-.")):
            a.plot(ns, out[k][:, 1], ls, color=fs.level_colour(9), lw=2, label=lab[k])
            a.fill_between(ns, out[k][:, 0], out[k][:, 2], color=fs.level_colour(9), alpha=0.12)
        a.axhline(sd, color=fs.NEUTRAL, lw=1, ls="--", label="guess the mean")
        a.set_title(QN[q].split(" (")[0]); a.set_xlabel("L9 runs used"); a.set_ylim(0, None); a.set_xticks([4, 8, 12, 16])
    ax[0].set_ylabel("held-out error (ln)"); ax[0].legend(frameon=False)
    f.savefig(OUT / "bo_level_lcurve.png", dpi=200); print("wrote", OUT / "bo_level_lcurve.png")


if __name__ == "__main__":
    if "--lcurve-only" not in sys.argv:
        fig_raw()
    fig_lcurve()
