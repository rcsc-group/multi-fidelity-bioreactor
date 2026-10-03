"""Falsifiable test of the h-kernel for grid extrapolation, one rpm at a time.

Model (Tuo-Wu-Yu decomposition, read in CONFIG 2209.13748 eq 6):
    z(h) = mu + delta(h) + e,   e ~ N(0, s^2),   flat prior on mu.
Two covariances for delta, both restated in Bect et al. 2103.14559 sec. 3:
    BM   (TWY1):  sigma^2 min(h, h')^(2p)          independent increments
    TWY2       :  sigma^2 (h h')^p  c(h - h')      Richardson-type, c = Matern-5/2
Prediction at h = 0 is the posterior of mu (PRE 2401.07562 eq 13).
Hyperparameters p, sigma, ell, s are marginalised on a grid (fully Bayesian):
    log p ~ N(log 1.5, 0.5^2), log sigma ~ U, log ell ~ U, log s ~ N(log 0.03, 0.7^2).
The noise prior s is an assumption (we have one replicate measurement only);
it is stated, not measured.

Output z = log(dtmix_chi) or -1/dtmix_chi (rate); h = cell size / L10 cell size.
Tests:
  T1 fit L6-L8, predict L9 (z-score, inside the 95% band?);
  T2 fit L6-L9, predict the L10 point at 32.5 rpm (one point; its protocol is a
     warm-start chain, sigma2_max = 0.20, so a miss is ambiguous).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs                       # noqa: E402
from scripts.pilot_truncated_chi import RPMS             # noqa: E402

FAM = {6: "kmix_l6_rpm{:g}", 7: "kmix_l7_rpm{:g}", 8: "kmix_l8_rpm{:g}", 9: "fig9_l9_rpm{:g}"}
H = {6: 16.0, 7: 8.0, 8: 4.0, 9: 2.0, 10: 1.0}
L10 = {32.5: "fig9_l10b_seg2"}

P = np.exp(np.linspace(np.log(0.3), np.log(4.0), 15))
SIG = np.exp(np.linspace(np.log(1e-3), np.log(1e2), 15))
ELL = np.exp(np.linspace(np.log(1.0), np.log(100.0), 7))
S = np.exp(np.linspace(np.log(0.005), np.log(0.2), 6))


def matern52(r, ell):
    a = np.sqrt(5) * np.abs(r) / ell
    return (1 + a + a * a / 3) * np.exp(-a)


def kdelta(kind, h1, h2, p, sig, ell):
    H1, H2 = np.meshgrid(h1, h2, indexing="ij")
    if kind == "BM":
        return sig ** 2 * np.minimum(H1, H2) ** (2 * p)
    return sig ** 2 * (H1 * H2) ** p * matern52(H1 - H2, ell)


def transform(y, tr):
    return np.log(y) if tr == "log" else -1.0 / y


def posterior(kind, h, z, hstar, s_scale, obs=True):
    """Mixture posterior of z at the cell sizes hstar (0 -> mu), moment-matched.

    obs=True adds the probe noise (predicts a new run); obs=False is the latent value.
    """
    hstar = np.atleast_1d(np.asarray(hstar, float))
    lws, ms, vs = [], [], []
    one = np.ones(h.size)
    ells = ELL if kind == "TWY2" else ELL[:1]
    for p in P:
        lp_p = norm.logpdf(np.log(p), np.log(1.5), 0.5)
        for sig in SIG:
            for ell in ells:
                Kd = kdelta(kind, h, h, p, sig, ell)
                k = kdelta(kind, hstar, h, p, sig, ell)                  # (m, n)
                kss = np.diag(kdelta(kind, hstar, hstar, p, sig, ell))
                for s in S:
                    lp_s = norm.logpdf(np.log(s), np.log(0.03), 0.7)
                    n2 = (s * s_scale) ** 2
                    K = Kd + n2 * np.eye(h.size)
                    try:
                        L = np.linalg.cholesky(K)
                    except np.linalg.LinAlgError:
                        continue
                    Ki = np.linalg.inv(K)
                    a = one @ Ki @ one
                    mu = one @ Ki @ z / a
                    r = z - mu
                    ll = -0.5 * r @ Ki @ r - np.log(np.diag(L)).sum() - 0.5 * np.log(a)
                    u = 1 - k @ Ki @ one
                    m = mu + k @ Ki @ r
                    v = kss - np.einsum("ij,jk,ik->i", k, Ki, k) + u * u / a + (n2 if obs else 0.0)
                    v = np.where(hstar == 0.0, 1.0 / a, v)       # delta(0) = 0 exactly
                    lws.append(ll + lp_p + lp_s); ms.append(m); vs.append(v)
    lw = np.array(lws)
    w = np.exp(lw - lw.max()); w /= w.sum()
    M, V = np.array(ms), np.array(vs)
    mean = w @ M
    var = w @ (V + M ** 2) - mean ** 2
    return mean, var


def data(chi):
    rows = {}
    for r in RPMS:
        d = {}
        for lv, t in FAM.items():
            f = ROOT / "runs" / t.format(r) / "results.json"
            if f.exists():
                y = json.load(open(f)).get(f"dtmix_{chi:.2f}", np.nan)
                if np.isfinite(y):
                    d[lv] = y
        rows[r] = d
    return rows


def main(chi=0.95):
    rows = data(chi)
    res = {}
    for tr in ("log", "rate"):
        for kind in ("BM", "TWY2"):
            zs, cover, lpd = [], [], []
            for r, d in rows.items():
                if not all(l in d for l in (6, 7, 8, 9)):
                    continue
                h = np.array([H[l] for l in (6, 7, 8)])
                y = np.array([d[l] for l in (6, 7, 8)])
                z = transform(y, tr)
                scale = 1.0 if tr == "log" else 1.0 / d[8]     # s is relative
                m, v = posterior(kind, h, z, H[9], scale)
                m, v = m[0], v[0]
                z9 = transform(d[9], tr)
                zz = (z9 - m) / np.sqrt(v)
                zs.append(zz)
                cover.append(abs(zz) < 1.96)
                # log density of dtmix itself (Jacobian of the transform)
                jac = 1.0 / d[9] if tr == "log" else 1.0 / d[9] ** 2
                lpd.append(norm.logpdf(z9, m, np.sqrt(v)) + np.log(jac))
            res[(tr, kind)] = (np.array(zs), np.mean(cover), np.sum(lpd))
            print(f"T1 chi={chi} {tr:4s} {kind:4s}: L9 inside 95% band {np.mean(cover):.0%} "
                  f"(n={len(zs)}), z = {np.round(zs, 1)}, sum log p(dtmix_L9) = {np.sum(lpd):.1f}")
    # T2 + figure: fit L6-L9 (log), show h -> 0 and L10
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(2, 5, figsize=(14, 5.6), sharex=True)
    hg = np.linspace(0.0, 17.0, 60)
    cols = {"BM": "#9467bd", "TWY2": "#d62728"}
    for ax, r in zip(axes.ravel(), RPMS):
        d = rows[r]
        lv = [l for l in (6, 7, 8, 9) if l in d]
        h = np.array([H[l] for l in lv])
        z = np.log([d[l] for l in lv])
        for kind in ("BM", "TWY2"):
            mm, vv = posterior(kind, h, z, hg, 1.0, obs=False)
            ss = np.sqrt(vv)
            ax.fill_between(hg, np.exp(mm - 1.96 * ss), np.exp(mm + 1.96 * ss),
                            color=cols[kind], alpha=0.15, lw=0)
            ax.plot(hg, np.exp(mm), color=cols[kind], lw=1.0, label=kind)
            if r in L10:
                y10 = json.load(open(ROOT / "runs" / L10[r] / "results.json"))[f"dtmix_{chi:.2f}"]
                m, v = posterior(kind, h, z, 1.0, 1.0)
                m, v = m[0], v[0]
                print(f"T2 {kind}: L10 at {r} rpm measured {y10:.1f} s, predicted "
                      f"{np.exp(m):.1f} s [{np.exp(m - 1.96 * np.sqrt(v)):.1f}, "
                      f"{np.exp(m + 1.96 * np.sqrt(v)):.1f}], z = {(np.log(y10) - m) / np.sqrt(v):+.1f}")
        ax.plot(h, np.exp(z), "ko", ms=3.5, label="L6-L9")
        if r in L10:
            ax.plot(1.0, json.load(open(ROOT / "runs" / L10[r] / "results.json"))[f"dtmix_{chi:.2f}"],
                    "kx", ms=7, mew=1.5, label="L10")
        ax.set_yscale("log")
        ax.set_ylim(np.exp(z).min() / 4, np.exp(z).max() * 12)
        ax.set_title(f"{r:g} rpm", fontsize=10, loc="left")
        ax.grid(**fs.GRID_KW)
    for ax in axes[-1]:
        ax.set_xlabel(r"cell size $h/h_{10}$")
    for ax in axes[:, 0]:
        ax.set_ylabel(rf"$\Delta t_{{{chi}}}$ (s)")
    hs_, ls_ = [], []
    for ax in axes.ravel():
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in ls_:
                hs_.append(hh); ls_.append(ll)
    axes[0, -1].legend(hs_, ls_, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    out = ROOT / f"experiments/multifidelity/h_kernel_test_chi{int(chi * 100)}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)


if __name__ == "__main__":
    main(float(sys.argv[1]) if len(sys.argv) > 1 else 0.95)
