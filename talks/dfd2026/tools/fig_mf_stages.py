"""Four click stages of the multi-fidelity fit (slide 4). One curve per stage, no ghost curves.
Stage k shows the L8 data and KRR fit f_L, the L10 data, ONE curve rho f_L + b + delta at a step of the NLL
optimisation, and the hyperparameter values at that step (the user asked to keep them).
Stages (one per click): 1 points only; 2 + L8 KRR f_L; 3 high fidelity enters (start, rho = 1, b = 0);\n4-5 iterations at 1/3 and 2/3 of the curve path; 6 final fit with 95 % band. Values panel from stage 3 on.
Derived from fig_mf_steps.py (same data, same optimiser). Writes public/mf_stage1..6.png.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from sklearn.kernel_ridge import KernelRidge

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO)]
from scripts import figstyle as fs  # noqa: E402
from scripts.fig_mf_grid import scalar_data  # noqa: E402

OUT = Path(__file__).resolve().parent.parent / "public"
Q, TH = "kLa_10", 7.0
C8, C10 = fs.level_colour(8), fs.level_colour(10)
MUTE = "#7a8797"
plt.rcParams.update({"mathtext.fontset": "cm", "font.family": "serif", "font.serif": ["CMU Serif", "DejaVu Serif"]})

SYN = "--real" not in sys.argv
if SYN:  # synthetic: truth log kLa(x); L8 = 1.3 x truth shifted down, plus a smooth bias, plus noise; L10 = truth plus small noise
    rng = np.random.default_rng(3)
    truth = lambda x: 2.3 + 0.9 * np.tanh((x - 27) / 4) + 0.15 * np.sin(x / 2.2)  # noqa: E731
    x8 = np.sort(rng.uniform(15, 37.5, 16))
    y8 = 0.55 * truth(x8) + 2.2 + 0.35 * np.sin(x8 / 3.0) + rng.normal(0, 0.04, x8.size)  # far from L10: other level, other slope, a bias
    x10 = np.array([16.5, 21.0, 26.0, 30.5, 36.0])
    y10 = truth(x10) + rng.normal(0, 0.02, x10.size)
else:
    (X8, y8), _, (X10, y10), _ = scalar_data(Q)
    x8, y8 = X8[X8[:, 1] == TH, 0], np.log(y8[X8[:, 1] == TH])
    x10, y10 = X10[X10[:, 1] == TH, 0], np.log(y10[X10[:, 1] == TH])
G = np.linspace(min(x8.min(), x10.min()) - 0.5, max(x8.max(), x10.max()) + 0.5, 300)
sc = lambda x: (x - 25.0) / 10.0  # noqa: E731

# KRR on L8, length scale by leave-one-out
best = None
for g in np.logspace(-1, 2, 40):
    for a in (1e-4, 1e-3, 1e-2):
        err = 0.0
        for i in range(len(x8)):
            m = np.arange(len(x8)) != i
            k = KernelRidge(kernel="rbf", gamma=g, alpha=a).fit(sc(x8[m])[:, None], y8[m] - y8.mean())
            err += (k.predict(sc(x8[i:i + 1])[:, None])[0] + y8.mean() - y8[i]) ** 2
        best = min(best or (np.inf,), (err, g, a))
krr = KernelRidge(kernel="rbf", gamma=best[1], alpha=best[2]).fit(sc(x8)[:, None], y8 - y8.mean())
fL = lambda x: krr.predict(sc(np.atleast_1d(x))[:, None]) + y8.mean()  # noqa: E731


def rbf(a, b, l):
    return np.exp(-0.5 * (sc(a)[:, None] - sc(b)[None, :]) ** 2 / l ** 2)


def parts(th):
    rho, b = th[:2]
    l, s, n = np.exp(th[2:])
    K = s ** 2 * rbf(x10, x10, l) + (n ** 2 + 1e-8) * np.eye(len(x10))
    L = np.linalg.cholesky(K)
    Ki = lambda v: np.linalg.solve(L.T, np.linalg.solve(L, v))  # noqa: E731
    beta = np.array([rho, b])  # rho and b are optimised with the GP hyperparameters, so the steps move the curve
    r = y10 - rho * fL(x10) - b
    return l, s, n, L, Ki, beta, r


def nll(th):
    *_, L, Ki, beta, r = parts(th)
    return 0.5 * r @ Ki(r) + np.log(np.diag(L)).sum()


def predict(th):
    l, s, n, L, Ki, beta, r = parts(th)
    ks = s ** 2 * rbf(G, x10, l)
    mean = beta[0] * fL(G) + beta[1] + ks @ Ki(r)
    var = s ** 2 - np.einsum("ij,ji->i", ks, Ki(ks.T))
    return mean, np.sqrt(np.clip(var, 0, None))


steps = [np.r_[1.0, 0.0, np.log([0.3, 0.3, 0.3])]]  # start: rho = 1, b = 0, i.e. the L8 curve itself
res = minimize(nll, steps[0], method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 4000}, callback=lambda t: steps.append(t.copy()))
print(f"L8 {len(x8)} pts, L10 {len(x10)} pts; KRR gamma {best[1]:.3g}; NLL {nll(steps[0]):.2f} -> {res.fun:.2f} in {len(steps) - 1} steps; rho, b = {parts(res.x)[5]}")


def base(ax):
    ax.set_yscale("log")
    if SYN:
        ax.set_yticks([]); ax.set_xticks([]); ax.minorticks_off()
    ax.set_xlim(G[0], G[-1])
    ax.set_ylim(np.exp(min(y8.min(), y10.min())) / 1.6, np.exp(max(y8.max(), y10.max())) * 1.6)
    ax.set_xlabel(r"$x$" if SYN else r"rocking speed (rpm)", color=MUTE, fontsize=14)
    ax.set_ylabel(r"$f$" if SYN else r"$k_La$ (1/h)", color=MUTE, fontsize=14)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTE)
    ax.tick_params(colors=MUTE, labelsize=12)
    ax.plot(x8, np.exp(y8), "o", color=C8, ms=6, mew=0, zorder=4)
    ax.plot(x10, np.exp(y10), "o", color=C10, ms=11, mew=0, zorder=5)


PCOL = {"rho": "#D55E00", "b": "#117733", "l": "#882255", "s": "#8c6d31", "n": "#555555"}
PTEX = {"rho": r"$\rho$", "b": r"$b$", "l": r"$\ell$", "s": r"$s$", "n": r"$\sigma_n$"}


def vals(th):
    return {"rho": th[0], "b": th[1], "l": np.exp(th[2]), "s": np.exp(th[3]), "n": np.exp(th[4])}


for frame in "ab":
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    base(ax)
    if frame == "b":
        ax.plot(G, np.exp(fL(G)), color=C8, lw=2.2, zorder=3)
    fig.subplots_adjust(left=0.08, right=0.74, bottom=0.12, top=0.97)
    fig.savefig(OUT / f"mf_{frame}.png", dpi=200, transparent=True)
    plt.close(fig)

f = np.array([nll(t) for t in steps])
means = np.array([predict(t)[0] for t in steps])
path = np.r_[0, np.cumsum(np.sqrt(((np.diff(means, axis=0)) ** 2).mean(1)))]
it = [0, int(np.argmax(path >= path[-1] / 3)), int(np.argmax(path >= 2 * path[-1] / 3)), len(steps) - 1]
# stage 1: points only; 2: + KRR f_L; 3: high fidelity enters (start of the optimisation); 4-5: iterations; 6: final fit + band
stages = [("pts", None), ("krr", None)] + [("hf", i) for i in it]
for k, (kind, idx) in enumerate(stages, 1):
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    base(ax)
    if kind != "pts":
        ax.plot(G, np.exp(fL(G)), color=C8, lw=2.2, zorder=3)
    if kind == "hf":
        m, sd = predict(steps[idx])
        ax.plot(G, np.exp(m), color=C10, lw=3.0, zorder=2)
        if idx == it[-1]:
            ax.fill_between(G, np.exp(m - 1.96 * sd), np.exp(m + 1.96 * sd), color=C10, alpha=0.15, lw=0, zorder=1)
        v = vals(steps[idx])
        fig.text(0.77, 0.85, r"$-\log p(\mathbf{y}_H\,|\,\theta)$" + f" = {f[idx]:.1f}", color=C10, fontsize=11)
        for r, key in enumerate(PCOL):
            fig.text(0.77, 0.72 - 0.09 * r, f"{PTEX[key]} = {v[key]:.2f}", color=PCOL[key], fontsize=14)
    fig.subplots_adjust(left=0.08, right=0.74, bottom=0.12, top=0.97)
    fig.savefig(OUT / f"mf_stage{k}.png", dpi=200, facecolor="white")
    plt.close(fig)
print("iteration stages at steps", it, "of", len(steps) - 1)
