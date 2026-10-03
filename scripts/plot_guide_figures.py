"""Teaching figures for the Yi-h guide. Toy problems with known answers (not our data).

G1  a QoI computed at four grid levels converges to the exact curve as h -> 0
G2  Richardson extrapolation at one input: four grids, a power-law fit, the limit
G3  a Gaussian process: prior samples, then the posterior after five observations
G4  Yi et al.'s two-fidelity idea: many cheap points, few expensive points, a linear transfer
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs   # noqa: E402

OUT = ROOT / "experiments/multifidelity/guide"
LEVELS = {"L6": 1.0, "L7": 0.5, "L8": 0.25, "L9": 0.125}
COLS = ["#c6dbef", "#6baed6", "#2171b5", "#08306b"]


def truth(x):
    return 1.0 + 0.5 * np.sin(2 * np.pi * x)


def level(x, h):
    """The QoI at cell size h: exact value plus an error that shrinks like h^1.5."""
    return truth(x) + h ** 1.5 * (0.6 + 0.3 * np.cos(3 * x))


def g1():
    x = np.linspace(0, 1, 200)
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for (name, h), c in zip(LEVELS.items(), COLS):
        ax.plot(x, level(x, h), color=c, lw=1.6, label=f"{name} (h = {h:g})")
    ax.plot(x, truth(x), "k--", lw=1.4, label="exact (h = 0)")
    ax.set_xlabel("input x (for example rpm, scaled)")
    ax.set_ylabel("QoI")
    ax.grid(**fs.GRID_KW)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "g1_levels.png", dpi=150, bbox_inches="tight")


def g2():
    x0 = 0.3
    h = np.array(list(LEVELS.values()))
    y = level(x0, h)
    # Richardson with three finest levels: ratio of increments gives p, then the limit
    r = (y[1] - y[2]) / (y[2] - y[3])
    p = np.log2(r)
    lim = y[3] - (y[2] - y[3]) / (2 ** p - 1)
    hg = np.linspace(0, 1.05, 100)
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.plot(h, y, "ko", ms=6, label="runs at four grids")
    ax.plot(hg, lim + (y[3] - lim) * (hg / h[3]) ** p, color="#d62728", lw=1.4,
            label=rf"fit $f_0 + a\,h^{{p}}$, p = {p:.2f}")
    ax.plot(0, lim, "*", color="#d62728", ms=14, label=f"extrapolated limit {lim:.3f}")
    ax.axhline(truth(x0), color="0.5", lw=0.8, ls="--", label=f"exact {truth(x0):.3f}")
    for (name, hh), yy in zip(LEVELS.items(), y):
        ax.annotate(name, (hh, yy), textcoords="offset points", xytext=(6, -12), fontsize=9)
    ax.set_xlabel("cell size h (coarsest = 1)")
    ax.set_ylabel(f"QoI at x = {x0}")
    ax.grid(**fs.GRID_KW)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "g2_richardson.png", dpi=150, bbox_inches="tight")
    print(f"G2: values {np.round(y, 4)}, R = {r:.3f}, p = {p:.3f}, limit = {lim:.4f}, exact = {truth(x0):.4f}")


def k(a, b, ell=0.2, s=0.6):
    return s ** 2 * np.exp(-0.5 * (a[:, None] - b[None, :]) ** 2 / ell ** 2)


def g3():
    rng = np.random.default_rng(3)
    x = np.linspace(0, 1, 200)
    xd = np.array([0.08, 0.3, 0.45, 0.7, 0.92])
    yd = truth(xd) + 0.03 * rng.standard_normal(xd.size)
    m0 = 1.0
    K = k(xd, xd) + 0.03 ** 2 * np.eye(xd.size)
    Ks = k(x, xd)
    mean = m0 + Ks @ np.linalg.solve(K, yd - m0)
    cov = k(x, x) - Ks @ np.linalg.solve(K, Ks.T)
    sd = np.sqrt(np.clip(np.diag(cov), 0, None))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
    prior = rng.multivariate_normal(np.full(x.size, m0), k(x, x) + 1e-9 * np.eye(x.size), 4)
    for s_ in prior:
        a1.plot(x, s_, lw=1.0, color="#6baed6")
    a1.fill_between(x, m0 - 1.96 * 0.6, m0 + 1.96 * 0.6, color="#6baed6", alpha=0.15, lw=0)
    a1.set_title("before data: prior samples", loc="left", fontsize=10)
    post = rng.multivariate_normal(mean, cov + 1e-9 * np.eye(x.size), 4)
    a2.fill_between(x, mean - 1.96 * sd, mean + 1.96 * sd, color="#2171b5", alpha=0.18, lw=0)
    for s_ in post:
        a2.plot(x, s_, lw=0.8, color="#6baed6")
    a2.plot(x, mean, color="#08306b", lw=1.6, label="posterior mean")
    a2.plot(xd, yd, "ko", ms=5, label="observations")
    a2.set_title("after five observations: posterior", loc="left", fontsize=10)
    for ax in (a1, a2):
        ax.set_xlabel("input x")
        ax.grid(**fs.GRID_KW)
    a1.set_ylabel("QoI")
    a2.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "g3_gp.png", dpi=150, bbox_inches="tight")


def g4():
    """Yi et al.'s three steps on a toy: KRR on cheap data; rho by GLS; GP on the residual."""
    rng = np.random.default_rng(7)
    x = np.linspace(0, 1, 200)
    lo = lambda z: level(z, 1.0)          # cheap: the coarsest grid
    hi = lambda z: level(z, 0.125)        # expensive: the finest grid
    xl = np.linspace(0, 1, 40)
    yl = lo(xl) + 0.02 * rng.standard_normal(xl.size)
    xh = np.array([0.05, 0.35, 0.65, 0.95])
    yh = hi(xh) + 0.02 * rng.standard_normal(xh.size)
    # step 1: kernel ridge regression of the cheap data (RBF kernel, ridge 1e-3)
    kr = lambda a, b: np.exp(-0.5 * (a[:, None] - b[None, :]) ** 2 / 0.15 ** 2)
    alpha = np.linalg.solve(kr(xl, xl) + 1e-3 * np.eye(xl.size), yl)
    fl = lambda z: kr(z, xl) @ alpha
    # steps 2-3: residual GP r ~ GP(0, k), noise 0.02; rho by generalised least squares
    kres = lambda a, b: k(a, b, ell=0.3, s=0.3)
    Kh = kres(xh, xh) + 0.02 ** 2 * np.eye(xh.size)
    F = np.column_stack([np.ones(xh.size), fl(xh)])
    Ki = np.linalg.inv(Kh)
    rho = np.linalg.solve(F.T @ Ki @ F, F.T @ Ki @ yh)
    res = yh - F @ rho
    ks = kres(x, xh)
    mean = rho[0] + rho[1] * fl(x) + ks @ Ki @ res
    sd = np.sqrt(np.clip(0.3 ** 2 - np.einsum("ij,jk,ik->i", ks, Ki, ks), 0, None))
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.plot(xl, yl, ".", color="#9ecae1", ms=6, label="cheap runs (40)")
    ax.plot(x, fl(x), color="#6baed6", lw=1.2, label="step 1: KRR fit of cheap runs")
    ax.plot(x, rho[0] + rho[1] * fl(x), color="#fdae6b", lw=1.2, ls="-.",
            label="step 2: transfer only")
    ax.fill_between(x, mean - 1.96 * sd, mean + 1.96 * sd, color="#d62728", alpha=0.15, lw=0)
    ax.plot(x, mean, color="#d62728", lw=1.5, label="step 3: transfer + residual GP, 95%")
    ax.plot(xh, yh, "o", color="#08306b", ms=7, label="expensive runs (4)")
    ax.plot(x, hi(x), "k--", lw=1.0, label="expensive level (true)")
    ax.set_xlabel("input x")
    ax.set_ylabel("QoI")
    ax.grid(**fs.GRID_KW)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "g4_yi.png", dpi=150, bbox_inches="tight")
    err = np.abs(mean - hi(x)).max()
    print(f"G4: rho0 = {rho[0]:.3f}, rho1 = {rho[1]:.3f}, max |error| of step 3 = {err:.3f}, "
          f"true within band: {np.mean(np.abs(mean - hi(x)) <= 1.96 * sd):.0%}")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(fs.rcparams())
    g1(); g2(); g3(); g4()
    print("saved to", OUT)
