"""Which uncertainty decreases as h -> 0? Synthetic 1-D illustration (not our data).

Model: y(h) = mu + delta(h) + e, flat-ish prior on mu (variance 100),
delta ~ GP(0, s_d^2 (h h')^p c_52(h - h'; ell)) (TWY2), e ~ N(0, s^2).
Truth: f(h) = 1 + 0.4 h^1.5. Data at h = 1, 1/2, 1/4, 1/8 (four levels).

Quantities against h:
  (1) prior sd of the discretisation error, s_d h^p               -> 0, monotone
  (2) posterior RMS error of level h, sqrt(E[delta(h)^2 | D])     -> 0 ("fidelity")
  (3) posterior sd of f(h), sqrt(Var[mu + delta(h) | D])          small at data, larger at h = 0
  (4) sd of mu after one more probe at h                          smaller for finer h
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

SD, P, ELL, S, V0 = 0.8, 1.5, 1.0, 0.02, 100.0
HD = np.array([1.0, 0.5, 0.25, 0.125])
rng = np.random.default_rng(0)
Y = 1.0 + 0.4 * HD ** 1.5 + S * rng.standard_normal(HD.size)


def c52(r):
    a = np.sqrt(5) * np.abs(r) / ELL
    return (1 + a + a * a / 3) * np.exp(-a)


def kd(a, b):
    A, B = np.meshgrid(a, b, indexing="ij")
    return SD ** 2 * (A * B) ** P * c52(A - B)


def posterior(hd, y, hg):
    """Posterior mean/var of mu, delta(hg), f(hg) given data at hd."""
    K = V0 + kd(hd, hd) + S ** 2 * np.eye(hd.size)
    Ki = np.linalg.inv(K)
    kmu = np.full(hd.size, V0)
    kdel = kd(hg, hd)                       # cov(delta(hg), y)
    kf = V0 + kdel                          # cov(f(hg), y)
    m_mu, v_mu = kmu @ Ki @ y, V0 - kmu @ Ki @ kmu
    m_del = kdel @ Ki @ y
    v_del = np.diag(kd(hg, hg)) - np.einsum("ij,jk,ik->i", kdel, Ki, kdel)
    m_f = kf @ Ki @ y
    v_f = V0 + np.diag(kd(hg, hg)) - np.einsum("ij,jk,ik->i", kf, Ki, kf)
    return m_mu, v_mu, m_del, v_del, m_f, v_f


def main():
    hg = np.linspace(0.0, 1.05, 211)
    m_mu, v_mu, m_del, v_del, m_f, v_f = posterior(HD, Y, hg)
    rms_err = np.sqrt(m_del ** 2 + np.maximum(v_del, 0))
    sd_mu_next = np.array([np.sqrt(posterior(np.r_[HD, h], np.r_[Y, 0.0], hg[:1])[1])
                           for h in hg[1:]])
    prior_sd = SD * hg ** P

    plt.rcParams.update(fs.rcparams())
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    a1.fill_between(hg, m_f - 1.96 * np.sqrt(v_f), m_f + 1.96 * np.sqrt(v_f), color="#1f77b4",
                    alpha=0.2, lw=0)
    a1.plot(hg, m_f, color="#1f77b4", lw=1.2, label="posterior of f(h), 95%")
    a1.plot(hg, 1 + 0.4 * hg ** 1.5, "k--", lw=0.8, label="truth")
    a1.plot(HD, Y, "ko", ms=4, label="probes (4 levels)")
    a1.set_xlabel(r"cell size $\bar h = h/h_c$"); a1.set_ylabel("QoI")
    a1.grid(**fs.GRID_KW)
    a1.legend(loc="upper left", frameon=False, fontsize=8)

    a2.plot(hg, prior_sd, color="0.5", lw=1.2, label="(1) prior sd of the error of level h")
    a2.plot(hg, rms_err, color="#2ca02c", lw=1.6, label="(2) posterior RMS error of level h")
    a2.plot(hg, np.sqrt(v_f), color="#1f77b4", lw=1.6, label="(3) posterior sd of f(h)")
    a2.plot(hg[1:], sd_mu_next, color="#d62728", lw=1.6,
            label=r"(4) sd of $f(0)$ after one more probe at h")
    for h in HD:
        a2.axvline(h, color="0.85", lw=0.6, zorder=0)
    a2.set_yscale("log")
    a2.set_xlabel(r"cell size $\bar h = h/h_c$"); a2.set_ylabel("standard deviation")
    a2.grid(**fs.GRID_KW)
    a2.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=8)
    fig.tight_layout()
    out = ROOT / "experiments/multifidelity/sigma_concepts.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print("saved", out)
    print(f"sd of f(0) now {np.sqrt(v_mu):.4f}; sd of f(h=1/8) {np.sqrt(v_f[np.argmin(abs(hg - 0.125))]):.4f}")
    for h in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03):
        k = np.argmin(abs(hg[1:] - h))
        print(f"  one more probe at h={h:.4f}: sd of f(0) -> {sd_mu_next[k]:.4f}")


if __name__ == "__main__":
    main()
