"""gcbml (Yi-h) fit of tau_mean_t and tau_95_t (hydro_dataset_v2.json), target h -> 0, with two h-kernels.

Question (2026-10-05): adding the pre-asymptotic levels L6-L7 widened the h -> 0 band (L6-L10 15% vs L8-L10 7%).
Hypothesis: TWY2 is a smooth (non-Markov) kernel in h, so every coarse level informs the extrapolation. The
lifted Brownian kernel with gamma = 0.5 (Boutelet & Sung 2503.23158 Sec 2.2) is k = min(h, h')^{2p}: Brownian
motion in t = h^{2p}, which has independent increments. Given the hyperparameters, a coarse level then adds an
independent increment and tells nothing about h = 0 beyond the finer level at the same x (the screening that
the recursive multi-fidelity form of Yi et al. has). It still couples through the shared hyperparameters.
Falsifiable: with LB(0.5), the L6-L10 band must move toward the L8-L10 band; if it does not, the coupling is in
the hyperparameters (p, sigma) and needs another remedy.

Model: log transform; input rpm scaled to [0, 1] over 15-37.5; hbar = 2^-(L - Lmin); constant mean, flat beta.
Prior scales (log units, estimates): S_mu 0.8, S_c 0.5, S_delta 0.5, S_noise 0.02.
Order prior (assumption A3 made explicit): log p ~ N(log 1.2, 0.45^2), 95% in [0.5, 2.9] -- the asymptotic range
of a scheme of formal order 1-2; Eca & Hoekstra treat p < 0.5 as anomalous.

Bands at h = 0: epistemic = quantiles of f(x, 0); total = f(x, 0) + e, e ~ N(0, s^2(x, 0)). There are no
replicate runs, so s is prior-driven. log s^2(x, 0) = m_s + zeta with zeta drawn from its marginal
N(0, sigma_z^2), without conditioning on the site values: this over-states the aleatoric part.

Run: KERNEL=twy2|lb uv run --project /oscar/data/dharri15/eaguerov/Github/gcbml \
         python scripts/gcbml_hydro_fit_v3.py N_WARMUP N_SAMPLES LMIN LMAX
Writes experiments/multifidelity/gcbml_hydro_v3_<kernel>_L<min>_<max>.{json,png}.
"""
import json
import math
import os
import sys
import time
from pathlib import Path

import jax
import jax.numpy as jnp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import gcbml  # noqa: E402,F401
from gcbml import inference, model, predict, transforms  # noqa: E402
from gcbml._config import bucket  # noqa: E402
from gcbml.data import PaddedData  # noqa: E402
from gcbml.priors import PriorScales  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import figstyle as fs  # noqa: E402

RPM_LO, RPM_HI = 15.0, 37.5
# (key, axis label, statistic for the fill channel): tau_mean is a time mean of the spatial mean -> hollow;
# tau_95 is a tail percentile (an extremum-type statistic) -> filled
QOIS = [("tau_mean_t", r"$\tau_{mean}$ (Pa)", "mean"), ("tau_95_t", r"$\tau_{95}$ (Pa)", "max")]
KERNEL = os.environ.get("KERNEL", "twy2")
# sensitivity knobs (G7-type prior-sensitivity runs): S_C, S_DELTA, S_NOISE, P_SD; TAG suffixes the output name
SCALES = PriorScales(S_mu=0.8, S_c=float(os.environ.get("S_C", 0.5)), S_delta=float(os.environ.get("S_DELTA", 0.5)),
                     S_noise=float(os.environ.get("S_NOISE", 0.02)), log_p_mean=math.log(float(os.environ.get("P_MED", 1.2))),
                     log_p_sd=float(os.environ.get("P_SD", 0.45)))
TAG = os.environ.get("TAG", "")
SHAPE = os.environ.get("SHAPE", "power")  # "saturating" needs the gcbml branch with err_shape (PYTHONPATH)
_extra = {} if SHAPE == "power" else {"shape": SHAPE}
CFG = model.ModelConfig(h_kernel=KERNEL, mean_basis="constant", beta_prior=None, increasing=True,
                        gamma_fixed=0.5 if KERNEL == "lb" else None, **_extra)


def padded(rows, key):
    n = len(rows)
    npad = bucket(n)
    X, H, y = np.zeros((npad, 1)), np.ones((npad, 1)), np.ones(npad)
    lmin = min(r["level"] for r in rows)
    X[:n, 0] = [(r["rpm"] - RPM_LO) / (RPM_HI - RPM_LO) for r in rows]
    H[:n, 0] = [2.0 ** -(r["level"] - lmin) for r in rows]
    y[:n] = [r[key] for r in rows]
    X[n:], H[n:], y[n:] = X[0], H[0], y[0]
    mask = np.arange(npad) < n
    return PaddedData(X=X, H=H, y=y, censored=np.zeros(npad, bool), run=np.where(mask, np.arange(npad), -1),
                      mask=mask)


def fit_one(key_name, n_warmup, n_samples, lmin, lmax, seed=0):
    rows = [r for r in json.load(open(ROOT / "experiments/multifidelity/hydro_dataset_v2.json"))
            if lmin <= r["level"] <= lmax]
    data = padded(rows, key_name)
    tf = transforms.get("log")
    z = jnp.asarray(tf.forward(jnp.asarray(data.y)))
    t0 = time.time()
    post = inference.fit(jax.random.key(seed), data, z, z, CFG, SCALES, n_controls=1, n_warmup=n_warmup,
                         n_samples=n_samples)
    t_fit = time.time() - t0
    params, zs = inference.flatten(post)
    xs = np.linspace(0, 1, 91)[:, None]
    mean, var = jax.vmap(lambda p, zz: model.predict_mu(p, data, zz, CFG, jnp.asarray(xs)))(params, zs)
    k1, k2, k3 = jax.random.split(jax.random.key(seed + 1), 3)
    draws, w = predict.sample_gaussian_draws(k1, mean, var, 4)
    epi = predict.summarize(draws, w, tf)
    # aleatoric sd at h = 0 per posterior draw (marginal zeta, see the module docstring)
    m_s = np.asarray(post.theta["m_s"]).reshape(-1)
    sig_z = np.exp(np.asarray(post.theta["log_sigma_z"]).reshape(-1))
    zeta = sig_z * np.asarray(jax.random.normal(k2, m_s.shape))
    s0 = np.exp(0.5 * (m_s + zeta))  # (S,)
    s0_rep = np.repeat(s0, 4)[:, None] * np.ones((1, xs.shape[0]))
    tot_draws = draws + s0_rep * np.asarray(jax.random.normal(k3, draws.shape))
    tot = predict.summarize(tot_draws, w, tf)
    lp = np.asarray(post.theta["log_p0"]).reshape(-1)
    diag = {k: [float(x) for x in v] for k, v in post.diagnostics.items()}
    worst = max(v[0] for v in diag.values())
    out = dict(qoi=key_name, kernel=KERNEL, levels=[lmin, lmax],
               rpm=(RPM_LO + xs[:, 0] * (RPM_HI - RPM_LO)).tolist(),
               m=np.asarray(epi["m"]).tolist(),
               epi_q025=np.asarray(epi["q025"]).tolist(), epi_q975=np.asarray(epi["q975"]).tolist(),
               tot_q025=np.asarray(tot["q025"]).tolist(), tot_q975=np.asarray(tot["q975"]).tolist(),
               sigma_epi=np.asarray(epi["sigma"]).tolist(),
               p0_q=np.quantile(np.exp(lp), [0.05, 0.5, 0.95]).tolist(),
               s0_q=np.quantile(s0, [0.05, 0.5, 0.95]).tolist(),
               fit_seconds=t_fit, max_rhat=worst, n_extensions=getattr(post.diagnostics, "n_extensions", 0),
               n_warmup=n_warmup, n_samples=n_samples)
    if SHAPE == "saturating":  # h_s in hbar units (coarsest fitted level = 1); also as an equivalent level
        hs = np.exp(np.asarray(post.theta["aux"]).reshape(-1))
        out["hs_q"] = np.quantile(hs, [0.05, 0.5, 0.95]).tolist()
        out["hs_level_q"] = (lmin - np.log2(np.quantile(hs, [0.95, 0.5, 0.05]))).tolist()
        print(f"  h_s 5/50/95% {np.round(out['hs_q'], 3)} = level {np.round(out['hs_level_q'], 1)}", flush=True)
    rel = np.asarray(out["sigma_epi"]) / np.asarray(out["m"])
    out["rel_sigma_epi_median"] = float(np.median(rel))
    print(f"{KERNEL} L{lmin}-{lmax} {key_name}: fit {t_fit:.0f} s; p0 5/50/95% {np.round(out['p0_q'], 2)}; "
          f"s0 5/50/95% {np.round(out['s0_q'], 3)}; sigma_epi/m median {np.median(rel):.3f}; "
          f"max rhat {worst:.3f}; extensions {out['n_extensions']}", flush=True)
    return rows, out


def main():
    nw, ns, lmin, lmax = (int(a) for a in sys.argv[1:5])
    plt.rcParams.update(fs.rcparams())
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    res = []
    for ax, (key, lab, stat) in zip(axes, QOIS):
        rows, out = fit_one(key, nw, ns, lmin, lmax)
        res.append(out)
        ax.fill_between(out["rpm"], out["tot_q025"], out["tot_q975"], color=fs.H0_COLOUR, alpha=0.18, lw=0,
                        label=r"$h \to 0$: total 95%")
        ax.fill_between(out["rpm"], out["epi_q025"], out["epi_q975"], color=fs.H0_COLOUR, alpha=0.45, lw=0,
                        label=r"$h \to 0$: epistemic 95%")
        ax.plot(out["rpm"], out["m"], color=fs.H0_COLOUR, ls=":", lw=2, label=r"$h \to 0$: median")
        for lev in range(lmin, lmax + 1):
            pts = sorted((r["rpm"], r[key]) for r in rows if r["level"] == lev)
            if pts:
                a = np.array(pts)
                ax.plot(a[:, 0], a[:, 1], **fs.series_kw(fs.level_colour(lev), fs.MK["tau"], stat=stat,
                                                         ours=True), label=f"L{lev}")
        ax.set_xlabel("rocking speed (rpm)")
        ax.set_ylabel(lab)
        ax.grid(**fs.GRID_KW)
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    stem = ROOT / f"experiments/multifidelity/gcbml_hydro_v3_{KERNEL}{TAG}_L{lmin}_{lmax}"
    fig.savefig(f"{stem}.png", dpi=150, bbox_inches="tight")
    json.dump(res, open(f"{stem}.json", "w"), indent=1)
    print("saved", stem)


if __name__ == "__main__":
    main()
