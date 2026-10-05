"""gcbml (Yi-h) fit of tau_mean_max and tau_95_qss on L6-L10, input = rpm (theta = 7 deg), target h -> 0.

Data: experiments/multifidelity/hydro_dataset.json (scripts/hydro_dataset.py). Transform log (positive QoI with
relative grid errors). One resolution component, hbar = 2^-(L - 6). TWY2 kernel, flat beta, constant mean.
Prior scales (log units, estimates): S_mu = 0.8 (tau varies ~5x over rpm), S_c = 0.5, S_delta = 0.5,
S_noise = 0.02 (replicates at L6-L8 agree within ~1%).
Outputs experiments/multifidelity/gcbml_hydro_<qoi>.json and gcbml_hydro_fit.png.

Run: uv run --project /oscar/data/dharri15/eaguerov/Github/gcbml python scripts/gcbml_hydro_fit.py
"""
import json
import sys
import time
from pathlib import Path

import jax
import jax.numpy as jnp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import gcbml  # noqa: F401
from gcbml import inference, model, predict, transforms
from gcbml._config import bucket
from gcbml.data import PaddedData
from gcbml.priors import PriorScales

ROOT = Path(__file__).resolve().parents[1]
RPM_LO, RPM_HI = 15.0, 37.5
QOIS = [("tau_mean_max", r"$\tau_{mean}$ (Pa)"), ("tau_95_qss", r"$\tau_{95}$ (Pa)")]
COLS = {6: "#c6dbef", 7: "#9ecae1", 8: "#4292c6", 9: "#2171b5", 10: "#08306b"}
SCALES = PriorScales(S_mu=0.8, S_c=0.5, S_delta=0.5, S_noise=0.02)
CFG = model.ModelConfig(h_kernel="twy2", mean_basis="constant", beta_prior=None, increasing=True)


def padded(rows, key):
    n = len(rows)
    npad = bucket(n)
    X = np.zeros((npad, 1))
    H = np.ones((npad, 1))
    y = np.ones(npad)
    X[:n, 0] = [(r["rpm"] - RPM_LO) / (RPM_HI - RPM_LO) for r in rows]
    H[:n, 0] = [2.0 ** -(r["level"] - min(rr["level"] for rr in rows)) for r in rows]
    y[:n] = [r[key] for r in rows]
    X[n:], H[n:], y[n:] = X[0], H[0], y[0]
    mask = np.arange(npad) < n
    return PaddedData(X=X, H=H, y=y, censored=np.zeros(npad, bool), run=np.where(mask, np.arange(npad), -1),
                      mask=mask)


def fit_one(key_name, n_warmup, n_samples, seed=0, min_level=6):
    rows = [r for r in json.load(open(ROOT / "experiments/multifidelity/hydro_dataset.json")) if r["level"] >= min_level]
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
    draws, w = predict.sample_gaussian_draws(jax.random.key(seed + 1), mean, var, 4)
    summ = predict.summarize(draws, w, tf)
    lp = np.asarray(post.theta["log_p0"]).reshape(-1)
    diag = {k: [float(x) for x in v] for k, v in post.diagnostics.items() if k.startswith(("log_p0", "log_sigma_mu"))}
    out = dict(qoi=key_name, rpm=(RPM_LO + xs[:, 0] * (RPM_HI - RPM_LO)).tolist(),
               m=np.asarray(summ["m"]).tolist(), q025=np.asarray(summ["q025"]).tolist(),
               q975=np.asarray(summ["q975"]).tolist(), q16=np.asarray(summ["q16"]).tolist(),
               q84=np.asarray(summ["q84"]).tolist(), sigma=np.asarray(summ["sigma"]).tolist(),
               p0_q=np.quantile(np.exp(lp), [0.05, 0.5, 0.95]).tolist(), fit_seconds=t_fit, diagnostics=diag,
               n_warmup=n_warmup, n_samples=n_samples)
    out["min_level"] = min_level
    json.dump(out, open(ROOT / f"experiments/multifidelity/gcbml_hydro_{key_name}_L{min_level}.json", "w"), indent=1)
    rel = np.asarray(out["sigma"]) / np.asarray(out["m"])
    print(f"{key_name}: fit {t_fit:.0f} s; p0 5/50/95% = {np.round(out['p0_q'], 2)}; "
          f"sigma_epi/m median {np.median(rel):.3f} max {rel.max():.3f}; diag {diag}")
    return rows, out


def main():
    nw, ns = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (500, 500)
    lmin = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, (key, lab) in zip(axes, QOIS):
        rows, out = fit_one(key, nw, ns, min_level=lmin)
        for lev in sorted(COLS):
            pts = sorted((r["rpm"], r[key]) for r in rows if r["level"] == lev)
            if pts:
                a = np.array(pts)
                ax.plot(a[:, 0], a[:, 1], "o", color=COLS[lev], ms=4, label=f"L{lev}")
        ax.fill_between(out["rpm"], out["q025"], out["q975"], color="#d95f02", alpha=0.15, lw=0,
                        label=r"$h \to 0$, 95% band")
        ax.fill_between(out["rpm"], out["q16"], out["q84"], color="#d95f02", alpha=0.35, lw=0,
                        label=r"$h \to 0$, 68% band")
        ax.plot(out["rpm"], out["m"], color="#d95f02", lw=2, label=r"$h \to 0$, median")
        ax.set_xlabel("rocking speed (rpm)")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3)
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False)
    fig.tight_layout()
    fig.savefig(ROOT / f"experiments/multifidelity/gcbml_hydro_fit_L{lmin}.png", dpi=150, bbox_inches="tight")
    print("saved figure")


if __name__ == "__main__":
    main()
