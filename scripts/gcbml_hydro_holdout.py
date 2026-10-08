"""Level hold-out on the real tau data: fit L{lmin}-L9, predict the observed L10 values (9 rpm), per model.

An out-of-sample test with real data: the finest level is never used in the fit. Models: error shape "power"
(current) and "saturating" with sat_lo_factor 1 (satA3), on L6-L9 and L7-L9. The prediction of an L10 run is
the noise-free level f(x, h10) plus the run noise (m_s, per posterior draw, no site-specific part).
Caveat: the L10 values come from chained runs (window 3 cycles after an rpm change), except 32.5 rpm.
Metrics per model: median |relative error| of the median, pooled 95% coverage, mean relative 95% width.

Run: PYTHONPATH=<gcbml worktree with err_shape>/src uv run --project <gcbml> python scripts/gcbml_hydro_holdout.py
"""
import json
import math
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

import gcbml  # noqa: F401
from gcbml import inference, model, predict, transforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.gcbml_hydro_fit_v3 import RPM_HI, RPM_LO, SCALES, padded  # noqa: E402

ALL = json.load(open(ROOT / "experiments/multifidelity/hydro_dataset_v2.json"))
KEYS = ("tau_mean_t", "tau_95_t")
SPECS = [("power", {}), ("satA3", {"shape": "saturating", "sat_lo_factor": 1.0})]


def holdout(key_name, lmin, extra, seed=0):
    cfg = model.ModelConfig(h_kernel="twy2", mean_basis="constant", beta_prior=None, increasing=True, **extra)
    rows = [r for r in ALL if lmin <= r["level"] <= 9]
    test = sorted((r for r in ALL if r["level"] == 10), key=lambda r: r["rpm"])
    data = padded(rows, key_name)
    tf = transforms.get("log")
    z = jnp.asarray(tf.forward(jnp.asarray(data.y)))
    post = inference.fit(jax.random.key(seed), data, z, z, cfg, SCALES, n_controls=1, n_warmup=1000, n_samples=1000)
    params, zs = inference.flatten(post)
    Xs = np.array([[(r["rpm"] - RPM_LO) / (RPM_HI - RPM_LO)] for r in test])
    Hs = np.full((len(test), 1), 2.0 ** -(10 - lmin))
    Ps = np.ones((len(test), 1))

    def one(p, zz):
        Pj = jnp.broadcast_to(p.P[0], Ps.shape)
        return model.predict_level(p, data, zz, cfg, Xs, Hs, Pj)

    mean, var = jax.vmap(one)(params, zs)
    msn = jnp.exp(jnp.asarray(post.theta["m_s"]).reshape(-1))  # run-noise variance per draw (pooled chains)
    var = var + msn[:, None]
    draws, w = predict.sample_gaussian_draws(jax.random.key(seed + 1), mean, var, 4)
    summ = predict.summarize(draws, w, tf)
    obs = np.array([r[key_name] for r in test])
    m, lo, hi = (np.asarray(summ[k]) for k in ("m", "q025", "q975"))
    rel = np.abs(m - obs) / obs
    return dict(rel_err=rel.tolist(), covered=((obs >= lo) & (obs <= hi)).tolist(), width=((hi - lo) / m).tolist())


def main():
    out = {}
    for key in KEYS:
        for name, extra in SPECS:
            for lmin in (6, 7):
                r = holdout(key, lmin, extra)
                out[f"{key}|{name}|L{lmin}-9"] = r
                print(f"{key:11s} {name:6s} L{lmin}-L9 -> L10: median |err| {np.median(r['rel_err']):.3f}, "
                      f"coverage {np.mean(r['covered']):.2f} ({sum(r['covered'])}/{len(r['covered'])}), "
                      f"mean 95% width {np.mean(r['width']):.2f}", flush=True)
    json.dump(out, open(ROOT / "experiments/multifidelity/gcbml_hydro_holdout.json", "w"), indent=1)


if __name__ == "__main__":
    main()
