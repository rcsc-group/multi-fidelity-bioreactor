"""Gate G4 (posterior-predictive check of the coarsest level, spec Step 3) on the hydro data, tau_mean_t / tau_95_t.

Starts from L{LMIN}-L{LMAX}; while G4 fails and at least 3 levels remain, removes the coarsest level and refits.
Uses gcbml.gates.g4_coarsest_level (campaign branch w6d; run with PYTHONPATH=<worktree>/src ahead of the installed
package). Same model and priors as gcbml_hydro_fit_v3.py (TWY2, scheme order prior).

Run: PYTHONPATH=/oscar/data/dharri15/eaguerov/Github/gcbml-wt/w6b/src uv run --project /oscar/data/dharri15/eaguerov/Github/gcbml \
         python scripts/gcbml_hydro_g4.py [LMIN] [LMAX]
"""
import json
import sys
from pathlib import Path

import jax
import jax.numpy as jnp
import numpy as np

import gcbml  # noqa: F401
from gcbml import gates, inference, transforms

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.gcbml_hydro_fit_v3 import CFG, SCALES, padded  # noqa: E402

LMIN = int(sys.argv[1]) if len(sys.argv) > 1 else 6
LMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 10
ALL = json.load(open(ROOT / "experiments/multifidelity/hydro_dataset_v2.json"))


def run(key_name, lmin):
    rows = [r for r in ALL if lmin <= r["level"] <= LMAX]
    data = padded(rows, key_name)
    tf = transforms.get("log")
    z = jnp.asarray(tf.forward(jnp.asarray(data.y)))
    post = inference.fit(jax.random.key(0), data, z, z, CFG, SCALES, n_controls=1, n_warmup=1000, n_samples=1000)
    params, zs = inference.flatten(post)
    block = np.array([i for i, r in enumerate(rows) if r["level"] == lmin])
    res = gates.g4_coarsest_level(params, data, zs, CFG, block)
    return res, np.quantile(np.exp(np.asarray(post.theta["log_p0"]).reshape(-1)), [0.05, 0.5, 0.95])


def main():
    for key in ("tau_mean_t", "tau_95_t"):
        lmin = LMIN
        while True:
            res, pq = run(key, lmin)
            st = res.stats
            print(f"{key} L{lmin}-{LMAX}: G4 {res.status}; p-value {st.get('p')}; stat {st.get('statistic')} "
                  f"dof {st.get('dof')}; order p 5/50/95% {np.round(pq, 2)}", flush=True)
            if res.status == "pass" or LMAX - lmin < 2:
                break
            lmin += 1
        print(f"{key}: levels kept L{lmin}-{LMAX}\n", flush=True)


if __name__ == "__main__":
    main()
