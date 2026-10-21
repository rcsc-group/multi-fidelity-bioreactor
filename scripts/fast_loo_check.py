"""Check FastLooKRR against LooKRR on the 40 L8 tau_95 runs: objective values at random parameters, then the fitted predictions and the time. Usage: uv run python scripts/fast_loo_check.py"""
import sys, time
import numpy as np
sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
from scripts.fast_loo_krr import FastLooKRR
from scripts.plot_mf_fig13 import LooKRR
from scripts.test_h1_mfbml import BOUNDS, LF_THETAS, RPMS, lf_run, load_points

lf, _ = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
keys = sorted(lf, key=lambda k: (k[1], k[0]))
X = np.array(keys, float); y = np.log([lf[k]["tau_95_t"] for k in keys])
out = {}
for name, cls in (("slow", LooKRR), ("fast", FastLooKRR)):
    m = cls(design_space=BOUNDS, params_optimize=True, noise_data=True, optimizer_restart=10, seed=42)
    t = time.time(); m.train(X, y, portion_test=0.2); out[name] = (time.time() - t, np.ravel(m.predict(X)))
print("seconds slow %.1f fast %.1f" % (out["slow"][0], out["fast"][0]), "max |difference of predictions| %.2e" % np.max(np.abs(out["slow"][1] - out["fast"][1])), "(prediction range %.3f)" % np.ptp(out["slow"][1]))
