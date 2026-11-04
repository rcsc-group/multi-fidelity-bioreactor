"""Learning curve for the two-level model H1 (direct L8 -> L10, raw scale), same protocol as h3_baselines_learning_curve.py.

Copied from the repo's test_h1_mfbml.py (model) and h3_baselines_learning_curve.py (protocol): same SIZES (3,4,6,8,10,12),
same DRAWS (15), same rng seed (20261021) and draw order, so the held-out subsets are the SAME as for the chain curve.
Model H1 (raw scale, as in test_h1_mfbml.loo): LF = KRR on all 40 L8 points; HF = FullBayesKRRLRGP on the n training L10 runs
with lf_model=LF (basis linear). Scored on the held-out L10 runs by rel RMS (same rel_rmse).
Usage: uv run python tools/h1_lc_copy.py tau_95_t,edr_max_t
Output: /oscar/data/dharri15/eaguerov/talks/dfd2026/data/lc_h1_<qs>.json
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
QS = tuple(sys.argv[1].split(",")) if len(sys.argv) > 1 else ("tau_95_t", "edr_max_t")
SIZES = (3, 4, 6, 8, 10, 12)
DRAWS = 15

from scripts.test_h1_mfbml import BOUNDS, LF_THETAS, RPMS, hf_points, hf_run, lf_run, load_points, train_lf2d  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402


def rel_rmse(pred, obs):
    return float(np.sqrt(np.mean((np.asarray(pred) / np.asarray(obs) - 1) ** 2)))


def main():
    warnings.filterwarnings("ignore")
    lf_pts, _ = load_points([(r, t) for r in RPMS for t in LF_THETAS], lf_run)
    hf_pts, _ = load_points(hf_points(), hf_run)
    hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
    X = np.array(hk, float)
    rng = np.random.default_rng(20261021)
    out = {}
    for q in QS:
        k8 = sorted(lf_pts, key=lambda k: (k[1], k[0]))
        lf = train_lf2d(np.array(k8), [lf_pts[k][q] for k in k8])
        y10 = np.array([hf_pts[k][q] for k in hk])
        lc = {}
        for n in SIZES:
            errs = []
            for _ in range(DRAWS):
                idx = rng.choice(len(hk), n, replace=False)
                tr = np.zeros(len(hk), bool)
                tr[idx] = True
                te = ~tr
                m = FullBayesKRRLRGP(BOUNDS, lf_model=lf, basis="linear").train(X[tr], y10[tr])
                pred = np.ravel(m.predict_mean(X[te]))
                errs.append(rel_rmse(pred, y10[te]))
            lc[n] = float(np.mean(errs))
            print(f"{q} n={n}: H1 mean rel RMS {100 * lc[n]:.1f} %  ({DRAWS} draws)", flush=True)
        out[q] = dict(learning_curve={str(n): {"h1": v} for n, v in lc.items()})
    dst = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026/data") / ("lc_h1_" + "_".join(QS) + ".json")
    dst.write_text(json.dumps(out, indent=1))
    print("wrote", dst)


if __name__ == "__main__":
    main()
