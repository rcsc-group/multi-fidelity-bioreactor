"""Dense L7 screen (200 maximin-LHS runs, scr_l7_*) as a landscape map (2026-11-11 07:45, user: "screen very finely,
then a few L9 to map out the landscape"). Per QoI (tau_95, dtmix_0.25 as rate, kLa_1T_10): fit a GP (FullBayes,
ordinary) to the L7 screen, predict at every L9 site (off-line and theta = 7, gate ignored), report Spearman rho
(L7 map, L9) and LOO-free calibrated error: L9 ~ a + b * L7map (LOO over L9 sites) vs L9 mean. Writes
state/l7_screen_map.json (grid predictions) for figures.
  uv run python scripts/bo2d/diag_l7_screen_map.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.bo2d import loop as LP  # noqa: E402
from scripts.bo2d import bo  # noqa: E402
from scripts.mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

C = LP.C
C.warm_run_gate = lambda d: ("PASS", {})
_row = C.campaign_row


def _safe_row(*a, **k):   # a few screen runs are short (empty QoI window): skip them
    try:
        return _row(*a, **k)
    except (ValueError, IndexError, KeyError, FileNotFoundError):
        return None


C.campaign_row = _safe_row


def clean(d):
    y = d.y.astype(float)
    return d[np.isfinite(y) & (y > 0)]


out = {}
G = bo.grid()
for q in ("tau_95", "dtmix_0.25", "kLa_1T_10"):
    d7 = clean(LP.as_objective(LP.rows(q, ["scr_l7_"], {7}), q))
    d9 = clean(LP.as_objective(LP.rows(q, ["bo1_", "bo4_", "fig9_l9_rpm"], {9}), q))
    X7, y7 = d7[["rpm", "theta"]].to_numpy(float), np.log(d7.y.to_numpy(float))
    m = FullBayesKRRLRGP(C.DOMAIN, basis="ordinary", grid=C.hf_grid(C.noise_floor(q), float(np.std(y7)))).train(X7, y7)
    X9, y9 = d9[["rpm", "theta"]].to_numpy(float), np.log(d9.y.to_numpy(float))
    p9 = np.ravel(m.predict_mean(X9))
    for sel, mk in (("all", np.ones(len(y9), bool)), ("off-line", ~np.isclose(X9[:, 1], 7)), ("theta=7", np.isclose(X9[:, 1], 7))):
        f, y = p9[mk], y9[mk]
        A = np.column_stack([np.ones_like(f), f]); e, e0 = [], []
        for i in range(len(y)):
            k = np.arange(len(y)) != i
            c, *_ = np.linalg.lstsq(A[k], y[k], rcond=None); e.append(y[i] - A[i] @ c); e0.append(y[i] - y[k].mean())
        r, r0 = np.sqrt(np.mean(np.square(e))), np.sqrt(np.mean(np.square(e0)))
        print(f"{q:11s} {sel:8s} n7={len(y7)} n9={len(y):2d} rho {spearmanr(f, y).correlation:+.2f}  "
              f"LOO calibrated {r:.3f} vs mean {r0:.3f} ({r / r0:.2f})", flush=True)
    out[q] = {"grid": G.tolist(), "mean": np.ravel(m.predict_mean(G)).tolist(), "X7": X7.tolist(), "X9": X9.tolist(),
              "y9": y9.tolist()}
(ROOT / "scripts/bo2d/state/l7_screen_map.json").write_text(json.dumps(out))
print("wrote state/l7_screen_map.json")
