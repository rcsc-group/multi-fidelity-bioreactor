"""Physics-informed LF for kLa (2026-11-10 22:50). L8 hydrodynamics track L9 (rho ~1), L8 kLa does not (rho 0.08,
L8/L9 ~3.7: under-resolved concentration boundary layer). Question: does an L8 HYDRO feature, or an eddy-cell form
kL ~ (eps nu)^(1/4) (Lamont & Scott 1970), rank L9 kLa_1T_10 better than L8 kLa itself?
Paired sets (gate ignored): pair12 (pr8_l8 vs bo1_l9), bo4 chain (bo4_l8 vs bo4_l9), theta = 7 line (truth snapshot medians).
Prints Spearman rho per feature per set; and LOO RMSE (ln) of y9 = a + b ln(feature) vs y9 = mean.
  uv run python scripts/bo2d/diag_hydro_kla.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
FEATS = ["kLa_1T_10", "vel_rms_qss", "edr_95", "ediss_mean_qss", "vor_mean", "vor_streaming", "tau_95_qss", "tau_98_max"]


def load(d):
    p = d / "results.json"
    return json.loads(p.read_text()) if p.exists() else None


def pairs(g8, g9):
    out = []
    for d9 in sorted((ROOT / "runs").glob(g9)):
        if "failed" in d9.name:
            continue
        key = d9.name.split("_l9_")[1]
        c = [d for d in (ROOT / "runs").glob(g8 + key) if "failed" not in d.name]
        if not c:
            continue
        r8, r9 = load(c[0]), load(d9)
        if r8 and r9 and np.isfinite(r9.get("kLa_1T_10", np.nan) or np.nan):
            out.append((r8, r9))
    return out


def loo(f, y):
    A = np.column_stack([np.ones_like(f), f]); e = []
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        c, *_ = np.linalg.lstsq(A[m], y[m], rcond=None); e.append(y[i] - A[i] @ c)
    e0 = [y[i] - y[np.arange(len(y)) != i].mean() for i in range(len(y))]
    return np.sqrt(np.mean(np.square(e))), np.sqrt(np.mean(np.square(e0)))


for name, g8, g9 in (("pair12", "pr8_l8_", "bo1_l9_*"), ("bo4", "bo4_l8_", "bo4_l9_*"), ("fig9_th7", "t6_l8_", "fig9_l9_*")):
    P = pairs(g8, g9)
    if len(P) < 5:
        print(name, "n", len(P), "too few"); continue
    y = np.log([p[1]["kLa_1T_10"] for p in P])
    print(f"\n{name}: n={len(P)} sd ln kLa9 {np.std(y):.2f}")
    feats = {k: np.array([p[0].get(k, np.nan) for p in P], float) for k in FEATS}
    nu_eps = feats["ediss_mean_qss"]
    feats["eddycell (eps)^1/4 * vel_rms^0"] = nu_eps ** 0.25
    feats["eddycell (edr_95)^1/4"] = feats["edr_95"] ** 0.25
    for k, f in feats.items():
        ok = np.isfinite(f) & (f > 0)
        if ok.sum() < 5:
            continue
        rho = spearmanr(f[ok], y[ok]).correlation
        r, r0 = loo(np.log(f[ok]), y[ok])
        print(f"  {k:34s} rho {rho:+.2f}  LOO {r:.3f} vs mean {r0:.3f} ({r / r0:.2f})")
