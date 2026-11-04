"""A1: chain learning curves with smaller, measured-only training sets.

Protocol identical to scripts/h3_baselines_learning_curve.py (chain_pred): L8 KRR (train_lf2d, log values) -> L9 GP
(FullBayesKRRLRGP, lf_model=L8 KRR, log values) -> L10 GP (lf_model=Offset(L9 model), log values). Held-out: the
same random subsets of the 13 L10 runs (rng seed 20261021, draws in the same order as h3, SIZES 3..12, DRAWS 15).
Variants:
  full     : 40 L8, 16 L9 (approved chain; must reproduce h3 chain values)
  warm9    : 40 L8, 6 L9 (warm x8_l9 only; no cold fig9_l9 runs)
  warm9_L8_20 : 20 L8 (the 6 sources of the warm L9 runs + 14 random others, seed 7), 6 L9
  warm9_L8_10 : 10 L8 (the 6 sources + 4 random others, seed 7), 6 L9
Cost (core-h) = sum of corehc_all per run id from cost_table.tsv (sacct, 2026-11-04 09:22), L10 own cost per run
(same convention as the approved chain curve). Warm L9 runs need their L8 source (checked in the L8 set).
Usage: uv run python tools/chain_variants_lc.py
Output: data/chain_variants_lc.json
"""
import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
warnings.filterwarnings("ignore")
import scripts.test_h1_mfbml as H1  # noqa: E402
import scripts.test_h2_mfbml as H2  # noqa: E402
from scripts.h3_band_check import Offset  # noqa: E402
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

QS = ("tau_95_t", "edr_max_t")
SIZES = (3, 4, 6, 8, 10, 12)
DRAWS = 15
TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")

COST = {}
with open("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        COST[r["run"]] = float(r["corehc_all"])


def c(name):
    return COST[name]


def rel_rmse(pred, obs):
    return float(np.sqrt(np.mean((np.asarray(pred) / np.asarray(obs) - 1) ** 2)))


lf_pts, _ = H1.load_points([(r, t) for r in H1.RPMS for t in H1.LF_THETAS], H1.lf_run)
l9_pts, _ = H1.load_points(H2.l9_points(), H2.l9_run)
hf_pts, _ = H1.load_points(H1.hf_points(), H1.hf_run)
hk = sorted(hf_pts, key=lambda k: (k[1], k[0]))
X = np.array(hk, float)
keys8_all = sorted(lf_pts, key=lambda k: (k[1], k[0]))
keys9_all = sorted(l9_pts, key=lambda k: (k[1], k[0]))
warm9 = [k for k in keys9_all if k[1] != 7.0]          # x8_l9 runs at theta 4 and 9
sources = [k for k in warm9]                           # source of x8_l9 at (r, t) is the L8 run at (r, t)
assert all(s in lf_pts for s in sources), "a warm L9 source is missing from the L8 set"
rng_sub = np.random.default_rng(7)
others = [k for k in keys8_all if k not in sources]


def pick(n_total):
    extra = rng_sub.choice(len(others), n_total - len(sources), replace=False)
    return sources + [others[i] for i in sorted(extra)]


VARIANTS = {
    "full": dict(keys8=keys8_all, keys9=keys9_all),
    "warm9": dict(keys8=keys8_all, keys9=warm9),
    "warm9_L8_20": dict(keys8=pick(20), keys9=warm9),
    "warm9_L8_10": dict(keys8=pick(10), keys9=warm9),
}


def costs_of(v):
    k8 = [H1.lf_run(r, t) for r, t in v["keys8"]]
    k9 = [H2.l9_run(r, t) for r, t in v["keys9"]]
    base = sum(c(n) for n in k8) + sum(c(n) for n in k9)
    inc = float(np.mean([c(H1.hf_run(r, t)) for r, t in hk]))
    return dict(n_L8=len(k8), n_L9=len(k9), base=round(base, 1), inc_per_L10=round(inc, 1))


rng = np.random.default_rng(20261021)
out = {}
for name, v in VARIANTS.items():
    cst = costs_of(v)
    out[name] = dict(costs=cst, keys8=[list(k) for k in v["keys8"]], keys9=[list(k) for k in v["keys9"]], err={})
    for q in QS:
        k8 = v["keys8"]
        lf = H1.train_lf2d(np.array(k8, float), np.log([lf_pts[k][q] for k in k8]))
        k9 = v["keys9"]
        m1 = FullBayesKRRLRGP(H1.BOUNDS, lf_model=lf, basis="linear").train(np.array(k9, float), np.log([l9_pts[k][q] for k in k9]))
        ad = Offset(m1)
        y10 = np.array([hf_pts[k][q] for k in hk])
        lc = {}
        for n in SIZES:
            errs = []
            for _ in range(DRAWS):
                idx = rng.choice(len(hk), n, replace=False)
                tr = np.zeros(len(hk), bool)
                tr[idx] = True
                te = ~tr
                m2 = FullBayesKRRLRGP(H1.BOUNDS, lf_model=ad, basis="linear").train(X[tr], np.log(y10[tr]))
                pred = np.exp(np.ravel(m2.predict_mean(X[te])))
                errs.append(rel_rmse(pred, y10[te]))
            lc[n] = float(np.mean(errs))
        out[name]["err"][q] = lc
        print(f"{name:12s} {q:10s} n_L8={cst['n_L8']:2d} n_L9={cst['n_L9']:2d} base={cst['base']:8.1f} "
              + " ".join(f"n{n}:{100 * lc[n]:.1f}%" for n in SIZES), flush=True)
    print(f"  costs: {cst}", flush=True)
    # reset rng per variant so each variant draws the same held-out subsets as the first (rng reseeded below)
    rng = np.random.default_rng(20261021)
(TALK / "data" / "chain_variants_lc.json").write_text(json.dumps(out, indent=1))
print("wrote", TALK / "data" / "chain_variants_lc.json")
