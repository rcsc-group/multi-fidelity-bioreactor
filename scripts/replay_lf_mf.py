"""PROPER mfbml with the repeated flow period as the low fidelity (pre-registered, diary 2026-10-21).
LF = looped-period replay QoI at 40 L8 conditions (experiments/multifidelity/replay_lf_l8.csv, scripts/collect_replay_lf.py); HF = the cold coupled L8 run at the same conditions (hf_ columns).
Learning curve: n = 4, 6, 8, 12, 16 HF runs, 20 random draws, scored on the other conditions (runs within 0.7 rpm at the same theta removed from the test set).
Methods: mfbml (FullBayesKRRLRGP, LF = KRR on the replay values), GPR with a linear mean prior, GPR with the replay value as prior mean plus one offset, replay x one constant.
RULE: REAL ADVANTAGE if mfbml or the replay-prior GPR has error <= 0.8 x the best HF-only method (GPR linear) at n <= 8 for >= 2 of the 3 QoIs. LOW FIDELITY INFORMATIVE if the
correlation of the residuals (after a quadratic trend in rpm and theta) of ln replay and ln HF is >= 0.5.
Usage: uv run python scripts/replay_lf_mf.py
"""
from __future__ import annotations

import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "scripts")]
from scripts.fig_threshold_mf import gpr_good_prior, gpr_lf_prior, quad  # noqa: E402
from scripts.test_h1_mfbml import BOUNDS, train_lf2d  # noqa: E402
from mfbml_local.krr_lr_gpr_fullbayes import FullBayesKRRLRGP  # noqa: E402

QS = ("kLa_1T_10", "kLa_1T_25", "dtmix_0.50")
SIZES = (4, 6, 8, 12, 16)
DRAWS = 20


def main():
    warnings.filterwarnings("ignore")
    rows = list(csv.DictReader((ROOT / "experiments/multifidelity/replay_lf_l8.csv").open()))
    out = {}
    for q in QS:
        R = [r for r in rows if r[q] not in ("", "nan", "None") and r[f"hf_{q}"] not in ("", "nan", "None") and float(r[q]) > 0 and float(r[f"hf_{q}"]) > 0]
        X = np.array([[float(r["rpm"]), float(r["theta"])] for r in R])
        lf, hf = np.log([float(r[q]) for r in R]), np.log([float(r[f"hf_{q}"]) for r in R])
        N = len(R)
        res = lambda y: y - quad(X) @ np.linalg.lstsq(quad(X), y, rcond=None)[0]   # noqa: E731
        corr_res = float(np.corrcoef(res(lf), res(hf))[0, 1])
        corr_raw = float(np.corrcoef(lf, hf)[0, 1])
        krr = train_lf2d(X, lf)
        rng = np.random.default_rng(11)
        curve = {}
        for n in SIZES:
            acc = {"MF": [], "GPR": [], "GPR_LF": [], "ratio": []}
            for _ in range(DRAWS):
                tr = np.zeros(N, bool)
                tr[rng.choice(N, n, replace=False)] = True
                near = np.zeros(N, bool)
                for i in np.where(tr)[0]:
                    near |= (np.abs(X[:, 0] - X[i, 0]) <= 0.7) & (X[:, 1] == X[i, 1])
                te = ~near
                if te.sum() < 5:
                    continue
                m = FullBayesKRRLRGP(BOUNDS, lf_model=krr, basis="linear").train(X[tr], hf[tr])
                acc["MF"].append(np.sqrt(np.mean((np.ravel(m.predict_mean(X[te])) - hf[te]) ** 2)))
                acc["GPR"].append(np.sqrt(np.mean((gpr_good_prior(X[tr], hf[tr], X[te]) - hf[te]) ** 2)))
                acc["GPR_LF"].append(np.sqrt(np.mean((gpr_lf_prior(X[tr], hf[tr], lf[tr], X[te], lf[te]) - hf[te]) ** 2)))
                acc["ratio"].append(np.sqrt(np.mean((lf[te] + np.mean(hf[tr] - lf[tr]) - hf[te]) ** 2)))
            curve[n] = {k: float(np.mean(v)) for k, v in acc.items()}
        out[q] = dict(n=N, corr_residual=corr_res, corr_raw=corr_raw, curve=curve, mean_ln_offset=float(np.mean(hf - lf)))
        print(f"{q}: {N} conditions; corr(ln replay, ln HF) {corr_raw:+.2f}, of the residuals after the (rpm, theta) trend {corr_res:+.2f}; mean ln(HF/replay) {np.mean(hf - lf):+.3f}")
        for n, c in curve.items():
            print(f"   n={n:2d}: MF {c['MF']:.2f}  GPR(lin) {c['GPR']:.2f}  GPR(replay prior) {c['GPR_LF']:.2f}  replay x const {c['ratio']:.2f}")
    adv = sum(1 for q in QS if all(min(out[q]["curve"][n]["MF"], out[q]["curve"][n]["GPR_LF"]) <= 0.8 * out[q]["curve"][n]["GPR"] for n in (4, 6, 8)))
    inform = all(out[q]["corr_residual"] >= 0.5 for q in QS)
    out["verdict"] = dict(advantage=("REAL ADVANTAGE" if adv >= 2 else "NO ADVANTAGE over the best HF-only method"), n_qois_meeting_rule=adv, low_fidelity_informative=bool(inform))
    print(out["verdict"])
    (ROOT / "experiments/multifidelity/replay_lf_mf.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
