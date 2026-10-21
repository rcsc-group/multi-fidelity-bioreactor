"""PI idea 1 (pre-registered, diary 2026-10-21): predict kLa at a later C* threshold from kLa at C* = 0.10.
Target ln kLa_1T_50 and ln kLa_1T_25 over the 89 L8 grid runs (theta 4, 7, 9). Leave-one-out, runs within 0.25 rpm at the same theta removed from training.
M0 constant; M1 GP on (rpm, theta); M2 linear in ln kLa_1T_10; M3 GP on (rpm, theta, ln kLa_1T_10).
USEFUL if the best kLa_1T_10 model has RMS ln error <= 0.8 x M1; else NO INFORMATION BEYOND (rpm, theta).
Usage: uv run python scripts/kla_threshold_extrapolation.py
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

ROOT = Path(__file__).resolve().parents[1]


def load():
    rows = []
    for pat in ("t6_l8_*", "t5_l8_*", "e8_l8_*"):
        for d in sorted((ROOT / "runs").glob(pat)):
            f = d / "results.json"
            if not f.exists():
                continue
            r = json.loads(f.read_text())
            v = [r.get(q) for q in ("kLa_1T_10", "kLa_1T_25", "kLa_1T_50")]
            if all(x and x > 0 for x in v):
                p = json.loads((d / "params.json").read_text())
                import math
                rows.append((p["omega_b"] * 60 / (2 * math.pi), p["theta_max"][0], *v))
    return np.array(rows)


def gp(X, y):
    k = ConstantKernel(1.0, (1e-2, 1e2)) * RBF(np.ones(X.shape[1]), (1e-1, 1e2)) + WhiteKernel(0.05, (1e-4, 1.0))
    return GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=2, random_state=0).fit(X, y)


def main():
    warnings.filterwarnings("ignore")
    D = load()
    rpm, th = D[:, 0], D[:, 1]
    lk10, lk25, lk50 = np.log(D[:, 2]), np.log(D[:, 3]), np.log(D[:, 4])
    Z = np.column_stack([(rpm - 25) / 7.5, (th - 6.5) / 2.5])
    out = {}
    for name, y in (("kLa_1T_50", lk50), ("kLa_1T_25", lk25)):
        err = {m: [] for m in ("M0", "M1", "M2", "M3")}
        for i in range(len(y)):
            tr = ~((np.abs(rpm - rpm[i]) <= 0.25) & (th == th[i]))
            err["M0"].append(y[tr].mean() - y[i])
            err["M1"].append(gp(Z[tr], y[tr]).predict(Z[i:i + 1])[0] - y[i])
            A = np.column_stack([np.ones(tr.sum()), lk10[tr]])
            c, *_ = np.linalg.lstsq(A, y[tr], rcond=None)
            err["M2"].append(c[0] + c[1] * lk10[i] - y[i])
            Z3 = np.column_stack([Z, (lk10 - lk10.mean()) / lk10.std()])
            err["M3"].append(gp(Z3[tr], y[tr]).predict(Z3[i:i + 1])[0] - y[i])
        rms = {m: float(np.sqrt(np.mean(np.square(e)))) for m, e in err.items()}
        best = min(rms["M2"], rms["M3"])
        verdict = "USEFUL" if best <= 0.8 * rms["M1"] else "NO INFORMATION BEYOND (rpm, theta)"
        out[name] = dict(rms=rms, verdict=verdict, n=len(y))
        print(f"{name}: n {len(y)}  RMS ln error  M0 constant {rms['M0']:.3f}  M1 GP(rpm,theta) {rms['M1']:.3f}  M2 linear(kLa_10) {rms['M2']:.3f}  M3 GP(rpm,theta,kLa_10) {rms['M3']:.3f}  -> {verdict}")
    (ROOT / "experiments/multifidelity/kla_threshold_extrapolation.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
