"""Physics-informed law for kLa (zero compute): does ln kLa(C* = 0.10) follow the small-eddy scaling kL a ~ a (eps nu)^(1/4) Sc^(-1/2), i.e. a slope of 0.25 on ln EDR?
Per level: regression ln kLa = c + s ln EDR (+ optional ln tau_95), leave-one-out RMS ln error, slope s, correlation. Levels: L8 (89 grid runs), L9 (theta 7), L10 (13+ runs).
EDR and tau_95 are the 3-cycle window QoIs (hydro_dataset_v2.qois) of the same flow state (before the oxygen release; L10: the training window run).
Usage: uv run python scripts/physics_law_kla_edr.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois  # noqa: E402
from scripts.test_h1_mfbml import hf_run  # noqa: E402

Q = "kLa_1T_10"


def val(run, q=Q):
    f = ROOT / "runs" / run / "results.json"
    v = json.loads(f.read_text()).get(q) if f.exists() else None
    return float(v) if v and v > 0 else np.nan


def rpm_th(run):
    p = json.loads((ROOT / "runs" / run / "params.json").read_text())
    return p["omega_b"] * 60 / (2 * math.pi), p["theta_max"][0]


def table():
    rows = {8: [], 9: [], 10: []}
    for pat in ("t6_l8_*", "t5_l8_*", "e8_l8_*"):
        for d in sorted((ROOT / "runs").glob(pat)):
            k = val(d.name)
            try:
                h = qois(d.name)
            except Exception:
                continue
            if np.isfinite(k) and h:
                rows[8].append((*rpm_th(d.name), k, h["edr_mean_t"], h["tau_95_t"], h["tau_mean_t"]))
    for d in sorted((ROOT / "runs").glob("fig9_l9_rpm*")):
        if "_ext" in d.name:
            continue
        k = val(d.name)
        try:
            h = qois(d.name)
        except Exception:
            continue
        if np.isfinite(k) and h:
            rows[9].append((*rpm_th(d.name), k, h["edr_mean_t"], h["tau_95_t"], h["tau_mean_t"]))
    for d in sorted((ROOT / "runs").glob("t10_l10_*")) + [ROOT / "runs" / "v_l10_rpm32.5_th7"]:
        k = val(d.name)
        r, t = rpm_th(d.name)
        src = f"x8_l10_rpm{r:g}_th{t:g}" if t == 5.5 else hf_run(r, t)
        try:
            h = qois(src)
        except Exception:
            continue
        if np.isfinite(k) and h:
            rows[10].append((r, t, k, h["edr_mean_t"], h["tau_95_t"], h["tau_mean_t"]))
    return {L: np.array(v) for L, v in rows.items() if len(v)}


def loo(X, y):
    err = []
    for i in range(len(y)):
        tr = np.arange(len(y)) != i
        A = np.column_stack([np.ones(tr.sum()), X[tr]])
        c, *_ = np.linalg.lstsq(A, y[tr], rcond=None)
        err.append(float(np.r_[1, X[i]] @ c - y[i]))
    return float(np.sqrt(np.mean(np.square(err))))


def main():
    T = table()
    for L, a in T.items():
        y = np.log(a[:, 2])
        e, t95, tm = np.log(a[:, 3]), np.log(a[:, 4]), np.log(a[:, 5])
        A = np.column_stack([np.ones(len(y)), e])
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = np.corrcoef(e, y)[0, 1]
        print(f"L{L}: n {len(y)}  ln kLa = {c[0]:+.2f} + {c[1]:.2f} ln EDR  (r {r:+.2f}); LOO RMS ln error: EDR only {loo(e[:, None], y):.3f}; "
              f"EDR + tau_95 {loo(np.column_stack([e, t95]), y):.3f}; constant {np.sqrt(np.mean((y - y.mean()) ** 2)):.3f}; "
              f"EDR + rpm + theta {loo(np.column_stack([e, a[:, 0], a[:, 1]]), y):.3f}; rpm + theta only {loo(np.column_stack([a[:, 0], a[:, 1]]), y):.3f}")
    np.save(ROOT / "experiments/multifidelity/physics_law_table.npy", np.array([[L, *row] for L, a in T.items() for row in a]))


if __name__ == "__main__":
    main()
