"""MF test: predict the C*=50% quantity from C*=10% runs (L9, rpm sweep, 7 deg).

The advisor's proposal one step down (no run reaches 90%): LF = the cheap
early-saturation value at every rpm, HF = the late value at 3 rpm. Three
targets:
  kLa : kLa_1T_10 -> kLa_1T_50          (whole-period fits, h^-1)
  time: cycles to C*=10% -> cycles to C*=50%  (after release)
  mix : dtmix at chi=0.50 -> dtmix at chi=0.95 (s)
Training sets: both endpoints (15, 37.5 rpm) + one interior rpm, all 7
choices; the other 6 rpm are scored. Same model as the Fig 13 MF panels
(plot_mf_fig13_unc: LOO-KRR on LF, then universal kriging with known noise, by REML).
HF noise: sigma_i = 5% of y_i. That is the 1-period kLa round-off spread
(6.8% between two repeats -> ~4.8% per run, diag_kla_phase_lock.py).
Baselines: GP on the 3 HF points alone; LF times the mean HF/LF ratio of the training points.
Falsifier for the proposal: MF does not beat both baselines in mean held-out error.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.plot_mf_fig13 import LooKRR                        # noqa: E402
from scripts.plot_mf_fig13_unc import HeteroUK                  # noqa: E402
from scripts.diag_early_late_transfer import cycles_to          # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]
REL_NOISE = 0.05


def data():
    rows = []
    for r in RPMS:
        run = f"fig9_l9_rpm{r:g}"
        d = json.load(open(ROOT / "runs" / run / "results.json"))
        t10, _ = cycles_to(run, 0.10)
        t50, _ = cycles_to(run, 0.50)
        rows.append((r, d.get("kLa_1T_10", math.nan), d.get("kLa_1T_50", math.nan), t10, t50,
                     d.get("dtmix_0.50", math.nan), d.get("dtmix_0.95", math.nan)))
    a = np.array(rows, float)
    return a[np.all(np.isfinite(a), axis=1)]


def one_case(x, lf, hf, train):
    lo, hi = x.min(), x.max()
    krr = LooKRR(design_space=np.array([[lo, hi]]), params_optimize=True,
                 noise_data=True, optimizer_restart=10, seed=42)
    krr.train(x.reshape(-1, 1), lf)
    f_lf = lambda z: np.ravel(krr.predict(np.asarray(z, float).reshape(-1, 1)))
    test = ~train
    sig = REL_NOISE * hf[train]
    mf = HeteroUK(lo, hi, lambda z: np.column_stack([np.ones_like(z), f_lf(z)])).fit(x[train], hf[train], sig)
    sf = HeteroUK(lo, hi, lambda z: np.ones((len(z), 1))).fit(x[train], hf[train], sig)
    out = {}
    for name, m in (("MF", mf), ("HF-only GP", sf)):
        mu, sd = m.predict(x[test])
        z = (mu - hf[test]) / np.sqrt(sd ** 2 + (REL_NOISE * hf[test]) ** 2)
        out[name] = (np.sqrt(np.mean(((mu - hf[test]) / hf[test]) ** 2)), np.mean(np.abs(z) <= 1.96))
    ratio = np.mean(hf[train] / lf[train])
    out["LF x ratio"] = (np.sqrt(np.mean(((ratio * lf[test] - hf[test]) / hf[test]) ** 2)), math.nan)
    return out


def main() -> None:
    a = data()
    x = a[:, 0]
    print(f"rpm with both ends: {list(x)}")
    for label, lf, hf in (("kLa 10% -> 50%", a[:, 1], a[:, 2]),
                          ("cycles to 10% -> to 50%", a[:, 3], a[:, 4]),
                          ("dtmix chi 0.50 -> 0.95", a[:, 5], a[:, 6])):
        print(f"\n== {label}   corr(log LF, log HF) = {np.corrcoef(np.log(lf), np.log(hf))[0, 1]:+.2f}")
        res = {}
        for mid in x[1:-1]:
            train = np.isin(x, [x[0], mid, x[-1]])
            for k, v in one_case(x, lf, hf, train).items():
                res.setdefault(k, []).append(v)
            print(f"  train 15/{mid:g}/37.5: " + "  ".join(
                f"{k} {v[0]:.0%}" for k, v in one_case(x, lf, hf, train).items()))
        for k, v in res.items():
            e = np.array([q[0] for q in v]); c = np.array([q[1] for q in v])
            cov = "" if np.all(np.isnan(c)) else f"  mean 95% coverage {np.nanmean(c):.0%}"
            print(f"  {k:11s} mean held-out rel. RMSE {e.mean():.0%} (range {e.min():.0%}-{e.max():.0%}){cov}")


if __name__ == "__main__":
    main()
