"""Pilot: predict dtmix at chi = 0.75 and 0.95 from the chi(t) curve up to 0.50 only.

Model: late mixing is dominated by the slowest mode, so
    1 - chi(t) ~ A exp(-lam t).
Fit ln(1 - chi) vs t by least squares on the samples with chi in [0.25, 0.50]
(the part of the curve that a run stopped at chi = 0.50 has). Then
    t(chi*) = (ln A - ln(1 - chi*)) / lam.
Compared with the measured dtmix_0.75 and dtmix_0.95, L9 rpm sweep and L8.
Falsifier: median |error| > 30% at chi = 0.95.
"""
import json, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _joined_cols, _t_scales, _COL_C2_LIQ_SUM, _COL_C2_LIQ_SUM2, _COL_F_LIQ  # noqa: E402

RPMS = [15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5]


def chi_curve(run):
    d = ROOT / "runs" / run
    t, s, s2 = _joined_cols(d, "tr_oxy.dat", [_COL_C2_LIQ_SUM, _COL_C2_LIQ_SUM2])
    _, f = _joined_cols(d, "vol_frac_interf.dat", [_COL_F_LIQ])
    fm = f.mean()
    var = s2 / fm - (s / fm) ** 2
    i0 = int(np.argmax(s > 1e-10 * fm))
    chi = np.clip(1 - var / var[i0], 0, 1)
    T_bio, _ = _t_scales(json.load(open(d / "params.json")))
    return (t[i0:] - t[i0]) * T_bio, chi[i0:]


def main():
    for tmpl in ("fig9_l9_rpm{:g}", "kmix_l8_rpm{:g}"):
        print(f"\n== {tmpl}")
        errs = {0.75: [], 0.95: []}
        for r in RPMS:
            run = tmpl.format(r)
            res = json.load(open(ROOT / "runs" / run / "results.json"))
            t, chi = chi_curve(run)
            k50 = int(np.argmax(chi >= 0.50))
            if chi[k50] < 0.50:
                continue
            m = (chi >= 0.25) & (np.arange(len(t)) <= k50)
            b, a = np.polyfit(t[m], np.log(1 - chi[m]), 1)      # ln(1-chi) = a + b t
            line = f"  {r:5g} rpm  lam={-b:.4f}/s"
            for c in (0.75, 0.95):
                pred = (math.log(1 - c) - a) / b
                meas = res.get(f"dtmix_{c:.2f}", math.nan)
                e = pred / meas - 1 if meas == meas else math.nan
                if e == e:
                    errs[c].append(abs(e))
                line += f" | chi {c}: pred {pred:7.1f} s, meas {meas:7.1f} s ({e:+.0%})"
            print(line)
        for c, v in errs.items():
            print(f"  median |error| at chi {c}: {np.median(v):.0%}  (n={len(v)})")


if __name__ == "__main__":
    main()
