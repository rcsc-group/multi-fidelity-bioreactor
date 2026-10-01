"""Is the 5-sample kLa estimator phase-locked to the rocking cycle?

The estimator (postprocess._kla_5pt_at_threshold) fits ln(1-C*) against t over the
5 output samples around the first C* >= theta crossing. With t_out = 0.02 and
T_per = 0.607, those 5 samples span 0.13 of a rocking period.

Hypothesis H: oxygen transfer is modulated within each cycle (interface area and
stirring change with phase), so the 5-sample slope measures the rate at the
PHASE where the crossing lands, not a cycle-averaged rate.

Predictions under H, each falsifiable:
 P1  the sliding 5-sample kLa(t) oscillates at f_b or 2 f_b with an amplitude
     >~ 20% of its one-period running mean, around C* = 0.25;
 P2  a fit over exactly one (or two) whole periods centred on the crossing
     differs from the 5-sample value by a phase-dependent amount;
 P3  for two runs that differ only by round-off (OMP=4 repeats, L4), the
     whole-period kLa agrees much better than the 5-sample kLa.
H is rejected if the P1 amplitude is < 10% AND P3 shows no improvement.

Usage: uv run python scripts/diag_kla_phase_lock.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.postprocess import _compute_c_star, _kla_5pt_at_threshold, _t_scales  # noqa: E402

RUNS = ([f"fig9_l9_rpm{r:g}" for r in (15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5)]
        + [f"fig10_l9_th{a:g}" for a in (2, 3, 4, 5, 6, 7)]
        + ["fig11_l10_rpm35", "fig11_l10_rpm37.5"]
        + ["diag_omp/omp4_a", "diag_omp/omp4_b", "diag_omp/omp1_a"])
THR = 0.25


def slope_fit(t, c):
    y = np.log(1.0 - np.clip(c, 0.0, 1.0 - 1e-10))
    return -np.polyfit(t, y, 1)[0]


def analyse(run):
    d = ROOT / "runs" / run
    p = json.load(open(d / "params.json"))
    _, T = _t_scales(p)
    t, c = _compute_c_star(d)
    if np.nanmax(c) < THR:
        return None
    i0 = int(np.argmax(c > 1e-6))                      # release
    k = int(np.argmax(c >= THR))
    dt = float(np.median(np.diff(t)))
    n_per = int(round(T / dt))

    # P1: sliding 5-sample kLa around the crossing, vs its one-period running mean
    lo, hi = max(i0 + 2, k - 2 * n_per), min(len(t) - 3, k + 2 * n_per)
    idx = np.arange(lo, hi)
    k5 = np.array([slope_fit(t[j - 2:j + 3], c[j - 2:j + 3]) for j in idx])
    ker = np.ones(n_per) / n_per
    smooth = np.convolve(k5, ker, mode="same")
    core = slice(n_per // 2, len(k5) - n_per // 2)
    resid = (k5 - smooth)[core]
    amp = float(np.std(resid) * math.sqrt(2) / np.mean(smooth[core]))   # sine amplitude / mean
    spec = np.abs(np.fft.rfft(resid - resid.mean()))
    freqs = np.fft.rfftfreq(len(resid), dt) * T                         # in units of f_b
    f_peak = float(freqs[1:][np.argmax(spec[1:])]) if len(spec) > 1 else float("nan")

    # P2: whole-period fits centred on the crossing
    def window_fit(n_periods):
        h = int(round(n_periods * n_per / 2))
        a, b = max(i0, k - h), min(len(t), k + h + 1)
        return slope_fit(t[a:b], c[a:b])
    k_5 = _kla_5pt_at_threshold(t, c, THR)
    phase = ((t[k] - t[i0]) / T) % 1.0
    return dict(run=run, k5=k_5, k1=window_fit(1), k2=window_fit(2), amp=amp,
                f_peak=f_peak, phase=phase, n_per=n_per)


def main() -> None:
    rows = [r for r in (analyse(x) for x in RUNS) if r]
    print(f"{'run':22s} {'kLa_5s':>8} {'kLa_1T':>8} {'kLa_2T':>8} {'5s/1T':>6} "
          f"{'osc amp':>8} {'f_peak/f_b':>10} {'phase':>6}")
    for r in rows:
        print(f"{r['run']:22s} {r['k5']:8.4f} {r['k1']:8.4f} {r['k2']:8.4f} "
              f"{r['k5'] / r['k1']:6.3f} {r['amp']:8.1%} {r['f_peak']:10.2f} {r['phase']:6.2f}")
    a = {r["run"]: r for r in rows}
    if "diag_omp/omp4_a" in a and "diag_omp/omp4_b" in a:
        x, y = a["diag_omp/omp4_a"], a["diag_omp/omp4_b"]
        for key, name in (("k5", "5-sample"), ("k1", "1-period"), ("k2", "2-period")):
            print(f"P3 OMP=4 repeat spread, {name:9s}: "
                  f"{abs(x[key] - y[key]) / (0.5 * (x[key] + y[key])):.1%}")


if __name__ == "__main__":
    main()
