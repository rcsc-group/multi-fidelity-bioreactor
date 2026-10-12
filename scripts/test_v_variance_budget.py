"""Test V (pre-registered, diary 2026-10-12): tracer variance budget, r = D_eff / D_phys.

In a closed incompressible liquid, advection conserves tracer variance; only diffusion destroys it:
    d sigma^2/dt = -2 D_eff <|grad c|^2>.
The solver applies D_phys (c2.D1 = 1/Pe) to the gradients it resolves, so
    r(t) = -(d sigma^2/dt) / (2 D_phys <|grad c|^2>)
is the factor by which the variance actually destroyed exceeds what physical diffusion on the resolved field
destroys. Nondimensional units: x in L_bio (domain L0 = 1, dx = 1/n), t in T_bio, D_phys = D_tracer_1 T_bio / L_bio^2.

sigma^2: the solver's liquid variance (chi_scale.sigma2_native, same as chi). <|grad c|^2>: mean of squared FACE
differences (dc * n)^2 over x-faces plus over y-faces, each over faces whose two cells are pure liquid
(f, cs >= 0.999). d sigma^2/dt: linear fit of sigma^2 over one period centred on each snapshot.
Rule: NUMERICAL-DOMINATED if median r (chi 0.5-0.95) >= 10 in every scal run; PHYSICAL if <= 2 in every scal run;
else MIXED. Usage: uv run python scripts/test_v_variance_budget.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.chi_scale import G_LIQ, _solver_reference, iter_snapshots, locate, sigma2_native  # noqa: E402

D_TRACER_1, L_BIO = 0.44e-9, 0.25  # src/BioReactor.c: D_tracer_1 (m^2/s), L_bio = geometry_a (m)
SCAL = [f"scal_l{L}_rpm{r}_th7" for r in ("32.5", "37.5") for L in (8, 9)]
REC = [f"rec_l8_rpm{r}_th7" for r in ("22.4", "22.5", "22.6", "25", "29.9", "30", "30.1", "32.5")]


def grad2(f, cs, c):
    n = c.shape[0]
    liq = (f >= G_LIQ) & (cs >= G_LIQ)
    mx = liq[:, 1:] & liq[:, :-1]
    my = liq[1:, :] & liq[:-1, :]
    gx = ((c[:, 1:] - c[:, :-1]) * n)[mx]
    gy = ((c[1:, :] - c[:-1, :]) * n)[my]
    if gx.size == 0 or gy.size == 0:
        return np.nan
    return float((gx * gx).mean() + (gy * gy).mean())


def series(run_id):
    loc = locate(run_id)
    if loc["fields_dir"] is None or loc["data_dir"] is None:
        return None
    ref = _solver_reference(loc)
    T_bio, T_per = ref["T_bio"], ref["T_per_nd"]
    d_phys = D_TRACER_1 * float(ref["params"].get("tracer_D_scale", 1.0)) * T_bio / L_BIO ** 2
    t, s2, g2 = [], [], []
    started = False
    for snap in iter_snapshots(loc["fields_dir"], ["f", "cs", "c2"]):
        f, cs, c = snap["f"], snap["cs"], snap["c2"]
        if not started:
            if float((cs * f * c).sum()) <= 1e-10 * float((cs * f).sum()):
                continue
            started = True
        t.append(snap["t"])
        s2.append(sigma2_native(f, cs, c, ref["f_mean"]))
        g2.append(grad2(f, cs, c))
    t, s2, g2 = map(np.array, (t, s2, g2))
    chi = 1 - s2 / ref["sigma0"]
    rate = np.full_like(t, np.nan)
    for i in range(len(t)):
        m = np.abs(t - t[i]) <= T_per / 2
        if m.sum() >= 5 and t[i] - T_per / 2 >= t[0] and t[i] + T_per / 2 <= t[-1]:
            rate[i] = -np.polyfit(t[m], s2[m], 1)[0]
    r = rate / (2 * d_phys * g2)
    return dict(run_id=run_id, t_s=(t - ref["t0"]) * T_bio, chi=chi, r=r, d_phys_nd=d_phys,
                n=int(np.sqrt(f.size)))


def main():
    out = {}
    for run in SCAL + REC:
        s = series(run)
        if s is None:
            print(f"{run}: missing")
            continue
        w = (s["chi"] >= 0.5) & (s["chi"] <= 0.95) & np.isfinite(s["r"])
        rw = s["r"][w]
        q = np.percentile(rw, [10, 50, 90]) if rw.size else [np.nan] * 3
        out[run] = dict(n_window=int(w.sum()), r_p10=float(q[0]), r_median=float(q[1]), r_p90=float(q[2]),
                        d_phys_nd=s["d_phys_nd"], grid_n=s["n"],
                        chi=s["chi"].tolist(), r=s["r"].tolist(), t_s=s["t_s"].tolist())
        print(f"{run:24s} n={s['n']:4d}  window snapshots {w.sum():4d}  r = D_eff/D_phys: "
              f"p10 {q[0]:8.1f}  median {q[1]:8.1f}  p90 {q[2]:8.1f}")
    med = [out[r]["r_median"] for r in SCAL if r in out]
    if len(med) < len(SCAL):
        verdict = "PENDING (missing scal runs)"
    elif all(m >= 10 for m in med):
        verdict = "NUMERICAL-DOMINATED"
    elif all(m <= 2 for m in med):
        verdict = "PHYSICAL"
    else:
        verdict = "MIXED"
    print("VERDICT:", verdict)
    (ROOT / "experiments/multifidelity/test_v_variance_budget.json").write_text(
        json.dumps(dict(verdict=verdict, runs=out), indent=1))


if __name__ == "__main__":
    main()
