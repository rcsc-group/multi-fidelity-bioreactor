"""Test: L8 -> L10 two-level warm start vs L9 -> L10 warm start (and cold L10 at 32.5 rpm).

PRE-REGISTERED entry (diary.md, verbatim):
- PRE-REGISTERED (compute, 3 L10 runs, ~1,500 core-h): L8 -> L10 two-level warm start (scripts/submit_l10_from_l8.py; same mechanics and N_CYCLES = 6 as submit_l10_fig13a_xlevel.py, which starts from L9). theta 7, rpm 17.5/25/32.5, source kmix_l8_rpm<r>. Reference: l10_fig13a_xlevel_rpm<r> (L10 from L9). QoIs tau_mean_t, tau_95_t, edr_mean_t over the hydro_dataset_v2 window (last 3 cycles before release).
  PASS if, for every QoI and rpm, |x8 - xlevel| / xlevel <= max(5%, 2 x the cycle-to-cycle relative sd of the xlevel run over the same 3 cycles). Also reported: per-cycle tau/EDR after the restart (does it settle within 3 cycles?).
  If FAIL: use L8 -> L9 -> L10 chains for HF points (cost about 2x) or more L10 cycles.

AMENDMENT (diary.md, "POWER of the pending tests", verbatim):
- L8 -> L10 warm-start validation power, amended: the reference l10_fig13a_xlevel (warm start from L9, 6 cycles) may carry the SAME transient bias as x8 (both 6 cycles after a coarse start); agreement would then not show that either is settled. Added reference: l10c_rpm32.5_seg1 (cold L10, 80 cycles) tau/EDR over the same window, at 32.5 rpm only. Reading at 32.5: x8 must also satisfy the same tolerance against l10c; if x8 matches xlevel but both miss l10c, the 6-cycle warm start is rejected for both.

Window: qois() from scripts/hydro_dataset_v2.py (imported, unmodified). It takes c_end = floor(t_checkpoint/T + n_mix_cycles),
the tracer release cycle (restart runs: checkpoint cycle + 6; cold run: t_checkpoint=0 -> 80), and the window is the last
3 cycles [c_end-3, c_end). Both run types therefore use "last 3 cycles before release". If a run has not reached c_end, qois()
silently clips c_end to the last whole cycle; this script flags that as INCOMPLETE (and does not use it for a verdict).

Per-cycle values (same columns/scales as qois): cycle k = steps with c_end-3+k <= t/T < c_end-3+k+1.
  tau_mean = mean of col 5 * rho U^2; edr = mean of col 9 * rho U^3/L;
  tau_95 per cycle = MEDIAN over the steps of that cycle of the per-step 95th percentile (col 2) * rho U^2, i.e. the
  qois() statistic restricted to one cycle. Cycle-to-cycle relative sd = sample sd (ddof=1) of the 3 per-cycle values / their mean.
Post-restart cycles: 6 cycles starting at t_checkpoint/T (steps with c0+k <= t/T < c0+k+1, k=0..5).

Choices not fixed by the spec (listed): (1) sample sd, ddof=1; (2) tolerance vs l10c uses the xlevel cycle sd at 32.5 rpm
("same tolerance"); (3) rel diff denominator is the reference value; (4) "x8 matches xlevel but both miss l10c" is evaluated
per all-QoI: xlevel "misses" l10c if any QoI fails the tolerance x8-vs-l10c style (|xlevel-l10c|/l10c > tol);
(5) overall PASS requires all 9 QoI/rpm comparisons vs xlevel and all 3 QoI vs l10c to pass; (6) settle check is descriptive:
reports each cycle's deviation of tau_mean/EDR from the mean of the last 3 of the 6 cycles.
Writes experiments/multifidelity/test_warmstart_x8.json.
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.hydro_dataset_v2 import qois, RHO, COL_EDISS  # noqa: E402
from scripts.postprocess import _load_dat, _t_scales  # noqa: E402

RPMS = ["17.5", "25", "32.5"]
QS = ["tau_mean_t", "tau_95_t", "edr_mean_t"]
X8 = "x8_l10_rpm{}_th7"
XL = "l10_fig13a_xlevel_rpm{}"
L10C = "l10c_rpm32.5_seg1"


def load(run):
    d = ROOT / "runs" / run
    if not (d / "params.json").exists() or not (d / "shear_stress.dat").exists():
        return None
    p = json.load(open(d / "params.json"))
    arr = _load_dat(d / "shear_stress.dat")
    if arr is None or len(arr) == 0:
        return None
    T_bio, T = _t_scales(p)
    a = p.get("geometry", {}).get("a", 0.25)
    U = a / T_bio
    return dict(p=p, arr=arr, T=T, scale=RHO * U**2, escale=RHO * U**3 / a, cyc=arr[:, 1] / T)


def cyc_vals(r, lo, n):
    out = []
    for k in range(n):
        m = (r["cyc"] >= lo + k) & (r["cyc"] < lo + k + 1)
        if m.sum() == 0:
            out.append(dict(tau_mean_t=None, tau_95_t=None, edr_mean_t=None, n_steps=0))
            continue
        a = r["arr"][m]
        out.append(dict(tau_mean_t=float(a[:, 5].mean() * r["scale"]), tau_95_t=float(np.median(a[:, 2]) * r["scale"]),
                        edr_mean_t=float(a[:, COL_EDISS].mean() * r["escale"]), n_steps=int(m.sum())))
    return out


def window_info(run):
    """QoIs plus completeness flag (reached release cycle?)."""
    r = load(run)
    if r is None:
        return None
    q = qois(run)
    c_rel = r["p"].get("t_checkpoint", 0.0) / r["T"] + r["p"]["n_mix_cycles"]
    complete = bool(r["cyc"].max() + 1e-6 >= math.floor(c_rel + 1e-6))
    return dict(r=r, q=q, complete=complete, c_rel=float(c_rel), cyc_max=float(r["cyc"].max()))


def main():
    res = {}
    for rpm in RPMS:
        x8, xl = window_info(X8.format(rpm)), window_info(XL.format(rpm))
        e = dict(x8_run=X8.format(rpm), xlevel_run=XL.format(rpm))
        for key, w in (("x8", x8), ("xlevel", xl)):
            if w is None:
                e[key] = "MISSING (no shear_stress.dat)"
            else:
                e[key] = dict(qois=w["q"], complete=w["complete"], release_cycle=w["c_rel"], last_cycle_reached=w["cyc_max"])
        if xl is not None:
            lo = xl["q"]["window"][0]
            pc = cyc_vals(xl["r"], lo, 3)
            sd = {}
            for q in QS:
                v = np.array([c[q] for c in pc], float)
                sd[q] = float(v.std(ddof=1) / v.mean()) if not np.isnan(v).any() else None
            e["xlevel_per_cycle_window"] = pc
            e["xlevel_cycle_rel_sd"] = sd
        if x8 is not None:
            c0 = x8["r"]["p"]["t_checkpoint"] / x8["r"]["T"]
            pcs = cyc_vals(x8["r"], c0, 6)
            e["x8_post_restart_cycles"] = pcs
            e["x8_post_restart_start_cycle"] = float(c0)
            if all(c["n_steps"] > 0 for c in pcs):
                ref = {q: np.mean([c[q] for c in pcs[3:]]) for q in ("tau_mean_t", "edr_mean_t")}
                e["x8_settle_dev_from_mean_cycles4_6"] = [
                    {q: float(c[q] / ref[q] - 1) for q in ref} for c in pcs]
        e["_x8"], e["_xl"] = x8, xl
        res[rpm] = e

    l10c = window_info(L10C)
    verdict_rows, ok_all, any_missing = [], True, False
    for rpm in RPMS:
        e = res[rpm]
        x8, xl = e.pop("_x8"), e.pop("_xl")
        if x8 is None or xl is None or not x8["complete"] or not xl["complete"]:
            any_missing = True
            e["comparison_vs_xlevel"] = "INCOMPLETE/MISSING"
            continue
        comp = {}
        for q in QS:
            sd = e["xlevel_cycle_rel_sd"][q]
            tol = max(0.05, 2 * sd)
            rd = abs(x8["q"][q] - xl["q"][q]) / xl["q"][q]
            comp[q] = dict(x8=x8["q"][q], xlevel=xl["q"][q], rel_diff=float(rd), tol=float(tol), verdict="PASS" if rd <= tol else "FAIL")
            ok_all &= rd <= tol
        e["comparison_vs_xlevel"] = comp
        if rpm == "32.5":
            if l10c is None or not l10c["complete"]:
                any_missing = True
                e["comparison_vs_l10c"] = "INCOMPLETE/MISSING" + ("" if l10c is None else f" (cold run at cycle {l10c['cyc_max']:.1f} of {l10c['c_rel']:.0f})")
            else:
                c2, c3 = {}, {}
                for q in QS:
                    tol = comp[q]["tol"]
                    rd = abs(x8["q"][q] - l10c["q"][q]) / l10c["q"][q]
                    rx = abs(xl["q"][q] - l10c["q"][q]) / l10c["q"][q]
                    c2[q] = dict(x8=x8["q"][q], l10c=l10c["q"][q], rel_diff=float(rd), tol=tol, verdict="PASS" if rd <= tol else "FAIL")
                    c3[q] = dict(xlevel=xl["q"][q], l10c=l10c["q"][q], rel_diff=float(rx), tol=tol, verdict="PASS" if rx <= tol else "FAIL")
                    ok_all &= rd <= tol
                e["comparison_vs_l10c"] = c2
                e["xlevel_vs_l10c"] = c3
                e["l10c_window"] = l10c["q"]["window"]
    for rpm in RPMS:
        pass
    if any_missing:
        overall = "INCOMPLETE (runs missing/unfinished; see per-rpm entries)"
    elif ok_all:
        overall = "PASS"
    else:
        e = res["32.5"]
        x8_ok_xl = all(v["verdict"] == "PASS" for v in e["comparison_vs_xlevel"].values())
        x8_miss = any(v["verdict"] == "FAIL" for v in e["comparison_vs_l10c"].values())
        xl_miss = any(v["verdict"] == "FAIL" for v in e["xlevel_vs_l10c"].values())
        overall = ("FAIL: x8 matches xlevel but both miss l10c -> 6-cycle warm start rejected for both"
                   if x8_ok_xl and x8_miss and xl_miss else "FAIL")
    out = dict(overall=overall, per_rpm=res, l10c_run=L10C, l10c_status=None if l10c is None else dict(
        complete=l10c["complete"], last_cycle_reached=l10c["cyc_max"], release_cycle=l10c["c_rel"], qois=l10c["q"]))
    json.dump(out, open(ROOT / "experiments/multifidelity/test_warmstart_x8.json", "w"), indent=1)

    print("OVERALL:", overall)
    for rpm in RPMS:
        e = res[rpm]
        print(f"--- {rpm} rpm: x8={e['x8'] if isinstance(e['x8'], str) else 'present, complete=%s' % e['x8']['complete']}; "
              f"xlevel complete={e['xlevel']['complete'] if isinstance(e['xlevel'], dict) else e['xlevel']}")
        if "xlevel_cycle_rel_sd" in e:
            print("  xlevel cycle rel sd:", {k: round(v, 4) for k, v in e["xlevel_cycle_rel_sd"].items()})
        for k in ("comparison_vs_xlevel", "comparison_vs_l10c", "xlevel_vs_l10c"):
            v = e.get(k)
            if isinstance(v, dict):
                for q, d in v.items():
                    print(f"  {k} {q}: rel_diff={d['rel_diff']:.4f} tol={d['tol']:.4f} {d['verdict']}")
            elif v:
                print(" ", k, v)
        if "x8_post_restart_cycles" in e:
            for i, c in enumerate(e["x8_post_restart_cycles"]):
                print(f"  x8 cycle {i+1}: tau_mean={c['tau_mean_t']} edr={c['edr_mean_t']} n={c['n_steps']}")
    print("l10c:", out["l10c_status"] if l10c is None else {k: v for k, v in out["l10c_status"].items() if k != "qois"})


if __name__ == "__main__":
    main()
