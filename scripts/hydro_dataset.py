"""Hydrodynamic QoIs (tau_mean_max, tau_95_qss) at theta = 7 deg on L6-L10, from ONE consistent code family.

Consistency audit (diary 2026-10-05): shear values changed with the 2026-09-15 mask fix; every run before it
is excluded. The families used are all post-fix and agree where they overlap (32.5 rpm):
  L6-L8: production "lean" binaries (kmix_/mf_/fig9 ladder runs, continuous from rest)
  L9:    fig9_l9_rpm* (prod-4b3a435; identical to the lean run fig10_l9_th7 at 32.5 rpm)
  L10:   l10_fig13a_xlevel_rpm* (cross-level restart, measured over 6 cycles); at 32.5 rpm it agrees with the
         cold L10 run kimcheck_l10_rpm32.5 within 0.5% (tau_mean) and 2% (tau_95)
Writes experiments/multifidelity/hydro_dataset.json: list of {level, rpm, run, tau_mean_max, tau_95_qss}.
"""
import glob
import json
import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "runs"
FIX_DATE = 1757908800  # 2026-09-15 00:00 UTC: shear mask fix
LEAN_PREFIXES = ("kmix_l", "mf_l", "fig9_ladder_l", "fig9_lean_l")


def rpm_of(p):
    return round(p["omega_b"] * 60 / (2 * math.pi), 1)


def load(run):
    d = RUNS / run
    return json.load(open(d / "params.json")), json.load(open(d / "results.json")), os.path.getmtime(d / "results.json")


def main():
    rows = {}
    for f in glob.glob(str(RUNS / "*/results.json")):
        run = Path(f).parent.name
        try:
            p, r, mt = load(run)
        except Exception:
            continue
        if (p.get("theta_max") or [None])[0] != 7.0 or mt < FIX_DATE:
            continue
        lev, b = p.get("fidelity"), os.path.basename(p.get("_binary", ""))
        ok = ((lev in (6, 7, 8) and "lean" in b and run.startswith(LEAN_PREFIXES) and p.get("t_checkpoint", 0) == 0)
              or (lev == 9 and run.startswith("fig9_l9_rpm") and "_ext" not in run)
              or (lev == 10 and run.startswith("l10_fig13a_xlevel_rpm")))
        if not ok:
            continue
        tm, t95 = r.get("tau_mean_max"), r.get("tau_95_qss")
        if tm is None or t95 is None or not (math.isfinite(tm) and math.isfinite(t95)):
            continue
        key = (lev, rpm_of(p))
        if key not in rows or mt > rows[key]["mtime"]:  # latest run per (level, rpm)
            rows[key] = dict(level=lev, rpm=rpm_of(p), run=run, tau_mean_max=tm, tau_95_qss=t95, mtime=mt)
    out = sorted(rows.values(), key=lambda r: (r["level"], r["rpm"]))
    for r in out:
        r.pop("mtime")
        print(f"L{r['level']} {r['rpm']:5} tau_mean={r['tau_mean_max']:.5f} tau95={r['tau_95_qss']:.5f}  {r['run']}")
    path = ROOT / "experiments/multifidelity/hydro_dataset.json"
    json.dump(out, open(path, "w"), indent=1)
    print("wrote", path, len(out), "rows")


if __name__ == "__main__":
    main()
