"""edr_95 backfill LIST (lead, 2026-11-04). Lists the database runs the fit uses for the EDR QoI that lack edr_95,
checks each has a checkpoint to restart from, and costs 3 cycles per run from sacct per-level rates.
Writes state/edr_backfill_list.csv. Submits NOTHING.

Rates: from the cost audit table (sacct, core-h per unit nondim time, per fidelity), median over runs with a span.
Cost per run = rate[level] * 3 * T_run, T_run = nondim period of the run (postprocess._t_scales).
Usage (compute job): uv run python scripts/campaign2d/backfill_list.py
"""
import csv
import json
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402

COST_TABLE = Path("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv")
OUT = ROOT / "scripts" / "campaign2d" / "state" / "edr_backfill_list.csv"
N_CYC = 3
RPM = (15.0, 37.5)
THETA = (2.0, 7.0)


def rates():
    """core-h per unit nondim time, median per fidelity, from sacct rows with a positive span."""
    per = {}
    for ln in COST_TABLE.read_text().splitlines()[1:]:
        f = ln.split("\t")
        try:
            fid, tf, tl, cc = int(f[1]), float(f[4]), float(f[5]), float(f[8])
        except ValueError:          # 'None' fields: run without a sacct span, skipped
            continue
        if tl > tf and cc > 0:
            per.setdefault(fid, []).append(cc / (tl - tf))
    return {k: statistics.median(v) for k, v in per.items()}


def main():
    kla = {r["run"]: r for r in __import__("csv").DictReader(open(ROOT / "experiments/multifidelity/kla_mix_dataset.csv"))}
    hyd = json.loads((ROOT / "experiments/multifidelity/hydro_dataset_v2.json").read_text())
    R = rates()
    rows = []
    for h in hyd:
        lev, run = int(h["level"]), h["run"]
        k = kla.get(run)
        if k is None:
            continue
        rpm, th = float(k["rpm"]), float(k["theta"])
        if not (RPM[0] <= rpm <= RPM[1] and THETA[0] <= th <= THETA[1]):
            continue
        d = ROOT / "runs" / run
        p = json.loads((d / "params.json").read_text())
        _, T = pp._t_scales(p)
        ck = (d / "checkpoint.dump").exists()
        cost = R.get(lev, math.nan) * N_CYC * T
        rows.append(dict(level=lev, run=run, rpm=rpm, theta=th, checkpoint=ck, core_h=round(cost, 2)))
    rows.sort(key=lambda r: (r["level"], r["rpm"], r["theta"]))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for lev in (8, 9, 10):
        sub = [r for r in rows if r["level"] == lev]
        tot = sum(r["core_h"] for r in sub if not math.isnan(r["core_h"]))
        print(f"L{lev}: {len(sub)} runs, no checkpoint {sum(not r['checkpoint'] for r in sub)}, "
              f"rate {R.get(lev, math.nan):.4g} core-h per nondim time, cost {tot:.0f} core-h")
    print(f"TOTAL {sum(r['core_h'] for r in rows if not math.isnan(r['core_h'])):.0f} core-h over {len(rows)} runs; "
          f"list written to {OUT} (not submitted)")


if __name__ == "__main__":
    main()
