"""Round 10: simulated cycles of warm runs, QoI-window check, and core-h per cycle with and without start-up.

Source files (mtime at writing):
  runs/<run>/logstats.dat   (solver time t and wall clock per logged step; restart runs start at t_checkpoint)
  runs/<run>/params.json    (T per run via scripts.postprocess._t_scales)
  experiments/cost_table.tsv in /oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv (corehc_all, ncpu; 2026-11-04 09:22)
  sacct node list: /oscar/scratch/eaguerov/tmp/cost_audit/sacct_jobs.psv (2026-11-04 09:21), jobmap_raw.txt
Per run:
  cycles_simulated = (t_last - t_first) / T   (t_first = first logged t)  [step A1]
  core-h per cycle (full)   = (wall_last - wall_first) / cycles x ncpu / 3600
  core-h per cycle (steady) = (wall_last - wall(t_first + 2T)) / (cycles - 2) x ncpu / 3600    [excludes the first 2 cycles]
  QoI window check (hydro_dataset_v2.qois): window = [c_end-3, c_end] in absolute cycles; must lie inside
  [t_first/T, t_last/T] for the run to contain it.
Output: data/warm_rate_r10.json
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
from scripts.postprocess import _t_scales  # noqa: E402
from scripts.hydro_dataset_v2 import qois  # noqa: E402
import scripts.test_h1_mfbml as H1  # noqa: E402
import scripts.test_h2_mfbml as H2  # noqa: E402

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
COST = {}
with open("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        COST[r["run"]] = r


def logstats(run):
    rows = []
    for line in open(REPO / "runs" / run / "logstats.dat"):
        if "t:" not in line or "Wall" not in line:
            continue
        parts = line.replace(":", " ").split()
        t = float(parts[parts.index("t") + 1])
        w = float(line.split("Wall clock time (s):")[1].split()[0])
        rows.append((t, w))
    rows = sorted(set(rows))
    return np.array(rows)


def wall_at(arr, t):
    return float(np.interp(t, arr[:, 0], arr[:, 1]))


def measure(run):
    p = json.load(open(REPO / "runs" / run / "params.json"))
    _, T = _t_scales(p)
    arr = logstats(run)
    t0, t1 = arr[0, 0], arr[-1, 0]
    w0, w1 = arr[0, 1], arr[-1, 1]
    cyc = (t1 - t0) / T
    ncpu = int(COST[run]["ncpu"])
    full = (w1 - w0) / cyc * ncpu / 3600.0 if cyc > 0 else None
    steady = None
    if cyc > 2.5:
        w2 = wall_at(arr, t0 + 2 * T)
        steady = (w1 - w2) / (cyc - 2) * ncpu / 3600.0
    q = qois(run)
    c_end = q["window"][1]
    c_lo = q["window"][0]
    inside = (c_lo >= t0 / T - 1e-6) and (c_end <= t1 / T + 1e-6)
    return dict(run=run, T=T, t_first=t0, t_last=t1, cycles=cyc, ncpu=ncpu, cost_sacct=float(COST[run]["corehc_all"]),
                wall_cost=(w1 - w0) * ncpu / 3600.0, full_rate=full, steady_rate=steady,
                window=[c_lo, c_end], window_inside_simulated=bool(inside),
                t_checkpoint=float(p.get("t_checkpoint", 0.0)), n_mix_cycles=int(p.get("n_mix_cycles", 0)))


groups = {
    "warm_L9_x8": [H2.l9_run(r, t) for r, t in H2.l9_points() if t != 7],
    "warm_L9_w9_th7": ["w9_l9_rpm17.5_th7", "w9_l9_rpm25_th7", "w9_l9_rpm32.5_th7"],
    "cold_L9_fig9": [H2.l9_run(r, t) for r, t in H2.l9_points() if t == 7],
    "warm_L10": [H1.hf_run(r, t) for r, t in H1.hf_points()],
}
out = {}
for g, names in groups.items():
    rows = []
    for n in names:
        if n not in COST or not (REPO / "runs" / n / "logstats.dat").exists():
            print("skip", g, n)
            continue
        try:
            m = measure(n)
        except Exception as e:  # noqa: BLE001
            print("fail", g, n, type(e).__name__, e)
            continue
        rows.append(m)
        print(f"{g:16s} {n:26s} cycles {m['cycles']:7.1f} ncpu {m['ncpu']:3d} full {m['full_rate'] if m['full_rate'] is None else round(m['full_rate'], 2)} "
              f"steady {None if m['steady_rate'] is None else round(m['steady_rate'], 2)} "
              f"window {m['window'][0]:.1f}-{m['window'][1]:.1f} inside {m['window_inside_simulated']}")
    cyc = [m["cycles"] for m in rows]
    full = [m["full_rate"] for m in rows if m["full_rate"] is not None]
    st = [m["steady_rate"] for m in rows if m["steady_rate"] is not None]
    out[g] = dict(n=len(rows), median_cycles=float(np.median(cyc)) if cyc else None,
                  median_full_rate=float(np.median(full)) if full else None,
                  median_steady_rate=float(np.median(st)) if st else None,
                  n_steady=len(st), rows=rows)
    print(f"== {g}: n={len(rows)} median cycles {out[g]['median_cycles']}, median full {out[g]['median_full_rate']}, "
          f"median steady {out[g]['median_steady_rate']} (n_steady {len(st)})")
(TALK / "data" / "warm_rate_r10.json").write_text(json.dumps(out, indent=1, default=float))
print("wrote data/warm_rate_r10.json")
