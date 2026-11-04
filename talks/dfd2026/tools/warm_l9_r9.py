"""Round 9 addendum: warm-L9 cost estimates and warm/cold QoI evidence.

Step 1  core-h per simulated cycle, per row of cost_table.tsv (2026-11-04 09:22):
        corehc_all / ((tlast - tfirst) / T_cycle), T_cycle = _t_scales(params of that run) (per run, depends on theta).
        Groups: fid 9 warm (cold=False), fid 9 cold (cold=True), fid 10 warm (cold=False). Median, IQR, n.
Step 2  warm-L9 cost estimate for each fig9_l9 theta-7 point:
        est = warm_median_per_cycle x (tail_cycles + WARMUP), with
        tail_cycles = post-release cycles of the cold run = t_end / T - n_mix_cycles (release at cycle n_mix_cycles = 80),
        WARMUP = 12 (warm start from L8, diary 2026-10-21 protocol). tail-only value also reported.
Step 3  paired evidence, warm vs cold at the same point (theta 7 unless noted):
        w9_l9_rpm<r>_th7 vs fig9_l9_rpm<r>; l9_warm_rpm<r> vs l9_cold_rpm<r>.
        QoIs: tau_95_t and edr_max_t (hydro_dataset_v2.qois), kLa_10 (results.json).
Output: data/warm_l9_r9.json
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

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
WARMUP = 12
rows = []
with open("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rows.append(r)


def params(run):
    p = REPO / "runs" / run / "params.json"
    return json.load(open(p)) if p.exists() else None


def per_cycle(r):
    if r["tlast"] in ("None", "") or r["tfirst"] in ("None", ""):
        return None
    p = params(r["run"])
    if p is None:
        return None
    _, T = _t_scales(p)
    cyc = (float(r["tlast"]) - float(r["tfirst"])) / T
    if cyc <= 0 or float(r["tfirst"]) < 0:
        return None
    return float(r["corehc_all"]) / cyc, cyc, T


def stats(vals):
    v = np.array(vals)
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    return dict(n=int(len(v)), median=float(med), q1=float(q1), q3=float(q3), iqr=float(q3 - q1))


groups = {"fid9_warm": [], "fid9_cold": [], "fid10_warm": []}
for r in rows:
    f = int(r["fid"])
    cold = r["cold"] == "True"
    key = None
    if f == 9 and not cold:
        key = "fid9_warm"
    elif f == 9 and cold:
        key = "fid9_cold"
    elif f == 10 and not cold:
        key = "fid10_warm"
    if key is None:
        continue
    pc = per_cycle(r)
    if pc is not None:
        groups[key].append((r["run"], pc[0], pc[1], pc[2]))
out = {"T_note": "T per run from scripts.postprocess._t_scales(params.json)", "step1": {}}
for k, v in groups.items():
    s = stats([x[1] for x in v])
    out["step1"][k] = s
    print(f"{k}: n={s['n']} median {s['median']:.2f} core-h/cycle, IQR {s['q1']:.2f}-{s['q3']:.2f} (IQR {s['iqr']:.2f})")
warm_med = out["step1"]["fid9_warm"]["median"]

# ---- step 2
fig9 = [r for r in rows if r["run"].startswith("fig9_l9_rpm") and "_ext" not in r["run"] and r["run"].count("_") <= 2]
est = {}
for r in fig9:
    p = params(r["run"])
    _, T = _t_scales(p)
    t_end = float(p["t_end"])
    n_mix = int(p["n_mix_cycles"])
    tail = t_end / T - n_mix
    est[r["run"]] = dict(t_end=t_end, T=T, n_mix=n_mix, tail_cycles=round(tail, 2),
                         cold_cost=float(r["corehc_all"]),
                         warm_est_tail_only=round(warm_med * tail, 1),
                         warm_est=round(warm_med * (tail + WARMUP), 1))
    print(r["run"], est[r["run"]])
out["step2"] = dict(warm_median_per_cycle=warm_med, warmup_cycles=WARMUP, runs=est)
vals = [e["warm_est"] for e in est.values()]
print(f"warm-L9 estimate over {len(vals)} fig9 runs: median {np.median(vals):.1f}, sum {np.sum(vals):.1f} core-h")
out["step2"]["sum_est"] = float(np.sum(vals))

# ---- step 3 (paired QoIs)
def qoi(run):
    d = {}
    q = qois(run)
    d["tau_95_t"] = q["tau_95_t"]
    d["edr_max_t"] = q["edr_max_t"]
    rj = REPO / "runs" / run / "results.json"
    d["kLa_10"] = json.load(open(rj)).get("kLa_10") if rj.exists() else None
    d["kLa_1T_10"] = json.load(open(rj)).get("kLa_1T_10") if rj.exists() else None
    return d


pairs = []
for rpm in ("17.5", "25", "32.5"):
    pairs.append((f"w9_l9_rpm{rpm}_th7", f"fig9_l9_rpm{rpm}"))
for rpm in ("22.5", "30", "37.5"):
    pairs.append((f"l9_warm_rpm{rpm}", f"l9_cold_rpm{rpm}"))
out["step3"] = []
for w, c in pairs:
    try:
        qw, qc = qoi(w), qoi(c)
    except Exception as e:  # missing shear_stress or results
        print("skip", w, c, e)
        continue
    rel = {k: (qw[k] / qc[k] - 1) if (qw[k] is not None and qc[k] not in (None, 0)) else None for k in ("tau_95_t", "edr_max_t", "kLa_10")}
    ln_kla = float(np.log(qw["kLa_10"] / qc["kLa_10"])) if (qw["kLa_10"] and qc["kLa_10"]) else None
    out["step3"].append(dict(warm=w, cold=c, rel_diff=rel, ln_kla_10=ln_kla, warm_vals=qw, cold_vals=qc))
    print(f"{w} vs {c}: tau95 {100*rel['tau_95_t']:+.1f} %  EDRmax {100*rel['edr_max_t']:+.1f} %  kLa10 {100*rel['kLa_10']:+.1f} % (ln {ln_kla:+.3f})"
          if rel["kLa_10"] is not None else f"{w} vs {c}: kLa10 missing")
(TALK / "data" / "warm_l9_r9.json").write_text(json.dumps(out, indent=1, default=float))
print("wrote data/warm_l9_r9.json")
