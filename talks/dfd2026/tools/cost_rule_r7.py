"""Round 7 accounting, ONE rule for every method.

Rule: a method is charged the sacct core-h (corehc_all, cost_table.tsv 2026-11-04 09:22) of every run id it needs, each
run id counted once. A warm-started run is charged its own cost plus its source run's cost, unless the source is already
charged to that method (then only its own cost).
  L8 run          : own cost (cold).
  L9 warm (x8_l9) : own cost; source = L8 run at the same (rpm, theta).
  L9 cold (fig9)  : own cost (cold, no source).
  L10 x8_l10      : own; source = L8 run at the same (rpm, theta).
  L10 xlevel      : own; source = fig9_l9_rpm<r> (cold L9).
Method sets:
  chain (full)  : 40 L8 + 16 L9 (incl. 10 fig9) ; each L10 run: own + source only if source not in set.
  chain (warm9*): L8 subset (incl. the 6 sources of x8_l9) + 6 warm L9 ; L10 run: own + source if missing (xlevel -> fig9).
  L10 GP (ours) : no L8/L9 set; each L10 run: own + source.
  H1            : 40 L8 ; each L10 run: own + source if missing.
  Kim code      : own 120 core-h per run (user's assumption), no sources.
Output: data/cost_rule_r7.json and printed table.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path("/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor")
sys.path[:0] = [str(REPO), str(REPO / "scripts")]
import scripts.test_h1_mfbml as H1  # noqa: E402
import scripts.test_h2_mfbml as H2  # noqa: E402

TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
COST = {}
with open("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        COST[r["run"]] = float(r["corehc_all"])
KIM = 120.0

hf = H1.hf_points()
hname = [H1.hf_run(r, t) for r, t in hf]
src = {}
for r, t in hf:
    src[H1.hf_run(r, t)] = f"fig9_l9_rpm{r:g}" if t == 7 else H1.lf_run(r, t)
# source names are checked against the cost table
for k, v in src.items():
    assert k in COST and v in COST, (k, v)

l8_40 = [H1.lf_run(r, t) for r in H1.RPMS for t in H1.LF_THETAS]
l8_40 = [n for n in l8_40 if n in COST]
l9_16 = [H2.l9_run(r, t) for r, t in H2.l9_points()]
warm6 = [H2.l9_run(r, t) for r, t in H2.l9_points() if t != 7]
srcs6 = [H1.lf_run(float(n.split("rpm")[1].split("_th")[0]), float(n.split("_th")[1])) for n in warm6]
assert all(s in COST for s in srcs6), srcs6


def inc_rule(set_names):
    """Mean over the 13 L10 runs of own + (source if the source is not in set_names)."""
    vals = []
    for n in hname:
        own = COST[n]
        s = src[n]
        vals.append(own + (0.0 if s in set_names else COST[s]))
    return float(np.mean(vals)), vals


def base_of(names):
    return float(sum(COST[n] for n in names))


out = {}
# chain, full
set_full = set(l8_40) | set(l9_16)
out["chain_full"] = dict(base=base_of(set_full), inc=inc_rule(set_full)[0], n_L8=len(l8_40), n_L9=len(l9_16))
# H1
out["H1"] = dict(base=base_of(l8_40), inc=inc_rule(set(l8_40))[0], n_L8=len(l8_40), n_L9=0)
# L10 GP (ours): base 0, every L10 run with its source
out["L10_GP_ours"] = dict(base=0.0, inc=inc_rule(set())[0], n_L8=0, n_L9=0)
# Kim code
out["Kim_code"] = dict(base=0.0, inc=KIM, n_L8=0, n_L9=0, note="assumed 120 core-h per run, no sources")

# cheap chain variants: keys from the run-time variant definition in tools/chain_variants_lc.py
variants = json.load(open(TALK / "data" / "chain_variants_lc.json"))
for name, v in variants.items():
    k8 = [H1.lf_run(r, t) for r, t in v["keys8"]]
    k9 = [H2.l9_run(r, t) for r, t in v["keys9"]]
    s8 = set(k8)
    s_all = s8 | set(k9)
    base = base_of(s_all)
    inc, _ = inc_rule(s_all)
    out[f"chain_{name}"] = dict(base=base, inc=inc, n_L8=len(k8), n_L9=len(k9))

print(f"{'method':24s} {'n_L8':>4s} {'n_L9':>4s} {'base core-h':>12s} {'per L10 run':>12s}")
for k, v in out.items():
    print(f"{k:24s} {v['n_L8']:4d} {v['n_L9']:4d} {v['base']:12.1f} {v['inc']:12.1f}")
print("per-L10 charges (own + source if missing) under chain_warm9_L8_10:")
_, vals = inc_rule(set(H1.lf_run(r, t) for r, t in variants["warm9_L8_10"]["keys8"]) | set(H2.l9_run(r, t) for r, t in variants["warm9_L8_10"]["keys9"]))
for n, vv in zip(hname, vals):
    print(f"   {n:26s} src {src[n]:24s} charged {vv:9.1f}")
(TALK / "data" / "cost_rule_r7.json").write_text(json.dumps(out, indent=1))
print("wrote data/cost_rule_r7.json")
