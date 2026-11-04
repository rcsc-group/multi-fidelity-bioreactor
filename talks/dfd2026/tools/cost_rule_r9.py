"""Round 9 addendum: one accounting rule, with every chain level warm-started.

Costs (core-h per run id):
  measured sacct corehc_all (cost_table.tsv, 2026-11-04 09:22), EXCEPT the 10 cold fig9_l9 theta-7 runs, which are
  replaced by the warm-L9 estimate (tools/warm_l9_r9.py step 2: warm median core-h/cycle x (tail cycles + 12)).
Rule (as round 7): a method is charged each run id it needs once; a warm run is charged its own cost plus its source's
cost unless the source is already in the method's set.
Methods: chain_full (40 L8, 16 L9 warm, 13 L10); chain variants (warm9 subsets, 6 warm L9 only);
H1 (40 L8, L10 runs); L10 GP ours (L10 runs with sources); Kim code (120 per run, assumed, no sources).
Output: data/cost_rule_r9.json
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
W9 = json.load(open(TALK / "data" / "warm_l9_r9.json"))["step2"]["runs"]
for k, v in W9.items():
    COST_W = None
# cost used for each fig9 theta-7 run: warm estimate
CW = dict(COST)
for k, v in W9.items():
    CW[k] = v["warm_est"]
KIM = 120.0
hf = H1.hf_points()
l10 = [H1.hf_run(r, t) for r, t in hf]
src = {}
for r, t in hf:
    src[H1.hf_run(r, t)] = f"fig9_l9_rpm{r:g}" if t == 7 else H1.lf_run(r, t)
for k, v in src.items():
    assert k in CW and v in CW, (k, v)
l8_40 = [H1.lf_run(r, t) for r in H1.RPMS for t in H1.LF_THETAS]
l8_40 = [n for n in l8_40 if n in COST]
l9_16 = [H2.l9_run(r, t) for r, t in H2.l9_points()]
warm6 = [n for n in l9_16 if not n.startswith("fig9")]
src9 = {}
for n in warm6:
    rr, tt = n.split("rpm")[1].split("_th")
    src9[n] = H1.lf_run(float(rr), float(tt))


def base(names):
    return float(sum(CW[n] for n in names))


def inc(set_names):
    vals = []
    for n in l10:
        s = src[n]
        vals.append(CW[n] + (0.0 if s in set_names else CW[s]))
    return float(np.mean(vals)), vals


out = {}
set_full = set(l8_40) | set(l9_16)
out["chain_full"] = dict(n_L8=len(l8_40), n_L9=len(l9_16), base=base(set_full), inc=inc(set_full)[0])
out["H1"] = dict(n_L8=len(l8_40), n_L9=0, base=base(set(l8_40)), inc=inc(set(l8_40))[0])
out["L10_GP_ours"] = dict(n_L8=0, n_L9=0, base=0.0, inc=inc(set())[0])
out["Kim_code"] = dict(n_L8=0, n_L9=0, base=0.0, inc=KIM)
for name, k8, k9 in [
    ("warm9", l8_40, warm6),
    ("warm9_L8_20", None, warm6),
    ("warm9_L8_10", None, warm6),
]:
    v = json.load(open(TALK / "data" / "chain_variants_lc.json"))[name.replace("warm9_L8_20", "warm9_L8_20")]
    kk8 = [H1.lf_run(r, t) for r, t in v["keys8"]]
    kk9 = [H2.l9_run(r, t) for r, t in v["keys9"]]
    s_all = set(kk8) | set(kk9)
    out[f"chain_{name}"] = dict(n_L8=len(kk8), n_L9=len(kk9), base=base(s_all), inc=inc(s_all)[0])

print(f"{'method':22s} {'n_L8':>4s} {'n_L9':>4s} {'base core-h':>12s} {'per L10 run':>12s} {'per L10 (L10-eq)':>16s}")
unit = out["chain_full"]["inc"]
for k, v in out.items():
    print(f"{k:22s} {v['n_L8']:4d} {v['n_L9']:4d} {v['base']:12.1f} {v['inc']:12.1f} {v['inc'] / unit:16.3f}")
print("L10 unit (own mean per L10 run, same code) = %.1f core-h" % unit)
out["L10_unit_core_h"] = unit
(TALK / "data" / "cost_rule_r9.json").write_text(json.dumps(out, indent=1))
print("wrote data/cost_rule_r9.json")
