"""Round 10 A2: warm-L9 cost estimate for the 10 fig9 theta-7 points.
est = median full core-h/cycle of the six measured x8_l9 warm runs (logstats wall clock x ncpu)
      x median simulated cycles of the same six runs (warm_rate_r10.json, group warm_L9_x8).
Output: data/warm_l9_r10.json (same keys as warm_l9_r9.json step2, for cost_rule).
"""
import json
from pathlib import Path
TALK = Path("/oscar/data/dharri15/eaguerov/talks/dfd2026")
W = json.load(open(TALK / "data" / "warm_rate_r10.json"))["warm_L9_x8"]
rate = W["median_full_rate"]
cyc = W["median_cycles"]
est = rate * cyc
r9 = json.load(open(TALK / "data" / "warm_l9_r9.json"))["step2"]["runs"]
runs = {k: dict(v, warm_est=round(est, 1), cycles_used=round(cyc, 2), rate_used=round(rate, 3)) for k, v in r9.items()}
out = dict(rate_core_h_per_cycle=rate, cycles=cyc, estimate_per_run=est, sum_est=est * len(runs), step2=dict(runs=runs, sum_est=est * len(runs)))
(TALK / "data" / "warm_l9_r10.json").write_text(json.dumps(out, indent=1))
print(f"rate {rate:.3f} core-h/cycle x cycles {cyc:.2f} = {est:.1f} core-h per fig9 point; sum over {len(runs)} = {est*len(runs):.1f}")
