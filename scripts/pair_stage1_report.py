"""Stage-1 hydro of the L10 roughness pairs (pre-registered 2026-10-19/21): EDR (Kim statistic edr_max_t and the time mean edr_mean_t) at 22.4 / 22.6 rpm and 29.9 rpm, theta 7.
RULE: REGIME EDGE AT L10 if |ln EDR(22.4) - ln EDR(22.6)| >= 0.15. Usage: uv run python scripts/pair_stage1_report.py"""
import math
import sys

sys.path[:0] = ["."]
from scripts.hydro_dataset_v2 import qois  # noqa: E402

q = {r: qois(f"x8_l10_rpm{r}_th7") for r in ("22.4", "22.6", "29.9")}
for r, v in q.items():
    print(r, {k: round(v[k], 5) for k in ("edr_mean_t", "edr_max_t", "tau_95_t", "tau_amp_t")}, "window", v["window"])
for k in ("edr_mean_t", "edr_max_t", "tau_95_t", "tau_amp_t"):
    print(k, "d ln (22.4 vs 22.6) = %+.3f" % math.log(q["22.4"][k] / q["22.6"][k]))
