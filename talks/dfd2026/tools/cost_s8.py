"""Actual training cost (core-h) of each method's training runs, from the sacct audit, per run id.

Per-run cost = corehc_all in /oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv (2026-11-04 09:22). That table is made by
cost_table.py: for each run id, sum over ALL its slurm jobs of Elapsed x AllocCPUS (completed or not).
Run-id sets come from the code that the learning curve uses (scripts/test_h1_mfbml.py, scripts/test_h2_mfbml.py):
  L8 set  = lf_run(r, t) for r in RPMS, t in LF_THETAS            (40 runs, cold from rest)
  L9 set  = l9_points() -> fig9_l9_rpm<r> (theta 7, cold) and x8_l9_rpm<r>_th<t> (theta 4, 9 at 17.5/25/32.5, warm from L8)
  L10 set = hf_points() -> l10_fig13a_xlevel_rpm<r> (theta 7, warm from fig9_l9_rpm<r>) and x8_l10_rpm<r>_th<t> (warm from L8)
Warm-start sources (rule from scripts/submit_l10_from_l8.py source_run and submit_l10_fig13a_xlevel.py):
  x8_l9_rpm<r>_th<t>   <- lf_run(r, t)                  (ASSUMED: same rpm/theta L8 checkpoint; not verified in code)
  x8_l10_rpm<r>_th<t>  <- lf_run(r, t)
  l10_fig13a_xlevel_rpm<r> <- fig9_l9_rpm<r>
Each method's training cost = union of run ids it needs, each counted once.
Writes data/costs.json. Data mtimes: cost_table.tsv 2026-11-04 09:22; repo test_h1/test_h2 modules (2026-10-xx).
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

COST = {}
with open("/oscar/scratch/eaguerov/tmp/cost_audit/cost_table.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        COST[r["run"]] = (float(r["corehc_all"]), r["cold"] == "True")


def c(name):
    if name not in COST:
        raise KeyError(f"no sacct row for {name}")
    return COST[name][0]


def src_l9(name):
    # the warm L9 point at (rpm, theta) starts from the L8 run at the same (rpm, theta)
    return None


L8 = [H1.lf_run(r, t) for r in H1.RPMS for t in H1.LF_THETAS]
L8 = [n for n in L8 if n in COST]
L9 = [H2.l9_run(r, t) for (r, t) in H2.l9_points()]
L10 = [H1.hf_run(r, t) for (r, t) in H1.hf_points()]
print(f"sets: L8 {len(L8)}, L9 {len(L9)}, L10 {len(L10)}")
missing = [n for n in L9 + L10 if n not in COST]
print("missing from cost table:", missing)

# warm-start sources
src = {}
for r, t in H2.l9_points():
    if t != 7:
        src[H2.l9_run(r, t)] = H1.lf_run(r, t)
for r, t in H1.hf_points():
    if t == 7:
        src[H1.hf_run(r, t)] = f"fig9_l9_rpm{r:g}"
    else:
        src[H1.hf_run(r, t)] = H1.lf_run(r, t)
print("sources (run <- source):")
for k in sorted(src):
    print(f"  {k} <- {src[k]}  [cold flag of run: {COST.get(k, (0, None))[1]}]")

# per-level medians inside the chain
cost_L8 = [c(n) for n in L8]
cost_L9 = [c(n) for n in L9]
cost_L10 = [c(n) for n in L10]
cold_L9 = [n for n in L9 if COST[n][1]]
warm_L9 = [n for n in L9 if not COST[n][1]]
print(f"L8 chain runs: median {np.median(cost_L8):.1f} (n={len(cost_L8)}), range {min(cost_L8):.1f}-{max(cost_L8):.1f}")
print(f"L9 chain runs: median {np.median(cost_L9):.1f} (n={len(cost_L9)}); cold {len(cold_L9)}, warm {len(warm_L9)}; "
      f"median warm {np.median([c(n) for n in warm_L9]):.1f}, cold {np.median([c(n) for n in cold_L9]):.1f}")
print(f"L10 chain runs: median {np.median(cost_L10):.1f} (n={len(cost_L10)}), range {min(cost_L10):.1f}-{max(cost_L10):.1f}")

# method base sets (union, each run counted once)
base_chain = set(L8) | set(L9)
base_ratio = set(L9) | {src[n] for n in L9 if n in src}
base_ratio = {s for s in base_ratio if s in COST}
base_ratio |= set(L9)
base_l10 = set()

def tot(s):
    return sum(c(n) for n in s)

out = {}
# chain: L10 runs add own cost (their sources are in base_chain)
inc_chain = np.mean([c(n) + (0 if (src.get(n) in base_chain or src.get(n) is None) else c(src[n])) for n in L10])
# ratio: L10 runs add own cost (sources already in base_ratio for these sets)
inc_ratio = np.mean([c(n) + (0 if (src.get(n) in base_ratio or src.get(n) is None) else c(src[n])) for n in L10])
# L10 only: each run brings its own source
inc_l10 = np.mean([c(n) + (0 if src.get(n) is None else c(src[n])) for n in L10])
incl = lambda n: c(n) + (c(src[n]) if (n in src and src[n] in COST) else 0.0)
lvl = dict(
    L8_own=float(np.median(cost_L8)), L8_n=len(cost_L8),
    L9_own=float(np.median(cost_L9)), L9_incl=float(np.median([incl(n) for n in L9])),
    L9_cold_incl=float(np.median([incl(n) for n in cold_L9])), L9_warm_incl=float(np.median([incl(n) for n in warm_L9])),
    L10_own=float(np.median(cost_L10)), L10_incl=float(np.median([incl(n) for n in L10])),
)
print("per-level medians (core-h):", {k: round(v, 1) for k, v in lvl.items()})
base_h1 = tot(set(L8))
inc_h1 = float(np.mean([c(n) + (0 if (src.get(n) in set(L8) or src.get(n) is None) else c(src[n])) for n in L10]))
print(f"H1 (L8 + n L10): base {base_h1:.1f}, inc per L10 run {inc_h1:.1f}")
out = dict(
    base_chain=tot(base_chain), inc_chain=float(inc_chain),
    base_ratio=tot(base_ratio), inc_ratio=float(inc_ratio),
    base_l10=0.0, inc_l10=float(inc_l10), base_h1=base_h1, inc_h1=inc_h1,
    n_base_chain=len(base_chain), n_base_ratio=len(base_ratio),
    median_L8=float(np.median(cost_L8)), median_L9=float(np.median(cost_L9)), median_L10=float(np.median(cost_L10)),
    median_L9_warm=float(np.median([c(n) for n in warm_L9])), median_L9_cold=float(np.median([c(n) for n in cold_L9])),
    L9_cold_runs=cold_L9, L9_warm_runs=warm_L9,
    levels=lvl, per_L10_own={n: c(n) for n in L10}, per_L10_source={n: (src.get(n), c(src[n]) if src.get(n) in COST else None) for n in L10},
)
for k, v in out.items():
    if not isinstance(v, (list, dict)):
        print(f"{k}: {v:.1f}" if isinstance(v, float) else f"{k}: {v}")
Path("/oscar/data/dharri15/eaguerov/talks/dfd2026/data").mkdir(parents=True, exist_ok=True)
json.dump(out, open("/oscar/data/dharri15/eaguerov/talks/dfd2026/data/costs.json", "w"), indent=1)
print("wrote data/costs.json")
