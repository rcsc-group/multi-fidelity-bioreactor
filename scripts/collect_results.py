"""Aggregate results.json files from a sweep into a single CSV.

Matches the column layout of Kim et al. (2025) Tables / Figs. 9–13.

Usage
-----
    python scripts/collect_results.py [--sweep config/sweep_fb_theta.json] [--out results.csv]

If --sweep is omitted, all runs/ subdirectories with a results.json are collected.
If --out is omitted, the CSV is written to experiments/<sweep_name>_results.csv
(or to stdout if no sweep config is given).

Output columns
--------------
  run_id, omega_b, fill_level, theta_max_0, fidelity
  kLa_10, kLa_25, kLa_50                   — 5-pt log-linear fit
  kLa_inst_10, kLa_inst_25, kLa_inst_50    — instantaneous dC*/dt / (1-C*)
  dtmix_0.50, dtmix_0.75, dtmix_0.95       — dimensional mixing time (s)
  vor_mean                                  — period-averaged |ξ| (1/s)
  rpm                                       — rocking frequency (rpm, derived from omega_b)
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).parents[1]

_RESULT_KEYS = [
    "kLa_10", "kLa_25", "kLa_50",
    "kLa_inst_10", "kLa_inst_25", "kLa_inst_50",
    "dtmix_0.50", "dtmix_0.75", "dtmix_0.95",
    "vor_mean",
]

_FIELDNAMES = [
    "run_id", "rpm", "omega_b", "fill_level", "theta_max_0", "fidelity",
    *_RESULT_KEYS,
]


def _load_run(run_dir: Path) -> dict | None:
    params_path  = run_dir / "params.json"
    results_path = run_dir / "results.json"
    if not params_path.exists() or not results_path.exists():
        return None
    try:
        params  = json.loads(params_path.read_text())
        results = json.loads(results_path.read_text())
    except (json.JSONDecodeError, OSError):
        return None

    omega_b = params.get("omega_b", float("nan"))
    row = {
        "run_id":      params.get("run_id", run_dir.name),
        "rpm":         round(omega_b * 60 / (2 * math.pi), 4),
        "omega_b":     omega_b,
        "fill_level":  params.get("fill_level", float("nan")),
        "theta_max_0": (params.get("theta_max") or [float("nan")])[0],
        "fidelity":    params.get("fidelity", float("nan")),
    }
    for key in _RESULT_KEYS:
        val = results.get(key, float("nan"))
        row[key] = val if val is not None else float("nan")
    return row


def collect(run_dirs: list[Path]) -> list[dict]:
    rows = []
    for d in sorted(run_dirs):
        row = _load_run(d)
        if row is not None:
            rows.append(row)
    # sort by (theta_max_0, fill_level, rpm) for readability
    rows.sort(key=lambda r: (r["theta_max_0"], r["fill_level"], r["rpm"]))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sweep", metavar="JSON",
                        help="Sweep config JSON (used to filter by run_id prefix; optional)")
    parser.add_argument("--out", metavar="CSV",
                        help="Output CSV path (default: stdout or experiments/<name>_results.csv)")
    parser.add_argument("--runs-root", metavar="DIR", default=str(_PROJECT_ROOT / "runs"),
                        help="Root directory containing run subdirs (default: runs/)")
    args = parser.parse_args()

    runs_root = Path(args.runs_root)
    run_dirs  = [d for d in runs_root.iterdir() if d.is_dir()]

    # If a sweep config is given, optionally filter by known run_ids
    # (Currently collects all — filtering by run_id prefix is future work)
    sweep_name = None
    if args.sweep:
        sweep_name = Path(args.sweep).stem

    rows = collect(run_dirs)
    if not rows:
        print("No completed runs found.", file=sys.stderr)
        sys.exit(1)

    # Determine output
    if args.out:
        out_path = Path(args.out)
    elif sweep_name:
        exp_dir = _PROJECT_ROOT / "experiments"
        exp_dir.mkdir(parents=True, exist_ok=True)
        out_path = exp_dir / f"{sweep_name}_results.csv"
    else:
        out_path = None  # stdout

    def _fmt(v):
        if isinstance(v, float) and math.isnan(v):
            return ""
        return str(v)

    if out_path:
        with open(out_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=_FIELDNAMES)
            w.writeheader()
            for row in rows:
                w.writerow({k: _fmt(row.get(k, "")) for k in _FIELDNAMES})
        print(f"Wrote {len(rows)} rows to {out_path}")
    else:
        w = csv.DictWriter(sys.stdout, fieldnames=_FIELDNAMES)
        w.writeheader()
        for row in rows:
            w.writerow({k: _fmt(row.get(k, "")) for k in _FIELDNAMES})


if __name__ == "__main__":
    main()
