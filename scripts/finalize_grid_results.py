"""Update hypothesis_ledger.json and run_provenance.csv from completed grid runs.

Usage:
    uv run python scripts/finalize_grid_results.py [H18 H19 ...]  # specific cells
    uv run python scripts/finalize_grid_results.py                  # all grid_H* runs
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
RUNS_ROOT    = PROJECT_ROOT / "runs"
SCRATCH_BASE = Path("/oscar/scratch/eaguerov/mpi_runs")
LEDGER       = PROJECT_ROOT / "experiments" / "hypothesis_ledger.json"
PROVENANCE   = PROJECT_ROOT / "experiments" / "run_provenance.csv"

COMMIT  = "15f9f03"
BINARY  = "BioReactor-mpi-stripped"
DATE    = "2026-06-18"

# Metadata that isn't in params.json directly
CELL_META = {
    "grid_H18": {"g_angle": 7, "c_angle": 90, "sigma": 0.5, "ref_kla": None},
    "grid_H19": {"g_angle": 7, "c_angle": 90, "sigma": 0.5, "ref_kla": None},
    "grid_H20": {"g_angle": 7, "c_angle": 90, "sigma": 0.5, "ref_kla": None},
    "grid_H21": {"g_angle": 7, "c_angle": 90, "sigma": 0.5, "ref_kla": None},
    "grid_H22": {"g_angle": 5, "c_angle": 90, "sigma": 0.5, "ref_kla": None},
    "grid_H23": {"g_angle": 4, "c_angle": 90, "sigma": 0.5, "ref_kla": None},
    "grid_H24": {"g_angle": 7, "c_angle": 90, "sigma": 0.5, "ref_kla": None},
    "grid_H25": {"g_angle": 7, "c_angle": 90, "sigma": 0.5, "ref_kla": 0.139},
}

CKPT_SOURCE = {
    "grid_H18": "001e47bf/checkpoint.dump",
    "grid_H19": "6439dcba/checkpoint.dump",
    "grid_H20": "08832477/checkpoint.dump",
    "grid_H21": "3438372e/checkpoint.dump",
    "grid_H22": "3c940365/checkpoint.dump",
    "grid_H23": "29dc4215/checkpoint.dump",
    "grid_H24": "63062937/checkpoint.dump",
    "grid_H25": "1f8c7d4b/checkpoint.dump",
}


def check_run(run_id: str) -> dict:
    run_dir   = RUNS_ROOT / run_id
    scratch   = SCRATCH_BASE / run_id
    res_file  = run_dir / "results.json"
    oxy_file  = (run_dir / "tr_oxy.dat") if (run_dir / "tr_oxy.dat").exists() \
                else (scratch / "tr_oxy.dat")

    result = {"run_id": run_id, "pass": False, "kLa_25": "NA", "kLa_50": "NA",
              "oxy_last": "NA", "t_last": "NA", "notes": []}

    if res_file.exists():
        r = json.loads(res_file.read_text())
        k25 = r.get("kLa_25", None)
        k50 = r.get("kLa_50", None)
        if k25 is not None and math.isfinite(float(k25)):
            result["kLa_25"] = f"{float(k25):.4f}"
        if k50 is not None and math.isfinite(float(k50)):
            result["kLa_50"] = f"{float(k50):.4f}"
    else:
        result["notes"].append("results.json missing")

    if oxy_file.exists():
        lines = [l for l in oxy_file.read_text().splitlines() if l.strip()]
        if lines:
            parts = lines[-1].split()
            try:
                result["t_last"]   = parts[1]
                result["oxy_last"] = parts[2]
            except IndexError:
                result["notes"].append("tr_oxy parse error")

    # Pass: kLa_25 is finite AND oxy_last > 0
    try:
        kla_ok = math.isfinite(float(result["kLa_25"]))
    except (ValueError, TypeError):
        kla_ok = False
    try:
        oxy_ok = float(result["oxy_last"]) > 0
    except (ValueError, TypeError):
        oxy_ok = False

    result["pass"] = kla_ok and oxy_ok
    return result


def update_ledger(run_id: str, info: dict):
    entries = json.loads(LEDGER.read_text())
    h_label = run_id.replace("grid_", "")  # "H18", "H19", ...
    for e in entries:
        if e.get("hypothesis", "").startswith(h_label + ":") and e["status"] == "pending":
            ref = CELL_META.get(run_id, {}).get("ref_kla")
            if info["pass"]:
                obs = (f"CONFIRMED: kLa_25={info['kLa_25']}, oxy_liq_sum={info['oxy_last']} "
                       f"at t={info['t_last']}; zero crashes.")
                if ref is not None:
                    try:
                        ratio = abs(float(info["kLa_25"]) - ref) / ref
                        obs += f" Agreement with ref kLa_25={ref}: |Δ|/ref={ratio:.1%}."
                    except Exception:
                        pass
                e["observation"] = obs
                e["status"] = "confirmed"
            else:
                e["observation"] = (f"FAILED or INCOMPLETE: kLa_25={info['kLa_25']}, "
                                    f"oxy={info['oxy_last']} at t={info['t_last']}. "
                                    + "; ".join(info["notes"]))
                e["status"] = "falsified"
    LEDGER.write_text(json.dumps(entries, indent=2))


def append_provenance(run_id: str, info: dict):
    p_file = RUNS_ROOT / run_id / "params.json"
    if not p_file.exists():
        return
    p = json.loads(p_file.read_text())
    g = p.get("geometry", {})
    jid = p.get("_slurm_job_id", "see-ledger")
    row = (
        f"{run_id},"
        f"{jid},"
        f"{DATE},"
        f"{BINARY},"
        f"{COMMIT},"
        f"{p.get('fidelity','?')},"
        f"{p.get('_ntasks','?')},"
        f"{p.get('omega_b','?')},"
        f"{CELL_META.get(run_id,{}).get('g_angle','?')},"
        f"110,"
        f"{CELL_META.get(run_id,{}).get('sigma','?')},"
        f"{CKPT_SOURCE.get(run_id,'?')},"
        f"{p.get('t_checkpoint','?')},"
        f"{info['t_last']},"
        f"NA,"
        f"{info['kLa_25']},"
        f"{info['kLa_50']},"
        f"NA,"
        f"{'pass' if info['pass'] else 'fail'},"
        f"geom_n={g.get('n','?')} fill={p.get('fill_level','?')} "
        f"omega_prev={p.get('omega_b_prev',0):.4f}"
    )
    with open(PROVENANCE, "a") as f:
        f.write(row + "\n")


def main():
    targets = sys.argv[1:] if sys.argv[1:] else None
    if targets:
        run_ids = [f"grid_{h}" if not h.startswith("grid_") else h for h in targets]
    else:
        run_ids = [d.name for d in RUNS_ROOT.iterdir()
                   if d.is_dir() and d.name.startswith("grid_H")]
        run_ids.sort()

    print(f"Checking {len(run_ids)} runs...")
    all_pass = True
    for run_id in run_ids:
        info = check_run(run_id)
        status = "PASS" if info["pass"] else "FAIL/INCOMPLETE"
        print(f"  {run_id:20s}  {status:12s}  kLa_25={info['kLa_25']:8s}  "
              f"oxy={info['oxy_last']:10s}  t={info['t_last']}")
        if info["pass"]:
            update_ledger(run_id, info)
            append_provenance(run_id, info)
        else:
            all_pass = False

    print()
    print("Ledger and provenance updated for completed runs.")
    return 0 if all_pass else 1


if __name__ == "__main__":
    raise SystemExit(main())
