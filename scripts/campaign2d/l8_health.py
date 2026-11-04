"""Health table for the campaign 2D L8 jobs (read-only). One row per job in state/jobs.jsonl:
sacct state, last logstats time t, cycles done (t / period), progress to t_end, NaN lines, error lines in the .err.
Run inside a compute job: uv run python scripts/campaign2d/l8_health.py
"""
import json
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import postprocess as pp  # noqa: E402

JOBS = ROOT / "scripts" / "campaign2d" / "state" / "jobs.jsonl"
SCRATCH = Path("/oscar/scratch/eaguerov/mpi_runs")
BAD = ("error", "abort", "segmentation", "fatal", "nan")
BENIGN = ("mca_", "libcuda")


def state_of(jobs):
    out = subprocess.run(["sacct", "-X", "-n", "-P", "-j", ",".join(jobs), "-o", "JobID,State,Elapsed"],
                         capture_output=True, text=True).stdout
    res = {}
    for line in out.splitlines():
        jid, st, el = (line.split("|") + ["", "", ""])[:3]
        res[jid] = (st, el)
    return res


def main():
    rows = [json.loads(x) for x in JOBS.read_text().splitlines()]
    st = state_of([r["job"] for r in rows])
    print(f"{'job':>8} {'run':30s} {'state':10s} {'elapsed':>10s} {'cycles':>7s} {'of':>4s} {'t':>7s} {'nan':>4s} {'errs':>5s}")
    for r in rows:
        run = r["run"]
        params = json.loads((ROOT / "runs" / run / "params.json").read_text())
        _, T = pp._t_scales(params)
        n_total = params["n_mix_cycles"] + 150.0
        state, el = st.get(r["job"], ("?", ""))
        logf = SCRATCH / run / "logstats.dat"
        t_last, nan, errs = float("nan"), 0, 0
        if logf.exists():
            lines = [ln for ln in logf.read_text(errors="replace").splitlines() if ln.startswith("i:")]
            if lines:
                t_last = float(lines[-1].split()[3])
            nan = sum("nan" in ln.lower() for ln in lines)
        errf = ROOT / "logs" / f"slurm_{r['job']}.err"
        if errf.exists():
            for ln in errf.read_text(errors="replace").splitlines():
                low = ln.lower()
                if any(b in low for b in BENIGN):
                    continue
                if any(b in low for b in BAD):
                    errs += 1
        cyc = t_last / T if not math.isnan(t_last) else float("nan")
        print(f"{r['job']:>8} {run:30s} {state:10s} {el:>10s} {cyc:7.1f} {n_total:4.0f} {t_last:7.2f} {nan:4d} {errs:5d}")


if __name__ == "__main__":
    main()
