"""Submit apples-to-apples benchmark: serial+checkpoint vs MPI+independent.

Both arms produce the same two kLa data points (omega_b[0] and omega_b[1])
at fidelities 4, 5, and 6.

Serial arm  — 2 SLURM jobs per fidelity:
  seg-0: BioReactor-video (OpenMP 4-thread), fresh start, omega_b[0]
  seg-1: BioReactor-video, checkpoint restart with omega_b[1]  ← self-submitted

MPI arm  — 2 SLURM jobs per fidelity (submitted upfront, both fresh):
  run-a: BioReactor-mpi-video (16 ranks), fresh start, omega_b[0]
  run-b: BioReactor-mpi-video (16 ranks), fresh start, omega_b[1]

Timing is read from the final logstats.dat Wall clock time entry.
"""
import json, math, subprocess, re, sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
RUNS_ROOT    = PROJECT_ROOT / "runs"
LOGS_DIR     = PROJECT_ROOT / "logs"
SCRATCH_BASE = Path("/oscar/scratch/eaguerov/mpi_runs")
SERIAL_TPL   = PROJECT_ROOT / "config" / "slurm_template.sh"
MPI_TPL      = PROJECT_ROOT / "config" / "slurm_mpi_template.sh"

OMEGA_B      = [1.5708, 1.8326]   # two omega_b values per arm
FIDELITIES   = [4, 5, 6]

BASE = {
    "n_harmonics": 1,
    "theta_max":      [7.0, 0.0, 0.0],
    "phi_angular":    [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h":    [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
    "fill_level": 0.5,
    "n_mix_cycles": 80,
}

WALLTIME = {4: "00:15:00", 5: "00:30:00", 6: "02:00:00"}


def t_mix_nd(params):
    """Non-dim t_mix = T_per_st * n_mix_cycles."""
    ob = params["omega_b"]
    L, H = params["geometry"]["a"], params["geometry"]["b"]
    th   = math.radians(params["theta_max"][0])
    T    = 2 * math.pi / ob
    V    = L/4 * (H + 0.5*L*math.tan(th))
    U    = V / (H*0.5) / T
    Tb   = L / U
    Tst  = T / Tb
    return Tst * params["n_mix_cycles"]


def t_end_nd(params):
    """t_end = t_mix + T_per_st * n_mix_cycles (kLa window = same length as warmup)."""
    ob = params["omega_b"]
    L, H = params["geometry"]["a"], params["geometry"]["b"]
    th   = math.radians(params["theta_max"][0])
    T    = 2 * math.pi / ob
    V    = L/4 * (H + 0.5*L*math.tan(th))
    U    = V / (H*0.5) / T
    Tb   = L / U
    Tst  = T / Tb
    tm   = Tst * params["n_mix_cycles"]
    return round(tm + Tst * params["n_mix_cycles"] + 50.0, 3)


def t_end_restart(params, t_ck):
    """t_end for restart seg: t_checkpoint + T_per_st * n_mix_cycles + kLa window."""
    ob = params["omega_b"]
    L, H = params["geometry"]["a"], params["geometry"]["b"]
    th   = math.radians(params["theta_max"][0])
    T    = 2 * math.pi / ob
    V    = L/4 * (H + 0.5*L*math.tan(th))
    U    = V / (H*0.5) / T
    Tb   = L / U
    Tst  = T / Tb
    return round(t_ck + Tst * 10 + Tst * 80 + 50.0, 3)


def write_run(run_id, params):
    rd = RUNS_ROOT / run_id
    rd.mkdir(parents=True, exist_ok=True)
    params["_canonical_run_dir"] = str(rd.resolve())
    (rd / "params.json").write_text(json.dumps(params, indent=2))
    return rd


def stage_mpi(run_id, params):
    sd = SCRATCH_BASE / run_id
    sd.mkdir(parents=True, exist_ok=True)
    (sd / "params.json").write_text(json.dumps(params, indent=2))
    return sd


def sbatch(args):
    r = subprocess.run(["sbatch"] + args, capture_output=True, text=True, check=True)
    return re.search(r"(\d+)", r.stdout).group(1)


jobs = []   # (label, run_id, job_id)

for fi in FIDELITIES:
    wt = WALLTIME[fi]

    # ── Serial arm ─────────────────────────────────────────────────────────
    # seg-0: fresh start, omega_b[0]
    s0_id = f"bench_f{fi}_serial_s0"
    s0p   = dict(BASE, omega_b=OMEGA_B[0], fidelity=fi)
    s0p["t_end"]        = t_end_nd(s0p)
    s0p["next_run_id"]  = f"bench_f{fi}_serial_s1"
    s0p["_walltime"]    = wt
    s0p["_mem"]         = "4G"
    s0p["run_id"]       = s0_id
    s0_rd = write_run(s0_id, s0p)

    # seg-1 params (pre-written so self-submit can find it)
    s1_id = f"bench_f{fi}_serial_s1"
    s1p   = dict(BASE, omega_b=OMEGA_B[1], fidelity=fi,
                 n_mix_cycles=10)
    t_ck  = s0p["t_end"]
    s1p["t_checkpoint"]   = t_ck
    s1p["omega_b_prev"]   = OMEGA_B[0]
    s1p["t_end"]          = t_end_restart(s1p, t_ck)
    s1p["_walltime"]      = wt
    s1p["_mem"]           = "4G"
    s1p["run_id"]         = s1_id
    s1_rd = write_run(s1_id, s1p)

    jid = sbatch([
        "--no-requeue",
        f"--job-name=bench-f{fi}-serial",
        "--ntasks=1", "--cpus-per-task=4",
        f"--time={wt}", "--mem=4G",
        f"--export=NONE,PARAMS={s0_rd}/params.json",
        str(SERIAL_TPL.resolve()),
    ])
    jobs.append((f"f{fi}_serial_s0", s0_id, jid))
    print(f"[f{fi}] serial  seg-0  → job {jid}  (self-submits seg-1 on completion)")

    # ── MPI arm ─────────────────────────────────────────────────────────────
    for idx, ob in enumerate(OMEGA_B):
        m_id = f"bench_f{fi}_mpi_{chr(ord('a')+idx)}"
        mp   = dict(BASE, omega_b=ob, fidelity=fi)
        mp["t_end"]      = t_end_nd(mp)
        mp["_walltime"]  = wt
        mp["_mem"]       = "2G"
        mp["_ntasks"]    = 16
        mp["run_id"]     = m_id
        m_rd = write_run(m_id, mp)
        sd   = stage_mpi(m_id, mp)

        jid = sbatch([
            "--no-requeue",
            f"--job-name=bench-f{fi}-mpi-{chr(ord('a')+idx)}",
            "--ntasks=16", "--cpus-per-task=1",
            f"--time={wt}", "--mem-per-cpu=2G",
            f"--export=NONE,PARAMS={sd}/params.json",
            str(MPI_TPL.resolve()),
        ])
        jobs.append((f"f{fi}_mpi_{chr(ord('a')+idx)}", m_id, jid))
        print(f"[f{fi}] mpi    run-{chr(ord('a')+idx)}  → job {jid}  (omega_b={ob})")

print()
print("Summary:")
print(f"{'Label':<25} {'run_id':<30} {'job_id'}")
for label, rid, jid in jobs:
    print(f"  {label:<23} {rid:<30} {jid}")
print()
print("Monitor: squeue -u eaguerov")
print("When done, run: uv run python scripts/benchmark_report.py")
