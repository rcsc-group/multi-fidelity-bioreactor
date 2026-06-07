"""2-segment L5 MPI smoke test.

Exercises the full production stack end-to-end:
  seg-0: fresh start → videos → checkpoint.dump → rsync → postprocess → self-submit
  seg-1: checkpoint restart → videos past t_mix → rsync → postprocess → results.json

Pass criteria (check after both segments complete):
  - runs/smoke_l5_seg0/results.json exists with finite kLa_25
  - runs/smoke_l5_seg1/results.json exists with finite kLa_25
  - frames/ directories in both scratch runs contain >1 frame
"""
import json, subprocess, re
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
MPI_TEMPLATE  = PROJECT_ROOT / "config" / "slurm_mpi_template.sh"
SCRATCH_BASE  = Path("/oscar/scratch/eaguerov/mpi_runs")
RUNS_ROOT     = PROJECT_ROOT / "runs"

def stage_and_submit(params: dict, ntasks: int = 16) -> str:
    run_id  = params["run_id"]
    run_dir = RUNS_ROOT / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    params["_canonical_run_dir"] = str(run_dir.resolve())
    (run_dir / "params.json").write_text(json.dumps(params, indent=2))

    scratch_dir = SCRATCH_BASE / run_id
    scratch_dir.mkdir(parents=True, exist_ok=True)
    (scratch_dir / "params.json").write_text(json.dumps(params, indent=2))

    cmd = [
        "sbatch", "--no-requeue",
        f"--job-name=smoke-l5",
        f"--ntasks={ntasks}",
        "--time=00:30:00",
        "--mem-per-cpu=2G",
        "--cpus-per-task=1",
        f"--export=NONE,PARAMS={scratch_dir}/params.json",
        str(MPI_TEMPLATE.resolve()),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    job_id = re.search(r"(\d+)", result.stdout).group(1)
    print(f"  Submitted {run_id} → job {job_id}")
    return job_id


# Canonical params matching production sweep (L5, omega_b sweep first two values)
BASE = {
    "fidelity": 5,
    "omega_b": 1.5708,         # seg-0: first omega_b
    "n_mix_cycles": 80,
    "n_harmonics": 1,
    "theta_max":      [7.0, 0.0, 0.0],
    "phi_angular":    [0.0, 0.0, 0.0],
    "omega_h": 0.0,
    "amplitude_h":    [0.0, 0.0, 0.0],
    "phi_horizontal": [0.0, 0.0, 0.0],
    "geometry": {"a": 0.25, "b": 0.071, "n": 8.0},
    "fill_level": 0.5,
    "_walltime": "00:30:00",
    "_mem": "2G",
    "_ntasks": 16,
    "_experiment_dir": str(PROJECT_ROOT / "experiments" / "smoke_l5"),
}

# Compute t_end for seg-0 (same formula as sweep.py)
import math
omega_b = BASE["omega_b"]
L, H    = BASE["geometry"]["a"], BASE["geometry"]["b"]
th      = math.radians(BASE["theta_max"][0])
T_per   = 2 * math.pi / omega_b
V       = L / 4 * (H + 0.5 * L * math.tan(th))
U       = V / (H * 0.5) / T_per
T_bio   = L / U
T_per_st = T_per / T_bio
t_mix_seg0 = T_per_st * 80
t_end_seg0 = t_mix_seg0 + T_per_st * 80   # 80 cycles warmup + 80 cycles kLa

# seg-1 restarts from seg-0 checkpoint; transition cycles only
t_end_seg1 = t_end_seg0 + T_per_st * 10 + 50.0   # 10 transition + buffer

seg0 = dict(BASE,
    run_id="smoke_l5_seg0",
    t_end=round(t_end_seg0, 3),
    next_run_id="smoke_l5_seg1",
)
seg1 = dict(BASE,
    run_id="smoke_l5_seg1",
    omega_b=1.8326,            # seg-1: second omega_b (restarts from seg-0 checkpoint)
    n_mix_cycles=10,           # transition only
    t_checkpoint=round(t_end_seg0, 3),
    t_end=round(t_end_seg1, 3),
)

print("Submitting L5 MPI smoke test (2 segments)...")
print(f"  seg-0: omega_b=1.5708  t_end={seg0['t_end']:.1f}  (fresh start)")
print(f"  seg-1: omega_b=1.8326  t_end={seg1['t_end']:.1f}  (checkpoint restart)")
print(f"  seg-1 will be self-submitted by seg-0 on completion\n")

# Write seg-1 params to runs dir now (self-submit reads it from there)
seg1_run_dir = RUNS_ROOT / "smoke_l5_seg1"
seg1_run_dir.mkdir(parents=True, exist_ok=True)
seg1_canon = dict(seg1, _canonical_run_dir=str(seg1_run_dir.resolve()))
(seg1_run_dir / "params.json").write_text(json.dumps(seg1_canon, indent=2))

j0 = stage_and_submit(seg0)
print(f"\nSeg-0 job: {j0}")
print(f"Seg-1 will be self-submitted by the SLURM script after seg-0 completes.")
print(f"\nMonitor:")
print(f"  squeue -u eaguerov")
print(f"\nPass criteria:")
print(f"  ls runs/smoke_l5_seg0/results.json runs/smoke_l5_seg1/results.json")
print(f"  python3 -c \"import json; print(json.load(open('runs/smoke_l5_seg0/results.json'))['kLa_25'])\"")
