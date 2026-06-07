"""Falsifiable experiment: same-omega_b MPI restart (su=1, no velocity rescaling).

Hypothesis A: rescaling (su != 1) is the bug. Prediction: this run passes.
Hypothesis B: MPI checkpoint restore itself is the bug. Prediction: this run crashes.

Restarts from smoke_l5_seg0 checkpoint (omega_b=1.5708) with the SAME omega_b.
t_checkpoint = seg0's t_end ≈ 97.3, n_mix_cycles=10, t_end ≈ 103.4.
"""
import json, subprocess, re, math
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
MPI_TEMPLATE  = PROJECT_ROOT / "config" / "slurm_mpi_template.sh"
SCRATCH_BASE  = Path("/oscar/scratch/eaguerov/mpi_runs")
RUNS_ROOT     = PROJECT_ROOT / "runs"

# Read seg-0 params to get exact omega_b and t_end
seg0_params = json.loads((RUNS_ROOT / "smoke_l5_seg0" / "params.json").read_text())
omega_b     = seg0_params["omega_b"]          # 1.5708 — SAME as seg-0
t_checkpoint = seg0_params["t_end"]           # ≈97.3

# Compute t_end: t_checkpoint + n_mix_cycles * T_per_st
L, H  = seg0_params["geometry"]["a"], seg0_params["geometry"]["b"]
th    = math.radians(seg0_params["theta_max"][0])
T_per = 2 * math.pi / omega_b
V     = L/4 * (H + 0.5*L*math.tan(th))
U     = V / (H*0.5) / T_per
T_bio = L / U
T_per_st = T_per / T_bio
n_mix = 10
t_end = round(t_checkpoint + T_per_st * n_mix + 50.0, 3)

params = {
    **{k: v for k, v in seg0_params.items() if not k.startswith("_")},
    "run_id": "mpi_same_omega_restart",
    "omega_b": omega_b,               # SAME omega_b → su = 1.0
    "omega_b_prev": omega_b,          # explicitly set so rescaling is no-op
    "n_mix_cycles": n_mix,
    "t_checkpoint": t_checkpoint,
    "t_end": t_end,
    "_ntasks": 16,
    "_walltime": "00:20:00",
    "_mem": "2G",
}

run_dir = RUNS_ROOT / params["run_id"]
run_dir.mkdir(parents=True, exist_ok=True)
params["_canonical_run_dir"] = str(run_dir.resolve())
(run_dir / "params.json").write_text(json.dumps(params, indent=2))

scratch_dir = SCRATCH_BASE / params["run_id"]
scratch_dir.mkdir(parents=True, exist_ok=True)
(scratch_dir / "params.json").write_text(json.dumps(params, indent=2))

# Stage checkpoint from seg-0
import shutil
ck = SCRATCH_BASE / "smoke_l5_seg0" / "checkpoint.dump"
shutil.copy2(ck, scratch_dir / "checkpoint.dump")

cmd = [
    "sbatch", "--no-requeue",
    "--job-name=mpi-same-omega",
    "--ntasks=16", "--time=00:20:00", "--mem-per-cpu=2G", "--cpus-per-task=1",
    f"--export=NONE,PARAMS={scratch_dir}/params.json,DUMP={scratch_dir}/checkpoint.dump",
    str(MPI_TEMPLATE.resolve()),
]
result = subprocess.run(cmd, capture_output=True, text=True, check=True)
job_id = re.search(r"(\d+)", result.stdout).group(1)

print(f"Experiment: same-omega_b MPI restart (su=1)")
print(f"omega_b    : {omega_b}  (same as seg-0 → rescaling is identity)")
print(f"t_checkpoint: {t_checkpoint:.2f}")
print(f"t_mix      : {t_checkpoint + T_per_st*n_mix:.2f}  (crash previously at ≈{t_checkpoint + T_per_st*n_mix - 0.1:.1f})")
print(f"t_end      : {t_end}")
print(f"Job        : {job_id}")
print()
print("Prediction A (rescaling is bug): run completes, finite kLa")
print("Prediction B (MPI restore is bug): run crashes at t≈t_mix, same as before")
