"""Results 1 (DFD talk) L8 fill-in at theta = 7 (user 2026-11-21: the ends of the rpm axis look empty).
Same protocol and binary source as the kmix_l8 sweep the figure uses (cold start, n_mix_cycles 80, lean MPI build of
f1c11e0, rebuilt because the scratch binary was purged). Control: ctl8_rpm35 (kmix_l8_rpm35 rerun) with the rebuilt binary must
match the original QoIs (checks that the rebuild is the same solver).
Stages: build | submit | check
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import replay_gate as rg  # noqa: E402
from scripts.simulate import submit_slurm  # noqa: E402

SRC = Path(sys.argv[2]) if len(sys.argv) > 2 else None
BIN = rg.BIN_DIR / "BioReactor-mpi-lean-f1c11e0"
RPM = [15.625, 16.25, 16.875, 18.125, 18.75, 19.375, 33.125, 33.75, 34.375, 35.625, 36.25, 36.875]
RUNS = ROOT / "runs"
TEMPLATE_RUN = RUNS / "kmix_l8_rpm35"


def params(rpm, run_id):
    p = json.loads((TEMPLATE_RUN / "params.json").read_text())
    p = {k: v for k, v in p.items() if not k.startswith("_") and not k.endswith("_prev")
         and k not in ("t_checkpoint", "restart_continue", "next_run_id")}
    p.update(run_id=run_id, omega_b=2 * 3.141592653589793 * rpm / 60.0, _binary=str(BIN))
    return p


if __name__ == "__main__":
    st = sys.argv[1]
    if st == "build":
        print(rg.build_binary("mpi-lean-f1c11e0", ["-DEXTRA_TRACERS=0"], src=SRC, mpi=True))
    elif st == "submit":
        assert BIN.exists(), BIN
        jobs = [(f"kmix_l8_rpm{r:g}", r) for r in RPM] + [("ctl8_rpm35", 35.0)]   # not a kmix_l8 name: the figure globs kmix_l8*
        for rid, r in jobs:
            if (RUNS / rid).exists():
                print("exists", rid); continue
            j = submit_slurm(params(r, rid), project_root=ROOT, runs_root=RUNS, walltime="03:00:00", ntasks=32,
                             template=ROOT / "config" / "slurm_mpi_template.sh", cpus=1, mem="4G",
                             exclude="node2356,node2410,node2411", partition="batch", allow_duplicate=True)   # 32 ranks: 5.4x the 8-rank step rate (measured 2026-11-21)
            print(rid, "job", j)
