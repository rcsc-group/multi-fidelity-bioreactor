"""F3b step 2: replay of the L8 record at scalar level LS on 16 ranks, flux projection on (Ls > 8).
Usage: drv_f3b_replay.py LS PERIODS RUN_ID WALLTIME [NTASKS]   (node2356 excluded)"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import replay_gate as rg

ls, periods, rid, wall = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3], sys.argv[4]
nt = int(sys.argv[5]) if len(sys.argv) > 5 else 16
REC = Path("/oscar/data/dharri15/eaguerov/replay_records/f3b_l8_rpm32.5")
cp = json.loads((rg.run_dir("f3b_coupled") / "params.json").read_text())
t_per = rg.t_period(cp)
p = rg.make_params(rid, rg.BIN_DIR / "BioReactor-replayL16", fidelity=8, n_mix_cycles=0, post_cycles=0,
                   t_checkpoint=80 * t_per, restart_continue=1, n_rec_periods=22, replay_loop=0, rec_dir=str(REC),
                   replay_scalar_level=ls, replay_diag=0, replay_project_flux=1 if ls > 8 else 0)
p["t_end"] = periods * t_per
job = rg.submit(p, checkpoint=rg.scratch_dir("f3b_coupled") / "replay_start.dump", walltime=wall, ntasks=nt, mem="4G",
                exclude="node2356")
print("job", job, rg.wait(job), flush=True)
