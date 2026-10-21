"""Periodic-flow replay at L8 (the PI idea: stitch a quasi-steady periodic flow): loop the first K recorded periods of the L8 record for 22 cycles, scalar level 8, 16 ranks.
Usage: drv_f3b_loop.py K RUN_ID WALLTIME   (submits only; node2356 excluded). Compare with the coupled run f3b_coupled (kLa_1T_10 46.93, kLa_1T_25 28.64, dtmix_0.50 11.794).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts import replay_gate as rg

k, rid, wall = int(sys.argv[1]), sys.argv[2], sys.argv[3]
REC = Path("/oscar/data/dharri15/eaguerov/replay_records/f3b_l8_rpm32.5")
cp = json.loads((rg.run_dir("f3b_coupled") / "params.json").read_text())
t_per = rg.t_period(cp)
p = rg.make_params(rid, rg.BIN_DIR / "BioReactor-replayL16", fidelity=8, n_mix_cycles=0, post_cycles=0,
                   t_checkpoint=80 * t_per, restart_continue=1, n_rec_periods=k, replay_loop=1, rec_dir=str(REC),
                   replay_scalar_level=8, replay_diag=0, replay_project_flux=0)
p["t_end"] = 22 * t_per
job = rg.submit(p, checkpoint=rg.scratch_dir("f3b_coupled") / "replay_start.dump", walltime=wall, ntasks=16, mem="4G", exclude="node2356")
print("job", job, flush=True)
