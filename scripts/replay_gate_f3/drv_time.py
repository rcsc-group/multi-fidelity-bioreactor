"""Cost of replay on fine scalar grids, driven by the existing 4-rank L6 record mM_coupled. Usage: drv_time.py PERIODS LS:PROJ ..."""
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/.claude/worktrees/agent-a0a362873e179658b")
from scripts import replay_gate as rg

per = float(sys.argv[1])
cases = [tuple(int(x) for x in a.split(":")) for a in sys.argv[2:]]
CID, SPIN, REC = "mM_coupled", 10, 3
import json
cp = json.loads((rg.run_dir(CID) / "params.json").read_text())
t_per = rg.t_period(cp)


def one(ls, proj):
    n = f"time_ls{ls}_p{proj}"
    p = rg.make_params(n, rg.BIN_DIR / "BioReactor-replayMc", fidelity=6, n_mix_cycles=0, post_cycles=0,
                       t_checkpoint=SPIN * t_per, restart_continue=1, n_rec_periods=REC, replay_loop=0,
                       rec_dir=str(rg.scratch_dir(CID)), replay_scalar_level=ls, replay_diag=0, replay_project_flux=proj)
    p["t_end"] = per * t_per
    jj = rg.submit(p, checkpoint=rg.scratch_dir(CID) / "replay_start.dump", walltime="08:00:00", ntasks=4, mem="4G")
    return n, jj, rg.wait(jj)


with ThreadPoolExecutor(len(cases)) as ex:
    for f in [ex.submit(one, *c) for c in cases]:
        print(f.result(), flush=True)
