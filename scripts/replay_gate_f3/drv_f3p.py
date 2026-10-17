"""F3a': L6 coupled record (checksum build) then replays at Ls 6,7,8 with and without flux projection, in parallel."""
import json
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/.claude/worktrees/agent-a0a362873e179658b")
from scripts import replay_gate as rg

CID, SPIN, PER = "f3p_coupled", 10, 5
rec = rg.BIN_DIR / "BioReactor-recP"
rep = rg.BIN_DIR / "BioReactor-replayP"
cp = rg.make_params(CID, rec, fidelity=6, n_mix_cycles=SPIN, post_cycles=PER, n_rec_periods=PER,
                    rec_dir=str(rg.scratch_dir(CID)))
t_per = rg.t_period(cp)
job = rg.submit(cp, walltime="01:00:00")
print("coupled", job, rg.wait(job), flush=True)


def one(ls, proj):
    rid = f"f3p_ls{ls}_p{proj}"
    p = rg.make_params(rid, rep, fidelity=6, n_mix_cycles=0, post_cycles=0, t_checkpoint=SPIN * t_per,
                       restart_continue=1, n_rec_periods=PER, replay_loop=0, rec_dir=str(rg.scratch_dir(CID)),
                       replay_scalar_level=ls, replay_diag=1, replay_project_flux=proj)
    p["t_end"] = PER * t_per
    j = rg.submit(p, checkpoint=rg.scratch_dir(CID) / "replay_start.dump", walltime="06:00:00")
    return rid, j, rg.wait(j)


with ThreadPoolExecutor(5) as ex:
    futs = [ex.submit(one, ls, pr) for ls, pr in ((6, 0), (7, 0), (7, 1), (8, 0), (8, 1))]
    for f in futs:
        print(f.result(), flush=True)
