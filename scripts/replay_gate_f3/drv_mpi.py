"""Usage: drv_mpi.py ID REC_TAG REPLAY_TAG NTASKS PERIODS LS:PROJ [LS:PROJ ...]   (L6 flow, release after 10 cycles)"""
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/.claude/worktrees/agent-a0a362873e179658b")
from scripts import replay_gate as rg

rid, rtag, ptag, nt, per = sys.argv[1:6]
nt, per = int(nt), int(per)
cases = [tuple(int(x) for x in a.split(":")) for a in sys.argv[6:]]
CID, SPIN = f"{rid}_coupled", 10
cp = rg.make_params(CID, rg.BIN_DIR / f"BioReactor-{rtag}", fidelity=6, n_mix_cycles=SPIN, post_cycles=per,
                    n_rec_periods=per, rec_dir=str(rg.scratch_dir(CID)))
t_per = rg.t_period(cp)
j = rg.submit(cp, walltime="02:00:00", ntasks=nt, mem="3G")
print("coupled", j, rg.wait(j), flush=True)


def one(ls, proj):
    n = f"{rid}_ls{ls}_p{proj}"
    p = rg.make_params(n, rg.BIN_DIR / f"BioReactor-{ptag}", fidelity=6, n_mix_cycles=0, post_cycles=0,
                       t_checkpoint=SPIN * t_per, restart_continue=1, n_rec_periods=per, replay_loop=0,
                       rec_dir=str(rg.scratch_dir(CID)), replay_scalar_level=ls, replay_diag=1, replay_project_flux=proj)
    p["t_end"] = per * t_per
    jj = rg.submit(p, checkpoint=rg.scratch_dir(CID) / "replay_start.dump", walltime="06:00:00", ntasks=nt, mem="3G")
    return n, jj, rg.wait(jj)


with ThreadPoolExecutor(max(1, len(cases))) as ex:
    for f in [ex.submit(one, *c) for c in cases]:
        print(f.result(), flush=True)
