import sys
import numpy as np
sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/.claude/worktrees/agent-a0a362873e179658b")
from scripts import replay_gate as rg

def tr(run):
    return rg.load_tr_oxy(rg.run_dir(run) / "tr_oxy.dat")

def cmp(a, b, label):
    d = rg.rel_diff_by_time(tr(a), tr(b), [2, 8, 9])
    print(f"{label:48s} common {d['n_common']:4d} of {d['n_a']}/{d['n_b']} max rel diff", {k: f"{v:.2e}" for k, v in d['max'].items()})

def chk(run):
    return np.loadtxt(rg.scratch_dir(run) / "rec_checksum.txt")

cmp("mM_coupled", "mM_ls6_p0", "MPI4 coupled vs MPI4 replay Ls=6")
cmp("mS_coupled", "mS_ls6_p0", "serial coupled vs serial replay Ls=6")
cmp("mS_coupled", "mM_coupled", "serial coupled vs MPI4 coupled")
cmp("mS_ls6_p0", "mM_ls6_p0", "serial replay Ls=6 vs MPI4 replay Ls=6")
cmp("mS_ls7_p1", "mM_ls7_p1", "serial replay Ls=7 proj vs MPI4 replay Ls=7 proj")
for a, b in (("mM_coupled", "mM_ls6_p0"), ("mS_coupled", "mS_ls6_p0")):
    ca, cb = chk(a), chk(b)
    n = min(len(ca), len(cb))
    print(a, b, "checksum rows", len(ca), len(cb), "max |diff| f, c2, oxy checksum cols:",
          [float(np.max(np.abs(ca[:n, c] - cb[:n, c]))) for c in (2, 4, 5)])
