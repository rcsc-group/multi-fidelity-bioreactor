import re
import sys
sys.path.insert(0, "/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor/.claude/worktrees/agent-a0a362873e179658b")
from scripts import replay_gate as rg

# L6 record: records per time unit and L8 flow: 1305 steps per period (kmix_l8 run), period = 0.60733 time
T = 0.6073291616284131
R6_PER_PERIOD = 306.7        # 920 records / 3 periods, measured
R8_PER_PERIOD = 1305.0       # kmix_l8_rpm32.5: 300141 iterations / 230 periods
rows = {}
for ls in (8, 9, 10):
    for pr in (0, 1):
        n = f"time_ls{ls}_p{pr}"
        pts = []
        for ln in (rg.scratch_dir(n) / "logstats.dat").read_text().splitlines():
            m = re.match(r"i: (\d+) t: (\S+) dt: (\S+) #Cells: (\d+) Wall clock time \(s\): (\S+)", ln)
            if m:
                pts.append((float(m.group(2)), float(m.group(5)), int(m.group(1))))
        (t0, w0, i0), (t1, w1, i1) = pts[0], pts[-1]
        rows[(ls, pr)] = {"dt_time": t1 - t0, "wall": w1 - w0, "iters": i1 - i0, "pts": len(pts)}
        print(n, rows[(ls, pr)])
print()
print("per record (L6 flow), wall s on 4 ranks: W = P + n*S, n = 2^(Ls-6)")
res = {}
for ls in (8, 9, 10):
    a, b = rows[(ls, 0)], rows[(ls, 1)]
    rec = a["dt_time"] / T * R6_PER_PERIOD
    W0, W1 = a["wall"] / rec, b["wall"] / (b["dt_time"] / T * R6_PER_PERIOD)
    n6 = 2 ** (ls - 6)
    S = W0 / n6
    P = W1 - W0
    res[ls] = (S, P)
    print(f"Ls={ls}: W(no proj)={W0:.3f} s/record  W(proj)={W1:.3f}  S={S:.4f} s/substep  P={P:.4f} s/projection")
print()
print("Extrapolation to the L8 flow, 22 post-release cycles, projection on, 4-rank wall x4 = core-seconds")
for ls in (8, 9, 10):
    S, P = res[ls]
    nsub = 2 ** (ls - 8)
    rec = 22 * R8_PER_PERIOD
    wall = rec * (P + nsub * S)
    print(f"Ls={ls}: records {rec:.0f}, substeps/record {nsub}, wall(4 ranks) {wall/3600:.1f} h, core-hours {4*wall/3600:.0f}")
