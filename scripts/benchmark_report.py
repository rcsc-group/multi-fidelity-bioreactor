"""Read timing + kLa from completed benchmark runs and print comparison table.

Run after benchmark_submit.py jobs complete.
Extrapolates wall-clock time to fidelity 10 assuming O(8^delta_f) scaling (2D CFD).
"""
import json, math
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[1]
RUNS_ROOT    = PROJECT_ROOT / "runs"
SCRATCH_BASE = Path("/oscar/scratch/eaguerov/mpi_runs")
FIDELITIES   = [4, 5, 6]
OMEGA_B      = [1.5708, 1.8326]


def wall_time_s(run_id, is_mpi=False):
    """Return wall-clock seconds from last logstats.dat entry, or None."""
    if is_mpi:
        p = SCRATCH_BASE / run_id / "logstats.dat"
    else:
        p = RUNS_ROOT / run_id / "logstats.dat"
    if not p.exists():
        return None
    last = p.read_text().strip().splitlines()[-1]
    for tok in last.split():
        if tok.replace('.','').isdigit():
            pass
    # parse: "i: N t: T dt: D #Cells: C Wall clock time (s): W CPU time (s): C"
    parts = last.split()
    for i, tok in enumerate(parts):
        if tok == "(s):" and i > 0 and parts[i-1] == "time":
            if "Wall" in parts[max(0,i-3):i]:
                return float(parts[i+1])
    return None


def kla(run_id):
    p = RUNS_ROOT / run_id / "results.json"
    if not p.exists():
        return None
    d = json.loads(p.read_text())
    return d.get("kLa_25")


def status(run_id, is_mpi=False):
    if (RUNS_ROOT / run_id / "results.json").exists():
        return "done"
    if is_mpi:
        p = SCRATCH_BASE / run_id / "logstats.dat"
    else:
        p = RUNS_ROOT / run_id / "logstats.dat"
    if p.exists():
        return "running/incomplete"
    return "not_started"


print("=" * 90)
print(f"{'Label':<28} {'Status':<18} {'Wall(s)':<12} {'kLa_25':<12}")
print("=" * 90)

serial_times = {}  # fidelity → (t_s0, t_s1)
mpi_times    = {}  # fidelity → (t_a, t_b)

for fi in FIDELITIES:
    s0_id = f"bench_f{fi}_serial_s0"
    s1_id = f"bench_f{fi}_serial_s1"
    ma_id = f"bench_f{fi}_mpi_a"
    mb_id = f"bench_f{fi}_mpi_b"

    t_s0 = wall_time_s(s0_id, is_mpi=False)
    t_s1 = wall_time_s(s1_id, is_mpi=False)
    t_ma = wall_time_s(ma_id, is_mpi=True)
    t_mb = wall_time_s(mb_id, is_mpi=True)

    serial_times[fi] = (t_s0, t_s1)
    mpi_times[fi]    = (t_ma, t_mb)

    for label, rid, t, is_mpi in [
        (f"f{fi}_serial_s0 (fresh, ob[0])", s0_id, t_s0, False),
        (f"f{fi}_serial_s1 (ckpt, ob[1])", s1_id, t_s1, False),
        (f"f{fi}_mpi_a     (fresh, ob[0])", ma_id, t_ma, True),
        (f"f{fi}_mpi_b     (fresh, ob[1])", mb_id, t_mb, True),
    ]:
        st  = status(rid, is_mpi)
        kv  = kla(rid)
        t_s = f"{t:.0f}" if t else "—"
        k_s = f"{kv:.4f}" if kv else "—"
        print(f"  {label:<26} {st:<18} {t_s:<12} {k_s}")
    print()

print("=" * 90)
print("Cost to produce 2 kLa points (wall-clock seconds):")
print(f"{'Fidelity':<12} {'Serial+ckpt':<18} {'MPI×2 indep':<18} {'Winner'}")
print("-" * 60)

serial_totals = {}
mpi_totals    = {}

for fi in FIDELITIES:
    t_s0, t_s1 = serial_times[fi]
    t_ma, t_mb = mpi_times[fi]
    st = (t_s0 or 0) + (t_s1 or 0) if (t_s0 and t_s1) else None
    mt = (t_ma or 0) + (t_mb or 0) if (t_ma and t_mb) else None
    serial_totals[fi] = st
    mpi_totals[fi]    = mt
    st_s = f"{st:.0f}s" if st else "—"
    mt_s = f"{mt:.0f}s" if mt else "—"
    winner = ("serial+ckpt" if (st and mt and st < mt)
              else "MPI×2" if (st and mt and mt < st)
              else "—")
    print(f"  {fi:<12} {st_s:<18} {mt_s:<18} {winner}")

# Extrapolate to fidelity 10
print()
print("Extrapolation to fidelity 10 (scaling ~8× per fidelity level in 2D):")
print(f"  (cost doubles in steps, quadruples in cells → 8× total per level)")
print()

for arm_name, totals in [("Serial+ckpt", serial_totals), ("MPI×2 indep", mpi_totals)]:
    # Fit log-linear: log(T) = a + b*f
    pts = [(f, t) for f, t in totals.items() if t]
    if len(pts) >= 2:
        fs = [p[0] for p in pts]
        ts = [math.log(p[1]) for p in pts]
        b = (ts[-1] - ts[0]) / (fs[-1] - fs[0])
        a = ts[0] - b * fs[0]
        t10 = math.exp(a + b * 10)
        t10_h = t10 / 3600
        print(f"  {arm_name}: scaling = {b/math.log(8):.2f}×/level (expected 1.0×)")
        print(f"  → Fidelity 10 estimate: {t10:.0f}s = {t10_h:.1f} hours (2 kLa points)")
    else:
        print(f"  {arm_name}: insufficient data to extrapolate")
    print()
