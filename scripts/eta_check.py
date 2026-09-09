import math
import numpy as np

def t_per_nd(rpm, theta):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)

T_ND = t_per_nd(22.5, 7.0)
print(f"T_per_nd(22.5rpm,7deg) = {T_ND:.5f}")

for label, path, elapsed_min in [
    ("cold_22.5", "/oscar/scratch/eaguerov/mpi_runs/l9_cold_rpm22.5/shear_stress.dat", 398.85),
    ("warm_22.5", "/oscar/scratch/eaguerov/mpi_runs/l9_warm_rpm22.5/shear_stress.dat", 398.85),
]:
    d = np.loadtxt(path, skiprows=1)
    t0, t1 = d[0, 1], d[-1, 1]
    cycles_done = (t1 - t0) / T_ND
    rate = elapsed_min / cycles_done
    print(f"{label}: t0={t0:.3f} t1={t1:.3f} cycles_done={cycles_done:.2f} "
          f"elapsed={elapsed_min:.1f}min -> {rate:.2f} min/cycle")
