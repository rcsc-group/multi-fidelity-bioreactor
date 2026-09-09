import math
import numpy as np

def t_per_nd(rpm, theta=7.0):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)

ELAPSED_H = 399.7 / 60.0
d_cold = np.loadtxt("/oscar/scratch/eaguerov/mpi_runs/l9_cold_rpm22.5/shear_stress.dat", skiprows=1)
d_warm = np.loadtxt("/oscar/scratch/eaguerov/mpi_runs/l9_warm_rpm22.5/shear_stress.dat", skiprows=1)
rate_h = max(ELAPSED_H / ((d_cold[-1,1]-d_cold[0,1]) / t_per_nd(22.5)),
            ELAPSED_H / ((d_warm[-1,1]-d_warm[0,1]) / t_per_nd(22.5)))

N_MIX = 45
cyc_left_cold = N_MIX - (d_cold[-1, 1] - d_cold[0, 1]) / t_per_nd(22.5)
cyc_left_warm = N_MIX - (d_warm[-1, 1] - d_warm[0, 1]) / t_per_nd(22.5)
h_first_free = min(cyc_left_cold, cyc_left_warm) * rate_h   # chain seg0 grabs this slot
h_both_free  = max(cyc_left_cold, cyc_left_warm) * rate_h   # cold_30 grabs this slot

N_SEG, CYC_PER_SEG = 8, 25
chain_h = N_SEG * CYC_PER_SEG * rate_h
cold_ref_h = N_MIX * rate_h

chain_done  = h_first_free + chain_h
cold30_done = h_both_free + cold_ref_h              # runs alongside the chain
cold375_done = cold30_done + cold_ref_h              # follows in the same slot after cold_30

print(f"rate: {rate_h*60:.2f} min/cycle")
print(f"chain starts ~{h_first_free:.1f}h, runs {chain_h:.1f}h -> done ~{chain_done:.1f}h")
print(f"cold_30 starts ~{h_both_free:.1f}h, runs {cold_ref_h:.1f}h -> done ~{cold30_done:.1f}h")
print(f"cold_37.5 starts ~{cold30_done:.1f}h, runs {cold_ref_h:.1f}h -> done ~{cold375_done:.1f}h")
print(f"\nALL DATA READY: ~{max(chain_done, cold375_done):.0f}h from now "
      f"({max(chain_done, cold375_done)/24:.1f} days) -- vs ~56h before reordering")
