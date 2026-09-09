import math
import numpy as np

def t_per_nd(rpm, theta=7.0):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)

ELAPSED_H = 399.7 / 60.0  # sacct: 6076948/6076949, 06:39:42

rates_min = []
for path in ["/oscar/scratch/eaguerov/mpi_runs/l9_cold_rpm22.5/shear_stress.dat",
             "/oscar/scratch/eaguerov/mpi_runs/l9_warm_rpm22.5/shear_stress.dat"]:
    d = np.loadtxt(path, skiprows=1)
    cycles = (d[-1, 1] - d[0, 1]) / t_per_nd(22.5)
    rates_min.append(ELAPSED_H * 60.0 / cycles)
rate_h = max(rates_min) / 60.0  # conservative (slower observed), in HOURS/cycle
print(f"measured rate: {rates_min[0]:.2f} / {rates_min[1]:.2f} min/cycle "
      f"-> using {rate_h*60:.2f} min/cycle (conservative)\n")

N_MIX = 45
d_cold = np.loadtxt("/oscar/scratch/eaguerov/mpi_runs/l9_cold_rpm22.5/shear_stress.dat", skiprows=1)
d_warm = np.loadtxt("/oscar/scratch/eaguerov/mpi_runs/l9_warm_rpm22.5/shear_stress.dat", skiprows=1)
cyc_left_cold = N_MIX - (d_cold[-1, 1] - d_cold[0, 1]) / t_per_nd(22.5)
cyc_left_warm = N_MIX - (d_warm[-1, 1] - d_warm[0, 1]) / t_per_nd(22.5)
h_first_free = min(cyc_left_cold, cyc_left_warm) * rate_h
h_both_free  = max(cyc_left_cold, cyc_left_warm) * rate_h
print(f"22.5rpm pair: {cyc_left_cold:.1f}/{cyc_left_warm:.1f} cycles left "
      f"-> first 32 CPUs free in {h_first_free:.1f}h, both free in {h_both_free:.1f}h\n")

# Best case: cold_30 and cold_37.5 slot into the two 32-CPU vacancies as
# they open (one at h_first_free, one at h_both_free), running in parallel
# with each other once both are in. Chain seg0 then needs 64 CPUs free
# again -> earliest after ONE of the two 45-cycle cold refs finishes.
cold_ref_h = N_MIX * rate_h
t_cold1_start, t_cold2_start = h_first_free, h_both_free
t_cold1_done = t_cold1_start + cold_ref_h
t_cold2_done = t_cold2_start + cold_ref_h
h_64_free_again = min(t_cold1_done, t_cold2_done)
print(f"cold_30/cold_37.5 (45 cycles each): start at {t_cold1_start:.1f}h / "
      f"{t_cold2_start:.1f}h, each takes {cold_ref_h:.1f}h -> "
      f"finish at {t_cold1_done:.1f}h / {t_cold2_done:.1f}h")
print(f"-> 64 CPUs free again (chain seg0 can start) at {h_64_free_again:.1f}h\n")

N_SEG, CYC_PER_SEG = 8, 25
chain_duration_h = N_SEG * CYC_PER_SEG * rate_h
print(f"chain: 8 segments x 25 cycles, sequential -> {chain_duration_h:.1f}h pure run time")
print(f"\nCHAIN COMPLETION ETA: ~{h_64_free_again + chain_duration_h:.0f}h from now "
      f"(best case: SLURM schedules the chain ahead of/alongside nothing else "
      f"competing for the freed slot)")
print(f"That's ~{(h_64_free_again + chain_duration_h)/24:.1f} days from now.")
