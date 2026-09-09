import math
import numpy as np

def t_per_nd(rpm, theta=7.0):
    w = rpm * 2 * math.pi / 60.0
    L, H = 0.25, 2 * 0.03575
    T = 2 * math.pi / w
    V = L / 4 * (H + 0.5 * L * math.tan(math.radians(theta)))
    U = V / (H * 0.5) / T
    return T / (L / U)

# seg2 (25.0rpm, run bce29aa5, job 6114128): started 09:59:06, now 3:28:45 elapsed
ELAPSED_MIN = 3 * 60 + 28.75
d = np.loadtxt("/oscar/scratch/eaguerov/mpi_runs/bce29aa5/shear_stress.dat", skiprows=1)
cycles_done = (d[-1, 1] - d[0, 1]) / t_per_nd(25.0)
rate = ELAPSED_MIN / cycles_done
print(f"seg2 (25rpm): {cycles_done:.2f} cycles done in {ELAPSED_MIN:.0f}min "
      f"-> {rate:.2f} min/cycle")

N_TRANSITION = 25
cycles_left_this_seg = N_TRANSITION - cycles_done
h_this_seg_left = cycles_left_this_seg * rate / 60
print(f"this segment: {cycles_left_this_seg:.1f} cycles left -> {h_this_seg_left:.1f}h")

# remaining segments after this one: 27.5, 30, 32.5, 35, 37.5 = 5 more, full 25 cycles each
remaining_full_segs = 5
h_remaining_segs = remaining_full_segs * N_TRANSITION * rate / 60
print(f"{remaining_full_segs} more full segments x {N_TRANSITION} cycles -> {h_remaining_segs:.1f}h")

total_h = h_this_seg_left + h_remaining_segs
print(f"\nCHAIN COMPLETION ETA: ~{total_h:.1f}h from now (~{total_h/24:.1f} days)")
print(f"  -> reaches 30rpm (next measured point) in ~{h_this_seg_left + N_TRANSITION*rate/60:.1f}h")
print(f"  -> reaches 37.5rpm (final, measured) in ~{total_h:.1f}h")
