"""Probe: 25 rpm theta 7 L10 runs (both names) against Kim's digitised values; diary 2026-10-22 quotes 0.00133 / 0.130.
Data: runs/x8_l10_rpm25_th7, runs/l10_fig13a_xlevel_rpm25 (shear_stress.dat via scripts/hydro_dataset_v2.qois)."""
import sys
sys.path[:0] = ["/oscar/data/dharri15/eaguerov/Github/multi-fidelity-bioreactor"]
from scripts.hydro_dataset_v2 import qois
for run in ["x8_l10_rpm25_th7", "l10_fig13a_xlevel_rpm25", "x8e_l10_rpm25_th7"]:
    try:
        q = qois(run)
        print(run, "tau_amp_t=%.6g edr_max_t=%.6g window=%s n_steps=%d" % (q["tau_amp_t"], q["edr_max_t"], q["window"], q["n_steps"]))
    except Exception as e:
        print(run, "ERROR", type(e).__name__, e)
print("Kim 25 rpm: tau_liq_mean 0.001247, Ediss_liq_mean 0.123088")
