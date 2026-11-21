"""Cheap check of the R9 sampling hypothesis: per-period oxygen uptake in the coupled record run, normalised by the
undersaturation (dC/dt / (C_sat - C) ~ per-period kLa). Its period-to-period spread is the size of error a replay
of K sampled periods can carry; compare to the -8.9 % R9 miss at (37.5, 7)."""
import numpy as np
S = "/oscar/scratch/eaguerov/mpi_runs"
T = 24.88 / 41
for s in ("rpm32.5_th7", "rpm37.5_th7"):
    d = np.loadtxt(f"{S}/ft9_rec_{s}/tr_oxy.dat", skiprows=1)
    t, c = d[:, 1], d[:, 2]
    edges = t[0] + T * (1 + np.arange(30))     # periods after release (release 1 cycle in)
    ce = np.interp(edges, t, c)
    csat = 0.1426 * 1.0                          # oxygen sum at saturation ~ liquid volume (b0 column)
    k = np.diff(ce) / (csat - 0.5 * (ce[1:] + ce[:-1]))
    r = k[1:] / k[:-1]
    print(s, "per-period kLa proxy: first 6", np.round(k[:6] / k[:6].mean(), 3),
          "| period-to-period change sd %.1f %%" % (100 * np.std(np.diff(k[2:]) / k[2:-1])))
