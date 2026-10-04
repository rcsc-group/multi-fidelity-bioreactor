"""Scaling of the dense GP algebra that Yi-h needs: Cholesky time and memory against n.

n is the number of scalar outputs in the likelihood (runs x outputs per run).
Measures, on this node, the wall time of one Cholesky (the cost of one MCMC
likelihood evaluation) and the peak resident memory of building K and factorising it.
Use: uv run python scripts/bench_gp_scaling.py
"""
import resource
import time

import numpy as np
import scipy.linalg as sla

def main():
    rng = np.random.default_rng(0)
    d = 10
    for n in (1000, 2000, 4000, 8000, 12000):
        x = rng.random((n, d))
        K = np.empty((n, n))
        for i in range(0, n, 500):                 # build in blocks to avoid an n*n*d temporary
            diff = np.abs(x[i:i + 500, None, :] - x[None, :, :]) / 0.5
            r = np.sqrt((diff ** 2).sum(-1))
            K[i:i + 500] = (1 + np.sqrt(3) * r) * np.exp(-np.sqrt(3) * r)
        K[np.diag_indices(n)] += 1e-2
        t = time.perf_counter()
        sla.cholesky(K, lower=True, overwrite_a=True, check_finite=False)
        dt = time.perf_counter() - t
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6   # kB -> GB
        print(f"n={n:6d}  matrix {8 * n * n / 1e9:6.2f} GB  peak RSS {rss:6.2f} GB  "
              f"cholesky {dt:7.2f} s  {n ** 3 / 3 / dt / 1e9:6.1f} GFLOP/s")
        del K

if __name__ == "__main__":
    main()
