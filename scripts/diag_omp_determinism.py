"""Is the L4 kLa of test_restart_continue deterministic under OpenMP?

CI (4-core runners, -fopenmp) has failed test_kla_survives_a_restart_only_when_
segments_are_stitched on every run since it was added. The continuous run's
kLa_10 was 1.057, 1.262, 1.321 and 1.394 on four runs with identical inputs,
while locally, on 1 core, the test passes at ratio 1.000.
Hypothesis: OpenMP reduction order makes the run nondeterministic, and the
5-point kLa_10 fit amplifies that past the test's 10% tolerance.
Falsifier: two OMP_NUM_THREADS=4 runs give the same kLa_10 (within 1e-6 relative),
or two OMP_NUM_THREADS=1 runs differ.

Usage (4 cores):  uv run python scripts/diag_omp_determinism.py OUTDIR
"""
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests.conftest import run_bioreactor                       # noqa: E402
from tests.verification.test_restart_continue import _params    # noqa: E402
from scripts.postprocess import _compute_c_star, _kla_5pt_at_threshold  # noqa: E402

out = pathlib.Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
for name, threads in (("omp4_a", "4"), ("omp4_b", "4"), ("omp1_a", "1"), ("omp1_b", "1")):
    d = run_bioreactor(_params(name, 20), out, timeout=7200,
                       env={"OMP_NUM_THREADS": threads})
    t, c = _compute_c_star(d)
    k = [_kla_5pt_at_threshold(t, c, thr) for thr in (0.10, 0.25, 0.50)]
    print(f"{name}: threads={threads} samples={len(t)} t_end={t[-1]:.4f} "
          f"kLa_10/25/50 = {k[0]:.6g} {k[1]:.6g} {k[2]:.6g}", flush=True)
