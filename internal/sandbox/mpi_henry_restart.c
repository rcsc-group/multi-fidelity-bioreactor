/**
Minimal reproducer: canonical henry.h + MPI + checkpoint restart.

Tests whether the stale non-leaf scalar pool issue (Bug 2) exists in
canonical henry.h, i.e. whether NaN appears in tracer diffusion specifically
on checkpoint restart under MPI.

EMBED is intentionally EXCLUDED to avoid the d ≤ 0 issue in canonical h_relax
(Bug 1). EMBED + large alpha is a separate upstream bug. Here we isolate Bug 2.

Bug 2 mechanism: automatic scalars (phi1, phi2, r, res) are allocated from a pool.
On a fresh run the pool is near-zero. On restart (after a segment that has filled
the pool with O(1) values), coarse non-leaf cells of these automatic scalars carry
stale pool values, because:
  - foreach() only visits LEAF cells
  - phi1/phi2 are never restricted in canonical henry.h (no restriction({phi1,phi2}))
  - r[] is never restricted before mg_solve
  - Non-leaf MPI ghost cells for r[], D.x[], beta.x[] may be stale if qcc's stencil
    analysis doesn't see the function-pointer accesses in h_relax/h_residual

Hypotheses:
  A: canonical henry.h + MPI + restart → NaN   (Bug 2 is upstream)
  B: canonical henry.h + MPI + restart → PASS  (Bug 2 was specific to henry_oxy2.h)

Usage (run via mpi_henry_restart.sh SLURM script):
  Run 1 (fresh):   srun -n 4 ./mpi_henry_restart
  Run 2 (restart): srun -n 4 ./mpi_henry_restart restart

PASS = no NaN detected throughout both runs.
FAIL = NaN/Inf detected on restart → Hypothesis A confirmed.
*/

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "henry.h"

// Two-phase properties: water (phase 1) and air (phase 2)
#define RHO1  1.0
#define RHO2  0.01
#define MU1   0.01
#define MU2   1e-4

// Henry's law tracer: small alpha so canonical h_relax stays well-conditioned
scalar oxygen[];
scalar * stracers = {oxygen};
#define D1_VAL    1e-4
#define ALPHA_VAL 2.         // small: beta << D/Delta, d > 0 guaranteed

// Simulation parameters
#define LEVEL       5        // 2^5 = 32 cells/side: fast, fits on 4 MPI ranks
#define T_DUMP      0.005    // checkpoint time (end of fresh run)
#define T_END       0.010    // end of restart run

int do_restart;  // 1 if argv[1] is provided

int main (int argc, char * argv[])
{
  do_restart = (argc > 1);

  rho1 = RHO1; rho2 = RHO2;
  mu1  = MU1;  mu2  = MU2;

  oxygen.D1    = D1_VAL;
  oxygen.D2    = D1_VAL * ALPHA_VAL;
  oxygen.alpha = ALPHA_VAL;

  L0 = 1.0;
  DT = 1e-3;
  init_grid (1 << LEVEL);
  run();
}

event init (t = 0)
{
  if (do_restart) {
    if (!restore (file = "checkpoint.dump")) {
      fprintf (stderr, "ERROR: restore() failed — run fresh first\n");
      exit (1);
    }
    fprintf (stderr, "[RESTART] restored at t=%g\n", t);
  } else {
    // Flat gas-liquid interface at y = 0.5
    fraction (f, y - 0.5);
    foreach ()
      oxygen[] = f[] > 0.5 ? 1.0 : 0.0;
    fprintf (stderr, "[FRESH] initialized at t=0\n");
  }
}

// Write checkpoint at T_DUMP (fresh run only)
event write_checkpoint (t = T_DUMP) {
  if (!do_restart) {
    dump (file = "checkpoint.dump");
    fprintf (stderr, "[DUMP] checkpoint at t=%g\n", t);
  }
}

// NaN guard: fires every step
event check_nan (i++)
{
  double omax = -1e300;
  int fail = 0;
  foreach (reduction(max:omax) reduction(||:fail)) {
    double v = oxygen[];
    if (v > omax) omax = v;
    if (!isfinite (v)) fail = 1;
  }
  if (fail) {
    fprintf (stderr,
      "FAIL: NaN/Inf in oxygen at t=%.6g i=%d\n", t, i);
    exit (1);
  }
  if (i % 10 == 0)
    fprintf (stderr,
      "OK: t=%.5g i=%d  oxygen_max=%.6g\n", t, i, omax);
}

// End of simulation
event end (t = do_restart ? T_END : T_DUMP)
{
  if (do_restart)
    fprintf (stderr,
      "PASS: restart run (t=%.3g..%.3g) — no NaN detected\n",
      T_DUMP, T_END);
  else
    fprintf (stderr,
      "PASS: fresh run (t=0..%.3g) — dump written\n", T_DUMP);
}
