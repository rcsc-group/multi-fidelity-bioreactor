// Primary author: Minki Kim
// Contributor: Radu Cimpeanu
// Contributor: Daniel Harris
// Date: 03/04/2025
//
// ============================================================================
// ANNOTATION LEGEND (this project's additions on top of the file above)
// ============================================================================
// This file began as Minki Kim et al.'s driver code, published at
// https://github.com/rcsc-group/BioReactor/tree/main/DriverCodes (predates
// the Basilisk version this fork compiles against). Everything below without
// a `[PROJECT ...]` tag is Kim et al.'s original code/comments, unmodified
// in substance even where reformatted.
//
//   [PROJECT ADDED: ...]    a block that does not exist in Kim's upstream file
//   [PROJECT CHANGED: ...]  upstream had different logic/values here; why we
//                           changed it
//   [PROJECT REMOVED: ...]  upstream had code here that this fork deletes
//                           outright (noted at the point of removal, or in
//                           a summary comment where nothing replaces it)
//
// Full account of what's changed and why is in diary.md (repo root) and the
// project's own commit history (`git log -- src/BioReactor.c`).
// ============================================================================

// [PROJECT ADDED, 2026-10-16] Passive-scalar record/replay (diary.md 2026-10-16). Oxygen `oxy` and
// tracer `c2` never enter rho or mu, so the flow does not depend on them. Two compile-time modes of
// THIS file (both off by default; the default build is the unchanged coupled solver):
//   -DRECORD_FLUX  coupled run that also writes, once per step from t_rec0 for n_rec_periods, the face
//                  flux uf exactly as the VOF step of f consumes it, plus a start checkpoint
//                  (replay_start.dump, written at the end of the step that lands on t_rec0).
//   -DREPLAY_FLUX  no Navier-Stokes: restore replay_start.dump (params: t_checkpoint = t_rec0,
//                  restart_continue = 1, n_mix_cycles = 0), then each step read (i,t,dt,uf) and run
//                  only what the coupled run does to f, phi1, phi2 and the scalar diffusion.
// What else touches f or the scalars between steps, and whether replay reproduces it:
//   vof_advection of f and f.tracers (phi1, phi2 of oxy and c2)   same code, same uf, dt, parity of i
//   conserving.h's momentum tracers q1, q2                         not needed: they never feed back on f
//   event adapt (centered.h): fractions_cleanup(cs,fs), uf=0 where fs=0   cleanup repeated after restore; the uf
//                                                                  zeroing is already in the recorded uf
//   event tracer (t=t_mix), event oxygen (t=t_mix; every step)    same bodies, triggered by t >= t_mix
//   iforce.h event acceleration: f = clamp(f,0,1) every step      same clamp (event clamp_f); this one mattered
//   remove_drop                                                    same event (a function of f only)
//   boundary conditions of f, oxy, c2, embed no-flux of the scalars  same declarations
//   NOT reproduced: anything the flow does to f through u itself (that is what the record replaces);
//   looped replay is exact only to the extent that the flow is periodic.
#ifndef RECORD_FLUX
#define RECORD_FLUX 0
#endif
#ifndef REPLAY_FLUX
#define REPLAY_FLUX 0
#endif
#if RECORD_FLUX && REPLAY_FLUX
#error "RECORD_FLUX and REPLAY_FLUX are mutually exclusive"
#endif

// Core Basilisk modules for embedding geometries and centered Navier-Stokes solver
#include "embed.h"
#if REPLAY_FLUX
// Frozen-flow stand-in for centered.h + two-phase.h + conserving.h: only what VOF advection of f, its
// tracers and the scalar diffusion read. `init` is registered here so that restore() runs before the
// replay events below it (events execute in registration order).
#include "run.h"                   // time loop, utils.h (DT, statsf, perf, minmod2)
#include "timestep.h"
#include "output.h"                // restore
event init (i = 0);
// Event order. qcc registers the events marked `last` after all the others, so in the coupled solver every
// user event (release, outputs, ...) runs at the START of a step, before the vof / diffusion / projection
// events of that step. The same declarations (as in centered.h) give the replay the same order.
event vof (i++,last);
event tracer_diffusion (i++,last);
face vector uf[];                  // read from the record, one step at a time
scalar f[], * interfaces = {f};
#include "vof.h"
#include "curvature.h"             // position(), used by the interface height outputs
#include "poisson.h"               // mg_solve, used by tracer_diffusion in henry_oxy2.h
#else
#include "navier-stokes/centered.h"
#endif

// Define dynamic viscosity as a function of the volume fraction (f)
// Uses harmonic averaging for two-phase flow
#define mu(f) ( 1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)  )

#if !REPLAY_FLUX
// Modules for two-phase flow (e.g., air and water) with surface tension and volume conservation
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"
#endif

// Modules for oxygen transfer + tracers
#include "henry_oxy2.h"

// Custom headers foroxy specific configurations and visualization
#if !REPLAY_FLUX
#include "view3.h"
#endif
#include "utils2.h"

// [PROJECT ADDED] Reinstates upstream's REMOVE_DROP feature (droplet/bubble
// removal via tag()), deleted during an earlier cleanup because the
// validated baseline (theta=7deg/32.5rpm) doesn't need it -- see the
// [PROJECT REMOVED] note below. Brought back as a runtime params.json
// toggle (params.remove_drop) rather than the original compile-time
// #define, so the sweep pipeline can test it on/off from the same binary.
// remove_droplets() signature verified against the canonical Basilisk
// install (tag.h) -- unchanged from what upstream called.
#include "tag.h"

// Mathematical constants (e.g., M_PI from math.h)
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

// [PROJECT ADDED] JSON parameter reader (jsmn-based). Upstream took L_bio,
// ANGLE, RPM as three positional argv[] floats; this project needed many
// more tunable parameters (fidelity, geometry shape, multi-harmonic forcing,
// checkpoint-restart timing, ...) for the sweep/optimization pipeline, so
// argv parsing was replaced with a params.json file. See main()'s
// [PROJECT CHANGED] block below for exactly which upstream argv reads this
// replaces.
#include "params_read.h"

// [PROJECT REMOVED] Upstream also defines these feature flags, along with
// the code blocks they gate, all deleted from this fork (unused for the
// validated baseline configuration; not required by Kim et al.'s published
// results at theta=7deg/32.5rpm, which use none of them):
//   CONTACT (contact-angle BC), OXYGEN_CIRCLE (circular oxygen init),
//   HORIZONTAL_MIXL/HORIZONTAL_MIXR (horizontal tracer mixing variants),
//   VERTICAL_MIXDOWN (downward tracer mixing variant), AMR (adaptive mesh
//   refinement — see MINLEVEL/MAXLEVEL note above; upstream itself ran with
//   AMR=0 for its own published results, so this is not a discrepancy vs.
//   Kim's actual runs), CFL_COND (custom CFL override), DUMP (upstream's
//   own dump mechanism — this fork has its own checkpoint-restart
//   dump/restore, added separately, see event dump_checkpoint below).
//   REMOVE_DROP (droplet/bubble removal) was also removed here, but later
//   reinstated as a runtime params.remove_drop toggle -- see the tag.h
//   include and remove_drop event below.
// Flags to control the inclusion of features (set to 1 = enable, 0 = disable)
#define EMBED            1   // Enable embedded boundary for solid geometry
#define OXYGEN           1   // Enable oxygen concentration simulation
#define OXYGEN_AIR       1   // Initial distribution (air side) of oxygen (if OXYGEN == 1)

// Mixing strategies for tracer release
#define TRACER           1   // Enable passive tracer simulation
#define VERTICAL_MIXUP   1   // Initial distribution (top side) of tracer: Vertical mixing (if TRACER == 1)

// Other simulation options
#define ACCELERATION     1   // Enable acceleration (rocking motion)
#define NORMCAL          1   // Calculate statistics (norms)
#if REPLAY_FLUX
#undef ACCELERATION
#define ACCELERATION     0   // no forcing: the recorded flux is the flow
#undef NORMCAL
#define NORMCAL          0   // replaced by the replay version of event normcal below
#endif
// [PROJECT ADDED] VIDEOS/DIAGNOSTICS have no upstream equivalent (upstream
// has no compile-time video toggle — its view3.h/draw3.h calls are
// unconditional — and no separate pressure-residual health-check output).
// `#ifndef` (rather than plain #define) so `make build-video`/`build-mpi-
// video` can override via `-DVIDEOS=1` without editing this file.
#ifndef VIDEOS
#define VIDEOS           0   // Videos: enable only for diagnostics; not needed in optimization loop
#endif
#ifndef DIAGNOSTICS
#define DIAGNOSTICS      0   // Write pressure_diag.dat (mgp.resa per step); health checks only, not production
#endif

// Output options
#define OUT_INTERFACE     1   // Save interface geometry



// ================================================================== //
//                       SIMULATION SETUP                             //
// ================================================================== //
int NN;  // Grid resolution: set from params.fidelity as 1<<fidelity in main()
// [PROJECT CHANGED] Upstream: `const double t_change = 30;` — a fixed 30
// PHYSICAL SECONDS before regular rocking motion is established (Main.tex
// line 430: "constant amplitude at the regular mode is attained after 30
// seconds"), i.e. Kim et al.'s own paper explicitly uses this criterion.
// Changed to a fixed cycle-count so the ramp scales with omega_b instead of
// being a fixed physical duration (commit 7ec98f9, "ramp duration = 3
// rocking cycles instead of hardcoded 30 s" -- rationale not recorded beyond
// the commit subject). NOTE (diary.md, 2026-07-28): because T_bio scales
// with 1/omega_b, a fixed-cycle-count ramp is proportionally MUCH shorter in
// real seconds at high RPM than at low RPM, whereas Kim's fixed-30s ramp is
// not. Tested whether skipping more analysis-window time reproduces the
// effect of a longer forcing ramp -- it does not (see diary.md) -- but that
// only tested the post-hoc window, not re-running with the actual 30s
// forcing profile itself. Still open.
// [PROJECT ADDED, 2026-09-07] Overridable via -DN_RAMP_CYCLES=<n> (mirrors
// the VIDEOS/DIAGNOSTICS pattern) so a slow-ramp test binary can be built
// without disturbing the default for every other binary -- direct test of
// the hysteresis/multiple-coexisting-branches hypothesis for the cold-
// start-beats-warm-start anomaly (diary.md 2026-09-07): if a quasi-static
// (long) ramp between conditions stays on the reference branch instead of
// jumping to the other one that a fast (3-cycle) ramp lands on, that's
// direct evidence for genuine bistability, not a residual bug.
#ifndef N_RAMP_CYCLES
#define N_RAMP_CYCLES 3            // Ramp duration in rocking cycles; t_change_st = N_RAMP_CYCLES * T_per_st
#endif
const double th_cont = 90;        // Contact angle for wetting conditions (degrees)
// Upstream's REMOVE_DROP constants, unchanged (see remove_drop event below).
const int    remove_minsize   = 20;      // minimum number of cells for a region to survive removal
const double remove_threshold = 1.0e-4;  // threshold for identifying disconnected fluid elements
double t_mix,t_dump;              // Time at which tracer is released, and dump file is saved (computed later)
const double nMix_cycle = 80;     // Number of cycles for tracer release (used to compute t_mix)
double t_end;                      // Final simulation time (simulation time unit); set from params.t_end in main()

// Output time intervals (derived from experimental timing)
// [PROJECT ADDED] `dt_video`/`t_out` have no upstream equivalent — Kim et
// al.'s driver has no shear-stress KPI pipeline at all (see the tau_95/98/
// 100/mean block in event normcal below, entirely new) and no video-frame
// dumping at a tunable interval.
// [PROJECT ADDED, 2026-07-30] frames_per_period wired from params.json so a
// one-off higher-cadence video render doesn't require editing/recompiling
// this file; default (5) preserves prior behavior exactly for existing
// production configs that don't set it.
double dt_video;  // Interval for video frames; set from T_per_st/params.frames_per_period in main()
// [PROJECT CHANGED, commit 1c3440c, 2026-07-14] was t_out=0.1 (~6.1 samples/
// rocking-period, T_per_nd~0.608 is RPM-independent by construction), raised
// to 0.02 (~30 samples/cycle) to at least match Kim et al. (2024)'s stated
// ~13 samples/cycle for their own extended-window turbulent-case sampling.
// CORRECTION (diary.md, 2026-07-28): the commit message claims this change
// alone moved tau_100_max from 0.37x to 1.03x Kim's value at 22.5 rpm/
// fidelity 9, citing experiments/hypothesis_ledger.json — that file was
// never actually updated with the entry, and the run pair that appears to
// be the actual A/B test
// (/oscar/scratch/eaguerov/mpi_runs/8994c04a.STALE_pre_tout_fix vs.
// _tout_test_22p5) shows NO meaningful difference when recomputed directly
// (0.371x vs. 0.372x Kim). That specific justification did not survive
// re-verification. t_out=0.02 is kept anyway (finer sampling is not
// expected to hurt, and it did fix a real aliasing problem in velocity
// time-series reconstruction — see diary.md), but do not cite the commit's
// stated tau_100_max improvement as evidence of anything.
const double t_out   = 0.02;      // Output interval for statistics [non-dim time]


// ================================================================== //
//                       REACTOR GEOMETRY                             //
// ================================================================== //
double LL = 1.0;        // Width of the reactor domain (non-dimensionalized by L_bio, so always 1)
double Ly = 0.286;      // Height of the reactor domain (non-dimensionalized by L_bio)
double y_init = 0.0;    // Initial liquid height (defines volume fraction f = 1 below this height)
double L_piv  = 0.143;  // Distance from center to pivot point (rocking axis), affects acceleration computation
double L_bio  = 0.25;   // Characteristic length (bag half-length, m); global so events can non-dimensionalize with it


// ================================================================== //
//                     MATERIAL PROPERTIES                            //
// ================================================================== //
const double rho_w = 1.0e3;   // Water density (kg/m^3)
const double rho_a = 1.225;   // Air density at 20°C (kg/m^3)
const double mu_w  = 1.0e-3;  // Water viscosity at 20°C (Pa·s)
const double mu_a  = 1.81e-5; // Air viscosity at 20°C (Pa·s)
const double grav  = 9.8;     // Gravitational acceleration (m/s^2)
const double sigma = 0.0728;  // Surface tension at the liquid-gas interface (N/m)

// Note: Phase 1 = liquid, Phase 0 = air
const double D_tracer_1 = 0.44e-9; // Diffusivity of tracers in water (m²/s)
const double D_tracer_2 = 1.0e-30; // Diffusivity of tracers in air (m²/s); near-zero in air
const double D_oxy_1 = 1.90e-9;    // Diffusivity of oxygen in water (m²/s)
const double D_oxy_2 = 1.98e-5;    // Diffusivity of oxygen in air (m²/s)
const double c_tracer_alpha = 1.0e30;    // Tracer solubility: c_water = c_tracer_alpha*c_air
const double c_oxy_alpha    = 1./30;     // Oxygen solubility: oxy_water = c_oxy_alpha*oxy_air (=1/30 tpiycal for O2 in water)

// ================================================================== //
//                 SCALARS AND FIELD POINTERS                         //
// ================================================================== //
// Parameters related to oscillation and tracer release
double Th_max, T_per, R_tr, x_tr, y_tr;

// Scalars to track tracer concentrations and oxygen
scalar c[], oxy[], c1[], c2[], c3[];   // for tracer and oxygen transfer
// Conditional stracers list: only include scalars that are active.
// Required for MPI correctness — unused scalars must not appear in stracers
// when their compile-time flags are off, or boundary reductions will
// reference uninitialised memory across MPI ranks.
// [PROJECT ADDED, 2026-09-15] c, c1 and c3 are the horizontal-left,
// horizontal-right and vertical-down tracer variants. Only ONE variant is
// ever initialised (VERTICAL_MIXUP -> c2), and postprocess reads only c2
// (dtmix) and oxy (kLa) -- so with the default flags c, c1 and c3 stay
// identically zero for the whole run while STILL costing, every timestep,
// two VOF-advected fields each (henry_oxy2.h event vof clones phi1/phi2 and
// appends both to f.tracers) plus a multigrid diffusion solve each
// (event tracer_diffusion). That is 3 of 5 tracers doing no work at full
// price. EXTRA_TRACERS=0 drops them. Default 1 preserves existing behaviour
// exactly; set 0 for mixing/kLa production runs.
#ifndef EXTRA_TRACERS
#define EXTRA_TRACERS 1
#endif
#if TRACER && OXYGEN && EXTRA_TRACERS
scalar * stracers = {c, oxy, c1, c2, c3};
#elif TRACER && OXYGEN
scalar * stracers = {oxy, c2};
#elif TRACER && EXTRA_TRACERS
scalar * stracers = {c, c1, c2, c3};
#elif TRACER
scalar * stracers = {c2};
#elif OXYGEN
scalar * stracers = {oxy};
#else
scalar * stracers = NULL;
#endif
double (* gradient) (double, double, double) = minmod2;   // Custom slope limiter function (used for scalar gradients)

// Buffers for file naming and output file pointers for statistics
char buf1[100], buf2[100], buf3[100], buf4[100];
FILE * fp_stats, * fp_norm, * fp_stats2, * fp_stats3, * fp_tau, * fp_edr;

// ── Steady-streaming vorticity, Kim's <|xi_bar_b'|> ───────────────────────
// [PROJECT ADDED, 2026-09-15] Kim's Fig 9/10/11/12 right-hand axis is the
// spatial mean of |vorticity OF THE TIME-AVERAGED FLOW| in the water. The
// overbar is INSIDE: average the velocity field over a rocking period FIRST,
// then take its curl. What this project reported as `vor_mean` was instead
// the time-average of the spatial mean of |instantaneous vorticity| -- a
// different quantity, dominated by the rocking oscillation rather than by the
// slow secondary circulation, and additionally normalised by the whole
// embedded fluid volume rather than the liquid. For an oscillatory flow the
// two can differ without limit (the mean flow can be ~0 while |omega| is
// large), so the old number was not comparable to Kim at all. See diary.md
// 2026-09-15 (12).
#ifndef STREAMING_VORTICITY
#define STREAMING_VORTICITY 1
#endif
#if REPLAY_FLUX                       // no velocity field in a replay
#undef STREAMING_VORTICITY
#define STREAMING_VORTICITY 0
#endif
#if STREAMING_VORTICITY
vector u_acc[];                       // running integral of u*f dt over one period
scalar f_acc[];                       // running integral of f dt (liquid residence)
double streaming_acc_dt = 0.;         // elapsed time in the current window
FILE * fp_streaming = NULL;
#endif

// Key physical and dimensionless parameters (computed in main)
double U0, Re_w, Re_a, We_w, Fr, rhor, mur, Pe_tracer_1, Pe_tracer_2, Pe_oxy_1, Pe_oxy_2, Th, Th_d, Th_2d, U_bio, w_bio, w_bio_st, T_per_st, T_bio, t_change_st;

// [PROJECT ADDED, 2026-09-07] Diagnostic-only mirror of the k=1 harmonic's
// ramped amplitude, updated every step inside event acceleration() so
// movies_output() (a separate event) can stamp the CURRENT (post-ramp-
// interpolation) envelope angle onto video frames -- not the same as Th,
// which is the instantaneous oscillating tilt, not the envelope. Used to
// visually verify checkpoint-restart chaining (diary.md 2026-09-07).
double _theta_env_deg = 0.0;
// [PROJECT ADDED] `params` holds everything read from params.json (see
// params_read.h); replaces the three scalar globals (L_bio/ANGLE/RPM) that
// upstream populated straight from argv.
BioreactorParams params;  // global so all events (acceleration, init) can access it
// [PROJECT ADDED] Checkpoint-restart bookkeeping. Upstream has no restart
// concept at all — every run is a fresh cold start. This project's sweeps
// are walltime-limited, so a run may need to continue across multiple SLURM
// jobs (or deliberately warm-start a new condition from a converged one);
// these three globals carry the state that decision needs. See main()'s
// [PROJECT ADDED] checkpoint-restart block below.
static double t_ramp_start      = 0.0;   // simulation time when the current ramp began
// [PROJECT ADDED, 2026-09-07] Added to the forcing's time argument so a
// restart resumes in phase with the state it restored. Nonzero only on a
// restart that changes omega_b or theta_max; see the derivation in main().
static double t_phase_offset    = 0.0;
static double t_dump_checkpoint = 0.0;   // simulation time to write checkpoint.dump

// [PROJECT ADDED, 2026-09-21] Graceful stop before the queueing system kills
// us. `event dump_checkpoint` fires exactly once, at t_end, so a run that
// runs out of clock before reaching it loses everything it has computed --
// measured 2026-09-21, when fig9_l9_rpm20 was 49% through 24 h of L9 with no
// checkpoint on disk and SLURM refused to raise TimeLimit on a RUNNING job.
//
// Modelled on Basilisk's own src/maxruntime.h: a deadline, not a signal.
// A deadline needs no async-signal-safe handler and does not depend on the
// queueing system delivering a signal to every MPI rank at a comparable
// timestep -- and the MPI_MAX reduction below makes all ranks agree on the
// same step for free. One rank entering the collective dump() alone would
// hang the job until the very kill this exists to pre-empt.
//
// The limit is declared through the environment rather than argv because
// argv[2] is already the restart file and the submit path is shared by every
// sweep script; the template exports it from SLURM's own accounting.
#define MAXRUNTIME_RESERVE_S 300   // time left for the dump: a 196 MB L10
                                   // checkpoint at 32 ranks writes in ~40 s,
                                   // and a truncated dump is worse than none
                                   // because it still looks resumable
static double _maxruntime_s   = HUGE;
static int    _emergency_dump = 0;   // set once, so we dump at most once

static const char * restart_file = NULL; // argv[2] if this is a restart run
// Upstream declares these for its `#if AMR` adaptive-refinement path (see
// upstream BioReactor.c:141,160-162,445-448). Upstream itself runs with
// `#define AMR 0` (uniform grid) for its published results, same as this
// fork — AMR was never actually exercised in either codebase's production
// runs. [PROJECT: the `#if AMR` refine() block itself was not carried over
// into this fork at all; these two variables are vestigial, kept only so
// nothing else referencing them fails to compile.]
int MINLEVEL, MAXLEVEL;    // Mesh refinement levels


// ================================================================== //
//                       MAIN FUNCTION                                //
// ================================================================== //

// [PROJECT FIX, 2026-09-30] The run stops at the SMALLEST period boundary at
// or after the requested end: n = max(1, ceil(x - PERIOD_TOL)), x = t/T. The
// old rule, (int)x + 1, gave k OR k+1 periods for a request of exactly k,
// depending on the last bit of x. Measured: 10 requested -> 11, 20 -> 20 on
// the same binary, which split test_restart_continue's two runs by 2 periods.
// scripts/periods.py implements the same expression, and chain.py/sweep.py use
// it to predict the dump time -- keep the two in step.
#define PERIOD_TOL 1e-5
static int next_period_count (double x) {
  int n = (int) ceil (x - PERIOD_TOL);
  return n < 1 ? 1 : n;
}

// ================================================================== //
#include "replay.h"

int main(int argc, char * argv[]){

  // [PROJECT CHANGED] Upstream: `double L_bio=atof(argv[1]); double ANGLE=
  // atof(argv[2]); double RPM=atof(argv[3]);` — three positional command-line
  // floats. Replaced with a single params.json (see params_read.h) because
  // this project's sweeps need many more knobs than upstream's three
  // (fidelity, geometry shape/fill, multi-harmonic forcing, checkpoint-restart
  // timing, MPI/video toggles, ...) and passing dozens of positional argv
  // floats does not scale.
  if (argc < 2) { fprintf(stderr, "Usage: BioReactor params.json\n"); return 1; }
  params = params_read(argv[1]);

  // Derive legacy scalars from params for the rest of main() unchanged.
  // omega_b (rad/s) and theta_max[0] (deg) replace the old ANGLE/RPM CLI args.
  // Multi-harmonic forcing and superellipse geometry will consume params directly
  // once the acceleration event and init event are extended (next steps).
  L_bio = params.geometry_a;  // bag half-length (m); sets all non-dim scales
  double ANGLE = params.theta_max[0];       // Fundamental rocking amplitude (degrees)
  double RPM   = params.omega_b * 60. / (2.*M_PI);  // Convert rad/s → RPM
  // [PROJECT CHANGED] Upstream: `const double NN = 64;` — a fixed resolution
  // (equivalent to fidelity 6 in this project's convention). Made runtime-
  // configurable because comparing accuracy/cost across resolutions is this
  // project's actual purpose (multi-fidelity optimization).
  NN    = 1 << params.fidelity;             // fidelity → grid cells per side (4→16, 7→128, 9→512)
  t_end = params.t_end;                    // non-dimensional sim time; 1 unit = T_bio seconds

  // [PROJECT CHANGED] Upstream: `L0 = LL;` with `LL` a plain (undimensioned)
  // constant =1.0. The `[0]` annotation is required by the current Basilisk
  // install's dimensional-analysis checking (see CLAUDE.md) — Kim et al.'s
  // upstream code predates that requirement, not a physics difference.
  L0 = 1. [0];  // [0] declares space dimensionless: simulation is fully non-dimensionalized (scaled by L_bio, T_bio, U_bio); Basilisk dimensional analysis requires annotations on literals, not variables
  // [PROJECT ADDED] Upstream never sets DT (uses Basilisk's default, unbounded).
  // DT: maximum allowed timestep (dimensionless).  HUGE is fine for serial/OpenMP
  // but MPI's global dtmax reduction can behave unexpectedly with infinity on
  // some implementations.  1.0 [0] (one non-dim time unit) is always larger than
  // the CFL-limited dt (~dx/U_max << 1) so it never artificially constrains the step.
  DT = 1. [0];
  origin(-L0/2., -L0/2.);   // Set coordinate origin to domain center

  init_grid(NN);

  // [PROJECT CHANGED] Upstream: `double Ly = 0.286;` — a fixed constant
  // (aspect ratio for the one geometry Kim et al. studied). Made runtime-
  // configurable (geometry.b/geometry.a) so the sweep/optimization pipeline
  // can vary bag shape as a parameter. At Kim's own values (a=0.25, b=0.071)
  // this recovers Ly=0.284, matching upstream's 0.286 to within rounding.
  // Bag geometry from params (must precede anything that uses Ly or H_bio)
  Ly = params.geometry_b / L_bio;  // dimensionless half-height; overwrites hardcoded 0.286

  // Rocking motion parameters
  R_tr  = 0.0084/L_bio;    // Tracer radius (scaled with domain)
  x_tr  = 0.;              // Tracer x-position (center)
  y_tr  = -Ly*0.5*0.5;     // y-position: halfway into the bottom half of domain
  Th_max = ANGLE*pi/180;   // Max angle in radians
  T_per  = 60./RPM;        // Period in seconds
  
  // Characteristic scales
  double H_bio,V_bio;
  // [PROJECT FIX, 2026-08-03, diary.md] `Ly` above is a HALF-height semi-axis
  // (matches its use in solid()/y_fill and docs_site/reference/params.md's
  // documented meaning of geometry.b), but H_bio must be the FULL bag height
  // to match upstream's convention (where the old fixed `Ly=0.286` constant
  // WAS the full height by construction). Missing the factor of 2 here made
  // every run's simulated bag exactly 2x taller than intended -- confirmed
  // by comparing ux_liq_vol (0.568, should be 0.286) against Kim's own
  // upstream driver recompiled and run at the same conditions, and by the
  // resulting u'_x,rms sitting at ~half of Kim's published Fig. A.16 value
  // while u'_y,rms matched. See diary.md 2026-08-03 for the full derivation.
  H_bio  = 2.*L_bio*Ly;
  V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
  U_bio  = V_bio/(H_bio*0.5)/T_per; // Characteristic velocity scale
  T_bio  = L_bio/U_bio;             // Characteristic time scale
  w_bio  = 2*pi/T_per;              // Angular velocity (rad/s)
  w_bio_st = w_bio*T_bio;           // Dimensionless angular velocity
  T_per_st = T_per/T_bio;           // Dimensionless period
  // [PROJECT FIXED, 2026-09-15, diary.md] Frame cadence must NOT divide the
  // rocking period exactly. It used to be dt_video = T_per_st/N, which makes
  // every frame land on the SAME N phases of the cycle forever -- measured on
  // real runs: a34fc4d4 (L9, 130 frames) and 57f68830 (L10, 50 frames) sample
  // only 5 distinct phases (0.0,0.2,0.4,0.6,0.8), so recording more frames
  // adds zero phase coverage and the true peak of any oscillating quantity is
  // never captured, only a max-over-5-samples lower bound. Every <tau> and
  // <eps> peak this project has reported from frames is therefore slightly
  // understated. (l8_coldstart_vid escaped only because its dt happened to
  // drift relative to its period.)
  //
  // Kim et al. hit the same trap and designed around it -- Main.tex: "a
  // constant time gap of 0.05 simulation time ... approximately 13 simulation
  // points per cycle but is intentionally misaligned with the period to
  // ensure convergence."
  //
  // Fix: offset the divisor by the golden ratio conjugate. phi = 0.618... is
  // the irrational that is hardest to approximate by rationals, so successive
  // frames land at maximally-spread phases and N frames over many periods
  // give near-uniform (low-discrepancy) phase coverage rather than N points.
  // Default raised 5 -> 13 to match Kim's cadence.
  {
    int nfp = params.frames_per_period > 0 ? params.frames_per_period : 13;
    dt_video = T_per_st / (nfp + 0.6180339887498949);
  }
  U0     = w_bio_st*Th_max;         // Initial rotational velocity
  t_change_st = N_RAMP_CYCLES * T_per_st;  // ramp over 3 cycles regardless of omega_b
  t_mix      = T_per_st*params.n_mix_cycles; // rocking cycles before tracer/oxygen start (wired from params.json)
  t_dump = t_mix;                   // Time to dump data (simulation time)

  // [PROJECT ADDED] Checkpoint restart — has no upstream equivalent at all;
  // Kim et al.'s driver is always a single cold-started run. Added because
  // this project's fidelity-9/10 sweeps need days per condition, longer than
  // one SLURM job's walltime allows, and because warm-starting a new
  // condition from an already-developed flow field (skipping its own
  // cold-start transient) makes chained sweeps far cheaper.
  // ── Checkpoint restart ─────────────────────────────────────────────────────
  // Detected via params.t_checkpoint > 0 (set by chain.py for restart segments).
  // params.t_end is a RELATIVE duration; the C code adds params.t_checkpoint.
  // All timing is computed HERE (before run()) so that the Basilisk event system
  // registers t_dump_checkpoint at the correct value — not as a placeholder that
  // gets overridden too late.  restore() itself is called in event init so that
  // the grid infrastructure is fully active.
  if (params.t_checkpoint > 0.0) {
    if (argc < 3) {
      fprintf (stderr, "ERROR: t_checkpoint > 0 in params but no dump file given (argv[2])\n");
      return 1;
    }
    restart_file     = argv[2];
    // [PROJECT FIXED, 2026-09-07] Restart forcing-phase error. Simulation
    // time is measured in units of T_bio = L_bio/U_bio, and U_bio depends on
    // BOTH omega_b and theta_max (see V_bio above), so the nondimensional
    // rocking period T_per_st differs between two segments that differ in
    // either. A checkpoint is written at an integer multiple of the SOURCE
    // segment's T_per_st (the n_per*T_per_st dump time below), i.e. at a
    // genuine zero-crossing of the run that wrote it -- but restore() hands
    // that t straight to a segment whose T_per_st is different, where it is
    // no longer a whole number of periods. The forcing
    // Th = Ak*sin(k*w_bio_st*t + phk) therefore resumes at an arbitrary
    // phase while the restored velocity and interface fields still
    // correspond to phase 0. This is the same wrong assumption the u/p/g
    // rescale in event init fixes ("U_bio is independent of theta_max"),
    // applied to the clock instead of to the fields.
    //
    // Measured phase errors into theta=7: source theta 6.9 gave -23.9deg,
    // 4.0 gave -53.7deg, 2.0 gave +116.6deg, and 7.0 gave exactly 0 --
    // which is the only reason trivial restarts ever looked clean. Because
    // the error is (period count)x(period ratio) mod 2pi it is pseudo-random
    // rather than proportional to Delta_theta_max, which is what made the
    // anomaly look like a step function: 0.1deg and 5deg steps both landed
    // ~25% high in quasi-steady tau_mean while 0deg was exact. Proven in
    // both directions (diary.md 2026-09-07 (4)): injecting a 30deg phase
    // error into an otherwise-clean restart with NO theta change reproduced
    // the escape (1.2596x, against 1.2547x for the real Delta_theta=0.1
    // case), and cancelling the phase error on the real Delta_theta=0.1
    // restart removed it entirely (0.9962x, against 0.9961x for a trivial
    // restart and 1.0000x for a cold start).
    //
    // The correction is a constant offset added to the forcing's time
    // argument, t_phase_offset, chosen so the forcing is at a zero-crossing
    // at the restart instant -- which is the phase the restored fields are
    // actually at. The scheduler's t is deliberately NOT touched: writing it
    // from event init desynchronizes the event system and dies with an FPE
    // in dtnext().
    //
    // The offset is derived from the restored t alone, NOT from the U_bio
    // ratio between the two segments. Deriving it from the ratio is the
    // obvious move and it is wrong for chains: a second segment that
    // changes nothing computes a ratio of 1 and therefore an offset of 0,
    // but its restored t is already displaced from a period boundary by the
    // FIRST segment's offset, so it resumes off-phase and escapes. Measured
    // on a theta 2 -> 7 -> 7 chain: segment 0 landed at 1.0003 and segment 1
    // at 1.2564, the full escape, one link down (diary.md 2026-09-07 (5)).
    //
    // Every checkpoint is written at a zero-crossing of ITS OWN writer's
    // forcing -- cold starts by n_per*T_per_st below, restarted segments by
    // the same expression shifted by their offset -- so the alignment
    // condition is simply w_bio_st*(t + t_phase_offset) == 2*pi*m, giving
    // the nearest-period-boundary form used here. That is provenance-free:
    // it needs nothing about the previous segment, self-heals over any
    // chain depth, and is bounded by T_per_st/2 rather than growing with t
    // the way the ratio form did (which reached 2.23 time units on the
    // theta 2 -> 7 case). For a single segment it reproduces the ratio
    // form's value exactly -- theta 6.9 -> 7 gives +23.86deg either way,
    // the value independently verified by the phi_angular cancellation
    // test.
    {
      // The exact dump time, read from the checkpoint itself. params
      // .t_checkpoint is chain.py's estimate from the source's last
      // shear_stress.dat sample and can be off by an output interval
      // (15.74 vs 15.750302 for the theta=6.9 source), which here would
      // leave a residual phase error of a few degrees.
      double t_ck = params.t_checkpoint;
      FILE *fp_ck = fopen (restart_file, "r");
      if (fp_ck) {
        double t_dumped;
        if (fread (&t_dumped, sizeof(double), 1, fp_ck) == 1 && t_dumped > 0.)
          t_ck = t_dumped;
        fclose (fp_ck);
      }
      // Nearest-boundary inference. Valid only for a checkpoint written
      // AT a zero-crossing, which every planned dump is; an emergency dump
      // is not, and says so in checkpoint.phase, read just below.
      t_phase_offset      = T_per_st*round (t_ck/T_per_st) - t_ck;
      params.t_checkpoint = t_ck;

      // [PROJECT ADDED, 2026-09-21] A recorded phase beats an inferred one.
      // Absent (every checkpoint written before 2026-09-21) we keep the
      // inference above, so old dumps restore exactly as they always did.
      {
        FILE *fp_ph = fopen ("checkpoint.phase", "r");
        if (fp_ph) {
          char buf[512];
          if (fgets (buf, sizeof buf, fp_ph)) {
            const char *k = strstr (buf, "\"t_phase_offset\"");
            const char *c = k ? strchr (k, ':') : NULL;
            if (c) {
              double rec = atof (c + 1);
              // The writer's forcing argument was w*(t + rec). Reproducing it
              // exactly is the whole point, so take it verbatim rather than
              // re-deriving anything from t_ck.
              t_phase_offset = rec;
              if (pid() == 0)
                fprintf (ferr, "restart: using recorded phase offset "
                         "%.6f from checkpoint.phase\n", rec);
            }
          }
          fclose (fp_ph);
        }
      }
    }
    // Smooth-step interpolation starts AT the checkpoint and runs N_RAMP_CYCLES forward.
    // alpha goes 0→1 over [t_checkpoint, t_checkpoint + N_RAMP_CYCLES*T_per_st].
    t_ramp_start     = params.t_checkpoint;
    // [OPEN, 2026-09-16] A split run does not exactly reproduce an unbroken
    // one: measured at fidelity 4, oxygen transfer over a 3-cycle continuation
    // segment lands ~7.6% below a continuous 6-cycle run. Two mechanisms were
    // hypothesised and BOTH were falsified by measurement:
    //   (a) oxygen replenishment pausing, because `event oxygen (t=t_mix; i++)`
    //       has t_mix recomputed relative to the new checkpoint. Firing it
    //       immediately (n_mix_cycles=0) gave 0.9238 vs 0.9236 -- no change.
    //       It accounts for only ~1.5% (n_mix=8, never firing, gives 0.9100).
    //   (b) the smooth-step ramp. A continuation sets no *_prev fields, so
    //       alpha interpolates from ZERO forcing over N_RAMP_CYCLES, which
    //       looks like it must underdrive the segment. Backdating t_ramp_start
    //       so alpha=1 from the first step made it WORSE, 0.9236 -> 0.9100.
    // Removing an under-driving ramp reducing oxygen transfer is not
    // understood, so no change is being shipped on the strength of it. The
    // guard in tests/verification/test_restart_continue.py is set at 10% so it
    // still catches total field loss and re-injection, which is its purpose;
    // this ~8% is a separate, smaller, and currently UNEXPLAINED restart
    // artifact. restart_continue=1 has not yet been used in production.
    t_mix            = params.t_checkpoint + T_per_st * params.n_mix_cycles;
    t_dump           = t_mix;
    {
      // The dump must land on a zero-crossing of THIS segment's forcing,
      // which t_phase_offset has shifted: the forcing argument is
      // w_bio_st*(t + t_phase_offset), so a whole number of periods means
      // t + t_phase_offset == n*T_per_st, not t == n*T_per_st. Without
      // this the fix would hand the next segment in a chain a checkpoint
      // that is itself off-phase, reintroducing the same bug one link down.
      double t_end_abs = params.t_checkpoint + params.t_end;
      int n_per        = next_period_count ((t_end_abs + t_phase_offset) / T_per_st);
      t_dump_checkpoint = n_per * T_per_st - t_phase_offset;
      t_end            = t_dump_checkpoint;
    }
  } else {
    // Fresh run: extend t_end to the next period boundary for clean phase alignment.
    int n_per = next_period_count (t_end / T_per_st);
    t_dump_checkpoint = n_per * T_per_st;
    t_end = t_dump_checkpoint;
  }

#if REPLAY_FLUX
  // A replay continues the recorded run from its start checkpoint; params.t_end is the duration to replay.
  // No forcing is evaluated, so there is no period boundary to round to (and no checkpoint to write).
  if (!restart_file) {
    fprintf (stderr, "ERROR: REPLAY_FLUX needs t_checkpoint > 0 in params and the start checkpoint as argv[2]\n");
    return 1;
  }
  t_end = params.t_checkpoint + params.t_end;
  t_dump_checkpoint = t_end;
#endif
#if RECORD_FLUX
  // Record window: from t_rec0 (default: the release time t_mix) for n_rec_periods periods (0 = to the end).
  t_rec0    = params.t_rec0 >= 0. ? params.t_rec0 : t_mix;
  t_rec_end = params.n_rec_periods > 0 ? t_rec0 + params.n_rec_periods*T_per_st : 1e30;
#endif

  // [PROJECT ADDED, 2026-09-21] Walltime deadline, if the queueing system
  // told us one. Absent (interactive runs, tests that do not set it) the
  // limit stays HUGE and the deadline event never fires, so behaviour is
  // exactly as before.
  {
    const char *s = getenv ("BIOREACTOR_MAXRUNTIME_S");
    if (s && *s) {
      double v = atof (s);
      if (v > 0.) {
        _maxruntime_s = v;
        if (pid() == 0)
          fprintf (ferr, "maxruntime: will checkpoint and stop %d s before "
                   "%.0f s of wall clock\n", MAXRUNTIME_RESERVE_S, v);
      }
    }
  }

  // Dimensionless numbers
  Re_w = rho_w*U_bio*L_bio/mu_w;            // Reynolds number of water
  Re_a = rho_a*U_bio*L_bio/mu_a;            // Reynolds number of air
  We_w = rho_w*U_bio*U_bio*L_bio/sigma;     // Weber number of water
  Fr   = U_bio/sqrt(grav*L_bio);            // Froude number
  rhor = rho_a/rho_w;                       // Density ratio
  mur  = mu_a/mu_w;                         // Viscosity ratio
  Pe_tracer_1 = U_bio*L_bio/(D_tracer_1*params.tracer_D_scale);    // Peclet number of tracer at phase 1 (water)
  Pe_tracer_2 = U_bio*L_bio/D_tracer_2;     // Peclet number of tracer at phase 2 (air)
  Pe_oxy_1    = U_bio*L_bio/(D_oxy_1*params.oxy_D_scale);       // Peclet number of oxygen at phase 1 (water)
  Pe_oxy_2    = U_bio*L_bio/D_oxy_2;        // Peclet number of oxygen at phase 1 (air) 

#if !REPLAY_FLUX
  rho1 = 1.0;             // Reference density of phase 1 (water) (scaled to 1)
  rho2 = rho1*rhor;       // Scaled density of phase 2 (air)
  mu1  = 1.0/Re_w;        // Dimensionless viscosity of phase1 (water)
  mu2  = mur*mu1;         // Scaled viscosity of phase 2 (air)
  f.sigma = 1.0 / We_w;   // Dimensionless surface tension
#endif
  
  // Tracer field configuration
  #if TRACER  
  // Define diffusivities in each phase (1 = liquid, 0 = air)
    c.D1 = 1./Pe_tracer_1;
    c.D2 = 1./Pe_tracer_2;
  
  // Tracer concentration ratio across interface: c_liquid = alpha*c_air; c1=alpha*c2
    c.alpha = c_tracer_alpha;

  // Apply slope limiter to gradients for numerical stability
    c.gradient = minmod2;

  // --- Optional extra tracers for comparisons --- //
    c1.D1 = 1./Pe_tracer_1;
    c1.D2 = 1./Pe_tracer_2;
    c1.alpha = c_tracer_alpha;
    c1.gradient = minmod2;
    c2.D1 = 1./Pe_tracer_1;
    c2.D2 = 1./Pe_tracer_2;
    c2.alpha = c_tracer_alpha;
    c2.gradient = minmod2;
    c3.D1 = 1./Pe_tracer_1;
    c3.D2 = 1./Pe_tracer_2;
    c3.alpha = c_tracer_alpha;
    c3.gradient = minmod2;
#endif

// Oxygen field configuration
#if OXYGEN
  oxy.D1 = 1./Pe_oxy_1;     // Oxygen diffusivities (water) 
  oxy.D2 = 1./Pe_oxy_2;     // Oxygen diffusivities (air) 
  oxy.alpha = c_oxy_alpha;  // Henry's law partition coefficient

  // Apply slope limiter for smooth gradients
  oxy.gradient = minmod2;
#endif

// Boundary conditions
// [VERIFIED UNCHANGED vs. upstream, diary.md 2026-07-28] Line-for-line
// identical to upstream's BC block, including the embed BCs below. Checked
// specifically because it looked, from the geometry formula alone, like this
// fork might be embedding all four walls where Kim only embeds two — it
// isn't; see the [PROJECT CHANGED] note at the superellipse solid() call.
#if !REPLAY_FLUX    // a replay has no velocity field
  u.n[left]  = dirichlet(0.);  // Set no-slip velocity boundary conditions (zero velocity)
  u.t[left]  = dirichlet(0.);
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);

// Solid geometry Boundary (if using embedded solids)
#if EMBED
  u.n[embed] = dirichlet(0.); // Set no-slip velocity boundary conditions (zero velocity)
  u.t[embed] = dirichlet(0.);
#endif
#endif // !REPLAY_FLUX

// Output Files
  char name[200],name2[200],name3[200],name4[200],name5[200];

  sprintf(name, "logstats.dat");          // Performance & time log
  sprintf(name2,"normf.dat");             // Norms and statistics of fields
  sprintf(name3,"vol_frac_interf.dat");   // Volume fraction & interface stats
  sprintf(name4,"tr_oxy.dat");            // Tracer and oxygen stats
  sprintf(name5,"shear_stress.dat");      // 98th-percentile absolute shear stress
  fp_stats = fopen(name, "w");
  fp_norm  = fopen(name2,"w");
  fp_stats2= fopen(name3,"w");
  fp_stats3= fopen(name4,"w");
  fp_tau   = fopen(name5,"w");
  // [PROJECT ADDED, 2026-11-04] edr_95 goes to its own file, not to
  // shear_stress.dat: that file must stay byte-identical to the pre-change build.
  fp_edr   = fopen("edr_pct.dat","w");
  fprintf(fp_edr, "i t edr_95 edr_max \n");
#if STREAMING_VORTICITY
  fp_streaming = fopen("streaming.dat","w");
  fprintf(fp_streaming, "i t vor_stream_mean vor_stream_max liq_vol window_dt \n");
#endif

  fprintf(fp_norm, "i t Omega_liq_avg Omega_liq_rms Omega_liq_vol Omega_liq_max ux_liq_avg ux_liq_rms ux_liq_vol ux_liq_max uy_liq_avg uy_liq_rms uy_liq_vol uy_liq_max ux_liq_savg uy_liq_savg uz_liq_savg uz_liq_rms \n");
  fprintf(fp_stats2, "i t f_liq_sum f_liq_interf posY_max posY_min posY_left posY_right \n");
  fprintf(fp_stats3, "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 c1_liq_sum c1_liq_sum2 c2_liq_sum c2_liq_sum2 c3_liq_sum c3_liq_sum2 b0 b1 b2 \n");
  fprintf(fp_tau,   "i t tau_95 tau_98 tau_100 tau_mean tau_100_strict tau_mean_strict tau_100_signed ediss_mean tau_mean_signed tau_kim_max tau_kim_mean ediss_kim_max ediss_kim_mean \n");

  NITERMAX = 1000;     // Max iterations per timestep
  TOLERANCE = 5.0e-4;  // // Solver tolerance (convergence criterion)
  
  // Run the simulation — Basilisk manages events from here
  run();
  
  // Close all output files
  fclose(fp_stats); fclose(fp_norm); fclose(fp_stats2); fclose(fp_stats3); fclose(fp_tau); fclose(fp_edr);
#if STREAMING_VORTICITY
  if (fp_streaming) fclose(fp_streaming);
#endif
#if RECORD_FLUX || REPLAY_FLUX
  if (fp_rec) fclose(fp_rec);
#endif
}

#include "restart.h"

#include "init.h"

#include "scalars.h"

#include "forcing.h"

#include "snapshots.h"

#include "checkpoint.h"

#include "diagnostics.h"

#include "replay_events.h"

#include "movies.h"