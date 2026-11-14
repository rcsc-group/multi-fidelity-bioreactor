// BioReactor.c -- two-phase flow (water + air) in a rocking cell-culture bag, with oxygen transfer and a
// passive tracer. Basilisk C (http://basilisk.fr), 2D by default.
//
// Origin: Minki Kim, Radu Cimpeanu and Daniel Harris, driver code of Kim et al. (2025),
// https://github.com/rcsc-group/BioReactor (DriverCodes). This project extends it: a params.json input,
// any grid level, checkpoint restart, more outputs (shear stress, dissipation, mixing), and a flux
// record/replay mode. The history of each change is in internal/docs/history/.
//
// How to read this code
//   1. This file: compile switches, global state, and main(), which only sets up the run.
//   2. The events do the work. They live in the headers included at the END of this file, one header
//      per concern. Basilisk runs the events of one step in the order it meets them in the source, so
//      the include order at the end of this file IS the event order. Do not reorder it.
//   3. src/README.md: event order table, compile switches, output files.
//   4. docs_site/developer/solver-pitfalls.md: Basilisk traps that cost us days. Read before changing
//      anything that touches restart, face fields or events.
//
// Units: everything is non-dimensional (length L_bio, time T_bio, velocity U_bio; see set_scales()).
//
// Build: `make build` (see the Makefile for the other variants). Run: `BioReactor params.json [restart.dump]`.

// ---------------------------------------------------------------------------------------------------
// Record/replay modes (both off by default). Oxygen and the tracer do not change density or viscosity,
// so they do not change the flow. That lets us store the flow once and re-use it:
//   -DRECORD_FLUX  normal coupled run; also writes the face flux uf of every step in the record window,
//                  and a start checkpoint replay_start.dump.
//   -DREPLAY_FLUX  no Navier-Stokes. Restores replay_start.dump, then reads uf step by step and advects
//                  only the volume fraction f and the scalars. Same code path for f and the scalars as the
//                  coupled run, so the result matches it (see replay.h for what is and is not reproduced).
// ---------------------------------------------------------------------------------------------------
#ifndef RECORD_FLUX
#define RECORD_FLUX 0
#endif
#ifndef REPLAY_FLUX
#define REPLAY_FLUX 0
#endif
#if RECORD_FLUX && REPLAY_FLUX
#error "RECORD_FLUX and REPLAY_FLUX are mutually exclusive"
#endif

// ---------------------------------------------------------------------------------------------------
// Basilisk modules
// ---------------------------------------------------------------------------------------------------
#include "embed.h"                 // the bag wall is an embedded boundary
#if REPLAY_FLUX
// Frozen-flow stand-in for centered.h + two-phase.h + conserving.h: only what the VOF advection of f,
// its tracers and the scalar diffusion need. `init` is declared first so that restore() runs before
// the replay events.
#include "run.h"                   // time loop, utils.h (DT, statsf, perf, minmod2)
#include "timestep.h"
#include "output.h"                // restore
event init (i = 0);
// Same declarations as centered.h: events marked `last` run after all user events of a step, so in both
// modes the user events run at the START of a step, before vof / diffusion.
event vof (i++,last);
event tracer_diffusion (i++,last);
face vector uf[];                  // read from the record, one step at a time
scalar f[], * interfaces = {f};
#include "vof.h"
#include "curvature.h"             // position(), for the interface height outputs
#include "poisson.h"               // mg_solve, for tracer_diffusion in henry_oxy2.h
#else
#include "navier-stokes/centered.h"
#endif

// Viscosity of a mixed cell: harmonic mean of the two phases, weighted by f.
#define mu(f) ( 1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)  )

#if !REPLAY_FLUX
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"   // momentum-conserving advection (large density ratio)
#endif

#include "henry_oxy2.h"            // two-phase scalars with a Henry's-law jump at the interface
#include "utils2.h"                // statsf2: statistics over the liquid only
#include "tag.h"                   // remove_droplets(), for params.remove_drop

#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

#include "params_read.h"           // params.json -> BioreactorParams

// ---------------------------------------------------------------------------------------------------
// Compile switches. Each #ifndef switch can be set with -DNAME=value (see the Makefile targets).
// ---------------------------------------------------------------------------------------------------
#define EMBED            1   // embedded boundary for the bag wall
#define OXYGEN           1   // oxygen field
#define OXYGEN_AIR       1   // oxygen starts in the air (the gas is the reservoir)
#define TRACER           1   // passive tracer fields
#define VERTICAL_MIXUP   1   // tracer c2 starts in the top half of the liquid
#define ACCELERATION     1   // rocking motion
#define NORMCAL          1   // per-step statistics (event normcal)
#if REPLAY_FLUX
#undef ACCELERATION
#define ACCELERATION     0   // no forcing: the recorded flux is the flow
#undef NORMCAL
#define NORMCAL          0   // replay.h has its own event normcal
#endif
#ifndef VIDEOS
#define VIDEOS           0   // write video frames (movies.h); diagnostics only
#endif
#ifndef DIAGNOSTICS
#define DIAGNOSTICS      0   // write pressure_diag.dat (multigrid residual per step); health checks only
#endif

// ---------------------------------------------------------------------------------------------------
// Run control
// ---------------------------------------------------------------------------------------------------
int NN;                           // grid cells per side = 2^params.fidelity

// The forcing amplitude ramps up smoothly over N_RAMP_CYCLES rocking periods (Kim et al. used a fixed
// 30 s instead). Override with -DN_RAMP_CYCLES=n.
#ifndef N_RAMP_CYCLES
#define N_RAMP_CYCLES 3
#endif
// Droplet removal (params.remove_drop): regions smaller than remove_minsize cells are removed.
const int    remove_minsize   = 20;
const double remove_threshold = 1.0e-4;
double t_mix;                     // release time of tracer and oxygen = n_mix_cycles periods
double t_end;                     // end time; rounded up to a period boundary in set_run_window()
double dt_video;                  // video frame interval; see set_scales() for why it is not T/N
const double t_out = 0.02;        // statistics output interval (about 30 samples per rocking period)

// ---------------------------------------------------------------------------------------------------
// Geometry (non-dimensional, length unit L_bio)
// ---------------------------------------------------------------------------------------------------
double LL = 1.0;        // domain width
double Ly = 0.286;      // bag half-height; overwritten from params.geometry
double y_init = 0.0;    // initial liquid height
double L_piv  = 0.143;  // distance from the centre to the rocking axis
double L_bio  = 0.25;   // bag half-length in metres (the length unit)

// ---------------------------------------------------------------------------------------------------
// Material properties (SI). Phase 1 = water, phase 2 = air.
// ---------------------------------------------------------------------------------------------------
const double rho_w = 1.0e3;   // water density (kg/m^3)
const double rho_a = 1.225;   // air density at 20 C (kg/m^3)
const double mu_w  = 1.0e-3;  // water viscosity at 20 C (Pa s)
const double mu_a  = 1.81e-5; // air viscosity at 20 C (Pa s)
const double grav  = 9.8;     // gravity (m/s^2)
const double sigma = 0.0728;  // surface tension (N/m)

const double D_tracer_1 = 0.44e-9; // tracer diffusivity in water (m^2/s)
const double D_tracer_2 = 1.0e-30; // tracer diffusivity in air: none
const double D_oxy_1 = 1.90e-9;    // oxygen diffusivity in water (m^2/s)
const double D_oxy_2 = 1.98e-5;    // oxygen diffusivity in air (m^2/s)
const double c_tracer_alpha = 1.0e30;    // tracer partition c_water/c_air: stays in the water
const double c_oxy_alpha    = 1./30;     // oxygen partition (Henry): oxy_water = oxy_air/30

// ---------------------------------------------------------------------------------------------------
// Fields
// ---------------------------------------------------------------------------------------------------
double Th_max, T_per, R_tr, x_tr, y_tr;   // max angle (rad), period (s), tracer patch radius and centre

// Tracers c, c1, c2, c3 are four initial shapes (line, left half, top half, circle) advected by the same
// flow; postprocessing reads only c2 (mixing time) and oxy (kLa). -DEXTRA_TRACERS=0 drops c, c1, c3,
// which saves about 40% of the run time. Every scalar in stracers must be initialised (MPI boundaries
// read all of them).
scalar c[], oxy[], c1[], c2[], c3[];
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
double (* gradient) (double, double, double) = minmod2;   // slope limiter for the scalars

FILE * fp_stats, * fp_norm, * fp_stats2, * fp_stats3, * fp_tau, * fp_edr;   // output files, see open_outputs()

// Steady-streaming vorticity (Kim et al., Figs. 9-12): average u over one rocking period FIRST, then take
// the curl, then average |curl| over the water. (The time average of |instantaneous vorticity| is a
// different quantity, dominated by the rocking itself.) Computed in forcing.h.
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

// Scales and dimensionless numbers, set in set_scales()
double U0, Re_w, Re_a, We_w, Fr, rhor, mur, Pe_tracer_1, Pe_tracer_2, Pe_oxy_1, Pe_oxy_2, Th, Th_d, Th_2d, U_bio, w_bio, w_bio_st, T_per_st, T_bio, t_change_st;

double _theta_env_deg = 0.0;   // current ramp envelope of the angle (deg), stamped on video frames
BioreactorParams params;       // everything read from params.json; global so events can read it

// Restart state. See set_run_window() and docs_site/explanation/checkpoint-restart.md.
static double t_ramp_start      = 0.0;   // time at which the current forcing ramp began
static double t_dump_checkpoint = 0.0;   // time at which checkpoint.dump is written

// Wall-clock deadline: if the queue says how long we have (env BIOREACTOR_MAXRUNTIME_S), write a
// checkpoint and stop MAXRUNTIME_RESERVE_S before the end. A deadline (as in Basilisk's maxruntime.h),
// not a signal: all MPI ranks agree on the step through one reduction. See checkpoint.h.
#define MAXRUNTIME_RESERVE_S 300   // a 196 MB L10 checkpoint on 32 ranks writes in ~40 s
static double _maxruntime_s   = HUGE;
static int    _emergency_dump = 0;   // set once: dump at most once

static const char * restart_file = NULL; // argv[2] if this is a restart run

// Number of whole rocking periods to run for a requested length of x periods: the smallest n >= x, at
// least 1. PERIOD_TOL absorbs round-off, so a request of exactly k periods gives k, never k+1.
// scripts/periods.py uses the same rule to predict the dump time: keep the two in step.
#define PERIOD_TOL 1e-5
static int next_period_count (double x) {
  int n = (int) ceil (x - PERIOD_TOL);
  return n < 1 ? 1 : n;
}

#include "forcing_state.h"   // forcing at any time, the checkpoint label, the source of a restart's ramp
#include "replay.h"    // record/replay state and functions (empty unless -DRECORD_FLUX or -DREPLAY_FLUX)

// ===================================================================================================
// Setup steps called by main(), in this order
// ===================================================================================================

// Grid, non-dimensional scales, ramp and release times.
static void set_scales (void)
{
  L_bio = params.geometry_a;                         // bag half-length (m): the length unit
  double ANGLE = params.theta_max[0];                // rocking amplitude (deg)
  double RPM   = params.omega_b * 60. / (2.*M_PI);   // rocking rate (rpm)
  NN    = 1 << params.fidelity;                      // level 7 -> 128 cells, 9 -> 512
  t_end = params.t_end;                              // requested length, non-dimensional

  L0 = 1. [0];  // [0]: space is dimensionless (dimensional analysis needs literals annotated)
  // Upper bound on the step. HUGE works in serial, but the MPI reduction of dtmax can misbehave with
  // infinity; one time unit is far above the CFL step, so it never limits it.
  DT = 1. [0];
  origin(-L0/2., -L0/2.);

  init_grid(NN);

  Ly = params.geometry_b / L_bio;  // bag half-height

  R_tr  = 0.0084/L_bio;    // tracer patch radius
  x_tr  = 0.;
  y_tr  = -Ly*0.5*0.5;     // halfway into the bottom half
  Th_max = ANGLE*pi/180;
  T_per  = 60./RPM;

  // Velocity scale: the liquid volume displaced per half period, over half the bag height.
  // H_bio is the FULL height (Ly is a half-height).
  double H_bio,V_bio;
  H_bio  = 2.*L_bio*Ly;
  V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
  U_bio  = V_bio/(H_bio*0.5)/T_per;
  T_bio  = L_bio/U_bio;
  w_bio  = 2*pi/T_per;              // rad/s
  w_bio_st = w_bio*T_bio;           // non-dimensional angular rate
  T_per_st = T_per/T_bio;           // non-dimensional period
  // Video frames must NOT divide the period exactly: with dt = T/N every frame falls on the same N phases
  // and peaks are missed. Dividing by N + 0.618 (golden ratio) spreads the phases evenly over many periods.
  {
    int nfp = params.frames_per_period > 0 ? params.frames_per_period : 13;
    dt_video = T_per_st / (nfp + 0.6180339887498949);
  }
  U0     = w_bio_st*Th_max;
  t_change_st = N_RAMP_CYCLES * T_per_st;    // ramp length
  t_mix      = T_per_st*params.n_mix_cycles; // release time of tracer and oxygen
}

// Start, release and end times. A fresh run ends on a period boundary (bag level). A restart (params
// .t_checkpoint > 0, dump file in argv[2]) continues from the dump; params.t_end is then a duration.
static void set_run_window (int argc, char * argv[])
{
  if (params.t_checkpoint > 0.0) {
    if (argc < 3) {
      fprintf (stderr, "ERROR: t_checkpoint > 0 in params but no dump file given (argv[2])\n");
      exit (1);
    }
    restart_file     = argv[2];
    // Exact dump time, from the first double of the dump file (params.t_checkpoint is an estimate).
    {
      FILE *fp_ck = fopen (restart_file, "r");
      if (fp_ck) {
        double t_dumped;
        if (fread (&t_dumped, sizeof(double), 1, fp_ck) == 1 && t_dumped > 0.)
          params.t_checkpoint = t_dumped;
        fclose (fp_ck);
      }
    }
    // The ramp starts at the checkpoint, from the forcing that made the restored state (the dump's label,
    // checkpoint.phase in the working directory) to the forcing in params.json: amplitude, phase and
    // frequency (forcing_state.h). Continuation = same forcing = no ramp. No label, or a label plus *_prev
    // keys, stops the run. Pitfall "restart-ramp" and pitfall "restart-phase" in docs_site/developer/solver-pitfalls.md.
    t_ramp_start = params.t_checkpoint;
#if REPLAY_FLUX
    // a replay has no forcing: keep the target, phase from the nearest period boundary
    fresh_source ();
    fsrc.cycles = round (params.t_checkpoint/T_per_st);
#else
    {
      int kind = restart_source ();
      if (kind < 0)
        return;      // refused: main() returns 1
      if (pid() == 0)
        fprintf (ferr, "restart: ramp from omega_b %.6g rad/s, theta_max %.6g deg (%s) to omega_b %.6g, "
                 "theta_max %.6g over %d periods; phase %.6f periods\n", fsrc.omega_b, fsrc.theta_max[0],
                 kind == 1 ? "label" : "old label + *_prev", params.omega_b, params.theta_max[0],
                 N_RAMP_CYCLES, fsrc.cycles);
      _restart_kind = kind;   // restart_forcing.json: first acceleration step (forcing.h)
    }
#endif
    // Release n_mix_cycles forcing periods after the restart (0: at the restart).
    t_mix = params.n_mix_cycles > 0 ?
      t_at_psi (2.*pi*(fsrc.cycles + params.n_mix_cycles)) : params.t_checkpoint;
    {
      // The dump must be on a zero crossing of THIS run's forcing phase psi; otherwise the next link of a
      // chain inherits an off-phase checkpoint.
      double t_end_abs = params.t_checkpoint + params.t_end;
      double rate;
      int n_per        = next_period_count (psi_at (t_end_abs, &rate)/(2.*pi));
      t_dump_checkpoint = t_at_psi (2.*pi*n_per);
      t_end            = t_dump_checkpoint;
    }
  } else {
    // Fresh run: end on the next period boundary.
    fresh_source ();
    int n_per = next_period_count (t_end / T_per_st);
    t_dump_checkpoint = n_per * T_per_st;
    t_end = t_dump_checkpoint;
  }

#if REPLAY_FLUX
  // A replay continues the recorded run from its start checkpoint; params.t_end is the duration to replay.
  // There is no forcing, so no period boundary to round to.
  if (!restart_file) {
    fprintf (stderr, "ERROR: REPLAY_FLUX needs t_checkpoint > 0 in params and the start checkpoint as argv[2]\n");
    exit (1);
  }
  t_end = params.t_checkpoint + params.t_end;
  t_dump_checkpoint = t_end;
#endif
#if RECORD_FLUX
  // Record window: from t_rec0 (default: the release time t_mix) for n_rec_periods periods (0 = to the end).
  t_rec0    = params.t_rec0 >= 0. ? params.t_rec0 : t_mix;
  t_rec_end = params.n_rec_periods > 0 ? t_rec0 + params.n_rec_periods*T_per_st : 1e30;
#endif

  // Wall-clock deadline, if the job script exported one.
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
}

// Dimensionless numbers, fluid properties, scalar diffusivities and boundary conditions.
static void set_physics (void)
{
  Re_w = rho_w*U_bio*L_bio/mu_w;            // Reynolds number, water
  Re_a = rho_a*U_bio*L_bio/mu_a;            // Reynolds number, air
  We_w = rho_w*U_bio*U_bio*L_bio/sigma;     // Weber number
  Fr   = U_bio/sqrt(grav*L_bio);            // Froude number
  rhor = rho_a/rho_w;                       // density ratio
  mur  = mu_a/mu_w;                         // viscosity ratio
  Pe_tracer_1 = U_bio*L_bio/(D_tracer_1*params.tracer_D_scale);    // Peclet numbers: tracer in water,
  Pe_tracer_2 = U_bio*L_bio/D_tracer_2;                            //   tracer in air,
  Pe_oxy_1    = U_bio*L_bio/(D_oxy_1*params.oxy_D_scale);          //   oxygen in water,
  Pe_oxy_2    = U_bio*L_bio/D_oxy_2;                               //   oxygen in air

#if !REPLAY_FLUX
  rho1 = 1.0;             // water density = 1
  rho2 = rho1*rhor;
  mu1  = 1.0/Re_w;
  mu2  = mur*mu1;
  f.sigma = 1.0 / We_w;
#endif

#if TRACER
  // D1, D2: diffusivity in water and air. alpha: partition c_water = alpha*c_air.
    c.D1 = 1./Pe_tracer_1;
    c.D2 = 1./Pe_tracer_2;
    c.alpha = c_tracer_alpha;
    c.gradient = minmod2;
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

#if OXYGEN
  oxy.D1 = 1./Pe_oxy_1;
  oxy.D2 = 1./Pe_oxy_2;
  oxy.alpha = c_oxy_alpha;  // Henry's law
  oxy.gradient = minmod2;
#endif

  // No slip on the box and on the bag wall (identical to Kim et al.).
#if !REPLAY_FLUX    // a replay has no velocity field
  u.n[left]  = dirichlet(0.);
  u.t[left]  = dirichlet(0.);
  u.n[right] = dirichlet(0.);
  u.t[right] = dirichlet(0.);
  u.n[top] = dirichlet(0.);
  u.t[top] = dirichlet(0.);
  u.n[bottom] = dirichlet(0.);
  u.t[bottom] = dirichlet(0.);
#if EMBED
  u.n[embed] = dirichlet(0.);
  u.t[embed] = dirichlet(0.);
#endif
#endif // !REPLAY_FLUX
}

// Time-series output files and their header lines. Columns: docs_site/reference/output-files.md.
static void open_outputs (void)
{
  fp_stats = fopen("logstats.dat", "w");          // step, time, wall clock
  fp_norm  = fopen("normf.dat","w");              // velocity and vorticity statistics in the liquid
  fp_stats2= fopen("vol_frac_interf.dat","w");    // liquid volume and interface heights
  fp_stats3= fopen("tr_oxy.dat","w");             // tracer and oxygen sums
  fp_tau   = fopen("shear_stress.dat","w");       // shear stress and dissipation
  fp_edr   = fopen("edr_pct.dat","w");            // dissipation percentiles
  fprintf(fp_edr, "i t edr_95 edr_max \n");
#if STREAMING_VORTICITY
  fp_streaming = fopen("streaming.dat","w");
  fprintf(fp_streaming, "i t vor_stream_mean vor_stream_max liq_vol window_dt \n");
#endif

  fprintf(fp_norm, "i t Omega_liq_avg Omega_liq_rms Omega_liq_vol Omega_liq_max ux_liq_avg ux_liq_rms ux_liq_vol ux_liq_max uy_liq_avg uy_liq_rms uy_liq_vol uy_liq_max ux_liq_savg uy_liq_savg uz_liq_savg uz_liq_rms \n");
  fprintf(fp_stats2, "i t f_liq_sum f_liq_interf posY_max posY_min posY_left posY_right \n");
  fprintf(fp_stats3, "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 c1_liq_sum c1_liq_sum2 c2_liq_sum c2_liq_sum2 c3_liq_sum c3_liq_sum2 b0 b1 b2 \n");
  fprintf(fp_tau,   "i t tau_95 tau_98 tau_100 tau_mean tau_100_strict tau_mean_strict tau_100_signed ediss_mean tau_mean_signed tau_kim_max tau_kim_mean ediss_kim_max ediss_kim_mean \n");
}

static void close_outputs (void)
{
  fclose(fp_stats); fclose(fp_norm); fclose(fp_stats2); fclose(fp_stats3); fclose(fp_tau); fclose(fp_edr);
#if STREAMING_VORTICITY
  if (fp_streaming) fclose(fp_streaming);
#endif
#if RECORD_FLUX || REPLAY_FLUX
  if (fp_rec) fclose(fp_rec);
#endif
}

int main (int argc, char * argv[])
{
  if (argc < 2) { fprintf(stderr, "Usage: BioReactor params.json [restart.dump]\n"); return 1; }
  params = params_read(argv[1]);

  set_scales ();
  set_run_window (argc, argv);
  if (_restart_refused) return 1;   // message printed by restart_refuse() (forcing_state.h)
  set_physics ();
  open_outputs ();

  NITERMAX = 1000;     // multigrid: max iterations per step
  TOLERANCE = 5.0e-4;  // multigrid: residual tolerance

  run();               // Basilisk runs the events below until t_end

  close_outputs ();
}

// ===================================================================================================
// Events, in execution order within a step. DO NOT REORDER (see the top of this file).
// ===================================================================================================
#include "restart.h"        // restart checks (functions only)
#include "init.h"           // event init: geometry, initial fields, or restore
#include "scalars.h"        // tracer and oxygen release
#include "forcing.h"        // rocking acceleration, streaming vorticity, droplet removal
#include "snapshots.h"      // uniform-grid field snapshots
#include "checkpoint.h"     // wall-clock stop and checkpoint dump
#include "diagnostics.h"    // time series: logstats, normcal, pressure diagnostics
#include "replay_events.h"  // record/replay hooks on the VOF step
#include "movies.h"         // video frames (-DVIDEOS=1)
