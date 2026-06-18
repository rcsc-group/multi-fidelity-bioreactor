// Primary author: Minki Kim
// Contributor: Radu Cimpeanu
// Contributor: Daniel Harris
// Date: 03/04/2025

// Core Basilisk modules for embedding geometries and centered Navier-Stokes solver
#include "embed.h"
#include "navier-stokes/centered.h"

// Define dynamic viscosity as a function of the volume fraction (f)
// Uses harmonic averaging for two-phase flow
#define mu(f) ( 1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2)  )

// Modules for two-phase flow (e.g., air and water) with surface tension and volume conservation
#include "two-phase.h"
#include "tension.h"
#include "navier-stokes/conserving.h"

// Modules for oxygen transfer + tracers
#include "henry_oxy2.h"

// Custom headers foroxy specific configurations and visualization
#include "view3.h"
#include "utils2.h"

// Mathematical constants (e.g., M_PI from math.h)
#define _USE_MATH_DEFINES
#include <math.h>

// JSON parameter reader (jsmn-based); replaces argv[1..3] with params.json
#include "params_read.h"

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
#define N_RAMP_CYCLES 3            // Ramp duration in rocking cycles; t_change_st = N_RAMP_CYCLES * T_per_st
const double th_cont = 90;        // Contact angle for wetting conditions (degrees)
double t_mix,t_dump;              // Time at which tracer is released, and dump file is saved (computed later)
const double nMix_cycle = 80;     // Number of cycles for tracer release (used to compute t_mix)
double t_end;                      // Final simulation time (simulation time unit); set from params.t_end in main()

// Output time intervals (derived from experimental timing)
const double dt_video= 0.6074/5;  // Interval for video frames
const double t_out   = 0.1;       // Output interval for statistics [non-dim time]


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
#if TRACER && OXYGEN
scalar * stracers = {c, oxy, c1, c2, c3};
#elif TRACER
scalar * stracers = {c, c1, c2, c3};
#elif OXYGEN
scalar * stracers = {oxy};
#else
scalar * stracers = NULL;
#endif
double (* gradient) (double, double, double) = minmod2;   // Custom slope limiter function (used for scalar gradients)

// Buffers for file naming and output file pointers for statistics
char buf1[100], buf2[100], buf3[100], buf4[100];
FILE * fp_stats, * fp_norm, * fp_stats2, * fp_stats3, * fp_tau;

// Key physical and dimensionless parameters (computed in main)
double U0, Re_w, Re_a, We_w, Fr, rhor, mur, Pe_tracer_1, Pe_tracer_2, Pe_oxy_1, Pe_oxy_2, Th, Th_d, Th_2d, U_bio, w_bio, w_bio_st, T_per_st, T_bio, t_change_st;
BioreactorParams params;  // global so all events (acceleration, init) can access it
// Restart / checkpoint support
static double t_ramp_start      = 0.0;   // simulation time when the current ramp began
static double t_dump_checkpoint = 0.0;   // simulation time to write checkpoint.dump
static const char * restart_file = NULL; // argv[2] if this is a restart run
int MINLEVEL, MAXLEVEL;    // Mesh refinement levels


// ================================================================== //
//                       MAIN FUNCTION                                //
// ================================================================== //
int main(int argc, char * argv[]){

  if (argc < 2) { fprintf(stderr, "Usage: BioReactor params.json\n"); return 1; }
  params = params_read(argv[1]);

  // Derive legacy scalars from params for the rest of main() unchanged.
  // omega_b (rad/s) and theta_max[0] (deg) replace the old ANGLE/RPM CLI args.
  // Multi-harmonic forcing and superellipse geometry will consume params directly
  // once the acceleration event and init event are extended (next steps).
  L_bio = params.geometry_a;  // bag half-length (m); sets all non-dim scales
  double ANGLE = params.theta_max[0];       // Fundamental rocking amplitude (degrees)
  double RPM   = params.omega_b * 60. / (2.*M_PI);  // Convert rad/s → RPM
  NN    = 1 << params.fidelity;             // fidelity → grid cells per side (4→16, 7→128, 9→512)
  t_end = params.t_end;                    // non-dimensional sim time; 1 unit = T_bio seconds

  L0 = 1. [0];  // [0] declares space dimensionless: simulation is fully non-dimensionalized (scaled by L_bio, T_bio, U_bio); Basilisk dimensional analysis requires annotations on literals, not variables
  // DT: maximum allowed timestep (dimensionless).  HUGE is fine for serial/OpenMP
  // but MPI's global dtmax reduction can behave unexpectedly with infinity on
  // some implementations.  1.0 [0] (one non-dim time unit) is always larger than
  // the CFL-limited dt (~dx/U_max << 1) so it never artificially constrains the step.
  DT = 1. [0];
  origin(-L0/2., -L0/2.);   // Set coordinate origin to domain center

  init_grid(NN);

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
  H_bio  = L_bio*Ly;
  V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
  U_bio  = V_bio/(H_bio*0.5)/T_per; // Characteristic velocity scale
  T_bio  = L_bio/U_bio;             // Characteristic time scale
  w_bio  = 2*pi/T_per;              // Angular velocity (rad/s)
  w_bio_st = w_bio*T_bio;           // Dimensionless angular velocity
  T_per_st = T_per/T_bio;           // Dimensionless period
  U0     = w_bio_st*Th_max;         // Initial rotational velocity
  t_change_st = N_RAMP_CYCLES * T_per_st;  // ramp over 3 cycles regardless of omega_b
  t_mix      = T_per_st*params.n_mix_cycles; // rocking cycles before tracer/oxygen start (wired from params.json)
  t_dump = t_mix;                   // Time to dump data (simulation time)

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
    // Smooth-step interpolation starts AT the checkpoint and runs N_RAMP_CYCLES forward.
    // alpha goes 0→1 over [t_checkpoint, t_checkpoint + N_RAMP_CYCLES*T_per_st].
    t_ramp_start     = params.t_checkpoint;
    t_mix            = params.t_checkpoint + T_per_st * params.n_mix_cycles;
    t_dump           = t_mix;
    {
      double t_end_abs = params.t_checkpoint + params.t_end;
      int n_per        = (int)(t_end_abs / T_per_st) + 1;
      t_dump_checkpoint = n_per * T_per_st;
      t_end            = t_dump_checkpoint;
    }
  } else {
    // Fresh run: extend t_end to the next period boundary for clean phase alignment.
    int n_per = (int)(t_end / T_per_st) + 1;
    t_dump_checkpoint = n_per * T_per_st;
    t_end = t_dump_checkpoint;
  }

  // Dimensionless numbers
  Re_w = rho_w*U_bio*L_bio/mu_w;            // Reynolds number of water
  Re_a = rho_a*U_bio*L_bio/mu_a;            // Reynolds number of air
  We_w = rho_w*U_bio*U_bio*L_bio/sigma;     // Weber number of water
  Fr   = U_bio/sqrt(grav*L_bio);            // Froude number
  rhor = rho_a/rho_w;                       // Density ratio
  mur  = mu_a/mu_w;                         // Viscosity ratio
  Pe_tracer_1 = U_bio*L_bio/D_tracer_1;     // Peclet number of tracer at phase 1 (water)
  Pe_tracer_2 = U_bio*L_bio/D_tracer_2;     // Peclet number of tracer at phase 2 (air)
  Pe_oxy_1    = U_bio*L_bio/D_oxy_1;        // Peclet number of oxygen at phase 1 (water)
  Pe_oxy_2    = U_bio*L_bio/D_oxy_2;        // Peclet number of oxygen at phase 1 (air) 

  rho1 = 1.0;             // Reference density of phase 1 (water) (scaled to 1)
  rho2 = rho1*rhor;       // Scaled density of phase 2 (air)
  mu1  = 1.0/Re_w;        // Dimensionless viscosity of phase1 (water)
  mu2  = mur*mu1;         // Scaled viscosity of phase 2 (air)
  f.sigma = 1.0 / We_w;   // Dimensionless surface tension
  
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

  fprintf(fp_norm, "i t Omega_liq_avg Omega_liq_rms Omega_liq_vol Omega_liq_max ux_liq_avg ux_liq_rms ux_liq_vol ux_liq_max uy_liq_avg uy_liq_rms uy_liq_vol uy_liq_max \n");
  fprintf(fp_stats2, "i t f_liq_sum f_liq_interf posY_max posY_min \n");
  fprintf(fp_stats3, "i t oxy_liq_sum oxy_liq_sum2 c_liq_sum c_liq_sum2 c1_liq_sum c1_liq_sum2 c2_liq_sum c2_liq_sum2 c3_liq_sum c3_liq_sum2 \n");
  fprintf(fp_tau,   "i t tau_95 tau_98 tau_100 tau_mean \n");

  NITERMAX = 1000;     // Max iterations per timestep
  TOLERANCE = 5.0e-4;  // // Solver tolerance (convergence criterion)
  
  // Run the simulation — Basilisk manages events from here
  run();
  
  // Close all output files
  fclose(fp_stats); fclose(fp_norm); fclose(fp_stats2); fclose(fp_stats3); fclose(fp_tau);
}


// ================================================================== //
//                      INITIAL CONDITIONS                            //
// ================================================================== //
event init (t = 0)
{
  if (restart_file) {
    // ── Checkpoint restart ────────────────────────────────────────────────
    // restore() is called here so that Basilisk's grid infrastructure is
    // fully active.  It sets t = t_checkpoint and restores all fields.
    // All timing was already computed in main() from params.t_checkpoint so
    // the event system has the correct t_dump_checkpoint before run().
    if (!restore (file = restart_file)) {
      fprintf (stderr, "ERROR: restore() failed to open '%s' — aborting\n", restart_file);
      exit (1);
    }
    // fs (embed face fractions) is a face field — excluded from Basilisk dumps.
    // After restore, fs=0 everywhere: the NS solver sees no solid walls and the
    // velocity collapses on the first timestep.  Re-compute fs from the same
    // static geometry used at fresh-start (rocking is body forces, not moving solid).
#if EMBED
    {
      double a_nd = params.geometry_a / L_bio;
      double b_nd = params.geometry_b / L_bio;
      if (params.geometry_n >= 8.)
        solid (cs, fs, intersection(a_nd - fabs(x), b_nd - fabs(y)));
      else
        solid (cs, fs, 1. - pow(fabs(x/a_nd), params.geometry_n)
                          - pow(fabs(y/b_nd), params.geometry_n));
    }
#endif
    // Rescale stored velocity and pressure to the new segment's non-dim frame.
    // U_bio ∝ omega_b (fixed geometry, theta_max) → scale = omega_b_prev / omega_b.
    // Without this, the restored velocity is 2× too large when frequency doubles.
    if (params.omega_b_prev > 0.) {
      double su = params.omega_b_prev / params.omega_b;
      foreach() {
        u.x[] *= su;
        u.y[] *= su;
        p[]   *= su * su;
        pf[]  *= su * su;   // half-step pressure (BCG tracer advection), same scale as p
      }
      // uf (face-centered velocity) must also be rescaled.  centered.h's
      // advection scheme (BCG), the VOF step, and the AMR adapt-at-t=0
      // all use uf directly before it is rebuilt from u.  Without this,
      // the first Poisson projection sees a factor-of-(1/su)^2 kinetic
      // energy mismatch and destroys the restored velocity field.
      foreach_face()
        uf.x[] *= su;
      boundary ({u.x, u.y, p, pf, uf.x, uf.y});
      // Scale g (combined pressure-gradient + acceleration term from BCG predictor).
      // g is stored in the checkpoint in seg0 non-dim units.  With ramp=1 immediately,
      // the acceleration part of g is correct (w_bio_st is constant across omega_b).
      // The pressure-gradient part scales as su² (same as p).  Scale by su² to keep
      // the pressure-gradient contribution accurate in the first BCG half-step.
      foreach() {
        g.x[] *= su * su;
        g.y[] *= su * su;
      }
      boundary ({g.x, g.y});

    }
    // Re-apply the prolongation/restriction setup from event defaults(i=0) in
    // henry_oxy2.h.  That event fires at i=0 on a fresh start but is skipped on
    // restart (i resumes from the checkpoint value, so i==0 is never seen again).
    // Basilisk's dump/restore does NOT serialize scalar function-pointer attributes
    // (restriction, prolongation, refine), so they revert to defaults.  Without
    // restriction_volume_average the multigrid h_relax propagates NaN from solid
    // cells (cm=0 → d=0 → c[]=n/0) into fluid cells → NaN kLa and tracer sums.
#if TREE
    for (scalar s in stracers) {
      s.refine = refine_embed_linear;
      set_prolongation (s, refine_embed_linear);
      set_restriction (s, restriction_volume_average);
    }
#endif

    // Re-apply centered.h's EMBED defaults (event defaults(i=0) is skipped on
    // restart because i resumes from the checkpoint value).  Without these,
    // p/pf/u/g use non-EMBED prolongation/restriction and uf uses
    // refine_face_solenoidal instead of refine_face.  The pressure gradient at
    // the embedded wall falls back to the generic gradient (no
    // pressure_embed_gradient), so the Poisson solve near the solid accumulates
    // small errors each timestep.  After ~300 steps (~1 T) the velocity
    // diverges → SIGFPE in tracer_diffusion/h_residual.
    // Confirmed experimentally: fresh run to t=15 is stable; restart crashes at
    // t=t_checkpoint+1.14T with both 4 and 16 MPI ranks.
#if TREE && EMBED
    uf.x.refine = refine_face;
    foreach_dimension()
      uf.x.prolongation = refine_embed_face_x;
    for (scalar s in {p, pf, u, g}) {
      s.refine = refine_embed_linear;
      set_prolongation (s, refine_embed_linear);
      set_restriction (s, restriction_embed_linear);
    }
    for (scalar s in {p, pf})
      s.embed_gradient = pressure_embed_gradient;
#endif // TREE && EMBED
    // Re-apply vof.h defaults: fraction_refine is set for f in vof.h's
    // event defaults(i=0), skipped on restart.  Without it AMR uses bilinear
    // prolongation for f, creating non-physical VOF fractions near the interface.
#if TREE
    f.refine = fraction_refine;
    set_prolongation (f, fraction_refine);
#endif // TREE

    // Reset ALL stracers to zero at EVERY multigrid level.  reset() zeroes
    // owned cells; boundary() communicates zeroed leaf ghost values.  But
    // coarse-level ghost cells (owned by a neighbouring MPI rank) are NOT
    // updated by boundary() — they retain stale checkpoint values.  When
    // h_relax runs at a coarse V-cycle level it reads a[1] which can be a
    // coarse ghost cell, producing a spurious non-zero restricted residual
    // (b[] >> 1) and driving a gradual multigrid divergence over ~5 periods.
    // restriction() propagates leaf=0 to coarse own cells and then runs
    // halo_restriction which communicates those zeroed coarse values to
    // neighbouring ranks' ghost slots, closing the ghost-cell gap.
    reset (stracers, 0.);
    boundary (stracers);
    restriction (stracers);
  } else {
    // ── Fresh start ────────────────────────────────────────────────────────
    // Parametric bag geometry (dimensionless semi-axes)
    double a_nd = params.geometry_a / L_bio;
    double b_nd = params.geometry_b / L_bio;  // == Ly

    // Fill level: liquid occupies fill_level fraction of bag height, measured from bottom
    double y_fill = b_nd * (2.*params.fill_level - 1.);
    fraction(f, y_fill - y);

    // Superellipse solid: |x/a|^n + |y/b|^n = 1
    // n >= 8 → perfect rectangle (avoids pow() singularities at sharp corners)
    #if EMBED
    if (params.geometry_n >= 8.)
      solid(cs, fs, intersection(a_nd - fabs(x), b_nd - fabs(y)));
    else
      solid(cs, fs, 1. - pow(fabs(x/a_nd), params.geometry_n)
                       - pow(fabs(y/b_nd), params.geometry_n));
    #endif
  }
}


// ================================================================== //
//                           TRACER SETUP                             //
// ================================================================== //
#if TRACER
event tracer(t = t_mix){

  double h_tr;
  h_tr = (M_PI*R_tr*R_tr);  // Area of circular tracer patch

  // circular shape tracer at the center of the liquid
  // fraction(c, -(sq(x-x_tr) + sq(y-y_tr) - sq(R_tr)) ); 

  // tracer released as a line (same area)
  // fraction(c, intersection( -(y-y_tr - 0.5*h_tr), -(-(y-y_tr + 0.5*h_tr)) ));

  // Vertical mixing-top side
  #if VERTICAL_MIXUP
  {
    // Midpoint of the actual liquid column — fill-level-aware.
    // y_fill = Ly*(2*fill-1), y_bot = -Ly → midpoint = Ly*(fill-1).
    // The old hardcoded -Ly/4 was only correct for fill=0.5 and left
    // fill<=0.3 with zero tracer coverage.
    double y_liq_mid = Ly * (params.fill_level - 1.0);
    foreach(){
      // f[]>0.5 (majority-liquid) instead of f[]==1 (exact) so that near-interface
      // cells are included — at coarse resolutions f never reaches exactly 1.0 in
      // the top liquid layer, causing c2 to remain 0 throughout the run.
      if ((f[] > 0.5) && (cs[]==1) && (y >= y_liq_mid))
        c2[] = 1.0;    // Upper half of liquid
    }
    // Synchronise MPI ghost cells after direct foreach() write.
    // Without this, ranks that border the upper liquid region see stale c2
    // values in neighbour halos at the next tracer_diffusion call, producing
    // a spurious large gradient (c_tracer_alpha=1e30 amplifies beta enormously)
    // that crashes h_relax under MPI FP trapping.  Same fix as boundary({oxy})
    // in the oxygen event.
    boundary ({c2});
  }
  #endif

}
#endif


// ================================================================== //
//                           OXYGEN SETUP                             //
// ================================================================== //
#if OXYGEN
event oxygen (t=t_mix; i++){

#if OXYGEN_AIR
  foreach(){
    if ((f[] == 0) && (cs[]==1))
      oxy[] = 1.;     // Oxygen in gas regions only
  }
  // Synchronise MPI ghost cells after direct foreach() write.
  // Without this, ranks that border gas regions see stale oxy=0 in
  // neighbour halos at the first tracer_diffusion call, producing a
  // spurious large gradient that crashes h_relax under MPI FP trapping.
  boundary ({oxy});
#endif
}
#endif


// ================================================================== //
//                           ACCELERATION                             //
// ================================================================== //
#if ACCELERATION
event acceleration(i++)
{
  // Smooth-step parameter interpolation: alpha: 0→1 over N_RAMP_CYCLES.
  // For fresh runs, *_prev fields are 0 → reproduces the original cold-start ramp.
  // For restarts, *_prev fields carry the previous segment's values → smooth transition
  // between two fully-forced steady states without ever underdriving the system.
  double elapsed  = t - t_ramp_start;
  double ramp_dur = N_RAMP_CYCLES * T_per_st;
  double x_ss     = (elapsed < ramp_dur) ? elapsed / ramp_dur : 1.0;
  double alpha    = 3.*x_ss*x_ss - 2.*x_ss*x_ss*x_ss;   // smooth-step ∈ [0,1]


  // Multi-harmonic angular forcing with smooth-step interpolation of amplitude and phase.
  // For each harmonic k: Ak and phk are interpolated from _prev → current over N_RAMP_CYCLES.
  Th = 0;  Th_d = 0;  Th_2d = 0;
  for (int k = 1; k <= params.n_harmonics; k++) {
    double wk  = k * w_bio_st;
    double Ak  = ((1.-alpha)*params.theta_max_prev[k-1]
                +     alpha *params.theta_max[k-1]) * pi / 180.;
    double phk =  (1.-alpha)*params.phi_angular_prev[k-1]
                +     alpha *params.phi_angular[k-1];
    Th    +=  Ak * sin(wk*t + phk);
    Th_d  +=  Ak * wk * cos(wk*t + phk);
    Th_2d += -Ak * wk*wk * sin(wk*t + phk);
  }

  // Multi-harmonic horizontal forcing with smooth-step interpolation.
  double x_acc = 0.;
  {
    double omega_h_now = (1.-alpha)*params.omega_h_prev + alpha*params.omega_h;
    if (omega_h_now > 0.) {
      double w_h_st = omega_h_now * T_bio;
      for (int k = 1; k <= params.n_harmonics; k++) {
        double wk_h = k * w_h_st;
        double Ah   = (1.-alpha)*params.amplitude_h_prev[k-1]
                    +     alpha *params.amplitude_h[k-1];
        double phh  = (1.-alpha)*params.phi_horizontal_prev[k-1]
                    +     alpha *params.phi_horizontal[k-1];
        x_acc += (Ah / L_bio) * sq(wk_h) * sin(wk_h*t + phh);
      }
    }
  }

  face vector av = a;
  // av.x: gravity + Coriolis + centrifugal + Euler (azimuthal) + horizontal translation
  // av.y: gravity + Coriolis + centrifugal + Euler (azimuthal)
  foreach_face(x)
    av.x[] = -sin(Th)/(Fr*Fr) + 2*Th_d*(u.y[] + u.y[-1,0])*0.5
    + Th_d*Th_d*(x+L_piv*sin(Th)) + Th_2d*(y+L_piv*cos(Th)) + x_acc;
  foreach_face(y)
    av.y[] = -cos(Th)/(Fr*Fr) - 2*Th_d*(u.x[] + u.x[0,-1])*0.5
    + Th_d*Th_d*(y+L_piv*cos(Th)) - Th_2d*(x+L_piv*sin(Th));
  a = av;
}
#endif


// Write a Basilisk checkpoint at the first complete period boundary after t_end.
// The checkpoint is always at θ=0 (zero-crossing) — clean phase alignment for
// the next segment's soft-start ramp.  Controlled by t_dump_checkpoint global.
event dump_checkpoint (t = t_dump_checkpoint) {
  if (pid() == 0)
    fprintf (ferr, "checkpoint: writing checkpoint.dump at t=%.4g\n", t);
  // p and pf have nodump=true by default in centered.h (they're reconstructed
  // from scratch after a restore, causing a first-step velocity crash).
  // Force them into the checkpoint so the restart Poisson solve starts from
  // the correct pressure and applies only a tiny correction.
  p.nodump = pf.nodump = false;
  dump (file = "checkpoint.dump");
  p.nodump = pf.nodump = true;
}

//  Log performance and runtime
event logstats (t+=0.1; t <= t_end) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    if (pid() == 0){
      fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
      fflush(fp_stats);
    }
}

// Pressure Poisson residuals — health-check build only (DIAGNOSTICS=1)
#if DIAGNOSTICS
event pressure_diagnostics (t+=0.1; t<=t_end) {
  static FILE *fp_pdiag = NULL;
  if (!fp_pdiag) {
    fp_pdiag = fopen("pressure_diag.dat", "w");
    fprintf(fp_pdiag, "i t mgp_resa mgu_resa mgp_i mgu_i\n");
  }
  if (pid() == 0) {
    fprintf(fp_pdiag, "%d %g %g %g %d %d\n", i, t, mgp.resa, mgu.resa, mgp.i, mgu.i);
    fflush(fp_pdiag);
  }
}
#endif

// Compute simulation statistics
#if NORMCAL
event normcal (t+=t_out; t<=t_end){
    //timing s = timer_timing (perf.gt, i, perf.tnc, NULL);

    scalar ux_liq[],uy_liq[],ux_liq_abs[],omega[],omega_liq[],oxy_liq[],f_liq[],posY[],c_liq[],c1_liq[],c2_liq[],c3_liq[];
    double omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max;
    double ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max;
    double f_liq_sum,f_liq_interf,posY_max,posY_min,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2;
    
    vorticity (u, omega); // vorticity

    // only liquid velocity
    foreach(){
      ux_liq[]  = u.x[]*f[];
      uy_liq[]  = u.y[]*f[];
      // BUG (upstream rcsc-group/BioReactor): oxy[]*f[] mixes gas-phase oxygen into
      // the integral at interface cells, inflating oxy_liq_sum above its physical bound.
      // FIX: extract only the liquid-phase contribution via Henry's law activity.
      // oxy[] = a*(f*alpha + (1-f))  where a = gas-equivalent activity.
      // Liquid-phase concentration (normalised by saturation = alpha) = f*a = f*oxy/(f*alpha+(1-f)).
      // Result: oxy_liq=0 in pure gas, oxy_liq=1 in pure liquid at saturation.
      oxy_liq[] = f[]*oxy[]/(f[]*c_oxy_alpha + (1.-f[]) + 1e-10);
      omega_liq[] = omega[]*f[];
      // BUG (upstream): (1-cs[])*f[] is zero everywhere inside the bag (cs=1 in fluid cells),
      // so statsf2(f_liq).sum ≈ 0.002 (only bag-wall cut cells) instead of ~0.28 (liquid volume).
      // Dividing oxy_liq_sum by this tiny number pushed C* >> 1 immediately.
      // FIX: use f[] directly — statsf2(f).sum = true liquid volume via embed-aware dv().
      f_liq[]   = f[];
      c_liq[]   = c[]*f[];
      c1_liq[]  = c1[]*f[];
      c2_liq[]  = c2[]*f[];
      c3_liq[]  = c3[]*f[];
    }

    position (f, posY, {0,1,0});  // (0,1,0) indicates the unit vector in the y-direction

    omega_liq_avg = normf(omega_liq).avg;
    omega_liq_rms = normf(omega_liq).rms;
    omega_liq_vol = normf(omega_liq).volume;
    omega_liq_max = normf(omega_liq).max;
    ux_liq_avg    = normf(ux_liq).avg;
    ux_liq_rms    = normf(ux_liq).rms;
    ux_liq_vol    = normf(ux_liq).volume;
    ux_liq_max    = normf(ux_liq).max;
    uy_liq_avg    = normf(uy_liq).avg;
    uy_liq_rms    = normf(uy_liq).rms;
    uy_liq_vol    = normf(uy_liq).volume;
    uy_liq_max    = normf(uy_liq).max;
    
    f_liq_sum     = statsf2(f_liq).sum;
    f_liq_interf  = interface_area(f);
    posY_max      = statsf(posY).max;
    posY_min      = statsf(posY).min;

    oxy_liq_sum   = statsf2(oxy_liq).sum;
    oxy_liq_sum2  = statsf2(oxy_liq).sum2;
    c_liq_sum     = statsf2(c_liq).sum;
    c_liq_sum2    = statsf2(c_liq).sum2;
    c1_liq_sum    = statsf2(c1_liq).sum;
    c1_liq_sum2   = statsf2(c1_liq).sum2;
    c2_liq_sum    = statsf2(c2_liq).sum;
    c2_liq_sum2   = statsf2(c2_liq).sum2;
    c3_liq_sum    = statsf2(c3_liq).sum;
    c3_liq_sum2   = statsf2(c3_liq).sum2;

    // ── 98th-percentile shear stress τ₉₈(t) ──────────────────────────────────
    // tau = |μ(∂u/∂y + ∂v/∂x)| in the liquid domain (f[] > 0.5).
    // Masked to bulk liquid to exclude interface artefacts.
    //
    // Implementation: two-pass 200-bin histogram — no dynamic memory, MPI-safe.
    //   Pass 1: global tau_max via statsf2 (Basilisk internally does MPI_Allreduce).
    //   Pass 2: local bin counts → MPI_Allreduce(SUM) → global percentile walk.
    //
    // Stencil mirrors vorticity() in basilisk/src/utils.h: uses face-centred
    // velocity values u.x[0,±1] and u.y[±1] for 2nd-order centred differences.
    // Pass 1: global max and mean via inline reduction (no scalar allocation).
    // Declaring scalar tau_liq[] inside an event leaks a Basilisk scalar on every
    // call and corrupts the global scalar list, causing segfaults at fidelity ≥7.
    double tau_max_val = 0., tau_sum = 0., tau_vol = 0.;
    foreach(reduction(max:tau_max_val) reduction(+:tau_sum) reduction(+:tau_vol)) {
      if (f[] > 0.5) {
        double du_dy = (u.x[0,1] - u.x[0,-1]) / (2.*Delta);
        double dv_dx = (u.y[1]   - u.y[-1])   / (2.*Delta);
        double tau   = mu(f[]) * fabs(du_dy + dv_dx);
        if (tau > tau_max_val) tau_max_val = tau;
        tau_sum += tau * (Delta*Delta);
        tau_vol += Delta*Delta;
      }
    }
    double tau_mean_val = (tau_vol > 0.) ? tau_sum / tau_vol : 0.;
    if (tau_max_val < 1e-14) tau_max_val = 1e-14;  // guard /0

    // Pass 2: 200-bin histogram — recomputes tau inline, no scalar needed
    #define TAU_BINS 200
    long bins[TAU_BINS];
    for (int k = 0; k < TAU_BINS; k++) bins[k] = 0;
    foreach() {
      if (f[] > 0.5) {
        double du_dy = (u.x[0,1] - u.x[0,-1]) / (2.*Delta);
        double dv_dx = (u.y[1]   - u.y[-1])   / (2.*Delta);
        double tau   = mu(f[]) * fabs(du_dy + dv_dx);
        int b = (int)(tau / tau_max_val * (TAU_BINS - 1));
        if (b < 0)         b = 0;
        if (b >= TAU_BINS) b = TAU_BINS - 1;
        bins[b]++;
      }
    }
    #if _MPI
    {
      long gbins[TAU_BINS];
      MPI_Allreduce(bins, gbins, TAU_BINS, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
      for (int k = 0; k < TAU_BINS; k++) bins[k] = gbins[k];
    }
    #endif

    // Walk bins once to find 95th and 98th percentiles; 100th == tau_max_val
    long total = 0, cumul = 0;
    for (int k = 0; k < TAU_BINS; k++) total += bins[k];
    double tau_95_val = tau_max_val;
    double tau_98_val = tau_max_val;
    int    found_95   = 0;
    for (int k = 0; k < TAU_BINS; k++) {
      cumul += bins[k];
      if (!found_95 && cumul >= (long)(0.95 * (double)total)) {
        tau_95_val = tau_max_val * (k + 1.0) / TAU_BINS;
        found_95   = 1;
      }
      if (cumul >= (long)(0.98 * (double)total)) {
        tau_98_val = tau_max_val * (k + 1.0) / TAU_BINS;
        break;
      }
    }
    #undef TAU_BINS

   // i, timestep, no of cells, real time elapsed, cpu time
   if (pid() == 0){

      fprintf(fp_norm, "%i %g %g %g %g %g %g %g %g %g %g %g %g %g \n",i,t,omega_liq_avg,omega_liq_rms,omega_liq_vol,omega_liq_max,ux_liq_avg,ux_liq_rms,ux_liq_vol,ux_liq_max,uy_liq_avg,uy_liq_rms,uy_liq_vol,uy_liq_max);
      fflush(fp_norm);

      fprintf(fp_stats2, "%i %g %g %g %g %g \n",i,t,f_liq_sum,f_liq_interf,posY_max,posY_min);
      fflush(fp_stats2);

      fprintf(fp_stats3, "%i %g %g %g %g %g %g %g %g %g %g %g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2,c1_liq_sum,c1_liq_sum2,c2_liq_sum,c2_liq_sum2,c3_liq_sum,c3_liq_sum2);
      //fprintf(fp_stats3, "%i %g %g %g %g %g \n",i,t,oxy_liq_sum,oxy_liq_sum2,c_liq_sum,c_liq_sum2);
      fflush(fp_stats3);

      fprintf(fp_tau, "%i %g %g %g %g %g \n", i, t, tau_95_val, tau_98_val, tau_max_val, tau_mean_val);
      fflush(fp_tau);
   }
}
#endif

// Make videos — field data exported as binary frames; rendered by scripts/render_videos.py
#if VIDEOS
int _vframe = 0;

event movies_output(t = t_mix; t += dt_video; t <= t_end)
{
  // MPI-safe video output:
  // interpolate() is collective (MPI_Allreduce per point) — ALL ranks must call it.
  // File I/O is guarded to rank 0 only so no two ranks open the same file.
  // _vframe increments on all ranks so the counter stays in sync.
#if _MPI
  if (pid() == 0)
#endif
  {
    if (_vframe == 0)
      system("mkdir -p frames");
  }

  int    n    = NN;  // NN = 1<<fidelity, set in main()
  double t_nd = t;
  double dx   = L0 / n;

  // Horizontal displacement (lab frame): X_lab = sum_k A_k * sin(k*w_h*t + phi_k)
  double xh_nd = 0.0;
  if (params.omega_h > 0.) {
    double w_h_nd = params.omega_h * T_bio;
    for (int k = 1; k <= params.n_harmonics; k++) {
      double wk = k * w_h_nd;
      xh_nd += (params.amplitude_h[k-1] / L_bio) * sin(wk*t + params.phi_horizontal[k-1]);
    }
  }

  // VOF field f interpolated onto n×n uniform grid, row-major (j=0 → y=Y0 = bottom).
  // Run on ALL ranks — interpolate() uses MPI_Allreduce internally and requires
  // every rank to participate for each (xi, yj) query.
  float *buf = (float *)malloc(n * n * sizeof(float));
  int idx = 0;
  for (int j = 0; j < n; j++) {
    double yj = Y0 + (j + 0.5) * dx;
    for (int i = 0; i < n; i++) {
      double xi = X0 + (i + 0.5) * dx;
      buf[idx++] = (float)interpolate(f, xi, yj);
    }
  }

  // Only rank 0 writes; all ranks already have the correct buf values.
#if _MPI
  if (pid() == 0)
#endif
  {
    char fpath[512];
    sprintf(fpath, "frames/frame_%06d.bin", _vframe);
    FILE *fp = fopen(fpath, "wb");
    if (fp) {
      fwrite(&n,     sizeof(int),    1, fp);
      fwrite(&t_nd,  sizeof(double), 1, fp);
      fwrite(&Th,    sizeof(double), 1, fp);
      fwrite(&xh_nd, sizeof(double), 1, fp);
      fwrite(buf,    sizeof(float), n * n, fp);
      fclose(fp);
    } else {
      fprintf(stderr, "movies_output: cannot open %s\n", fpath);
    }
  }

  free(buf);
  _vframe++;
}
#endif

