
//# Advection/diffusion of a soluble tracer

attribute {
  double D1, D2, alpha;  //D1 for f = 1, D2 for f = 0
  scalar phi1, phi2; // private
}

extern scalar * stracers;

/**
## Defaults

On trees we need to ensure conservation of the tracer when
refining/coarsening. */

#if TREE
event defaults (i = 0)
{
  for (scalar s in stracers) {
#if EMBED
    // Basilisk >= 2025: set_prolongation() replaces the combined assignment
    // s.refine = s.prolongation = ... used in the original henry_oxy2.h.
    s.refine = refine_embed_linear;
    set_prolongation (s, refine_embed_linear);
#else
    s.refine  = refine_linear;
#endif
    // Basilisk >= 2025: set_restriction() replaces direct attribute assignment;
    // s.dirty was removed from _Attributes (cache invalidation is now internal).
    set_restriction (s, restriction_volume_average);
  }
}
#endif // TREE


// Advection

static scalar * phi_tracers = NULL;

event vof (i++)
{

  phi_tracers = f.tracers;

  ///*
  for (scalar c in stracers) {

    scalar phi1 = new scalar, phi2 = new scalar;

    c.phi1 = phi1, c.phi2 = phi2;

    // scalar_clone: src/grid/cartesian-common.h
    scalar_clone (phi1, c);
    scalar_clone (phi2, c);
    phi2.inverse = true;

    f.tracers = list_append (f.tracers, phi1);
    f.tracers = list_append (f.tracers, phi2);

    foreach() {
      double a = c[]/(f[]*c.alpha + (1. - f[]));
      phi1[] = a*f[]*c.alpha;
      phi2[] = a*(1. - f[]);
    }

    // Restrict phi1/phi2 from leaf cells to all coarse levels before VOF
    // advection accesses them.  On restart, scalar pool memory may contain NaN
    // at coarse levels (the foreach() above only sets leaf cells).  VOF
    // advection reads coarse values at coarse-fine boundaries → NaN in oxy.
    restriction ({phi1, phi2});
  }
  //*/
}

// Diffusion
///*
struct HDiffusion {
  face vector D;     // alpha
  face vector beta;  // newly added for the second term of diffusion
#if EMBED
  double (* embed_flux) (Point, scalar, vector, double *);
#endif
};

// Harmonic restriction for face diffusivity D.
//
// Default restriction_face uses arithmetic averaging: D_coarse = (D_a + D_b)/2.
// For high-contrast D (D_gas/D_liquid ~ 6e4 here), arithmetic restriction
// produces D_coarse ~ D_gas at every coarse level, making the Gauss-Seidel
// smoother in h_relax overshoot by ~D_gas/D_liquid per level → catastrophic
// multi-level runaway.
//
// Harmonic mean 2ab/(a+b) ≈ 2*min(a,b) for a >> b.  This weights the slow
// (liquid) phase as the transport bottleneck, giving D_coarse ~ 2*D_liquid
// instead of ~D_gas/2, and stabilises the V-cycle for any contrast ratio.
//
// NOTE: harmonic restriction is physically correct in the NORMAL direction
// (series resistance across the interface) but slightly under-estimates the
// TANGENTIAL conductance (parallel paths of D_liq and D_gas).  The multigrid
// coarse solver is less accurate but numerically stable; fine-level accuracy
// is preserved because the finest-level D.x[] is computed from the harmonic
// VOF-weighted formula in tracer_diffusion and is never changed.
static inline void face_harmonic (Point point, vector v)
{
  foreach_dimension()
    foreach_blockf (v.x) {
#if dimension == 1
      v.x[] = fine(v.x,0);
      v.x[1] = fine(v.x,2);
#elif dimension == 2
      {
        double a = fine(v.x,0,0), b = fine(v.x,0,1);
        v.x[] = (a + b > 1e-300) ? 2.*a*b/(a + b) : 0.;
      }
      {
        double a = fine(v.x,2,0), b = fine(v.x,2,1);
        v.x[1] = (a + b > 1e-300) ? 2.*a*b/(a + b) : 0.;
      }
#else // dimension == 3
      {
        double a = fine(v.x,0,0,0), b = fine(v.x,0,1,0),
               c = fine(v.x,0,0,1), d = fine(v.x,0,1,1);
        double denom = (a > 0 ? 1./a : 0.) + (b > 0 ? 1./b : 0.) +
                       (c > 0 ? 1./c : 0.) + (d > 0 ? 1./d : 0.);
        v.x[] = (denom > 1e-300) ? 4./denom : 0.;
      }
      {
        double a = fine(v.x,2,0,0), b = fine(v.x,2,1,0),
               c = fine(v.x,2,0,1), d = fine(v.x,2,1,1);
        double denom = (a > 0 ? 1./a : 0.) + (b > 0 ? 1./b : 0.) +
                       (c > 0 ? 1./c : 0.) + (d > 0 ? 1./d : 0.);
        v.x[1] = (denom > 1e-300) ? 4./denom : 0.;
      }
#endif
    }
}

static inline void restriction_face_harmonic (Point point, scalar s)
{
  face_harmonic (point, s.v);
}

static void h_relax (scalar * al, scalar * bl, int l, void * data)
{
  // similar to relax in poisson.h
  // alpha = D; no lambda
  scalar a = al[0], b = bl[0];
  struct HDiffusion * p = (struct HDiffusion *) data;
  face vector D = p->D, beta = p->beta;

  scalar c = a;
  foreach_level_or_leaf (l) {
    // -lambda[] = cm[]/dt
    double n = - sq(Delta)*b[], d = cm[]/dt*sq(Delta);
    foreach_dimension() {
      // Pe-limit: after restriction, coarse-level beta can still drive Pe > 2.
      // Clamp |beta| at each level so the Gauss-Seidel update is monotone.
      double b1 = beta.x[1], b0 = beta.x[];
      {
        double bmax1 = 2.*D.x[1]/Delta, bmax0 = 2.*D.x[]/Delta;
        if (b1 >  bmax1) b1 =  bmax1; else if (b1 < -bmax1) b1 = -bmax1;
        if (b0 >  bmax0) b0 =  bmax0; else if (b0 < -bmax0) b0 = -bmax0;
      }
      n += D.x[1]*a[1] + D.x[]*a[-1] +
	Delta*(b1*a[1] - b0*a[-1])/2.; // added terms in henry.h
      d += D.x[1] + D.x[] -
	Delta*(b1 - b0)/2.; // added terms in henry.h
    }

    ///*
#if EMBED
    if (p->embed_flux){
      double c;
      double e = embed_flux (point, a, D, &c);
      n -= c*sq(Delta);
      d += e*sq(Delta);
    }
    if (d <= 0.)   // was: if (!d) — also catches d<0 from coarse-level beta terms
      c[] = b[] = 0.;
    else
#endif // EMBED
    //*/
    // Guard: beta terms can drive d <= 0 for large-alpha tracers at thin VOF
    // interface cells (ff→0 makes beta ≈ −D/ff → large).  Clamp d to a floor
    // proportional to the diagonal cm/dt term so the relaxation stays bounded.
    {
      double d_floor = cm[]/dt * sq(Delta) * 1e-12;
      // !(d >= d_floor) catches NaN as well as d < d_floor: ordered >= with NaN
      // returns false, so !false = true → clamp.  Plain (d < d_floor) leaves NaN
      // unchanged because NaN comparisons always return false.
      if (!(d >= d_floor)) d = d_floor;
      if (!isfinite(n)) n = 0.;
      c[] = n/d;
    }
  }
}
//*/

///*
static double h_residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct HDiffusion * p = (struct HDiffusion *) data;
  // similar to residual in poisson.h
  // alpha = D; lambda = beta;
  face vector D = p->D, beta = p->beta;
  double maxres = 0.;

#if TREE
  // Zero ALL res[] cells (leaf AND non-leaf) before the leaf foreach() writes
  // leaf residuals.  foreach() only visits OWNED leaf cells; ghost leaf cells
  // (owned by other MPI ranks) retain whatever was in the scalar pool from a
  // previous use of this slot.  mg_cycle calls restriction(res) which averages
  // over ALL children including ghost leaf ones — stale ghost leaf pool values
  // (which can be ~1e32 after a prior blow-up) propagate into the level-0
  // restricted residual and cascade through h_relax.
  // Zeroing ALSO ghost leaf cells via foreach_cell() ensures they start at 0
  // before mpi_boundary_restriction updates them from the owning rank's correct
  // computed values.  Non-zero pool values in any cell are thus eliminated.
  foreach_cell()
    res[] = 0.;
  // Compute residual inline using D/beta directly (no local face vector g[]).
  // The original g[] TREE version stored face fluxes in a local face vector and
  // then read g.x[1] in foreach().  On MPI restart, foreach_face() only writes
  // *owned* faces; ghost faces at rank boundaries retain pool-SNaN.  boundary({g})
  // does not fill them for locally-allocated face vectors in Basilisk MPI.
  // Instead, use D.x[1]/beta.x[1] directly — their ghost faces are properly
  // filled by boundary({D,beta}) in tracer_diffusion before mg_solve.
  // For a uniform grid (all cells same level) this is numerically identical to
  // the g[] conservative formulation; the difference only arises at coarse/fine
  // AMR interfaces, which are absent at fixed fidelity.
  foreach (reduction(max:maxres)) {
    // -lambda[] = cm[]/dt;
    res[] = b[] + cm[]/dt*a[];
    foreach_dimension() {
      // Pe-limit beta at the right (+1) and left (0) faces independently.
      double bmax_r = 2.*D.x[1]/Delta, b_lim_r = beta.x[1];
      double bmax_l = 2.*D.x[]/Delta,  b_lim_l = beta.x[];
      if (b_lim_r >  bmax_r) b_lim_r =  bmax_r;
      else if (b_lim_r < -bmax_r) b_lim_r = -bmax_r;
      if (b_lim_l >  bmax_l) b_lim_l =  bmax_l;
      else if (b_lim_l < -bmax_l) b_lim_l = -bmax_l;
      res[] -= (D.x[1]*face_gradient_x (a, 1) + b_lim_r*face_value (a, 1) -
                D.x[]*face_gradient_x (a, 0) - b_lim_l*face_value (a, 0))/Delta;
    }

    ///*
    //EMBED
#if EMBED
    if (p->embed_flux){
      double c, e = embed_flux (point, a, D, &c);
      res[] += c - e*a[];
    }
#endif // EMBED
    //*/

    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#else // !TREE
  // "naive" discretisation (only 1st order on trees) //
  foreach (reduction(max:maxres)) {
    res[] = b[] + cm[]/dt*a[];
    foreach_dimension()
      res[] -= (D.x[1]*face_gradient_x (a, 1) -
		D.x[0]*face_gradient_x (a, 0) +
		beta.x[1]*face_value (a, 1) -
		beta.x[0]*face_value (a, 0))/Delta;
    ///*
    //EMBED
#if EMBED
    if (p->embed_flux) {
      double c, e = embed_flux (point, a, D, &c);
      res[] += c - e*a[];
    }
#endif // EMBED
    //*/

    if (fabs (res[]) > maxres)
      maxres = fabs (res[]);
  }
#endif // !TREE
  return maxres;
}
//*/

event tracer_diffusion (i++)
{
  free (f.tracers);
  f.tracers = phi_tracers;

  for (scalar c in stracers) {

    /**
    The advected concentration is computed from $\phi_1$ and $\phi_2$ as
    $$
    c = \phi_1 + \phi_2
    $$
    and these fields are then discarded. */

    scalar phi1 = c.phi1, phi2 = c.phi2, r[];

    foreach() {
      c[] = phi1[] + phi2[];
      r[] = - cm[]*c[]/dt;
    }
    delete ({phi1, phi2});
    // Propagate the leaf-only phi1+phi2 update to all coarse tree levels so
    // mg_solve starts from a consistent state.  Without this, stale checkpoint
    // values at coarse levels cause the first V-cycle to diverge in restarts.
    restriction ({c});
    // Also restrict r so that coarse non-leaf r[] cells carry meaningful values
    // (averaged from leaf) before mg_solve.  halo_restriction (triggered by the
    // stencil extension below) communicates coarse r[] to MPI neighbour ghost
    // cells; without this, those ghost cells receive pool garbage → O(1e141)
    // b[] in h_relax at coarse levels → FE_OVERFLOW SIGFPE on restart.
    restriction ({r});

    /**
    The diffusion equation for $c$ is then solved using the multigrid
    solver and the residual and relaxation functions defined above. */
    ///*
    face vector D[], beta[];
    foreach_face() {
      double ff = (f[] + f[-1])/2.;
      D.x[] = fm.x[]*c.D1*c.D2/(c.D1*(1. - ff) + ff*c.D2);
      beta.x[] = - D.x[]*(c.alpha - 1.)/
	(ff*c.alpha + (1. - ff))*(f[] - f[-1])/Delta;
      // Limit |beta| to cell-Pe ≤ 2: the central-difference Gauss-Seidel
      // relaxation diverges when |beta|*Delta/D > 2.  Thin VOF interface
      // cells (ff→0) with large alpha drive Pe → ∞ and blow up h_relax.
      {
        double beta_max = 2.*D.x[]/Delta;
        if (beta.x[] > beta_max) beta.x[] = beta_max;
        if (beta.x[] < -beta_max) beta.x[] = -beta_max;
      }
    }

    // Sync leaf-level ghost faces across MPI ranks before restricting.
    // face vector D[], beta[] are freshly allocated each call: ghost faces at
    // MPI-rank boundaries retain pool-initialized NaN until boundary() fills them.
    // Without this, restriction() averages NaN ghost children into coarse values,
    // which halo_restriction then propagates to neighbouring ranks' ghost cells.
    // h_relax reads D.x[1]/D.y[0,1] at every multigrid level — if those ghost
    // faces are NaN, FE_INVALID fires. boundary() here overwrites the pool NaN
    // with the correct neighbour-rank values before restriction proceeds.
    boundary ({D, beta});
    // Use harmonic restriction for D at all coarse levels.
    // Arithmetic restriction (the default restriction_face) gives
    // D_coarse ≈ D_gas at coarse levels for D_gas/D_liquid ~ 6e4, causing
    // multigrid Gauss-Seidel overcorrection by ~10^4 per level → runaway.
    // Harmonic restriction gives D_coarse ≈ 2*D_liquid, stabilising the
    // V-cycle without affecting fine-level physics (D is recomputed each step).
    D.x.restriction = restriction_face_harmonic;
    D.y.restriction = restriction_face_harmonic;
    restriction ({D, beta, cm});
    struct HDiffusion q;
    q.D = D;
    q.beta = beta;

    // Extend qcc's stencil analysis to non-leaf levels for ALL fields that
    // h_relax accesses at coarse ghost cells via foreach_dimension().
    //
    // h_relax uses foreach_dimension(), which in 2D expands to x and y iterations.
    // In the y-direction it reads D.y[1] and beta.y[1] — i.e., D.y[] and beta.y[]
    // of the ghost coarse cell above the rank boundary.  Without D.y[], beta.x[],
    // and beta.y[] in the stencil, halo_restriction never communicates those face-
    // vector components for non-leaf ghost cells; they retain pool-garbage values.
    // h_relax then computes c[] = n/d with garbage ghost inputs → multigrid diverges
    // → oxy reaches ~1e275 → statsf2 overflow → 0*Inf = NaN → FE_INVALID crash.
    //
    // KEY DISTINCTION (confirmed by experiment):
    //   WRONG: (void)(D.x[1]) — adds OFFSET cell (x+1,y) to stencil → widens MPI
    //          halo for restriction at all levels → ~7000x slowdown (commit b3e4163).
    //   RIGHT: (void)(D.y[]) — adds y-COMPONENT of same ghost cell → same halo
    //          width, only more fields per ghost cell → no slowdown.
    //
    // c[], r[] — solution and RHS scalars (covers h_relax reads of a[1]=c[] and b[])
    // D.x[], D.y[] — x and y face-vector components of diffusivity
    // beta.x[], beta.y[] — x and y face-vector components of advective coefficient
    foreach_cell()
      if (!is_leaf(cell))
        (void)(c[] + r[] + D.x[] + D.y[] + beta.x[] + beta.y[]);

    scalar aa = c;
#if EMBED
    if (!q.embed_flux && aa.boundary[embed] != symmetry)
      q.embed_flux = embed_flux;
#endif //EMBED

    mg_solve ({c}, {r}, h_residual, h_relax, &q);
  }
}

/**
## References

~~~bib
@article{haroun2010,
  title = {Volume of fluid method for interfacial reactive mass transfer:
           application to stable liquid film},
  author = {Haroun, Y and Legendre, D and Raynal, L},
  journal = {Chemical Engineering Science},
  volume = {65},
  number = {10},
  pages = {2896--2909},
  year = {2010},
  doi = {10.1016/j.ces.2010.01.012}
}

@hal{farsoiya2021, hal-03227997}
~~~
*/
