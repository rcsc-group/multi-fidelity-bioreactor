
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
    if (!d)
      c[] = b[] = 0.;
    else
#endif // EMBED
    //*/
    // Guard: beta terms can drive d <= 0 for large-alpha tracers at thin VOF
    // interface cells (ff→0 makes beta ≈ −D/ff → large).  Clamp d to a floor
    // proportional to the diagonal cm/dt term so the relaxation stays bounded.
    {
      double d_floor = cm[]/dt * sq(Delta) * 1e-12;
      if (d < d_floor) d = d_floor;
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
  /// conservative coarse/fine discretisation (2nd order) //
  face vector g[];
  // Zero non-leaf res[] cells before the leaf foreach() writes leaf residuals.
  // foreach() only visits LEAF cells; coarse (non-leaf) cells retain whatever
  // was in the scalar pool from a previous use of this slot.  mg_cycle calls
  // restriction(res) which averages over ALL children including non-leaf ones —
  // stale coarse pool values (up to ~1e20) propagate into the level-0 restricted
  // residual and blow up h_relax.  This is harmless on fresh runs (pool near-zero)
  // but causes NaN kLa on every checkpoint restart.
  foreach_cell()
    if (!is_leaf(cell))
      res[] = 0.;
  foreach_face()
    g.x[] = 0.;
  foreach_face() {
    // Pe-limit beta at each face (coarse-level beta can exceed 2D/Delta after
    // restriction even when leaf-level values were Pe-limited in tracer_diffusion).
    double b_lim = beta.x[];
    double bmax  = 2.*D.x[]/Delta;
    if (b_lim >  bmax) b_lim =  bmax;
    if (b_lim < -bmax) b_lim = -bmax;
    g.x[] = D.x[]*face_gradient_x (a, 0) + b_lim*face_value (a, 0);
  }
  foreach (reduction(max:maxres)) {
    // -lambda[] = cm[]/dt;
    res[] = b[] + cm[]/dt*a[];
    foreach_dimension()
      res[] -= (g.x[1] - g.x[])/Delta;

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

    restriction ({D, beta, cm});
    struct HDiffusion q;
    q.D = D;
    q.beta = beta;

    // Touch non-leaf cells of c, r, D via foreach_cell() to extend qcc's stencil
    // analysis for these scalars to include non-leaf levels.  Without this, qcc
    // only sees r (and c, D) accessed via the EMBED-filtered leaf foreach() above,
    // giving an incomplete stencil.  tree_restriction → halo_restriction uses the
    // stencil to compute coarse halos; with the incomplete stencil, it produces
    // O(1e141) coarse residuals from O(0.07) leaf values, blowing up h_relax and
    // producing NaN kLa on every checkpoint restart.  These foreach_cell() reads
    // extend the stencil at qcc compile time — they are live (not dead code) but
    // are no-ops at runtime since they write nothing.
    foreach_cell()
      if (!is_leaf(cell)) (void)(c[] + r[] + D.x[]);

    // from mgstats poisson
    ///*
    scalar aa = c;
#if EMBED
    if (!q.embed_flux && aa.boundary[embed] != symmetry)
      q.embed_flux = embed_flux;
#endif //EMBED
    //*/

    mg_solve ({c}, {r}, h_residual, h_relax, &q);

    //*/
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
