# Geometry shape test (n=8 vs n=60) -- INVALID, superseded

**Do not use these results.** `n8_control` (run `46acc8f0`) and `n60_sharp`
(run `c7c20700`) produced bit-identical `results.json` (every field, full
double precision) because both values hit the same code branch in
`src/BioReactor.c`:

```c
// n >= 8 -> perfect rectangle (avoids pow() singularities at sharp corners)
if (params.geometry_n >= 8.)
  solid (cs, fs, intersection(a_nd - fabs(x), b_nd - fabs(y)));
else
  solid (cs, fs, 1. - pow(fabs(x/a_nd), params.geometry_n)
                    - pow(fabs(y/b_nd), params.geometry_n));
```

Any `n >= 8` collapses to the exact same sharp-rectangle implicit function,
independent of the actual exponent value. `n=8` and `n=60` are therefore
identical inputs to the solver -- this test compared a condition against
itself.

**What this does settle:** our default geometry (`n=8`) is already a
mathematically sharp rectangle (`intersection(a-|x|, b-|y|)`), not a rounded
superellipse as the pure math formula would suggest for finite n. The
"~8% corner chamfer" concern raised before running this test does not apply
here -- that arithmetic is for the `pow()` branch, which `n=8` never reaches.

**What is still a live, untested difference:** Kim et al.'s actual upstream
driver (`github.com/rcsc-group/BioReactor/DriverCodes/BioReactor.c`) does not
embed the x-direction walls at all -- they are plain grid-aligned domain
edges (`u.n[left]=dirichlet(0.)`, `u.n[right]=dirichlet(0.)`). Embedding
(`solid()`/cs/fs cut-cells) is used only for the y-direction (top/bottom
plates). Our fork routes all four walls through the embedded-boundary
machinery, even in the `n>=8` sharp-rectangle branch. This is the same
family of issue as the missing `fm`/`cm` metric-factor bug found in the
shear-stress stencil (embedded boundaries need those corrections; plain
box boundary conditions don't). Testing this properly requires an actual
source change (a genuine "plain box BC on x, embed only on y" code path),
not a parameter sweep -- pending user decision on whether to add it.
