# Why non-dimensionalization

Every `t`, `omega_b`, and `tau_*` value in a raw `.dat` file is
non-dimensional — scaled to the bag's own characteristic length, velocity,
and time, not to seconds or Pascals. `postprocess.py` converts back to
physical units when it writes `results.json`. This page is why that
conversion exists and how it's derived, so the constants in
[params.json reference](../reference/params.md) and
[Output files reference](../reference/output-files.md) aren't just numbers
to memorize.

## The problem it solves

A rocking bioreactor's behavior depends on its size, fill level, and rocking
speed all at once — a bigger bag at the same frequency sloshes differently
than a small one. Comparing two simulations' raw velocity or shear-stress
numbers directly only makes sense if their bags are identical. Scaling
every quantity by the bag's *own* characteristic velocity and length makes
the numbers comparable across different geometries and frequencies — one
non-dimensional time unit means "one characteristic sloshing timescale has
passed," regardless of whether that bag is small and fast or large and slow.

## The actual derivation (`BioReactor.c`)

```c
H_bio  = L_bio*Ly;
V_bio  = L_bio/4*(H_bio + 0.5*L_bio*tan(Th_max));
U_bio  = V_bio/(H_bio*0.5)/T_per;   // characteristic velocity scale
T_bio  = L_bio/U_bio;               // characteristic time scale
```

`L_bio` is `geometry.a` (the bag half-width) — everything else derives from
it. `V_bio` is an estimate of the characteristic sloshing volume swept per
rocking period, built from the bag's geometry and its maximum tilt angle
(`Th_max`). Dividing that by half the bag height and by the rocking period
`T_per` gives `U_bio`, a characteristic sloshing *velocity* — and once you
have a characteristic length and velocity, `T_bio = L_bio / U_bio` falls out
as the characteristic *time*.

Everything the solver reports in non-dimensional form is scaled against
these three: a non-dimensional time `t` is `t_physical / T_bio`; a
non-dimensional velocity is `u_physical / U_bio`.

## Converting the KPIs that matter

**kLa** (mass-transfer rate) is inherently a rate — physical units h⁻¹ — so
it converts via the time scale alone:

```
kLa_physical = kLa_nd × 3600 / T_bio
```

**Shear stress (τ)** converts differently, because `BioReactor.c` sets
`rho1=1, mu1=1/Re_w` internally, which makes the *dimensionless* group
`τ_nd = τ_dim / (ρ_w U_bio²)` — not `τ_dim × T_bio / μ_w`, which is the more
commonly assumed form and will silently give you the wrong answer by a
factor related to the Reynolds number if you use it here:

```
τ_physical [Pa] = τ_nd × ρ_w × U_bio²        (U_bio = geometry.a / T_bio)
```

This is the exact conversion `postprocess.py` uses for every `tau_*` key in
`results.json`, and it's what makes those values directly comparable to
Kim et al.'s dimensional reported values — see
[Validating against Kim et al. (2024)](kim-et-al-validation.md).

## What this buys you in practice

- `t_buffer` (the kLa measurement window) is sized in non-dimensional time,
  so the same config value works across different `omega_b` — see the
  [Glossary](../glossary.md) entry for the actual sizing rule of thumb.
- `n_mix_cycles` (rocking cycles before oxygen injection) is a cycle count,
  not a time — it's automatically consistent across frequencies because a
  "cycle" already encodes the period.
- Grid-convergence and fidelity comparisons ([Fidelity guide](../reference/fidelity-guide.md))
  are meaningful specifically *because* every fidelity level is solving the
  same non-dimensional problem — only the mesh resolution changes.
