# params.json reference

All fields for `runs/<run_id>/params.json`.

## Motion

| Field | Type | Units | Default | Description |
|-------|------|-------|---------|-------------|
| `omega_b` | float | rad/s | required | Fundamental rocking angular frequency (how fast the bag rocks back and forth). 1 Hz = 2π ≈ 6.28 rad/s |
| `n_harmonics` | int | — | 1 | Number of active harmonics (1–3). A harmonic is a frequency component; 1 means pure sinusoidal rocking at `omega_b`. Vectors always padded to length 3 |
| `theta_max` | float[3] | degrees | [7,0,0] | Maximum rocking angle per harmonic. Index 0 is the fundamental (dominant) harmonic. Typical range: 2–15 degrees |
| `phi_angular` | float[3] | rad | [0,0,0] | Phase offset of rocking per harmonic (delays or advances the timing). Index 0 is **always forced to 0** — it is the global time reference |
| `omega_h` | float | rad/s | 0.0 | Horizontal translation frequency. Set to 0.0 to disable horizontal motion |
| `amplitude_h` | float[3] | m | [0,0,0] | Horizontal translation amplitude per harmonic (how far the bag slides sideways) |
| `phi_horizontal` | float[3] | rad | [0,0,0] | Phase offset of horizontal translation per harmonic |

## Geometry

| Field | Type | Units | Default | Description |
|-------|------|-------|---------|-------------|
| `geometry.a` | float | m | 0.25 | Bag half-width (horizontal semi-axis; half the total bag width) |
| `geometry.b` | float | m | 0.071 | Bag half-height (vertical semi-axis; half the total bag height) |
| `geometry.n` | float | — | 8.0 | Superellipse exponent controlling bag shape: n=2 gives an ellipse, n≥8 gives a rounded rectangle |
| `fill_level` | float | fraction | 0.5 | Fraction of bag volume filled with liquid (0 = empty, 1 = full). Typical range: 0.3–0.7 |

## Simulation control

| Field | Type | Units | Default | Description |
|-------|------|-------|---------|-------------|
| `run_id` | string | — | required | Unique label for this run; output files go to `runs/{run_id}/` |
| `fidelity` | int | — | required | Basilisk grid level; the computational mesh is 2^fidelity × 2^fidelity cells. Higher = more accurate and slower. See the [Fidelity guide](fidelity-guide.md) |
| `n_mix_cycles` | int | — | 80 | Number of complete rocking cycles to run before injecting oxygen. Used to let the flow field reach a steady state before kLa measurement begins |
| `t_end` | float | non-dim | 250.0 | When to stop the simulation (in non-dimensional time). Computed automatically by `simulate.py`, `chain.py`, and `sweep.py`; only set manually for custom runs |

## Checkpoint restart (set automatically by chain.py and sweep.py)

Do not set these manually — they are populated by the sweep scripts.

| Field | Type | Description |
|-------|------|-------------|
| `t_checkpoint` | float | Absolute non-dim time at which the restored checkpoint was saved (0 for fresh runs) |
| `omega_b_prev` | float | Rocking frequency of the segment that wrote the checkpoint; used to smoothly ramp to the new frequency |
| `theta_max_prev` | float[3] | Rocking amplitude of the previous segment |
| `phi_angular_prev` | float[3] | Rocking phase of the previous segment |
| `amplitude_h_prev` | float[3] | Horizontal translation amplitude of the previous segment |
| `phi_horizontal_prev` | float[3] | Horizontal translation phase of the previous segment |
| `omega_h_prev` | float | Horizontal translation frequency of the previous segment |

## Notes

**`t_end` is non-dimensional.** One non-dim time unit ≈ T_bio = L_bio / U_bio seconds,
where U_bio is the characteristic sloshing velocity (function of geometry, fill, and omega_b).
`simulate.py` computes `t_end = t_mix + t_buffer` automatically from `n_mix_cycles` and the
config's `t_buffer`; you only need to set `t_end` manually for custom runs.

**`phi_angular[0]` is always 0.** It is the global time-origin reference for the
rocking phase; it is physically redundant and is overridden at parse time.

**`omega_b` and `omega_h` are independent in the model.** On a physical platform
they are driven by the same motor (ω_h = ω_b). For pure rocking with no horizontal
translation, set `omega_h: 0.0` and all `amplitude_h` to zero.
