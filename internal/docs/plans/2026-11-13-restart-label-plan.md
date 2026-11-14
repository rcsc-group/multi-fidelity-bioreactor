# Restart label and frequency ramp

Decision (user, 2026-11-13): "Label + freq ramp". A restart ramps from the forcing that made the restored
state to the forcing in params.json. The dump carries that source forcing in its label; params.json holds the
target only. Silent defaults are not allowed.

## Why

- Today the source forcing comes from `*_prev` keys in params.json. Every script that writes a restart must
  copy them; when a script forgets, the solver ramps from rest without a word. This happened to all 8 L9
  campaign warm starts and to the fig9 extension (+8..+49 % velocity at the seam).
- params.json + checkpoint.phase hold the same facts twice (user: "smells as duplicate info").
- The frequency is not ramped: a cross-case restart jumps omega_b in one step.

## Design

1. **Label** = `checkpoint.phase`, one JSON line, written at every dump (planned and early). The name stays:
   the staging scripts and the slurm template already move it with the dump. Keys:
   - `t` (dump time, run units), `cycles` (forcing phase at the dump, in periods since the run's origin,
     unit-free), `cycles_h` (same, horizontal shaking);
   - units of the dumped fields: `U_bio` [m/s], `T_bio` [s];
   - forcing at the dump instant (instantaneous, so a dump inside a ramp is right): `omega_b` [rad/s],
     `n_harmonics`, `theta_max` [deg], `phi_angular` [rad], `omega_h` [rad/s], `amplitude_h` [m],
     `phi_horizontal` [rad];
   - `Th` [rad], `Th_d` [rad/s] at the dump (for the continuity test);
   - `T_per_st`, `t_phase_offset` kept for old readers.
2. **Restart** (not REPLAY_FLUX): read the label from the working directory.
   - no label -> error, exit 1.
   - label with forcing + any `*_prev` in params -> error (one source of truth).
   - legacy label (no `omega_b`) + `*_prev` in params -> WARNING on stderr, use `*_prev` (old rule).
   - legacy label without `*_prev` -> error.
3. **Ramp**: amplitude, phase and frequency go from the source to the target with the smooth step
   alpha = 3x^2 - 2x^3 over N_RAMP_CYCLES target periods. Angle phase Psi(t) = 2 pi cycles0 + integral of w;
   integral of alpha = D (x^3 - x^4/2). Th_2d adds A k w' cos. Same source = target -> no ramp, bit-for-bit
   the old continuation.
4. **Units**: velocity/pressure rescale uses su = U_bio(label) / U_bio(run).
5. **Times on the phase**: stop time = first t with Psi(t) = 2 pi n after t_checkpoint + t_end; release
   time t_mix = t where Psi - Psi0 = 2 pi n_mix_cycles. Found by bisection on the monotone Psi.
   With no ramp in w they equal the old formulas.
6. **Restart record**: `restart_forcing.json` at restart: Th, Th_d (SI) at t_checkpoint from the new run's
   forcing. Guard: equal to the label's Th, Th_d (continuity of angle and rate).
7. Python: remove `*_prev` writers (chain.py, submit_continue.py, campaign, bundle cases, run scripts used by
   tests); simulate.py guard: a restart needs a label next to the dump and must not set `*_prev`.

## Tests (computational experiments)

- Bundle L6 rerun with the new src: all guards; `cont`, `kill_resume`, `xl6` with no `*_prev` must match
  cold as before; S4 `caseB_warm` re-scored (frequency now ramps).
- New guards: `test_restart_forcing_is_continuous[cont|kill_resume|xl6|caseB_warm]`;
  `test_restart_refuses_prev_keys`, `test_restart_refuses_no_label` (two short L4 runs, a few steps).
- Mutant `prev_from_rest` (label forcing ignored -> ramp from zero) must be caught.
- Python unit tests for the writers and the submit guard.
