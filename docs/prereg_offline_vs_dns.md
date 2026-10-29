# Pre-registration: do the hydrodynamics-only scalar solvers converge to their own DNS? (2026-10-29)

Written before any of these runs. The thresholds below do not change after results arrive.

## Question
Kim (meeting 2026-10-29): neither offline workflow (his rebuild, our stored flux) appears to converge to the full DNS. Our 2026-10-29 figures cannot settle this: the two solvers ran different conditions (ours 32.5 rpm, his 20 rpm). Retracted as evidence.

## Fixed setup
- Condition: 20 rpm, theta 7, a = 0.25 m, fill 0.5 (Kim's default case; same geometry as our params.json).
- Grid: ours L8 (adaptive, max 256), his N = 256 (uniform). Second grid only in T4.
- Protocol, both solvers: release tracer and oxygen at cycle 20 (after his amplitude ramp ends at cycle 10), then 3 cycles. The offline solver replays exactly those 3 cycles. No looping, except in factor W.
- Reference for each offline run = the DNS of the same solver, same run, same window.

## Metrics (over cycles 0-3 after release)
- M1 tracer mass in the liquid: max |m(t)/m(0) - m_DNS(t)/m_DNS(0)|
- M2 mixing chi(t) = 1 - var/var0 of the upper-half tracer: max |chi - chi_DNS|
- M3 liquid oxygen: max |O/O_DNS - 1|, in each solver's own definition
- M4 kLa at C* = 0.10 (only if the window reaches C* = 0.10; else M3 only)

## Noise floor F (measured, not assumed)
Two DNS of the same condition that differ only in the MPI rank count: the round-off is different and nothing else is. F = M1-M3 between them. Both solvers.

## Tests
- **T0 DNS agreement (ours vs his DNS).** Pass: kLa within 10 % (if available), max |chi_ours - chi_his| <= 0.05, and tracer mass drift <= 0.5 % in both. Fail: the two DNS solve different problems, so they must be explained before any offline result is compared across solvers.
- **T1 exact replay (both solvers).** The offline solver reads the DNS's own per-step projected face flux, with no hold and no loop. Pass: M1-M3 <= max(F, 1e-3). Fail: an implementation defect in that offline solver. For Kim, this needs a test build that dumps uf every step. This is a scratch copy, not his code.
- **T2 one-factor ladder (both solvers), each factor added alone to T1:**
  - S (space): face field rebuilt as `fm*face_value(u,0)`. Ours: REC_KIM_FLUX. His: line 815.
  - H (time): the field is held for T/100 between snapshots.
  - W (wrap): one recorded cycle looped 3 times.
  - Factor material iff M1 > max(2F, 0.5 %) or M2 > max(2F, 0.02) or M3 > max(2F, 2 %).
- **T3 his workflow as shipped** (S+H+W) vs his DNS, and ours as shipped (W only) vs our DNS. "Converges to the DNS" iff no metric exceeds the material threshold.
- **T4 grid (only if T2 finds a material factor):** repeat that factor at ours L9 / his N = 512. A spatial rebuild error should shrink with the grid. If it does not, the error has another cause.

## Power and limits
- The offline-vs-own-DNS comparisons are deterministic. Chaos enters only through F, which is measured.
- One condition. A factor that matters only at other rpm is not detected; if T3 is within a factor of 2 of a threshold, repeat at 32.5 rpm.
- Three cycles. Slow drifts (tracer leak about 1 % per cycle) are detectable; a smaller drift is not.
- The definitions of liquid oxygen differ between the solvers. M3 is relative within each solver; across solvers (T0) only chi and kLa are compared.

## Gates before submission (code, not docs)
- G0: ours, `n_rec_periods = 3` records 3 cycles and the replay runs them without wrap (src/BioReactor.c lines 691, 1153).
- G1: his test build dumps uf every step; a 6-step dry run is byte-identical in the scalars to the stock build.
- G2: smoke test at the coarsest grid for each run type before any 256 run.
