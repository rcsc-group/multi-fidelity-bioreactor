# Meeting 2026-10-29

## 1. Kim's offline scalar solver: what he does, what we do, the gap
Repo: `rcsc-group/BioReactor3D`, folder `dev/rocking-bioreactor-2d-offline/` (read from the code, not the README).

**What he does** (the face velocity that moves the tracer and oxygen is rebuilt from cell values):
1. Producer stores only cell-centred velocity: `codes/flowfield_out_MK.h:277-279` (`a_ux[k] = u.x[]`, `a_uy[k] = u.y[]`, `a_f[k] = f[]`), called at `codes_bio_run/BioReactor.c:595` (`ff_dump_binary("Data_dump", N, t, u, f, omega)`).
2. Offline solver rebuilds the face field: `codes/TracerTest_test.c:815`: `uf.x[] = fm.x[]*face_value(u.x,0);`
3. The "face velocity" in his ASCII dump is also an average of cell values: `codes_bio_run/BioReactor.c:569-573` (`uf_MK.x[] = 0.5*(u.x[] + u.x[-1])`).
4. Reading a stored `uf` is switched off: `TracerTest_test.c:545-546`, commented with "MK changed - 260715".
5. No projection anywhere in `TracerTest_test.c` (searched for project, div, poisson).

**What we do** (we store the flux the projection produced, and replay exactly that):
1. Record: `src/BioReactor.c:479-493` (`rec_gather` writes `uf.x[]` and `uf.x[1]` of every cell), called from `rec_write_step` (line 601) at line 3336.
2. That `uf` is divergence-free because Basilisk's `navier-stokes/centered.h:430` runs `project (uf, p, alpha, dt, ...)` before it is used.
3. Switch to emulate his choice: `src/BioReactor.c:444` `REC_KIM_FLUX`, same formula as his line 815, applied to the same run.

**The gap, same metric in all three** (`fig_divergence_compare.png`): max over fully open cells of the net face-flux imbalance, divided by max |uf|:

| | median | max |
|---|---|---|
| ours, projected flux (L8) | 1.4e-3 | 3.3e-3 |
| ours, rebuilt like him (L8, same flow) | 5.8e-2 | 1.4e-1 |
| his own solver (N=128) | 1.2e-1 | 3.3e-1 |

**Effect on the scalars**
- In our solver, L8, 32.5 rpm, 7 deg, 22 cycles, one run: kLa at 10 % saturation 56.3 against 45.7 per hour; tracer in the liquid 0.80 of its start after 20 cycles (0.99 with the true flux); mixing time about 4 % lower.
- In his solver against his own producer (N=128, 3 cycles): tracer in the liquid falls 2.3 % (cycle 1), 4.8 % (cycle 2), 7.2 % (about 2.9 cycles) while the producer stays flat. Oxygen in the liquid is 64 % (cycle 1) and 77 % (cycle 2) higher than the producer's. His README says tracer grows about 5 % and oxygen is 7-8 % slower: the code run shows the opposite sign.

**Caveats:** N=128 not 512; 3 cycles; one run each; his producer flow after the dump window is not the replayed cycle; the amplitude ramp in his case ends at cycle 10, so the stored cycle may not be periodic; his producer job was still running at the time. Not for the talk.

## 2. DFD talk deck (draft)
- `2026-10-29_dfd_deck_draft.pdf`: title, definitions, methods. Other slides are empty. The title slide has since been made light; this PDF is the earlier render.
- Source: `talks/dfd2026/`. Next: Motivation (three storyboards to choose from).
