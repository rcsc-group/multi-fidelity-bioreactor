# Meeting 2026-10-29

## DFD talk deck (draft)
- `2026-10-29_dfd_deck_draft.pdf`: title, definitions, methods. Other slides are empty. The title slide has since been made light; this PDF is the earlier render.
- Source: `talks/dfd2026/`. Next: Motivation (three storyboards to choose from).

## Kim's offline scalar solver (dev/rocking-bioreactor-2d-offline, read from the code)
- It stores only the cell-centred velocity of one cycle. It rebuilds the face velocity as `fm*face_value(u.x,0)` (`TracerTest_test.c`, line 815). It stores no true face flux and does no projection.
- Test in OUR solver (L8, 32.5 rpm, 7 deg, 22 cycles; one run, one condition): the true recorded flux has `d/max|uf|` = 1.4e-3 (median). The rebuilt field has 5.8e-2 (about 40 times more).
- Effect on the scalars: kLa at 10 % saturation 56.3 against 45.7 per hour. Tracer in the liquid falls to 0.80 of its start after 20 cycles (0.99 with the true flux). Mixing time changes by about 4 %.
- Open: his own code is being run now (producer N=128, then his offline solver) to check his reported errors. Not finished at the time of writing.
- Not for the talk yet.
