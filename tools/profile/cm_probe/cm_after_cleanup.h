/**
# Arm "once"/"every" of the cm probe: restrict cs (= cm) right after fractions_cleanup() of event adapt
(centered.h), at step 0 only (CM_ONCE) or every step. A second event adapt runs after the one of centered.h
(same-name events run in definition order). Diary 2026-11-23, pre-registration "cm restriction". */

#if CM_ONCE
event adapt (i = 0) {
#else
event adapt (i++) {
#endif
  restriction ({cm});
}
