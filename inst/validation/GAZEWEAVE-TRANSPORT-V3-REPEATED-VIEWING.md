# Transport v3 full-cohort repeated-viewing court

Protocol: `gazeweave-transport-v3/3.0.2`
Court seed: `20260823`
Shuffle seed: `20260824`
Status: **advance**
Retrieval outcomes at decision time: **sealed**

## Decision

Transport v3 passes the prospectively frozen study positive-control rule. On
the full 1,295-target P1-to-P4 comparison, mean held-out information is
`0.28464` bits (crossed participant/item bootstrap 95% interval `0.18206` to
`0.38573`). The paired intact-minus-mean(reversal, shuffle) contrast is
`0.05237` bits (`0.00720` to `0.09536`). Both estimates are positive and both
intervals exclude zero.

This is a sensitivity result, not evidence that Transport v3 is superior to
the fair comparators. Its P1-to-P4 information is effectively tied with
Transport v2 (`0.28869`), registered density (`0.28040`), density sigma 80
(`0.27069`), and Replay (`0.25714`) at the resolution of their crossed
intervals.

## Repetition trajectory

| Comparison | Trials | Mean information (bits) | Crossed 95% interval | Mean rank | Top-one credit | Alignment convergence |
|---|---:|---:|---:|---:|---:|---:|
| P1 to P2 | 1,295 | 0.46321 | 0.32347 to 0.59790 | 2.029 | 0.517 | 0.99938 |
| P1 to P3 | 1,295 | 0.34113 | 0.23076 to 0.45066 | 2.175 | 0.435 | 0.99954 |
| P1 to P4 | 1,295 | 0.28464 | 0.18206 to 0.38573 | 2.244 | 0.410 | 0.99985 |
| P2 to P3 | 1,295 | 0.30631 | 0.18375 to 0.41722 | 2.167 | 0.444 | 0.99985 |
| P3 to P4 | 1,295 | 0.30586 | 0.19616 to 0.41031 | 2.224 | 0.411 | 0.99969 |

The initial-to-later trajectory declines, as did the previously frozen Replay
trajectory. This is retained as a descriptive habituation-compatible pattern;
no monotonic increase was required or tested.

## Common P1-to-P4 panel

All rows below use the same 1,295 targets, K=5 candidate sets, outer
participant/item folds, registration opportunity, uniform prior, and crossed
bootstrap plan.

| Method | Mean information (bits) | Crossed 95% interval | Top-one credit |
|---|---:|---:|---:|
| Transport v2 | 0.28869 | 0.16929 to 0.39947 | 0.439 |
| Transport v3 | 0.28464 | 0.18206 to 0.38573 | 0.410 |
| Registered density ridge | 0.28040 | 0.16340 to 0.39211 | 0.425 |
| Density sigma 80 | 0.27069 | 0.17219 to 0.36541 | 0.425 |
| Replay | 0.25714 | 0.14877 to 0.35815 | 0.422 |
| MultiMatch learned ridge | 0.21000 | 0.11954 to 0.29745 | 0.398 |
| MultiMatch position EMD | 0.21003 | 0.12097 to 0.29516 | 0.405 |
| Density sigma 160 | 0.17258 | 0.09587 to 0.24215 | 0.383 |
| MultiMatch position | 0.16597 | 0.08221 to 0.24063 | 0.377 |
| MultiMatch direction | 0.01386 | -0.00859 to 0.03444 | 0.249 |
| MultiMatch vector | 0.01158 | -0.01528 to 0.03724 | 0.269 |
| MultiMatch length | 0.00385 | -0.00767 to 0.01515 | 0.230 |
| MultiMatch duration | -0.00232 | -0.00892 to 0.00333 | 0.197 |

The fixed raw MultiMatch dimensions are displayed through their response-blind
training calibration so their information scores share the same probability
scale. Their raw similarity values remain available as diagnostics in the
local comparator checkpoints.

## Calibration, isolation, and numerical evidence

- Every one of the 20 outer task/fold fits used inner out-of-fold
  episode-scale temperature calibration with response-blind effective-fixation
  shrinkage; training sizes were 320 to 328 and every inner match-key overlap
  count was zero.
- Held-out ECE was `0.01332`, `0.00690`, `0.00545`, `0.01289`, and `0.00987`
  for P1-P2, P1-P3, P1-P4, P2-P3, and P3-P4 respectively.
- Each outer fold had zero participant overlap and zero item overlap. Every
  candidate set contained exactly five candidates and one true item.
- Alignment convergence was 0.99938 to 0.99985 across intact tasks, above the
  frozen 0.99 floor. P1-P4 row-level all-five-candidates convergence was
  0.99923.
- Candidate scoring for all 9,065 intact/control rows took 1,710.84 seconds in
  aggregate. A separately profiled representative P1-P4 fold, including warp
  fitting and nested calibration, took 191.85 seconds: 42.99 seconds for fit
  and calibration and 148.81 seconds for 960 intact/control score rows
  (0.155 seconds per row).

## Frozen-control and solver audit

Reversal and shuffle reorder complete fixation tuples, preserving coordinates,
durations, and therefore spatial-duration density. Shuffle is deterministic by
participant/item key. Controls are scored with the intact fold's warp and
calibration; they never refit the model.

The first full pass was rejected before interpretation because candidate-set
strata incorrectly forced identity calibration. A second diagnostic pass
revealed only 97% alignment convergence: when a native failure entered the
pure-R oracle, a stalled coverage-continuation start discarded the canonical
structural start. The repair makes that structural start a fallback only after
continuation fails and prefers no new objective or threshold. Frozen
native/reference parity and the aggregate v3 test court passed after the
repair. All provisional checkpoints were deleted and the complete court was
rerun from empty state. The values in this report come only from that final
run.

## Privacy boundary

The script reads only study-presentation fixations at this gate. Retrieval
response, correctness, confidence, saliency, old/lure, and phase values were
not read or merged. The two restricted CSVs and all participant-linked
repeated-viewing checkpoints remain local under explicit `.gitignore` rules.
No trial-level score is included in the repository or package artifact.
