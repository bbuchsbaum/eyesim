# GazeWeave multi-presentation Replay measurement verdict

Status: complete and frozen for behavioral testing
Protocol version: 1
Seed: 20260828
Mote: `bd-01M0GTPJWX3TP6S65QW6RP2FJ8`

## Verdict

The measurement court passed every prospectively frozen gate. The behavioral
score is now fixed as `replay_all4_shrunk_exhaustive`: an equal-weight mixture
of four separately aligned study paths, exhaustive candidates within the
participant and held-out item fold, training-only temperature calibration,
and response-blind shrinkage toward the candidate prior from retrieval
effective fixation count.

Passing the measurement gates does not establish retrieval reinstatement. On
the full 3,000-ms recognition interval, item-identifying evidence was small and
uncertain. The important positive result is that the same Replay family showed
large repeated-viewing sensitivity and a reliable loss of information under
complete temporal reversal. The weak retrieval result therefore cannot be
explained simply by a wholly insensitive implementation.

## Full-cohort measurement results

All methods scored the same 2,055 trials from 45 participants against exact
candidate pools containing 5--33 items.

| Method | Mean bits | Crossed 95% interval | Mean log loss | Top-one credit |
|---|---:|---:|---:|---:|
| P4 Replay | 0.0054 | [-0.0077, 0.0175] | 3.1826 | 0.0516 |
| All-four Replay, unshrunk | 0.0185 | [-0.0035, 0.0388] | 3.1736 | 0.0710 |
| All-four Replay, shrunk | 0.0209 | [-0.0047, 0.0445] | 3.1719 | 0.0710 |
| Density sigma 80 | 0.0217 | [-0.0010, 0.0438] | 3.1713 | 0.0662 |
| Density sigma 160 | 0.0186 | [-0.0016, 0.0377] | 3.1735 | 0.0633 |

All-four unshrunk Replay improved over P4 by 0.0130 bits, but its crossed 95%
interval included zero, [-0.0034, 0.0273]. Reliability shrinkage added 0.0024
bits, also uncertain, [-0.0019, 0.0069]. Density sigma 80 had the best mean log
loss by a very small margin. GazeWeave did not outperform density on retrieval
measurement in this court.

## Sensitivity and temporal-order control

Presentation 4 was identified from the equal-weight mixture of presentations
1--3 with 0.7014 bits of held-out information, crossed 95% interval [0.5085,
0.8955]. Reversing presentation 4 while preserving every location and fixation
duration reduced information by 0.0312 bits, [0.0153, 0.0482]. Thus Replay is
sensitive to both episode identity and local chronology when they are present.

## Fixation-count assessment

Retrieval effective fixation count had only a weak association with unshrunk
Replay information: standardized slope 0.0645, crossed 95% interval [-0.0175,
0.1500], Spearman rho 0.0405. Within an episode, presentations with more
effective fixations received somewhat greater posterior responsibility
(Spearman rho 0.151), but equal presentation weights remain the declared prior.
No fixation-count weighting was introduced.

## Calibration and computation

Every all-four calibration fold selected a non-negative reliability parameter
and achieved training log loss no worse than its own temperature-only boundary
solution. All scores were finite, all forward--backward fits converged, and
every participant and item was disjoint across outer training and evaluation.

Across four folds, elapsed time was approximately:

- P4 Replay: 191 seconds;
- all-four Replay: 725 seconds;
- density sigma 80 and 160 together: 14 seconds;
- intact and reversed study controls: 798 seconds.

The private checkpoints and trial-level candidate tables remain in the
Git-ignored local results directory. Only the frozen algorithm may now enter
the prespecified behavioral court; no saliency subgroup, time-window,
nonlinear, or ordinal search is authorized by this verdict.
