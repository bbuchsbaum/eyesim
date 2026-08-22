# GazeWeave Transport v3 frozen simulation court

Protocol version: 3.0.2
Gate: GW-18.8
Seed: 20260822
Verdict: all mandatory gates passed

The court uses only generated public fixtures. It scores the same eight
exhaustive candidates, uniform design prior, two outer item folds, and common
bootstrap resamples for Transport v3, Transport v2, stabilized Replay,
registered density, density sigma 80 and 160, all raw MultiMatch dimensions,
and a nested ridge MultiMatch composite. Outer and inner item overlaps are zero
in every receipt. No retrieval field or private dataset is read.

## Mandatory Transport v3 ledger

| Gate | Observed | Threshold | Result |
|---|---:|---:|---|
| Representation invariance | all 7 checks | all | pass |
| Null mean information, absolute | 0.000886 bits | <= 0.02 | pass |
| Null top-1 minus prior | 0 | <= 0.02 | pass |
| Held-out ECE | 2.78e-17 | <= 0.10 | pass |
| Matched false-positive rate | 0.0625 | <= 0.075 | pass |
| Finite converged rate | 1.00 | >= 0.99 | pass |
| Global time-dilation change | 0 | <= 1e-8 | pass |
| Intact minus reversal information | 0.3595 bits | > 0 | pass |
| Partial coherent sequence intermediate | true | true | pass |
| Wrong-candidate registration guard | all rows | all | pass |

Coverage refinement changed candidate log scores by at most `9.65e-4` and
information by `4.60e-4` bits. Split/merge, unit conversion, candidate order,
batch size, and harmless timing dilation were exact within their frozen
tolerances. Reversal changed sequence evidence while the registered density
score was unchanged, demonstrating the intended fixed-density chronology
contrast.

## Predictive and resource evidence

Transport v3 achieved 2.675 mean signal information bits, 0.225 mean log loss,
rank 1 on every signal item, power 1.0 at the matched 0.0625 false-positive
rate, and convergence 1.0. The two complete shared folds took 8.421 and 7.937
seconds including every comparator, nested fitting, signal scoring, and null
scoring; serialized in-memory fold results were 2.56 MB each. Pair and
scope-matched backend performance remain governed by the separately frozen
GW-18.5 benchmarks.

The machine-readable court manifest, bounded ledger, predictive summaries,
operating characteristics, common-bootstrap intervals, perturbation table,
registration audit, resource table, fold receipts, serialized scientific
result, session information, and MD5 manifest are in
`gaze-weave-transport-v3-court-results/`.
