# Running the GazeWeave comparison

From the package checkout:

```r
devtools::load_all()
source("inst/validation/gaze-weave-comparison.R")

validation_result <- run_gaze_weave_validation(
  n_participants = 24,
  n_items = 6,
  sample_sizes = c(8, 16, 32),
  bootstrap_repetitions = 500,
  randomization_draws = 199,
  signal_probabilities = c(0.25, 0.5, 0.75, 1),
  seed = 20260815,
  output_dir = "inst/validation/results",
  verbose = TRUE
)
```

Read `PROTOCOL.md` before `REPORT.md`. The `results` directory contains the
participant-level retrieval measures, aggregate retrieval and power tables,
ground-truth diagnostic recovery, the exact run configuration, and session
information.

The density baseline uses the duration-weighted Gaussian time marginal with the
same three spatial scales as GazeWeave. MultiMatch dimensions are never averaged.
Both comparator families are scored raw and after the same independently fitted
registration used by GazeWeave.

## GazeWeave v2 evidence program

The current evidence program is documented by:

- `GAZEWEAVE-V2-CONTRACT.md`: frozen estimands and interpretation contract;
- `GAZEWEAVE-V2-COURT.md` and `gaze-weave-v2-results/`: synthetic and
  invariance court;
- `GAZEWEAVE-REAL-CONTROLS.md`, `GAZEWEAVE-REAL-VERDICT.md`, and
  `gaze-weave-real-results/`: original participant-by-item real court; and
- `GAZEWEAVE-REPLICATION-CONTROLS.md`,
  `GAZEWEAVE-REPLICATION-VERDICT.md`, and
  `gaze-weave-replication-results/`: prospectively frozen disjoint replication
  of resolution-invariant Replay calibration; and
- `GAZEWEAVE-PROBE-DELAY-CONTROLS.md` and
  `GAZEWEAVE-PROBE-DELAY-VERDICT.md`: phase-resolved, same-study probe and
  post-probe persistence court. Its raw inputs and provisional results remain
  local and Git-ignored; and
- `GAZEWEAVE-RECOGNITION-SENSITIVITY.md`,
  `GAZEWEAVE-RECOGNITION-SENSITIVITY-VERDICT.md`, and
  `gaze-weave-recognition-sensitivity.R`: frozen secondary analysis of whether
  full recognition-trial information varies with saliency or response
  correctness. Its trial-linked outputs remain local and Git-ignored; and
- `GAZEWEAVE-RECOGNITION-FULL-COHORT.md`,
  `GAZEWEAVE-RECOGNITION-FULL-COHORT-VERDICT.md`, and
  `gaze-weave-recognition-full-cohort.R`: maximal 36-participant, 1,295-trial
  combined-window court with fixed K=5 candidates, raw density at sigma 80 and
  160 pixels, all six MultiMatch dimensions, learned composites, and elastic
  matching. Participant-linked outputs remain local and Git-ignored.

The original real court is retained as negative calibration evidence. The
disjoint replication repairs the repeated-viewing failure but still fails the
blank-screen imagery gate. `GAZEWEAVE-DECISION.md` records the resulting
no-default, no-superiority decision. The probe-delay court also returns
`not_supported`: it recovers repeated-viewing rank structure but finds no
reliable item information during the full or late post-probe delay. The
secondary recognition-sensitivity court is also `not_supported`: neither
GazeWeave engine showed a reliable saliency, correctness, or interaction effect
in the complete probe-plus-delay interval. The maximal-cohort court supersedes
that eight-participant pilot for same-study sensitivity claims and reaches the
same `not_supported` conclusion. Transport was descriptively competitive but
did not establish better power than density or MultiMatch.
