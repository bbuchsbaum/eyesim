# GazeWeave resolution-invariant real-data replication

Status: frozen before outcome scoring
Protocol version: 2
Date frozen: 2026-08-15
Final seed: `20260818`
Mote gate: `bd-01M049BMXCY58VEPY612M48TSN`

## Purpose

This court evaluates one algorithm correction identified by a completion
audit. Replay previously exponentiated a total log likelihood from a fixed
normalized-duration grid and regularized its temperature toward one. Candidate
contrasts in that total grow approximately with `grid_size`, so the prior was
not invariant to the numerical resolution of the same gaze-time occupation.

The corrected candidate score is mean log predictive density per normalized
duration bin:

\[
\bar\ell(Y\mid X)=\ell(Y\mid X)/L,
\]

where `L` is `grid_size`. The raw total HMM log likelihood remains an alignment
diagnostic. Inner-fold calibration estimates a residual temperature on
`bar(l)` with its fixed log-temperature prior centred at one. Transition-model
selection is unchanged because every training trial in a specification has
the same `L`.

The original synthetic and real-control results remain immutable. This is a
prospectively frozen replication on participants and items absent from the
original court. It uses the same parent dataset and is therefore an internal,
not external-dataset, replication.

## Data-only cohort construction

Preprocessing and eligibility are exactly those in
`GAZEWEAVE-REAL-CONTROLS.md`. Before any model was fit or outcome inspected,
the selector:

1. removed the original eight participants and eight items;
2. retained participants with at least 100 eligible participant-item paths;
3. sampled four young and four older participants with seed `20260818`;
4. removed the original items from the complete-item pool; and
5. sampled eight items complete for every selected participant.

The frozen participants are `9`, `18`, `111`, `130`, `300`, `302`, `316`, and
`325`. The frozen items are `4`, `11`, `24`, `62`, `63`, `80`, `92`, and `101`.
The eligible unused pool contained 39 young and 21 older participants; 79
unused items were complete across the selected participants. These identities
were derived only from phase presence, fixation-count eligibility, and the
declared exclusions.

## Evaluation contract

- Tasks are first-to-fourth repeated viewing and fourth-study-to-blank-screen
  imagery.
- The four Cartesian participant-fold by item-fold cells exclude every
  evaluated participant and item from nuisance fitting and model selection.
- Each held-out path is scored among four held-out within-participant item
  templates under a uniform prior.
- Conditions are observed signal, fixation-order shuffle at fixed spatial and
  duration occupation, within-participant wrong item, and participant-generic
  gaze.
- Models, preprocessing, warps, calibration, and ridge weights use only
  unaltered signal rows in outer training data.
- Candidate pools, folds, warp policy, priors, and endpoints are identical for
  every method.

All model specifications are frozen to the original real court: isotropic
contraction plus translation; Transport v2 coverage, omission, spatial,
chronology, continuation, and iteration grids; Replay 48-bin duration grid,
maximum skip two, Student emissions, transition grid, and inner item-fold
calibration; registered multiscale density; all six MultiMatch dimensions;
elastic consensus; and the same nested-ridge grids. The sole algorithm change
is Replay candidate-score normalization described above.

## Endpoints and uncertainty

The primary endpoint is held-out `gaze_info_bits`. The court also reports log
loss, top-one fractional credit, mean rank, Brier score, five-bin ECE, a
1000-draw crossed participant-by-item bootstrap, power at the empirical 95th
percentile of wrong-item plus generic-gaze information, empirical null error,
order-ablation contrasts, registration diagnostics, posterior/coverage
diagnostics, convergence, runtime, memory, and the full-window Replay
sensitivity analysis.

## Frozen gates

The replication passes only if all seven machine gates pass:

1. Replay and Transport each have positive repeated-viewing mean information
   and top-one credit above `1/4`.
2. At least one GazeWeave engine has positive blank-screen mean information and
   top-one credit above `1/4`.
3. At least one GazeWeave engine loses repeated-viewing information after order
   shuffling, while registered density is invariant within numerical tolerance.
4. Each engine's pooled empirical false-positive rate is at most `0.075`.
5. Registration increases density null top-one credit by no more than `0.05`.
6. At least 99% of engine candidate sets converge.
7. Fold receipts show zero participant and item overlap between training and
   evaluation.

Replay earns a directional default only if all gates pass, its pooled signal
log loss is no more than `0.05` nats worse than Transport, its top-one credit
is no more than `0.05` worse, and it is not materially worse than the strongest
fair baseline (at most `0.10` log-loss and top-one deficits). Superiority
requires crossed-bootstrap 95% intervals favouring GazeWeave for both
information and log-loss advantage against the strongest baseline.

No threshold or model parameter may change after final scoring. Failure is a
bounded scientific result and leaves the epic open; it is not permission to
retune this replication.
