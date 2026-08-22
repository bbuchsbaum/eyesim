# GazeWeave comparative validation report

Run date: 2026-08-15. The canonical configuration is in
`results/run-config.csv`; participant-level and aggregate results are retained
beside this report.

## Verdict

This simulation supplies bounded evidence that GazeWeave works as intended and
can improve discrimination and power when item identity is carried by temporal
order or an incomplete replay. It does not support a universal superiority
claim.

- Registered geometry: GazeWeave, registered density, and several registered
  MultiMatch dimensions all achieved perfect retrieval. Under the 25% signal
  stress at n = 8, estimated power was 0.930 for GazeWeave, 0.952 for registered
  density, and 0.774 for Holm-corrected registered MultiMatch. GazeWeave is not
  the best method in this regime.
- Order at fixed density: GazeWeave and MultiMatch vector similarity both
  achieved perfect retrieval; density was exactly at chance. Under 25% signal
  at n = 8, GazeWeave power was 0.934 versus 0 for registered density and 0.738
  for Holm-corrected registered MultiMatch.
- Partial replay: GazeWeave achieved pairwise AUC 0.990 and top-1 credit 0.951.
  The strongest registered MultiMatch dimension, position, achieved AUC 0.935
  and top-1 0.750; registered density remained at chance. Under 25% signal at
  n = 8, GazeWeave power was 0.874 versus 0 for density and 0.510 for
  Holm-corrected MultiMatch.

Thus GazeWeave meets the frozen superiority rule for partial-replay
discrimination and for low-prevalence power in the two temporally informative
scenarios. It does not exceed all comparators for geometry or preserved-order
retrieval.

## Retrieval at full signal

| Scenario | Method | Pairwise AUC | Top-1 credit |
|---|---|---:|---:|
| Registered geometry | GazeWeave | 1.000 | 1.000 |
| Registered geometry | Registered density | 1.000 | 1.000 |
| Registered geometry | Registered MultiMatch position | 1.000 | 1.000 |
| Order at fixed density | GazeWeave | 1.000 | 1.000 |
| Order at fixed density | Registered density | 0.500 | 0.167 |
| Order at fixed density | Registered MultiMatch vector | 1.000 | 1.000 |
| Partial replay | GazeWeave | 0.990 | 0.951 |
| Partial replay | Registered density | 0.500 | 0.167 |
| Partial replay | Registered MultiMatch position | 0.935 | 0.750 |

Chance top-1 credit is 1/6 = 0.167. Ties receive fractional credit, with a
1e-10 relative tolerance so floating-point roundoff cannot create arbitrary
ranks.

## Low-prevalence power stress

The table reports power when item-specific labels are retained with probability
0.25 and otherwise drawn uniformly. Parentheses are normal-approximation 95%
Monte Carlo intervals for 500 simulated studies; these quantify simulation
error, not uncertainty about generalization to real data.

| Scenario | n | GazeWeave | Registered density | Registered MultiMatch, Holm |
|---|---:|---:|---:|---:|
| Registered geometry | 8 | 0.930 (0.908-0.952) | 0.952 (0.933-0.971) | 0.774 (0.737-0.811) |
| Registered geometry | 16 | 0.998 (0.994-1.000) | 0.998 (0.994-1.000) | 0.958 (0.940-0.976) |
| Registered geometry | 32 | 1.000 | 1.000 | 1.000 |
| Order at fixed density | 8 | 0.934 (0.912-0.956) | 0.000 | 0.738 (0.699-0.777) |
| Order at fixed density | 16 | 0.998 (0.994-1.000) | 0.000 | 0.962 (0.945-0.979) |
| Order at fixed density | 32 | 1.000 | 0.000 | 1.000 |
| Partial replay | 8 | 0.874 (0.845-0.903) | 0.000 | 0.510 (0.466-0.554) |
| Partial replay | 16 | 0.992 (0.984-1.000) | 0.000 | 0.778 (0.742-0.814) |
| Partial replay | 32 | 1.000 | 0.000 | 0.982 (0.970-0.994) |

Across all scenarios, sample sizes, and signal levels, estimated type-I error
ranged from 0.028 to 0.066 for GazeWeave, 0 to 0.058 for registered density, and
0.004 to 0.038 for the Holm-corrected registered MultiMatch family. All satisfy
the declared maximum of 0.075.

The stress grid was added only after the initial full-signal run put both
GazeWeave and MultiMatch at power 1.0. No model or comparator parameter changed,
but the stress results are explicitly an addendum rather than a preregistered
analysis.

## Diagnostic recovery

| Scenario | Known-link coupling mass | Omission AUC | Convergence |
|---|---:|---:|---:|
| Registered geometry | 1.0000 | not applicable | 1.000 |
| Order at fixed density | 1.0000 | not applicable | 1.000 |
| Partial replay | 0.9726 | 1.000 | 1.000 |

For registered geometry, recovered scale had bias -0.00010 and RMSE 0.00351.
These results support the fidelity of the warp, coupling, and missing-mass
diagnostics under known synthetic correspondences. They are evidence for
mechanistic inspectability, not evidence that human readers understand the gaze
braid more readily than MultiMatch output.

## Remaining limits

- The generator deliberately isolates conditions GazeWeave is designed to
  model; real gaze may violate these assumptions.
- The bank contains 24 simulated participants, six templates, and short paths.
- Hyperparameters were frozen before the comparative run but have not been
  estimated from independent real positive controls.
- `gaze_bits` remains contrastive compatibility, not a calibrated likelihood
  ratio.
- No repeated-viewing or looking-at-nothing dataset has been tested.
- No reader study establishes comparative human interpretability.

The next external-validity gate is therefore a fully held-out repeated-viewing
positive control followed by a looking-at-nothing dataset, with participant and
item cross-fitting and this same registered-baseline policy.
