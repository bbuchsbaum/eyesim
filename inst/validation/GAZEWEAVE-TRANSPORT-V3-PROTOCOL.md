# GazeWeave Transport v3 frozen protocol

Status: frozen before Transport v3 implementation
Protocol version: `3.0.2`
Date frozen: 2026-08-20
Mote gate: `bd-01M0H43RT37FKQ2TF2TF36A8QH`
Machine-readable contract: `gaze-weave-transport-v3-manifest.json`

## Purpose and blinding boundary

This protocol defines Transport v3 before changing the Transport estimator or
opening any additional retrieval result. It is an estimator and evaluation
contract, not a claim that Transport is a probability model or that it will
outperform Replay, density, or MultiMatch.

Model form, hyperparameters, thresholds, and no-go rules may use mathematical
invariants, public synthetic data, and the response-blind repeated-viewing
positive control. Retrieval response, correctness, confidence, saliency,
old/lure status, and all retrieval subgroups remain sealed until the v3 code,
configuration, folds, candidates, and score files have been hashed at GW-18.10.
No result observed after a gate is opened may relax that gate's threshold.

Private Wynn data, folds derived from it, checkpoints, and trial-level output
remain local and Git-ignored. The protocol, fixture bank, synthetic court, and
all examples committed with the package contain only generated public data.

## Scientific roles and primary endpoint

Transport is a symmetric, registered, coverage-conditioned alignment. Its
coupling is an optimized correspondence. Coupling diffuseness is alignment
ambiguity induced partly by correspondence regularization; it is not
psychological uncertainty and must never be labelled posterior replay
probability. Replay remains the directional encoding-to-recall model. Density
and MultiMatch remain distinct spatial and scanpath comparators.

For a held-out row with true candidate `k*`, declared candidate prior
`pi[k]`, and normalized candidate probability `p[k]`, the sole default
inferential endpoint is

```text
gaze_info_bits = log2(p[k*] / pi[k*]).
```

`odds_bits`, log loss, Brier score, calibration, rank, coverage, selection,
spatial fidelity, chronology, contraction, correspondence, and solver
diagnostics are secondary. They do not become additional primary outcomes.

## Inputs, coordinates, and duration mass

A path contains finite `x`, `y`, `onset`, and non-negative `duration`.
Zero-duration rows are discarded. Remaining duration is normalized within an
episode to unit mass. A path with no positive-duration fixation is invalid.
One-fixation paths are valid spatial objects with an empty chronology; their
chronology residual is zero by convention and their lack of order information
is reported.

The scientific courts use degrees of visual angle on a declared `24 x 18`
degree screen. Pixel inputs must be converted once using recorded screen
geometry before scoring. Spatial bandwidths, residuals, and warp translations
use the same declared coordinate unit. Equivalent pixel/degree conversions
must agree within the frozen unit-invariance tolerance.

Exact adjacent colocated fixations are coalesced before duration
normalization. A non-zero location tolerance is sensitivity-only and must be
declared in degrees. A return to the same location after an intervening
fixation is not coalesced. Proportional splitting or merging of a fixation is
therefore a representation change, not a scientific perturbation.

## Order clock and local edges

The default memory clock is ordinal. Each fixation has directed edges to the
next `k = 2` events when they exist. Physical-time edges are an explicit
sensitivity analysis and are never mixed silently with ordinal edges. Uniform
temporal dilation leaves duration mass, ordinal edges, and the default result
unchanged.

Let `R_X` and `R_Y` be the directed binary local-edge matrices for reference
and source paths. With normalized correspondence `Gamma`, define selected
unit marginals `u = Gamma 1` and `v = Gamma' 1`, and

```text
A = sum[ii'] R_X[ii'] u[i] u[i']
B = sum[jj'] R_Y[jj'] v[j] v[j']
C = sum[ii'jj'] R_X[ii'] R_Y[jj'] Gamma[ij] Gamma[i'j']
r_T = 1 - 2 C / (A + B), when A + B > 1e-15.
```

When both supported edge masses are zero, `r_T = 0`; when exactly one is
zero, `r_T = 1`. Numerical results are clipped to `[0, 1]` only within 64
machine epsilons of the interval. `C` is evaluated with matrix products, never
by allocating an explicit four-index tensor. This directed-Dice residual is
conditioned on supported edges, so complete reversal cannot become cheaper
merely because a path is longer.

## Coverage, selection, and correspondence

Let `a` and `b` be unit duration masses. The transported plan is separated as

```text
Pi = M Gamma,  0 <= M <= 1,  sum(Gamma) = 1,
M Gamma 1 <= a,  M Gamma' 1 <= b.
```

`M` is matched coverage. `Gamma` is normalized correspondence. `M u` and
`M v` are selected reference and source mass. These are separate reported
objects. At `M = 0`, candidate-specific fit is exactly neutral, no division by
`M` occurs, and conditional residuals are `NA` with an explicit zero-coverage
status.

For positive coverage, define:

```text
r_S(Gamma, W) = sum[ij] Gamma[ij] c(x[i], W(y[j]))
r_A(Gamma)    = JS(u, a)
r_B(Gamma)    = JS(v, b)
r_I(Gamma)    = KL(Gamma || u outer v)
r_W(W)        = declared distance from the outer-training warp centre
```

`JS` is Jensen-Shannon divergence in nats with zero terms handled by
continuity. It makes selective use of one convenient fixation explicit while
remaining invariant to proportional refinement. Dominated marginal
constraints remain mandatory. `r_I` is correspondence mutual information and
is used only as numerical smoothing.

The fixed-coverage scientific energy is

```text
E_0(M, Gamma, W) =
  M * {r_S + lambda_T r_T + lambda_A r_A + lambda_B r_B}
  + lambda_W r_W.
```

The optimization energy is

```text
E_eps = E_0 + epsilon * M * r_I.
```

Candidate evidence uses `E_0` evaluated at the regularized solution; `r_I`
is excluded. Coverage, normalized correspondence, selected source mass,
selected target mass, conditional spatial residual, conditional chronology
residual, correspondence smoothing, and warp penalty are always stored as
separate fields. No component is called uncertainty about memory.

The prospective default weights in degree coordinates are `lambda_T = 2`,
`lambda_A = 0.5`, `lambda_B = 0.5`, and `lambda_W = 1`. Spatial cost is the
frozen Gaussian-mixture cost at scales `0.75` and `1.5` degrees with weights
`0.7` and `0.3`. Any replacement value must be selected by inner out-of-fold
log loss before the corresponding outer fold is scored and recorded as a
protocol deviation; retrieval outcomes cannot select it.

## Coverage integration

Coverage is integrated under a prospectively fixed `Beta(2, 2)` prior using
twelve-point Gauss-Legendre quadrature on `[0, 1]`. The nodes are approximately

```text
0.0092197, 0.0479414, 0.1150487, 0.2063410, 0.3160843, 0.4373833,
0.5626167, 0.6839157, 0.7936590, 0.8849513, 0.9520586, 0.9907803.
```

Thus coherent replay below 50 percent is represented, including a node near
0.21. Prior density is incorporated into normalized quadrature weights, not
added a second time to energy. A twenty-four-point rule is the required
refinement check. `M = 0` is evaluated as a boundary diagnostic even though
it has zero mass under the default continuous prior. Posterior coverage is an
alignment-profile weight conditional on the chosen energy model; it is not a
memory posterior.

The candidate compatibility is the log integral of `exp(-E_0)` over coverage
using the normalized prior/quadrature weights. All candidates use identical
nodes and solver policy.

## Multiple study episodes

Each study presentation is prepared and aligned independently. Episodes are
never concatenated, and no chronology edge crosses a presentation boundary.
After applying the same calibrated energy scale, the candidate score is the
log of an equal-prior mixture of its valid episode likelihoods.

For missing episodes, the evaluation row uses the intersection of valid
presentation indices across every candidate in its common candidate pool.
Weights are renormalized equally over that shared set and the omitted indices
are recorded. If the intersection is empty, the row is missing and contributes
to the declared failure rate, not a neutral score. Reordering episode rows
cannot change evidence. The one-episode case is exactly the single-template
contract. Outcome-driven episode weighting is prohibited.

## Registration policy

The default warp is isotropic contraction about the known screen centre plus
a participant/session translation. Warp distributions and point estimates
are fitted only in the outer-training data, without the held-out
participant-item cell, then applied identically to every candidate in a row.
Pair-specific registration is prohibited. The identity warp is always
available. Affine registration is sensitivity-only.

Every comparator receives the same candidate-invariant registration
opportunity: its declared raw version and, where meaningful, the same
outer-trained registered coordinates. A comparator may not receive true-item
registration unavailable to its wrong candidates.

## Candidate pools, priors, folds, and calibration

The default pool is every permitted within-participant candidate on common
support. If exhaustive scoring is infeasible, averaged deterministic panels
are frozen before outcomes are opened, use every true item equally often, and
must reproduce exhaustive information within `0.01` bits on the public court.
Candidate priors are the actual design priors; they are uniform only when the
design is uniform. The same pools and priors are used by every comparator.

Outer folds hold out participant-item cells and are shared by all methods.
All preprocessing, warp fitting, hyperparameter selection, coverage-policy
selection, reliability fitting, and temperature fitting occur without the
held-out cell. Inner folds are grouped by item within the outer-training set,
with at least two candidate-bearing folds. Degenerate strata use the frozen
boundary solution and record the fallback.

For candidate log compatibilities `s[k]`, temperature `tau`, reliability
weight `rho`, and prior `pi[k]`,

```text
q[k] = softmax(log(pi[k]) + s[k] / tau)
p[k] = rho * q[k] + (1 - rho) * pi[k].
```

`tau` has a unit-scale Gaussian penalty on `log(tau)` centred at zero and is
bounded to `[0.05, 20]`. Reliability uses response-blind path support only and
includes `rho = 1` as the temperature-only boundary. Its inner out-of-fold log
loss must be no worse than that boundary within `1e-10`; otherwise the boundary
is used. Fixation count, effective count, duration concentration, entropy, and
generic gaze quality are diagnostics and cannot become item-specific evidence
without prospective inner-fold selection by proper score.

## Missing data and common support

Rows missing the true candidate, the declared prior, all valid episodes, or a
finite score from any mandatory comparator are excluded from the primary
common-support comparison and counted by method and reason. Optional
MultiMatch dependency failures are explicit and cannot silently shrink only
one method's support. Secondary method-specific support summaries are labelled
as such. Failed alignments never receive a favourable sentinel score.

## Baseline panel

Every scientific court contains, on identical outer folds, candidate pools,
priors, registration opportunities, and uncertainty resamples:

1. Transport v2;
2. stabilized directional Replay;
3. registered multiscale density;
4. fixed density with `sigma = 80` pixels;
5. fixed density with `sigma = 160` pixels;
6. raw MultiMatch position, direction, length, duration, and shape dimensions;
7. a ridge-regularized learned MultiMatch composite fitted in nested folds;
8. Transport v3 reference and optimized backends.

## Seeds and stopping rules

Seeds are `20260820` for the v2 golden fixture, `20260821` for the synthetic
court, `20260822` for calibration and null-prior recovery, `20260823` for the
repeated-viewing bootstrap, `20260824` for performance ordering, and
`20260825` for the sealed retrieval bootstrap. R's `L'Ecuyer-CMRG` generator
is used for parallel resampling, and candidate order is independently permuted
from the fixed performance seed.

The entropic solver stops only when dominated-marginal feasibility is at most
`1e-8` and relative scientific-energy change is at most `5e-5`, or after 1000
iterations. The optional scientific-polish solver stops at conditional-gradient
gap `<= 1e-7`, relative energy change `<= 1e-8`, or 200 iterations. A maximum
iteration exit is reported as non-converged even when its score is finite.
At least 99 percent of mandatory court alignments must be finite and converged.

## Frozen quantitative gates

The authoritative machine values are in the JSON manifest. In summary:

- identical local order has residual `<= 1e-10`; complete reversal has
  residual `>= 0.95`; reversal residuals for lengths 4, 8, 16, 32, and 64 span
  no more than `0.05`;
- ordinal-clock uniform dilation changes results by at most `1e-10`;
- proportional split/merge changes energy, coverage, and selected mass by at
  most `1e-8` after coalescing;
- coverage quadrature refinement changes pair energy by at most `1e-3` and
  calibrated candidate information by at most `0.01` bits;
- optimized and reference energies agree within `1e-8`, coverage within
  `1e-10`, barycentric coordinates within `1e-6` degrees, warp parameters
  within `1e-10`, and candidate rank exactly;
- median 15-by-15 pair runtime is at most 15 ms and at least four times faster
  than the frozen approximately 60 ms GW-14 reference; a representative full
  fold is at least twice as fast as GW-14 with bounded peak memory;
- null mean information is within `0.02` bits of zero, null top-one credit is
  within `0.02` of its prior expectation, held-out ECE is at most `0.10`, and
  mandatory finite/converged rate is at least `0.99`;
- repeated viewing advances only when intact mean held-out information is
  positive and the paired intact-minus-mean(reversal, shuffle) contrast is
  positive. At least one of those two crossed 95-percent bootstrap intervals
  must have a lower endpoint above zero.

Synthetic advancement additionally requires chronology destruction at fixed
density to reduce Transport v3 evidence, coherent partial replay to lie
between intact and destroyed controls, matched false-positive rate at most
`0.075`, and no mandatory representation, fold-isolation, calibration,
convergence, or performance failure.

## Gate decisions and no-go rules

The human-readable decision table is
`GAZEWEAVE-TRANSPORT-V3-DECISIONS.md`. A failed invariant, reference-parity,
fold-isolation, calibration, or convergence gate returns the algorithm to the
responsible pre-study gate. Thresholds remain frozen. Failure of the runtime
gate keeps the pure-R estimator experimental and prohibits promoting it as the
default exhaustive engine. Failure of the repeated-viewing advance rule seals
retrieval outcomes and leads to a bounded no-go/experimental-role report.

Transport v3 is promoted only if it passes invariants, calibration, runtime,
synthetic discrimination, and the study positive control. If it is stable and
scientifically distinctive but tied with fair baselines on real item
information, it remains the symmetric explanatory alignment engine. If it
fails stability or runtime, the readable reference may remain as an
experimental audit implementation. The later Wang result may revise the role
but cannot tune this frozen estimator.

## Response-blind implementation addendum A: coverage resolution

The first public pairwise refinement fixture, run before study or retrieval
scoring, found that the originally declared 8-node rule differed from its
16-node check by `0.00117`, exceeding the unchanged `0.001` energy threshold.
The default and refinement rules were therefore increased to 12 and 24 nodes.
On the same fixture their maximum pair-energy difference was `0.0000257`.
No scientific threshold, response-derived choice, coverage prior, objective
weight, or outcome was changed. Protocol version `3.0.1` records this
response-blind correction.

## Response-blind implementation addendum B: continuation ceiling

The first scope-matched 315-candidate performance fold found one candidate
whose high-entropy continuation stage reached 500 iterations while its final
stage converged. The unchanged feasibility and relative-energy tolerances
passed when the maximum was increased to 1000 iterations. Under that ceiling,
315 of 315 deterministic 15-fixation candidate alignments converged in 13.955
seconds, versus the recorded 830.849-second optimized GW-14 fold (59.54x),
with a 120,869,488-byte serialized result. Protocol version `3.0.2` records
this response-blind stopping correction; no objective, threshold, study path,
or retrieval outcome changed.
