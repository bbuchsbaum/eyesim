# Specify edge-normalized GazeWeave Transport

Transport separates matched coverage, normalized correspondence,
selected source and target mass, spatial fidelity, directed local order,
correspondence smoothing, and candidate-invariant registration. Coverage
is integrated by Gauss-Legendre quadrature under a fixed beta prior.

## Usage

``` r
gaze_transport_spec(
  spatial = gaze_gaussian_mixture(c(0.75, 1.5), weights = c(0.7, 0.3)),
  chronology = gaze_order_neighbours(neighbours = 2),
  coverage_prior = c(2, 2),
  coverage_nodes = 12L,
  temporal_weight = 2,
  selection_weights = c(0.5, 0.5),
  entropy_schedule = c(0.15, 0.05, 0.015),
  temperature_bounds = c(0.05, 20),
  warp = gaze_warp_none(),
  screen = NULL,
  maxit = 1000L,
  step_size = 2,
  tolerance = 1e-06,
  projection_maxit = 1000L,
  projection_tolerance = 1e-08,
  projection_method = c("auto", "standard", "log"),
  multistart = 1L,
  backend = c("auto", "optimized", "reference"),
  polish = c("none", "audit"),
  polish_maxit = 200L,
  polish_gap_tolerance = 1e-07,
  polish_relative_tolerance = 1e-08,
  reliability = c("none", "effective_fixations"),
  reliability_kappa_bounds = c(0, 100),
  calibration_folds = 2L,
  calibration_seed = 20260822L,
  log_temperature_prior_sd = 1,
  log1p_kappa_prior_sd = 1,
  revision = "2026.10",
  calibration_control = NULL
)
```

## Arguments

- spatial:

  A \[gaze_gaussian_mixture()\] specification.

- chronology:

  Ordinal chronology from \[gaze_order_neighbours()\].

- coverage_prior:

  Two positive beta-prior shape parameters.

- coverage_nodes:

  Number of default coverage quadrature nodes.

- temporal_weight:

  Non-negative directed local-order weight.

- selection_weights:

  Non-negative reference and source selection weights.

- entropy_schedule:

  Positive correspondence-smoothing continuation values; the last value
  defines the optimized objective. The smooth first stage of the default
  reduces, but does not remove, the dependence of the local optimum
  reached on \`step_size\`.

- temperature_bounds:

  Calibration temperature bounds.

- warp:

  Candidate-invariant cross-fitted warp specification.

- screen:

  Optional screen geometry.

- maxit, step_size:

  Mirror-descent iteration limit and largest step.

- tolerance:

  Stopping tolerance in units of the regularized objective (nats): a
  stage stops when a first-order model predicts at most \`tolerance\` of
  remaining decrease at the current plan (see \`revision\`). On a
  60-pair review set, a tight continuation of nodes stopped at the
  default \`1e-6\` found a further decrease of at most \`1e-6\` at the
  95th percentile; tighter values did not remove the rare slow nodes and
  multiplied projection-limited stalls, and looser ones (\`1e-5\` or
  more) stopped some nodes at worse local solutions.

- projection_maxit, projection_tolerance, projection_method:

  Fixed-mass projection controls.

- multistart:

  One or two common structural starts. Both backends always add the
  adjacent-coverage continuation, and \`2\` adds a spatial start
  (natively supported). On the review set, \`2\` lowered the
  95th-percentile gap to a multistart oracle only from 0.0056 to 0.0053
  nats at 1.48 times the runtime, and on a second review set it never
  moved a score by more than 2e-5, so it is opt-in.

- backend:

  Pair solver backend. \`"reference"\` is the readable R oracle;
  \`"optimized"\` uses the estimator-preserving batched implementation;
  \`"auto"\` uses the optimized backend with a reference fallback. A
  specification the native backend cannot solve (\`projection_method =
  "log"\` or \`"standard"\`) is routed to the reference backend for
  every pair, a node that reaches \`maxit\` is recorded as
  \`"not_converged"\` and scored rather than raising an error, and a
  per-pair numerical fallback raises a warning of class
  \`gaze_transport_backend_fallback\`. The backend used is recorded in
  every alignment's \`convergence\` element.

- polish:

  Optional scientific-objective Frank-Wolfe audit. The default keeps the
  entropic solution; \`"audit"\` polishes every coverage node.

- polish_maxit, polish_gap_tolerance, polish_relative_tolerance:

  Conditional-gradient stopping controls.

- reliability:

  Response-blind calibration policy. \`"effective_fixations"\` learns
  effective-fixation shrinkage (kappa) on inner out-of-fold predictions
  and includes the temperature-only solution as an exact boundary. It
  requires \`gaze_calibration_control(method = "global")\`. The default,
  \`"none"\`, leaves sparse rows to the evidence-scaled temperature (see
  \[gaze_calibration_control()\]).

- reliability_kappa_bounds:

  Non-negative bounds for the shrinkage scale.

- calibration_folds, calibration_seed:

  Inner calibration-fold policy.

- log_temperature_prior_sd, log1p_kappa_prior_sd:

  Penalties of the global calibration (\`gaze_calibration_control(method
  = "global")\`).

- revision:

  Solver revision, recorded in the specification and in every alignment
  and fit. The only supported revision is \`"2026.10"\`, the corrected
  solver and estimand. The chronology residual is \`1 - 2A / (R + S +
  1e-3)\`, continuous as the selected edge mass vanishes. The stopping
  rule is first-order: a stage stops when the plan is feasible, no cell
  outside the support is trapped, and a first-order model predicts at
  most \`tolerance\` of remaining decrease. It does not certify a local
  optimum: slow directions and saddles can remain (on the review set
  about 1 1.4e-3). Trapped cells are reseeded. Mirror steps are capped
  so the exponent never saturates. A line search limited by projection
  noise records the stage as \`"stalled_projection_limited"\` and
  \`maxit\` as \`"not_converged"\`; both are scored. Standard Sinkhorn
  is finished by a dual Newton projection and then log-domain Sinkhorn.
  Every coverage node is solved from the independent and the
  continuation start, replacing the reference backend's silent cold
  restart. The result is a local optimum of a non-convex objective.
  Scores fall short of a multistart oracle (one-sided, mostly in null
  and sparse-source pairs): by up to 0.12 nats with 95th percentile
  0.0056 on one 60-pair review set, and up to 0.117 with 95th percentile
  0.023 (null pairs 0.033, matched 0.006) on a second 88-pair, half-null
  set. \`backend = "auto"\` warns about and records every numerical
  fallback. The polish Frank-Wolfe gap is \`NA\` with a reason on a
  selection boundary. The earlier solver revision \`"2026.08"\` was
  retired in eyesim 0.3.0 because of documented defects (a discontinuous
  chronology residual, a relative-change stopping rule, and silent
  backend fallbacks and cold restarts); requesting it is an error.
  Specifications saved before revisions existed, which were solved as
  \`"2026.08"\`, are an error too: re-create them with
  \`gaze_transport_spec()\`. Install the git tag
  \`gazeweave-2026.08-frozen\` to reproduce 2026.08 results, including
  the frozen validation courts.

- calibration_control:

  A \[gaze_calibration_control()\]. \`NULL\` uses its defaults.

## Value

A frozen \`gaze_transport_spec\`.

## Experimental

GazeWeave is experimental. Its interface, defaults and result columns
may change in future releases; record the package version with any
results.
