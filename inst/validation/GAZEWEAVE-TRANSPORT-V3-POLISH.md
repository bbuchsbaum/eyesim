# GazeWeave Transport v3 scientific-polish audit

Protocol version: 3.0.2
Gate: GW-18.6
Status: audit backend accepted; default remains entropic-only

## Frozen method

The optional `polish = "audit"` path starts from each fixed-coverage entropic
solution and applies conditional-gradient updates to the frozen scientific
objective. It does not optimize the correspondence-smoothing term. Each linear
minimization problem is an exact dominated-marginal transport program solved by
`lpSolve`; deterministic grid-bracketed scalar minimization chooses the step.
The update is a convex combination of feasible couplings, so coverage and both
dominated marginals remain feasible.

The audit reports the scientific objective trace, conditional-gradient gap,
step, relative improvement, exact coverage residual, row/column dominance
residuals, non-negativity residual, stopping reason, and the entropic starting
objective. Warp fitting remains outside pair scoring, so polishing cannot alter
the frozen candidate-invariant warp.

`lpSolve` is a suggested package and is required only when the audit is
requested. POT was checked as the proposed optional differential oracle but was
not available in the project Python environment. It is not a runtime or test
dependency.

## Fixed fixture court

The deterministic four-fixation fixture in
`tests/testthat/test_gaze_weave_transport_v3_polish.R` produced:

| Quantity | Result |
|---|---:|
| Entropic optimized log score | -0.004006057 |
| Polished optimized log score | -0.003862281 |
| Entropic optimized/reference difference | 0 |
| Polished optimized/reference difference | 1.83e-13 |
| Total scientific-energy improvement across two coverage nodes | 2.87e-4 |
| Maximum final conditional-gradient gap after 30 updates | 3.02e-6 |
| Maximum feasibility residual | 2.22e-16 |
| Independent-versus-spatial-start energy spread | 8.54e-14 |
| Independent-versus-spatial-start barycentric RMSE | 3.14e-16 |

The court also permutes candidate order and batch size, compares the optimized
and pure-R reference backends, freezes coverage, verifies unchanged warp
parameters, and checks candidate-rank agreement across independent and
spatial-vertex starts. Elementwise coupling equality is deliberately not an
acceptance condition.

## Default decision

The scientific objective improved and the exact feasibility and invariance
checks passed. However, on the same small fixture the mean pair time over 20
runs rose from 0.00165 seconds for the entropic optimized backend to 0.21950
seconds with polishing, about 133 times slower. The 30-update fixture also did
not reach the strict 1e-7 gap threshold, although it remained monotone and
feasible.

Therefore Transport v3 keeps `polish = "none"` as its frozen default.
Conditional-gradient polishing is retained as an explicit diagnostic/audit
mode, not as the production scorer. This decision uses only invariant synthetic
fixtures and no sealed retrieval behavior or private Wynn outcome.
