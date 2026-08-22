# GazeWeave benchmark baseline

Run `gaze-weave.R` after loading or installing `eyesim`. The benchmark uses a
deterministic smooth path, one optimizer start, and 200 maximum iterations.

Initial local baseline on 2026-08-15 (macOS arm64, R 4.5.1):

| Reference fixations | Source fixations | Median seconds |
|---:|---:|---:|
| 20 | 20 | 0.089 |
| 50 | 50 | 0.273 |

These measurements document the initial performance envelope; they are not
portable hard thresholds. Future performance checks should compare relative
changes on a stable runner and keep numerical-correctness tests separate.

## Directional Replay

Run `gaze-weave-replay.R` after loading or installing `eyesim`. Replay uses a
fixed 64-point duration grid and a low-rank restart decomposition. Its
forward-backward recursion is `O(grid_size * encoding_fixations * max_skip)`;
the dense transition matrix is constructed only for inspection and tiny-oracle
tests.

Record scoring-only timings separately for 20- and 100-fixation templates. The
benchmark fits the training model once before timing repeated held-out pair
scores. Treat these as runner-specific baselines, not package pass/fail limits.

Initial local baseline on 2026-08-15 (macOS arm64, R 4.5.1, five repetitions):

| Encoding fixations | Recall fixations | Grid points | Median seconds |
|---:|---:|---:|---:|
| 20 | 20 | 64 | 0.002 |
| 100 | 100 | 64 | 0.006 |

## Coverage-conditioned Transport v2

Run `gaze-weave-transport-v2.R` after loading or installing `eyesim`. This
benchmark exercises two fixed coverage values and two entropy-continuation
values with one deterministic start. The inspectable R solver evaluates the
relational gradient with matrix products; for similarly sized paths its
leading per-iteration cost is cubic in fixation count. Runtime also scales
with the coverage profile, continuation schedule, starts, and projection
iterations.

Record convergence together with runtime. These timings are a performance
envelope, not a correctness threshold; algebraic, feasibility, and solver
tests remain the correctness court.

Initial local baseline on 2026-08-15 (macOS arm64, R 4.5.1, three
repetitions):

| Reference fixations | Source fixations | Coverage values | Continuation values | Median seconds | Convergence rate |
|---:|---:|---:|---:|---:|---:|
| 10 | 10 | 2 | 2 | 3.213 | 1.00 |
| 20 | 20 | 2 | 2 | 5.778 | 1.00 |

## Edge-normalized Transport v3

Run `gaze-weave-transport-v3.R` after loading `eyesim`. The benchmark compares
the readable R oracle with the RcppArmadillo backend on the frozen 15-by-15
pair and a deterministic exhaustive candidate fold. Both backends use the same
ascending within-candidate coverage continuation. No state crosses candidate
boundaries, and native plans are re-evaluated by the R scientific objective.

The native backend requires Rcpp and RcppArmadillo at build time. Runtime use
is optional through `backend = "reference"`; `backend = "auto"` falls back to
the R implementation when native policy or numerical support is unavailable.
The package registers native routines and uses no GPU or platform-specific
instruction set. CRAN and sanitizer checks remain release gates.

Frozen local baseline on 2026-08-20 (macOS arm64, R 4.5.1):

| Scope | Evaluations | R reference median | Optimized median | Speedup |
|---|---:|---:|---:|---:|
| 15 x 15 pair | 1 | 47 ms | 8 ms | 5.88x |
| Representative exhaustive fold | 4 | 5.798 s | 41 ms | 141.41x |
| GW-14-scale deterministic fold | 315 | 830.849 s | 13.955 s | 59.54x |

Pair values use five repetitions after one native warm-up; fold values use
three repetitions. The machine-readable record is
`gaze-weave-transport-v3-baseline.csv`. These runner-specific measurements
pass the frozen 15 ms, 4x pair and 2x fold gates. Correctness remains governed
by the separate native-versus-R differential court.

The 315-candidate row uses the authoritative optimized GW-14 elapsed time from
its Mote completion receipt and a scope-matched Transport-v3 run with 315
deterministic 15-fixation candidates. All 315 v3 alignments converged; the
serialized result occupied 120,869,488 bytes.
