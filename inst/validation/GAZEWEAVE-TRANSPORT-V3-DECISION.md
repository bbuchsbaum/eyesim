# Transport v3 benchmark and product-role decision

Protocol: `gazeweave-transport-v3/3.0.2`
Epic: GW-18
Decision gate: GW-18.11
Decision date: 2026-08-21
Retrieval behavior opened: yes, only after the frozen score bundle and hashes
were verified

## Decision

Retain Transport v3 as GazeWeave's symmetric explanatory alignment engine.
It is suitable for planned questions about matched coverage, optimized spatial
correspondence, bounded local order, and coherent partial replay when neither
path has a privileged generative direction.

Do not present Transport v3 as a default retrieval predictor or as generally
superior to Replay, registered density, Transport v2, or MultiMatch. It passed
the invariant, calibration, numerical, runtime, simulation, and full-cohort
study-sensitivity gates, but it tied the fair comparators on repeated-viewing
item information and did not pass the frozen retrieval or behavioral-prediction
uncertainty gates.

The independent Wang result remains confirmatory. It may revise this role, but
it cannot tune the estimator or alter the completed Transport-v2 or v3 courts.

## Gate ledger

| Gate | Evidence | Result |
|---|---|---|
| Estimator lock | Public Transport-v2 golden fixture and MD5 manifest | Pass; scientific outputs agree within `1e-8` |
| Bounded chronology | Analytic and generated paths of length 1--64 | Pass; reversal sensitivity span `<= 0.05` |
| Representation | Split/merge, units, duplicate locations, order, batch size | Pass at frozen tolerances |
| Coverage | 12- versus 24-node quadrature | Pass; maximum pair-energy change `9.65e-4`, information change `4.60e-4` bits |
| Native parity | Random, sparse, extreme, short, and public path fixtures | Pass; score tolerance `1e-8` |
| Scientific polishing | Conditional-gradient audit of the unregularized objective | Correct and feasible; retained as audit-only because it was about 133x slower on the fixed fixture |
| Nested calibration | Inner item-held-out temperature and response-blind reliability | Pass; no outer or inner match-key overlap, null/ECE gates pass |
| Synthetic court | Common exhaustive candidates and comparator panel | Pass; mean signal information `2.675` bits, matched false-positive rate `0.0625`, convergence `1.00` |
| Pair runtime | Frozen 15-by-15 fixture | Pass; 8 ms optimized versus 47 ms reference, `5.88x` |
| Scope-matched runtime | 315 deterministic candidate alignments | Pass; 13.955 s versus 830.849 s GW-14, `59.54x`; all converged |
| Study positive control | 1,295 targets, K=5, P1--P4 and controls | Pass; P1--P4 `0.28464` bits (`0.18206`, `0.38573`); intact-control `0.05237` (`0.00720`, `0.09536`) |
| Retrieval | Exact frozen 0--3000 ms, 1,295 targets, four equal-prior episodes | No claim; `0.02719` bits (`-0.00305`, `0.05618`) and principal comparator contrasts include zero |
| Old-item response prediction | Cross-fitted quality base plus item-specific v3 information | No claim; gain `-0.01508` bits (`-0.04322`, `0.00624`) |

Runtime measurements are local macOS arm64/R 4.5.1 evidence, not a
cross-platform performance guarantee. The native implementation has no GPU or
platform-specific instruction-set requirement; platform release evidence is
owned by package CI and R CMD check.

## Public contract

- `gaze_transport_v3_cv()` scores every permitted candidate in each declared
  contrast, with actual design priors and participant/item-disjoint outer and
  inner fitting.
- `broom::tidy()` exposes identifiers and the sole primary endpoint,
  `gaze_info_bits`, by default.
- `gaze_transport_v3_result()` nests matched coverage, spatial RMSE with its
  declared coordinate unit, local-order preservation, contraction scale,
  template rank, episode count and weights, and solver stability.
- One presentation is the exact identity case of the fixed equal-prior
  likelihood mixture. Several presentations remain separate and receive equal
  prior weight over the common valid presentation-index intersection. They are
  never concatenated.
- `autoplot()` provides the registered overlay, optimized-correspondence gaze
  braid, raw-versus-registered paths, bounded diagnostics, and a one-score
  evidence ledger. Selecting an episode changes only an explanatory alignment
  panel, not the all-episode evidence score.
- Every Transport label says optimized correspondence or alignment. No
  Transport coupling is called a posterior replay probability.

## Invariance and failure contract

Transport v3 is invariant, within frozen tolerances, to candidate order, batch
size, global time dilation, declared coordinate-unit conversion, and fixation
splitting that preserves mass and local-edge representation. It is deliberately
sensitive to reversal, local swaps, coherent partial coverage, and candidate
pool/prior changes.

Do not interpret a result as item evidence when participant/item folds leak,
candidate pools differ across compared conditions, coordinate units are
undeclared, a non-identity warp is fitted on the scored pair, inner data are too
thin for the declared calibration, common episode support differs without the
recorded intersection rule, or a required alignment fails to converge. These
are design or diagnostic failures rather than evidence of absent replay.

## API-friction decisions

| Workflow | Observed friction | Decision |
|---|---|---|
| Read one fit row | Raw result tables expose many calibrated and diagnostic columns beside the primary endpoint | Default `tidy()` is restricted to identifiers plus `gaze_info_bits`; the audit surface is explicit and nested |
| Inspect several presentations | A plotted pair could be mistaken for the evidence-producing mixture | Plot captions name the selected episode and state that the score retains all equal-prior episodes |
| Interpret spatial RMSE | Native coordinates previously left the unit implicit | The public audit record declares screen, spatial, or `native_coordinate_units` |
| Use the common dispatcher | Transport v3 needs episode and design-prior arguments absent from the legacy dispatcher | Retain the dedicated `gaze_transport_v3_cv()` entry point; do not hide estimand-defining arguments in metadata |
| Read the coupling | Transport and Replay braids can look superficially similar | Engine-specific titles state optimized alignment versus posterior replay semantics |

## Claims

Permitted:

- calibrated symmetric gaze alignment with bounded local-order sensitivity;
- study/repeated-viewing sensitivity under the frozen full-cohort design;
- refinement-stable coverage integration and native/reference score parity;
- the recorded runner-specific performance improvements; and
- a small positive retrieval item-information point estimate with uncertainty
  that includes zero.

Prohibited:

- reliable retrieval reinstatement or old-item response prediction;
- superiority over Replay, density, Transport v2, or MultiMatch;
- interpreting optimized correspondence mass as posterior uncertainty;
- treating alignment diagnostics as multiple inferential endpoints;
- saliency, confidence, subgroup, or nonlinear rescue claims; and
- changing the algorithm from the completed retrieval result.

## Privacy, provenance, and handoff

Public generated fixtures, aggregate reports, benchmark records, protocol
files, and the retrieval freeze receipt live under `inst/validation` and
`inst/benchmarks`. Restricted Wynn CSVs, participant-linked repeated-viewing
checkpoints, retrieval score checkpoints, and behavior merges are Git-ignored
and excluded by `.Rbuildignore`; they must not enter a source or binary package.

The pre-outcome retrieval bundle hashes are:

| Artifact | MD5 |
|---|---|
| `source-manifest.csv` | `06ad00cecc77acf079d02414819bd52c` |
| `design-manifest.csv` | `e2d6509e001e9edd35f74aff797fbe0a` |
| `thresholds.csv` | `e6febb599fda20b1648ed32539d6c191` |
| `freeze-receipt.csv` | `023bae869e39e9ea727a68aaea7073ff` |

GW-10 receives this role decision and frozen provenance without reopening its
existing Transport-v2 court. Any later confirmatory data must be evaluated
against the frozen v3 specification and manifest, never used to tune them.
