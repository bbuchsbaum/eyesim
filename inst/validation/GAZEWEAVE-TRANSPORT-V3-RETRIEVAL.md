# Transport v3 frozen 0--3000 ms retrieval court

Protocol: `gazeweave-transport-v3/3.0.2`
Outer/bootstrap seed: `20260825`
Measurement frozen: 2026-08-21 05:37:20 UTC
Behavior opened only after measurement hashes verified: **yes**

## Decision

Transport v3 does not earn a retrieval-prediction claim in this cohort. It has
the largest point estimate of item information in the exact 0--3000 ms window,
but its crossed interval includes zero, its advantages over the principal fair
comparators include zero, and it does not improve held-out old-item response
log loss.

This result is retained without changing the algorithm, time window,
calibration, saliency form, subgroup, or nonlinear specification. Transport v3
continues to have positive study sensitivity and distinctive chronology
behavior; the product-role decision therefore belongs to GW-18.11 rather than
to a post-hoc retrieval rescue.

## Pre-outcome freeze

Before behavioral values were merged, the court wrote and verified:

- 1,295 score-blind retrieval rows from the exact `[0,3000)` ms window;
- 5,180 study rows, exactly four presentations per participant/item;
- K=5 candidate sets with a uniform prior and four equal-prior episode
  likelihoods per candidate;
- the existing four crossed participant/item outer folds and exhaustive
  candidate plan from the full-cohort comparator court;
- a 12-file source/protocol manifest and six private design artifacts;
- the v3 specification, thresholds, session information, and freeze receipt.

Public manifest hashes are:

| Artifact | MD5 |
|---|---|
| `source-manifest.csv` | `06ad00cecc77acf079d02414819bd52c` |
| `design-manifest.csv` | `e2d6509e001e9edd35f74aff797fbe0a` |
| `thresholds.csv` | `e6febb599fda20b1648ed32539d6c191` |
| `freeze-receipt.csv` | `023bae869e39e9ea727a68aaea7073ff` |

The private design contents and trial-linked score checkpoints remain under
the Git-ignored retrieval-results directory. Source and checkpoint manifests
still verified after the behavioral court (`SOURCE_HASH_OK=TRUE`,
`CHECKPOINT_HASH_OK=TRUE`).

## Frozen measurement

Every method uses the same 1,295 trials, outer folds, K=5 candidate panels,
uniform prior, registration opportunity, and crossed bootstrap draws.

| Method | Mean item information (bits) | Crossed 95% interval |
|---|---:|---:|
| Transport v3, four episodes | 0.02719 | -0.00305 to 0.05618 |
| Transport v2 | 0.00728 | -0.01544 to 0.02900 |
| MultiMatch position | 0.00604 | -0.00876 to 0.01977 |
| Density sigma 80 | 0.00594 | -0.01284 to 0.02366 |
| MultiMatch position EMD | 0.00347 | -0.01216 to 0.01840 |
| Replay | 0.00316 | -0.02083 to 0.02585 |
| Density sigma 160 | 0.00314 | -0.01031 to 0.01640 |
| Registered density ridge | 0.00242 | -0.01763 to 0.02131 |
| MultiMatch vector | 0.00039 | -0.00052 to 0.00134 |
| MultiMatch direction | 0.00013 | -0.00321 to 0.00357 |
| MultiMatch length | -0.00010 | -0.00236 to 0.00186 |
| Learned MultiMatch ridge | -0.00083 | -0.01672 to 0.01382 |
| MultiMatch duration | -0.00208 | -0.00750 to 0.00269 |

The paired v3-minus-comparator point estimates are `0.01990` bits versus v2,
`0.02402` versus Replay, `0.02125` versus density sigma 80, `0.02476` versus
registered density, and `0.02802` versus learned MultiMatch. All principal
intervals include zero; only the contrast with learned MultiMatch excludes
zero (`0.00406` to `0.05287`). This isolated contrast does not override the
prespecified common-panel conclusion.

Transport v3's mean log loss is `1.59059`, mean rank `2.781`, and top-one
credit `0.2479`. The four folds have zero participant and item overlap, all
rows use all four episodes, and candidate-alignment convergence is `0.99954`.
All four calibrations are genuinely inner out-of-fold, use 315--332 training
rows, and have zero match-key overlap. Fold wall times are 524.1, 599.7, 521.8,
and 528.7 seconds (2,174.3 seconds total).

## Canonical old-item behavioral result

The behavioral endpoint is the out-of-fold log-loss gain from adding
item-specific `gaze_info_bits` to a base response model. The base model already
contains saliency, log effective fixations, and total path duration; generic
gaze quality and response-blind reliability therefore remain separate from
item-specific evidence. Training and evaluation participant/item cells do not
overlap.

For 630 old items with a non-missing 1--4 response, Transport v3 changes
held-out response information by `-0.01508` bits (crossed 95% interval
`-0.04322` to `0.00624`). Its full-model log loss is `0.31203`, compared with
`0.30157` for the quality-only base model. Thus v3 does not improve canonical
old-item response prediction.

No method provides a stable positive old-item predictive gain in this panel.
Several weak dimensions produce large negative held-out gains, illustrating
why in-sample coefficient significance is not the decision endpoint.

The explanatory crossed participant/item GLMM gives a standardized v3 slope
of `0.1639` (standard error `0.1525`, `p=0.283`). The fit is singular and is
reported only as a diagnostic, as prospectively specified.

## Declared secondary outcomes

- On 641 usable lure trials, Transport v3 changes false-alarm prediction by
  `+0.00279` bits (95% `-0.00710` to `0.01311`): no reliable gain.
- Across 1,271 non-missing old and lure responses, the descriptive Spearman
  association between the 1--4 response code and v3 information is `-0.0601`.
  This graded-confidence analysis is secondary and does not revise the
  primary result.

## Claims

Permitted: Transport v3 is calibrated, computationally viable, study-sensitive,
chronology-sensitive, and has a small positive retrieval item-information
point estimate under four-episode likelihood mixing.

Prohibited: reliable retrieval reinstatement, superiority over v2/Replay/
density, improved old-item response prediction, a confidence effect, or a
saliency-specific rescue.
