# GazeWeave probe-delay persistence verdict

Decision date: 2026-08-20
Protocol: `GAZEWEAVE-PROBE-DELAY-CONTROLS.md`, version 2
Mote: `bd-01M0FFZNNG9PGWQ50KPZEGP4XA`
Verdict: **not supported**

## Evidential boundary

This is a prospectively frozen, phase-resolved reanalysis of the Wynn
probe-delay memory study. It is useful development evidence because it tests
short probe paths, longer post-probe paths, partial replay, registration, and
all declared comparators on real fixation data. It is not an independent
replication, does not close the Wang-data gate, and does not justify a universal
GazeWeave default.

The fixation exports and all provisional court outputs are local-only and
Git-ignored. No participant-level data or derived trial scores are authorized
for publication by this court.

## Court receipt

- 46 participants were present in both exports; 22 had the complete intended
  study and retrieval trial counts and 18 met the two-fold fixation rule.
- The frozen cohort contains eight participants, two from each
  counterbalancing group, and 80 participant-item trials.
- Each outer evaluation cell contains four held-out participants with five
  held-out candidate images per participant. Every training cell excludes all
  evaluated participants and all identities in the evaluated item fold.
- All 9,600 trial-method-condition rows were scored, every Replay and Transport
  fit converged, every candidate set had size five, and all 24 fold audits had
  zero participant and item overlap.
- Empirical null error was 0.0531 for every calibrated method, as expected from
  the pooled 95th-percentile wrong-item and generic-gaze threshold.

## Primary evidence

Information is `log2(p_true / prior_true)` with chance top-one credit 0.20.
Intervals are 1,000-draw crossed participant-by-item bootstrap intervals.

| Task | Method | Information bits | 95% interval | Top-one credit |
|---|---|---:|---:|---:|
| Repeated viewing | Transport v2 | 0.341 | [-0.202, 0.783] | 0.513 |
| Repeated viewing | Replay | 0.142 | [-0.322, 0.476] | 0.475 |
| Repeated viewing | Registered density ridge | 0.326 | [-0.204, 0.797] | 0.513 |
| Repeated viewing | MultiMatch ridge | 0.135 | [-0.205, 0.388] | 0.388 |
| Probe | Transport v2 | -0.006 | [-0.163, 0.128] | 0.225 |
| Probe | Replay | 0.018 | [-0.067, 0.115] | 0.188 |
| Delay | Transport v2 | 0.003 | [-0.167, 0.137] | 0.275 |
| Delay | Replay | -0.068 | [-0.309, 0.089] | 0.225 |
| Delay | Registered density ridge | -0.004 | [-0.160, 0.109] | 0.213 |
| Delay | MultiMatch ridge | -0.012 | [-0.102, 0.069] | 0.225 |
| Late delay | Transport v2 | -0.004 | [-0.127, 0.095] | 0.263 |
| Late delay | Replay | -0.033 | [-0.196, 0.104] | 0.250 |

The repeated-viewing point estimates and ranks are a useful positive-control
signal, although their crossed intervals remain wide with only eight
participants. The retrieval result is much clearer in direction: no method has
a delay interval excluding zero, and neither GazeWeave engine has positive
mean information in late delay. At the frozen null threshold, delay power was
0.0875 for Transport, 0.10 for Replay, 0.05 for registered density, 0.0625 for
MultiMatch, and 0.0125 for elastic matching.

The predeclared persistence rule therefore returns `not_supported`. Probe and
early-delay point estimates cannot rescue that rule. The saliency and old/lure
descriptives are heterogeneous and are not evidence for a post-hoc subgroup.

## What the diagnostics say

- Reordering did not reduce held-out delay information: signal-minus-shuffle
  was -0.0046 bits for Transport, -0.0095 for Replay, and -0.0234 for
  MultiMatch. Density was invariant to numerical precision. This court offers
  no evidence that local chronology improves retrieval identification.
- Retrieval gaze was substantially more spatially constrained than study gaze.
  The cross-fitted source-to-reference correction scale was 1.76 for the probe,
  1.28 for the full delay, 1.17 for early delay, and 1.40 for late delay,
  compared with 0.98 for repeated viewing. These are correction parameters,
  not direct psychological contraction estimates.
- Replay assigned about 70% of retrieval duration to replay states and
  Transport selected about 69--70% coverage even while contrastive item
  information was near zero. Coverage is therefore an explanatory alignment
  diagnostic, not evidence of item-specific reinstatement by itself. The
  contrastive endpoint prevented a plausible-looking alignment from becoming
  a false positive.
- Mean spatial residual was about 86 px for repeated viewing, versus 115--143
  px across retrieval windows. This separation is consistent with the stronger
  repeated-viewing identification result.

## Parsimony and comparator conclusion

Replay completed the full fitted court in 25.9 seconds, compared with 3,776.1
seconds for the current R Transport v2 implementation and 222.8 seconds for all
baselines. Replay is therefore far more computationally parsimonious, but it
did not provide better predictive evidence here. Transport and registered
density were tied on repeated-viewing rank and were both near zero during
retrieval. MultiMatch was available for delay, combined, and repeated viewing;
only 7.5% of probe pairs had the three fixations needed for it, so excluding it
from the probe court was necessary rather than advantageous.

This court supports the existing no-default decision. It demonstrates that the
phase-resolved machinery works end to end and that `gaze_info_bits` remains
appropriately skeptical under central clustering and flexible registration. It
does not demonstrate that GazeWeave has more power than registered density or a
learned MultiMatch composite, nor does it establish post-probe item-specific
reinstatement in this frozen sample.

## Reproduction and local artifacts

From the package root:

```r
devtools::load_all()
source("inst/validation/gaze-weave-probe-delay.R")
run_gaze_weave_probe_delay_court(
  "inst/validation/gaze-weave-probe-delay-results"
)
```

The ignored local results directory contains the configuration, fold audit,
path coverage, method and bootstrap summaries, operating characteristics,
alignment and warp diagnostics, runtime receipts, session information, and an
MD5 manifest. The raw CSV checksums are verified separately against
`test_data/wynn_probe_delay/manifest.csv` before every non-smoke run.
