# GazeWeave maximal-cohort recognition verdict

Verdict date: 2026-08-20
Protocol version: 1
Mote: `bd-01M0FRA5QSTQZJ4GMNVH6MBN6R`
Status: `not_supported`

## Bottom line

The maximal-cohort court does **not** establish that gaze reinstatement varies
with probe saliency, response correctness, or their interaction during the
combined 500-ms probe plus 2500-ms delay interval. This is true for GazeWeave
Transport v2 and Replay, and no density, MultiMatch, or elastic comparator
survived its descriptive multiplicity adjustment.

Transport v2 had the largest mean calibrated item information in this court,
but the advantage was tiny and was not a prespecified superiority test. The
result therefore supports neither a power advantage nor a universal default.
It does show that all methods can be compared on the same held-out candidate
task, and that the Transport implementation can now complete that court much
more efficiently without changing its scores.

## Design receipt

- 46 participants were present in the local exports.
- 1,348 participant-item trials from 45 participants met the combined-window
  fixation criteria.
- The score-blind fixed-K support rule retained 1,295 trials from 36
  participants, or 96.1% of eligible trials.
- Each held-out trial used the same five candidates for every method.
- Four crossed participant-by-item folds evaluated 332, 316, 315, and 332
  trials with zero participant or item overlap.
- Both GazeWeave engines produced finite scores and converged for all 1,295
  evaluated candidate sets.
- All 2,000 shared crossed-bootstrap draws were valid. The fixed-effect design
  condition number was 4.23.
- Mixed models converged, although most random-intercept fits were singular;
  the crossed bootstrap remains the primary uncertainty procedure.

Participant-linked inputs, trial scores, candidate scores, checkpoints, and
fitted objects remain local and Git-ignored.

## Prespecified GazeWeave effects

Values are information bits. Family intervals are the Bonferroni 99.1667%
intervals for the six engine-by-contrast tests.

| Engine | Contrast | Estimate | 95% interval | Family interval | Bootstrap p |
|---|---|---:|---:|---:|---:|
| Transport v2 | Saliency 20 to 100 | +0.0227 | [-0.0436, +0.0878] | [-0.0697, +0.1099] | 0.509 |
| Transport v2 | Correct at saliency 60 | -0.0114 | [-0.0591, +0.0346] | [-0.0797, +0.0594] | 0.603 |
| Transport v2 | Interaction per 20 saliency points | -0.0072 | [-0.0361, +0.0239] | [-0.0463, +0.0326] | 0.634 |
| Replay | Saliency 20 to 100 | +0.0236 | [-0.0542, +0.0934] | [-0.0851, +0.1215] | 0.493 |
| Replay | Correct at saliency 60 | +0.0007 | [-0.0545, +0.0527] | [-0.0763, +0.0747] | 0.988 |
| Replay | Interaction per 20 saliency points | +0.0002 | [-0.0305, +0.0333] | [-0.0459, +0.0484] | 0.989 |

Every ordinary and familywise interval includes zero. All Transport signs
agreed with the mixed model and were retained in all 36 leave-one-participant-
out fits, so the null verdict is not caused by a stability or convergence gate.

## Method landscape

| Method | Mean information bits | Mean log loss | Mean rank | Top-one credit |
|---|---:|---:|---:|---:|
| Transport v2 | +0.00728 | 1.60439 | 2.851 | 0.241 |
| Raw density, sigma 80, calibrated | +0.00594 | 1.60532 | 2.844 | 0.235 |
| Raw MultiMatch position, calibrated | +0.00604 | 1.60525 | 2.842 | 0.242 |
| Raw MultiMatch position EMD, calibrated | +0.00347 | 1.60703 | 2.869 | 0.231 |
| Replay | +0.00316 | 1.60725 | 2.836 | 0.239 |
| Raw density, sigma 160, calibrated | +0.00314 | 1.60726 | 2.853 | 0.241 |
| Registered density ridge | +0.00242 | 1.60776 | 2.869 | 0.234 |
| Registered elastic ridge | +0.00043 | 1.60914 | 3.032 | 0.175 |
| Registered MultiMatch ridge | -0.00083 | 1.61002 | 2.889 | 0.224 |

The remaining separately calibrated MultiMatch dimensions ranged from
-0.00208 to +0.00039 mean bits. None of the 33 comparator conditional tests
had a Benjamini-Hochberg adjusted value below 0.90.

The five-candidate chance top-one rate is 0.20. Several methods were
descriptively near 0.24, but the preregistered court did not define a clustered
inferential test of overall top-one performance. These values are therefore
weak descriptive evidence of item-identifying structure, not a positive
sensitivity gate and not evidence that Transport is more powerful.

## Transport optimization receipt

Profiling a 15-by-15 production-style alignment attributed 83% of runtime to
the R log-domain masked Sinkhorn projection. Transport now uses standard-domain
BLAS matrix scaling and automatically falls back to the original log-domain
implementation when numerical stability or convergence requires it.

The equivalence court includes randomized projection comparisons, fixed-mass
marginal invariants, an extreme-range fallback, complete alignment comparisons,
symmetry, chronology sensitivity, and convergence. It passes 117 focused
assertions. On five real cohort-spanning trials, the maximum absolute fast-
versus-reference score difference was `4.44e-16`, and the maximum coupling
difference was `2.22e-16`.

| Fold | Trials | Backend | Seconds | Convergence |
|---|---:|---|---:|---:|
| 1 | 332 | Log-domain reference | 3485.0 | 1.00 |
| 2 | 316 | Log-domain reference | 2859.8 | 1.00 |
| 3 | 315 | Fast with fallback | 830.8 | 1.00 |
| 4 | 332 | Fast with fallback | 596.6 | 1.00 |

Like-sized comparisons imply a 3.4-fold improvement for folds 2 versus 3 and
a 5.8-fold improvement for folds 1 versus 4. Replay still required only 13.8
to 14.9 seconds per fold, so it remains much more computationally parsimonious.
Transport's remaining cost comes largely from exhaustive training-candidate
calibration and repeated candidate batches. Changing calibration to fixed K=5
may be defensible, but it changes the calibration estimand and was not used to
accelerate this frozen court.

## Interpretation

GazeWeave remains more parsimonious than MultiMatch at the inferential level:
it returns one proper-score endpoint with coverage, warp, and chronology as
diagnostics. It is more interpretable than density similarity when the question
is *how* a partial replay unfolded, because its coupling or replay posterior can
be shown as a gaze braid. Those are structural advantages, not demonstrated
power advantages.

For this dataset, the honest conclusion is:

> Full-trial gaze carried at most very weak item-identifying information, and
> no evaluated method showed reliable modulation by saliency or recognition
> correctness. GazeWeave was competitive but not superior.

This court supersedes the eight-participant pilot for same-study sensitivity
claims. It remains an internal testbed and does not close the independent
Wang-data validation gate.
