# GazeWeave comparative validation protocol

This protocol is frozen before inspecting the comparative results. Its purpose
is to test bounded scientific claims, not to tune GazeWeave until it wins.

## Claims under test

1. Cross-fitted registration preserves item identification under a shared
   contraction and translation.
2. Local temporal transport distinguishes item-specific sequence when spatial
   occupancy is identical.
3. Unbalanced transport retains useful identification and explanatory
   correspondence under dropped and inserted fixations.

The validation does not establish performance on real looking-at-nothing data,
likelihood calibration, or reader-level interpretability.

## Frozen methods

- GazeWeave uses spatial sigmas 2.5, 5, and 10 pixels, local-order horizon 0.3,
  temporal weight 4, unmatched-mass penalties 2 and 2, entropy 0.03, and two
  deterministic optimizer starts.
- The density baseline is the duration-weighted Gaussian time marginal at the
  same three sigmas. Cosine similarity is calculated per scale and averaged.
- MultiMatch contributes its six dimensions separately. No composite or
  post-result selection is allowed.
- Density and MultiMatch are evaluated both on raw source coordinates and after
  the exact independently learned transform used by GazeWeave.
- MultiMatch family-level detection uses Holm correction across its six
  dimensions, separately for raw and registered coordinates.

## Scenarios

- `registered_geometry`: six spatially distinct paths undergo participant-level
  scale, translation, timing dilation, and spatial noise. The warp is learned
  from a separate calibration set.
- `order_at_fixed_density`: all six templates contain exactly the same points
  and duration mass but traverse them in different orders. Source timing is
  dilated and coordinates are jittered.
- `partial_replay`: the fixed-density construction is followed by two dropped
  fixations and one inserted fixation.

Each scenario is generated independently for 24 simulated participants with six
candidate templates. The random seed is 20260815.

## Endpoints

Discrimination is summarized by pairwise AUC, fractional top-1 credit under
ties, and reciprocal rank. Power is estimated by participant bootstrap for
sample sizes 8, 16, and 32 using 500 repetitions and one-sided alpha 0.05.
Each simulated study is evaluated against 199 conditional random-label draws;
this avoids relying on a Gaussian approximation for tied or discrete metrics.
GazeWeave uses participant-mean `gaze_bits`. Scalar baselines use
participant-mean true-minus-average-nonmatch similarity.

Type-I behavior is estimated by treating one random-label draw as the observed
null study and the remaining exchangeable draws as its reference distribution.
This preserves the complete candidate-score structure while breaking item
correspondence.

### Power-stress addendum

The initial full-signal run reached power 1.0 for both GazeWeave and
Holm-corrected MultiMatch at the smallest sample size, so it could not resolve a
power difference. After observing that ceiling, but without changing any model
or comparator parameter, a declared stress grid was added. Item-specific labels
are retained with probability 0.25, 0.50, 0.75, or 1.00; otherwise a candidate
label is drawn uniformly. This models the prevalence of item-specific replay
and produces power curves below the ceiling. Because this addendum was motivated
by the first-stage result, it is identified explicitly rather than described as
preregistered.

## Interpretability diagnostics

- Registration: bias and RMSE of the independently recovered scale.
- Correspondence: proportion of transported mass assigned to known fixation
  links.
- Omissions: AUC for assigning more missing reference mass to deliberately
  dropped than retained fixations.
- Numerical audit: convergence rate of true-template alignments.

These are ground-truth diagnostic-recovery tests. They do not show that human
readers understand a gaze braid better than a vector of scores.

## Decision rules

GazeWeave can be described as showing superior simulated discrimination in a
scenario only if its mean pairwise AUC and top-1 credit exceed every registered
scalar baseline. A power advantage requires higher estimated power than both
registered density and Holm-corrected registered MultiMatch while its estimated
type-I error is no greater than 0.075. Claims remain scenario-specific.
