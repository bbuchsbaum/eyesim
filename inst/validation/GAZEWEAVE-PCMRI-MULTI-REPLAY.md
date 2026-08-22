# GazeWeave multi-presentation Replay reliability court

Status: frozen before implementation scoring
Protocol version: 1
Seed: 20260828
Mote: `bd-01M0GTPJWX3TP6S65QW6RP2FJ8`

## Scientific question

The previous recognition courts used study presentation 4 as a single encoding
template. This court asks whether retrieval measurement improves when all four
observed encoding paths inform the episode, the nonmatch denominator is exact,
and low-information retrieval paths are allowed to abstain.

This is a measurement-development court. Response correctness, confidence,
and saliency effects are not inspected until the score and its calibration are
frozen.

## Encoding episode

For participant `p`, item `k`, and retrieval path `Y`, retain the four study
paths `X_pkr` separately. The candidate score is the equal-weight mixture of
resolution-normalized Replay likelihoods:

```text
s_pk(Y) = log[(1/4) sum_r exp{s_pkr(Y)}]
```

where each `s_pkr` is the mean HMM log likelihood per normalized duration bin.
No path is concatenated, and no transition is introduced between study
presentations. The mixture posterior over presentations is retained as a
diagnostic. Equal presentation weights are primary. Fixation-count weighting
is not introduced unless a later training-only court establishes that it
improves proper scores.

The nuisance warp is fit from all four training presentations with equal
presentation weight and is then fixed for every candidate in the held-out
fold.

## Candidate evidence

For each evaluated participant and item fold, use every eligible studied item
for that participant in the fold. The candidate pool therefore varies with
available observations but contains the true item exactly once and at least
one nonmatch. Candidate prior is uniform over the declared exhaustive pool.

The exact exhaustive score is primary. Deterministic K=5 panels are reconstructed
only from already-computed candidate scores to quantify denominator sampling
variance; they cannot replace the exhaustive endpoint unless a prospectively
declared approximation achieves a prespecified reliability tolerance.

## Response-blind reliability shrinkage

Candidate evidence is first calibrated with a training-only temperature. A
candidate-invariant support statistic is computed from the retrieval path
after coalescing immediately adjacent identical fixation atoms:

```text
n_eff = 1 / sum_i a_i^2
q(Y)  = n_eff / (n_eff + kappa)
```

The final posterior is

```text
p_final(k | Y) = q(Y) p_replay(k | Y) + [1 - q(Y)] prior(k).
```

Temperature and non-negative `kappa` are estimated jointly from inner-held-out
candidate log loss, with weak penalties toward temperature one and `kappa = 0`.
Neither response, accuracy, saliency, probe type, candidate identity, nor the
true-candidate indicator may enter `q`. Splitting an identical fixation into
adjacent duration-preserving atoms must not change `n_eff` or the final score.

## Generic gaze-quality channel

The alignment retains, as diagnostics only:

- raw and coalesced fixation counts;
- effective fixation count;
- retained gaze duration;
- duration concentration;
- reliability weight `q`;
- Replay/background coverage and spatial residual;
- presentation-mixture weights and effective template count.

These quantities may predict behavior as generic engagement or measurement
quality. They are not included in `gaze_info_bits` except through the declared
response-blind shrinkage and are never labelled item reinstatement.

## Data and folds

- Use the complete recognition interval `[0, 3000)` ms.
- Keep old and lure trials with at least one retained recognition fixation and
  all four study presentations, each with at least one retained fixation.
- Use the existing 800 by 600 image rectangle and greater-than-80-ms retained
  duration rule.
- The frozen eligible pool is 2,056 trials from 46 participants. The exhaustive
  candidate rule retains 2,055 trials from 45 participants and 120 items:
  1,011 old and 1,044 lure trials.
- Keep the fixed item assignment seed `20260820` and four crossed outer folds.
  Score training and evaluation participants and items are disjoint.
- The 90 unique participant-by-item-fold pools contain 5--33 items (median 25,
  mean 22.83). Across target trials, where larger pools necessarily occur more
  often, the median candidate count is 27 and the mean is 25.19.

## Frozen comparisons

1. presentation-4 Replay, exhaustive candidates, temperature only;
2. all-four equal-mixture Replay, exhaustive candidates, temperature only;
3. all-four equal-mixture Replay, exhaustive candidates, joint temperature and
   reliability shrinkage (primary);
4. density sigma 80 and sigma 160 using the equal-weight mean of four
   separately normalized study density maps and the identical exhaustive
   candidate pools. Each raw cosine score receives a training-only temperature;
   neither density method is registered.

The existing K=5 P4 courts remain historical references, not co-primary tests.

## Measurement analyses

- held-out mean `gaze_info_bits`, log loss, rank, and top-one credit;
- paired crossed participant-by-item intervals for the three Replay variants;
- log-loss improvement from shrinkage overall and within prespecified
  effective-fixation bands `[1,2)`, `[2,3)`, `[3,5)`, and `5+`;
- monotonic association between `n_eff` and unshrunk information, with crossed
  uncertainty;
- presentation responsibility versus presentation number and presentation
  `n_eff`, as diagnostics of the equal-weight mixture;
- exact denominator versus deterministic K=5 panel variance;
- repeated-viewing positive control (presentation 4 identified from the
  equal-weight mixture of presentations 1--3) and a complete-reversal
  negative control that preserves locations and fixation durations;
- runtime and candidate-count scaling.

No recognition-response model is fit during this stage.

## Advancement gates

The improved score is frozen for behavioral testing only if:

1. analytic, differential, order-invariance, fixation-refinement, and
   one-template-equivalence tests pass;
2. every evaluated item has four separate template components, one true
   candidate, finite normalized probabilities, and zero train/evaluation
   participant or item overlap;
3. shrinkage never worsens inner calibration loss relative to its own fitted
   temperature-only boundary solution beyond numerical tolerance;
4. low-support paths move toward, not away from, the candidate prior;
5. repeated-viewing sensitivity remains positive and complete reversal loses
   information;
6. the exhaustive score and its candidate-pool provenance are reproducible;
7. retrieval measurement results and runtime are reported whether favorable or
   not.

Passing these gates freezes the algorithm; it does not establish a memory
effect. Only then may the old-trial behavioral court compare saliency plus
generic gaze quality against the incremental item-specific score with crossed
participant-by-item uncertainty.

Raw fixation data, participant-linked scores, candidate tables, calibration
fits, and bootstrap draws remain local-only and Git-ignored.
