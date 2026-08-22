# GazeWeave recognition sensitivity verdict

Verdict date: 2026-08-20
Protocol: `GAZEWEAVE-RECOGNITION-SENSITIVITY.md`, version 1
Mote: `bd-01M0FQ85Z8T3FW3Y5F1P9FE25V`
Verdict: `not_supported`

## Bottom line

Neither GazeWeave engine showed reliable sensitivity to probe saliency,
response correctness, or their interaction in the complete 0--3000 ms
recognition interval. All six ordinary 95% GazeWeave intervals included zero,
so the court failed before the stricter familywise and stability gates were
needed.

This result does not show that saliency or successful recognition is unrelated
to gaze reinstatement in general. It shows that, in this small secondary sample,
the held-out item-information scores did not resolve those conditional
associations. The result cannot rescue the failed post-probe persistence gate,
select a default engine, establish superiority over the baselines, or close the
independent Wang-data gate.

## Frozen full-trial effects

Effects are in held-out item-information bits. `saliency_20_to_100` is the
predicted change from saliency 20 to 100; `correct_at_60` is correct minus
incorrect at saliency 60; and `interaction_per_20` is the change in the
correctness contrast for a 20-point saliency increase. The family interval is
the prespecified Bonferroni 99.1667% interval over both engines and all three
contrasts.

| Engine | Contrast | Estimate | 95% interval | Family interval |
|---|---|---:|---:|---:|
| Transport v2 | Saliency 20 to 100 | +0.144 | [-0.425, +0.729] | [-0.706, +1.446] |
| Transport v2 | Correct at saliency 60 | -0.067 | [-0.498, +0.320] | [-0.818, +0.516] |
| Transport v2 | Saliency by correctness | -0.018 | [-0.276, +0.256] | [-0.585, +0.378] |
| Replay | Saliency 20 to 100 | +0.252 | [-0.306, +0.962] | [-0.485, +1.926] |
| Replay | Correct at saliency 60 | -0.171 | [-0.755, +0.222] | [-1.189, +0.395] |
| Replay | Saliency by correctness | -0.036 | [-0.330, +0.273] | [-0.793, +0.425] |

The saliency point estimates were positive for both engines, but the intervals
were broad and crossed zero substantially. Correct trials did not have higher
GazeWeave information at the centered saliency level; the estimates instead
were slightly negative and similarly uncertain.

## Probe and delay localization

The phase-localization analyses were also null. During the 0--500 ms probe,
Transport v2 estimated -0.068 bits from saliency 20 to 100 and +0.022 bits for
correctness at saliency 60; Replay estimated -0.014 and +0.032 bits. During the
500--3000 ms delay, the corresponding estimates were +0.134 and -0.084 bits
for Transport v2, and +0.323 and -0.172 bits for Replay. Every ordinary and
familywise interval included zero.

The complete-trial score is computed by fitting the entire scanpath and is not
an arithmetic average of separately calibrated probe and delay scores. The
phase analyses therefore localize the absence of a reliable effect; they are
not decomposition terms for the combined result.

## Baseline comparison

Registered density and the learned MultiMatch composite had no ordinary 95%
full-trial effect. Registered elastic matching produced the only ordinary
interval excluding zero: correct trials were +0.151 bits higher at saliency 60,
with a 95% interval of [+0.012, +0.365]. Its wider family-comparable interval,
[-0.021, +0.514], included zero, and neither its probe nor delay localization
showed an ordinary effect. This is a secondary comparator finding, not evidence
of GazeWeave sensitivity or elastic-matching superiority.

## Integrity and limitations

- The frozen analysis used 80 trials from eight participants: 59 correct and
  21 incorrect, spanning 61 items.
- The fixed-effects design condition number was 4.37, below the frozen
  instability threshold of 30.
- Exactly 1,995 of 2,000 shared crossed participant-by-item bootstrap draws
  were valid for every fitted contrast.
- All optimizers converged. Thirteen of 14 participant/item random-intercept
  audit fits were singular, including every full-trial fit. This is expected
  from the sparse item replication and is why the crossed bootstrap, not the
  mixed model, supplies primary uncertainty.
- The analysis was proposed after the primary persistence result and after
  some saliency-stratified delay descriptives had been seen. It is explicitly
  secondary and partially post hoc.
- Correctness was strongly associated with probe type: 43 of 48 old trials but
  16 of 32 lure trials were correct. The frozen model therefore adjusted for
  probe type rather than interpreting a raw correct/incorrect difference.

The evidence supports a practical conclusion: this dataset is suitable for
testing the complete recognition interval, but it does not establish the
desired sensitivity. A larger prospectively analyzed sample with more errors,
more repeated items, and a prespecified saliency-by-correctness hypothesis is
needed before using these associations as a validation claim.

## Reproduction and local-data policy

From the package root:

```r
source("inst/validation/gaze-weave-recognition-sensitivity.R")
run_gaze_weave_recognition_sensitivity(
  "inst/validation/gaze-weave-recognition-sensitivity-results"
)
```

The raw Wynn exports, GW-11 trial scores, and this court's participant-linked
outputs are Git-ignored and must not be published. The protocol, runner, tests,
and this aggregate human verdict are publishable because they contain no
participant-level gaze records.
