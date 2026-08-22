# Behavioral verdict for multi-presentation Replay

Status: complete
Protocol version: 1
Seed: 20260829
Measurement score frozen before response merge: yes

## Verdict

The confirmatory behavioral endpoint is **not supported**. On 989 old-item
trials, frozen all-four Replay improved held-out prediction by 0.00243 bits per
trial relative to saliency, generic gaze quality, and candidate-pool size. The
crossed participant-by-item 95% interval was [-0.01252, 0.01765].

Replay also failed to add predictive information beyond density sigma 80:
-0.00181 bits, 95% interval [-0.01505, 0.01106]. The result does not support a
claim that item-specific Replay information predicts recognition in this data.

## Prespecified predictive results

| Added gaze score | Mean behavioral bits | Crossed 95% interval |
|---|---:|---:|
| All-four shrunk Replay (primary) | 0.00243 | [-0.01252, 0.01765] |
| All-four unshrunk Replay | 0.00103 | [-0.01440, 0.01624] |
| Presentation-4 Replay | -0.01420 | [-0.05568, 0.01703] |
| Density sigma 80 | 0.00243 | [-0.00655, 0.01179] |
| Replay beyond density sigma 80 | -0.00181 | [-0.01505, 0.01106] |

All behavioral models were trained and evaluated on disjoint participants and
items. Their common base model contained linear saliency, log effective
fixation count, retained gaze duration, and log candidate count. Continuous
predictors were standardized inside each training fold.

## Why the earlier mixed-model result looked stronger

The descriptive crossed-random-intercept logistic model reproduced a nominally
positive Replay association:

```text
b = 0.300, SE = 0.120, z = 2.49, p = 0.0127
```

The fit converged and was not singular. However, the prespecified crossed
participant-by-item bootstrap gave a Replay coefficient of 0.289 with interval
[-0.094, 0.728], and the out-of-fold predictive endpoint was null. The ordinary
GLMM standard error is therefore too optimistic for the scientific conclusion
we need here. It is a useful descriptive specification check, not evidence that
survived the stronger inferential court.

Saliency remained the clear behavioral predictor: the per-10-point coefficient
was 0.287 with crossed interval [0.125, 0.500]. Generic gaze support was not
decisive after adjustment: log effective fixation count had interval [-0.337,
0.480], and retained duration had interval [-0.090, 0.659]. Candidate-pool size
also had no detectable association.

## Scientific interpretation

The post-mortem is now fairly sharp:

- GazeWeave is not simply broken. It strongly identifies repeated study viewing
  and detects the loss caused by complete reversal.
- Using all four study presentations materially improves the point estimate
  over presentation 4 alone, and response-blind reliability shrinkage behaves
  correctly.
- In the 3,000-ms probe-plus-delay retrieval interval, item-identifying gaze
  evidence is nevertheless small, uncertain, and no better than density sigma
  80.
- The nominal old-response GLMM association does not yield reliable crossed
  or out-of-fold evidence and is not specific beyond density.

No saliency subgroup, alternate time window, ordinal response, nonlinear term,
or phase-aware Replay search was performed. Those would be new exploratory or
future-study analyses, not rescues of this frozen court.
