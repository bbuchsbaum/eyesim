# Own-versus-other-participant recognition decomposition

This exploratory court separates two sources of study-to-recognition gaze
similarity. The canonical model uses the participant's own four study
episodes. The control model replaces each candidate with four study episodes
from other participants who viewed the same image version. Both models
use the frozen Transport v3 specification, warp, calibration, trials, and
five-candidate panels.

The other-participant model measures image-driven spatial organization without
using the retrieval participant's encoding gaze. Its donors are chosen
deterministically, without response, saliency, accuracy, confidence, or gaze
score. Every candidate receives the four distinct presentation numbers. Donor
selection maximizes the number of distinct other participants; only sparse
item/version cells reuse a donor across presentations. The fold audit reports
the observed donor support.

`newtest` trials are the no-own-memory control. Each has no study episode for
the retrieval participant, so it can be scored only against other-participant
same-image templates. Candidate items come from the frozen item fold, and all
study donors are selected before any response is inspected.

The measurement endpoints are true-minus-nonmatch candidate margin, top-one
credit relative to 20% chance, and the paired own-minus-group differences.
Uncertainty uses the crossed participant/item bootstrap. Only after the
measurement comparison is complete does the ordinal 1--4 newness analysis ask
whether own-template margin adds held-out information beyond group-template
margin, saliency, effective fixation count, and retained duration.

Exact visibility masks and visible-versus-hidden region analyses are excluded
from this court by scope. Therefore a positive group-template result remains a
perceptual or image-structure control, not a memory effect.

Run the four folds and then finalize:

```r
source("inst/validation/gaze-weave-own-group-decomposition.R")
run_gaze_weave_own_group_measurement()
result <- own_group_finalize()
```
