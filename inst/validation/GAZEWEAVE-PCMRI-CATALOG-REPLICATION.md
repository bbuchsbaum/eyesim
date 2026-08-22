# GazeWeave pcmri catalog replication

Status: frozen before reconstruction scoring
Protocol version: 1
Seed: 20260820
Mote: `bd-01M0GK8JRQNKXFC6P68WWGGZV7`

## Purpose and evidential boundary

This court reconstructs the public *pcmri eye-movement reinstatement --
behavioral catalog (PROBE-WINDOW)* from the private Wynn fixation exports and
compares its density and MultiMatch results with GazeWeave. It is a same-data
method comparison, not an independent replication.

The public report is an analysis of precomputed tables, not a complete metric
pipeline. Its embedded source loads three CSV files that were not deployed,
and it does not declare the density bandwidth or all permutation-generation
details. Consequently, two results must remain distinct:

1. **report audit**: reproduce values printed by the deployed report from its
   embedded tables and code;
2. **raw-data reconstruction**: recompute the closest documented metric from
   the private raw fixation exports.

The reconstruction uses duration-weighted density maps with `sigma = 80` px,
an 80 by 60 grid, and the 800 by 600 image rectangle. These settings come from
the repository's same-paradigm analysis scripts. They are strong local
provenance, but are not claimed to be the exact settings of the undeployed
public-report generator.

Raw fixation files, participant-linked metrics, candidate scores, model
checkpoints, and bootstrap draws remain local-only and Git-ignored. Only
aggregate protocols and verdicts may be published.

## Audit findings that constrain interpretation

The deployed report calls itself a probe-window analysis but describes and
loads the complete test-picture, mask, and delay interval of about 3 seconds.
For this court, `combined` means the half-open interval `[0, 3000)` ms: the
500-ms probe plus the 2500-ms delay. Probe-only and delay-only analyses are
separate sensitivity analyses and must not be substituted for this estimand.

The report's numerical tables and executable code take precedence over nearby
prose. Several narrative summaries appear to have survived from an older
render and contradict the current tables. In particular:

- the printed old-trial density model reports a positive response association
  (`z = 2.13`, `p = .034`), while later prose quotes weaker older values;
- the printed decomposition retains both item-specific density difference and
  nonspecific permutation similarity, while its prose describes the former as
  only a trend;
- all five MultiMatch chance differences are nominally positive in the printed
  table, but the summary says only position and duration reinstate;
- no MultiMatch dimension predicts old-trial responses in the printed
  dimension-wise models;
- the printed combined model gives no shape effect (`p = .537`), while the
  summary claims that shape predicts recognition;
- the subject random-slope summary reports 93 percent positive slopes, not
  every subject.

The public analysis runs many uncorrected fixed-effect tests, generally uses
random intercepts without corresponding within-subject random slopes, and
does not evaluate held-out calibration or item identification. Its findings
are therefore exploratory associations, not evidence that density has greater
predictive power than GazeWeave.

## Frozen reconstruction

### Geometry and paths

- Read `test_data/wynn_probe_delay/study_fix_input_new.csv` and
  `test_data/wynn_probe_delay/testdelay_fix_input_matched.csv` through the
  existing verified importers.
- Translate the image rectangle x = 112--912 and y = 84--684 to 800 by 600
  analysis coordinates.
- For the **catalog reconstruction**, use all valid study fixations from all
  four presentations to form each participant-by-image density template.
  Pool fixation mass across presentations.
- For the **fair method comparison**, use presentation 4 as the common
  reference for every method. This is the already-scored GW-13 estimand and
  avoids inventing transitions or giving density a four-presentation template
  while sequence methods receive only one path.
- A later multi-presentation sequence analysis may score the four study paths
  separately and marginalize their evidence, but may not concatenate them.
- Use the complete recognition interval `[0, 3000)` ms.
- Retain old and lure trials for own-template reinstatement. Foils have no
  participant-specific studied target and are reserved for group-template or
  generic-gaze analyses.

### Common trial support

The primary fair-comparison table includes every old/lure trial with:

- a corresponding four-presentation study template;
- at least three retained study fixations in every presentation;
- at least three retained recognition fixations in the combined interval;
- a nonzero response for behavioral models.

Report full-data density descriptives additionally on every trial for which
that metric is defined, but do not compare methods on unequal support. Every
method receives the same target trials and the same nonmatching candidate
identities.

### Candidate and permutation contract

The public report uses matched-minus-random-other-image similarity, whereas
GazeWeave's primary endpoint is held-out item information. Report both without
equating them:

1. `catalog_diff`: matched compatibility minus the mean compatibility of the
   same fixed nonmatching images;
2. `gaze_info_bits`: held-out information gain for the true image relative to
   its candidate prior.

Use five candidates per trial, including the true image exactly once, drawn
deterministically within participant and item fold. The same four nonmatches
must be used for density, MultiMatch, Replay, and Transport v2. Registration,
model fitting, and temperature calibration must exclude the evaluated
participant and item. Candidate order may not affect initialization or score.

`catalog_diff` is used only to reproduce the public report's estimand.
`gaze_info_bits` remains the primary fair method-comparison endpoint.

### Methods

- duration-weighted Pearson density similarity, `sigma = 80` px, 80 by 60
  output grid;
- the five standard MultiMatch dimensions, plus position EMD as an eyesim
  diagnostic;
- GazeWeave Replay as the primary GazeWeave engine;
- GazeWeave Transport v2 as a slower secondary engine;
- registered density and ridge-trained MultiMatch composites as fair learned
  baselines where already defined by the GazeWeave court.

The report's five separate MultiMatch dimensions are descriptive. For claims
of comparative predictive performance, the baseline is the nested-CV ridge
composite, not the best post hoc dimension.

## Analyses

### Metric sensitivity

For each method, report held-out mean information bits, crossed
participant-by-item 95 percent intervals, top-one credit, mean rank, and log
loss. Also report participant-weighted mean `catalog_diff` and its crossed
interval. Study P1-to-P4 results from GW-15 remain the positive-control ceiling
and are not recomputed under a different definition.

### Catalog reconstruction

Reconstruct these prespecified report questions:

- mean matched-minus-nonmatch reinstatement for own-template density;
- old versus lure difference;
- linear saliency slope per 10 percentage points and its lure interaction;
- old-trial `said_old` association, with saliency controlled;
- lure false-alarm association;
- correlations between density and MultiMatch differences;
- old-trial behavior models for each MultiMatch dimension and the nested-CV
  composite.

Code response 1 or 2 as `said_old = 1`, response 3 or 4 as `said_old = 0`, and
response 0 as missing. Code `sal10 = (saliency - 60) / 10`.

For direct numerical comparison with the report, fit its crossed random-
intercept model `(1 | participant) + (1 | item)`. Label those estimates
`catalog-compatible`. The inferential result uses a participant random slope
for the within-participant gaze predictor where identifiable and a crossed
participant-by-item bootstrap. Disagreement between these two specifications
is reported, not optimized away.

### Multiplicity and claims

- Treat the density old-trial response association as the single catalog
  replication target.
- Treat the five MultiMatch chance tests as one family and adjust them
  together.
- Treat dimension-wise behavior models as exploratory unless the learned
  composite is positive out of fold.
- Do not call a method more powerful from a smaller p-value. Compare methods
  by held-out log loss/information and by bootstrap contrasts on common trials.
- Do not give mechanistic labels such as familiarity, recollection, or fluency
  to a similarity association without an independently identifying design.

### Behavioral prediction follow-up

Because the initial catalog-compatible model gave Replay a nominal old-response
association, test that signal predictively. For each crossed outer fold, train
two logistic models on trials whose participants and items are both absent from
the evaluation fold:

```text
null: said_old ~ sal10
full: said_old ~ sal10 + standardized gaze_info_bits
```

Standardize gaze information using training data only. The endpoint is held-out
log-loss improvement `(loss_null - loss_full) / log(2)` in bits per trial.
Use identical old-trial support for all methods and paired participant-by-item
bootstrap weights. Report each method's crossed 95 percent interval and the
Replay-minus-comparator intervals; control the five Replay comparisons as one
family. Also bootstrap the all-trial standardized association with the same
crossed weights. Replay has behavioral evidence only if its own predictive
interval and the familywise Replay-minus-density interval are both above zero.

## Gates

The reconstruction passes only if:

1. all source checksums match the private manifest;
2. trial-window, geometry, response, condition, and repetition mappings pass
   explicit tests;
3. density sigma 80 reproduces the public density result's sign and broadly
   comparable scale, with any exact mismatch attributed rather than hidden;
4. candidate identities and support are identical across compared methods;
5. cross-fitting has zero evaluated-participant and evaluated-item leakage;
6. Replay and density are both complete before comparative conclusions;
7. a Transport comparison is described as pending until all prespecified
   folds finish;
8. aggregate results distinguish report replication, raw reconstruction,
   item-identification evidence, and behavioral association.
