# Evaluate replayable mismatched-template draws

Evaluate replayable mismatched-template draws

## Usage

``` r
mismatch_draws(
  strata,
  refs,
  n = 100,
  trial = "trial",
  seed = 1,
  assignments = NULL,
  replace_ref = c("self_item", "all"),
  restandardize = TRUE,
  standardize = c("within_subject", "none"),
  by = c("fixation", "trial"),
  control_policy = c("exclude_replacement", "fixed")
)
```

## Arguments

- strata:

  A case-control table from
  [`sample_controls()`](https://bbuchsbaum.github.io/eyesim/reference/sample_controls.md).

- refs:

  A reference set.

- n:

  Number of null draws.

- trial:

  Trial column, nested within the reference subject column.

- seed:

  Seed, or NULL to use the session RNG.

- assignments:

  Saved output of
  [`mismatch_assignments()`](https://bbuchsbaum.github.io/eyesim/reference/mismatch_assignments.md).
  If NULL, create `n` draws using `seed`; supplied assignments are
  evaluated without RNG use.

- replace_ref:

  Replace only self_item or the entire reference tuple.

- restandardize:

  Refit scaling in every null draw. FALSE applies the original observed
  fit. Each paired matched/mismatched comparison uses the same fit
  (estimated on retained mismatched scores when TRUE). This choice alone
  does not establish an exact randomization test: calibration also
  requires an appropriate exchangeability assumption.

- standardize:

  Score standardization, as in
  [`score_strata()`](https://bbuchsbaum.github.io/eyesim/reference/score_strata.md).

- by:

  Contrast resolution, as in
  [`case_control_contrast()`](https://bbuchsbaum.github.io/eyesim/reference/case_control_contrast.md).

- control_policy:

  By default, `"exclude_replacement"` removes observed controls donated
  by any trial of either the case's item or the replacement item, using
  the saved donor metadata. Cases are always retained. `"fixed"` retains
  all controls: this item-discrimination baseline can be negative under
  a genuine matched effect because it includes the replacement item's
  own test fixations.

## Value

A list with `observed` contrasts on the original controls,
`paired_observed` matched contrasts on each draw's retained controls and
scaling, `contrasts` for mismatched templates, assignments, original
observed scaling, and settings. Both paired tables carry `null_draw` and
identical row keys. Compare `paired_observed` with `contrasts` to use
the same controls and units. Evaluation proceeds one draw at a time.

## Details

Exclusion uses donor items, including repeated trials of either item,
rather than assuming trial labels equal item labels. It does not redraw
controls. A stratum with no remaining controls is an error: a donor pool
containing only the case and replacement items cannot support the paired
comparison, regardless of the number of controls drawn. In particular,
two test trials fail when they contain those two items. Generated
uniform/map controls have no donor item and are retained. For older
saved strata without donor metadata, complete case rows can supply the
donor keys; otherwise resample and save the strata or choose the fixed
policy. No automatic p-values or calibrated effect decomposition are
computed.
