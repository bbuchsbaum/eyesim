# Draw replayable case-control locations

Creates explicit assignments. Scoring this table never samples again.

## Usage

``` r
sample_controls(
  fix_tab,
  n = 20,
  subject = "subject",
  trial = "trial",
  stratify_by = NULL,
  replace = "auto",
  seed = 1,
  draws = 1,
  method = c("observed", "uniform", "map"),
  xbounds = NULL,
  ybounds = NULL,
  control_map = NULL,
  item = if ("item" %in% names(fix_tab)) "item" else NULL
)
```

## Arguments

- fix_tab:

  A flat fixation table with finite x/y and case metadata.

- n:

  Number of controls per case.

- subject, trial:

  Columns identifying subjects and trials (trial is nested within
  subject).

- stratify_by:

  Additional columns restricting observed donor pools.

- replace:

  `"auto"` uses replacement only if an observed pool is too small;
  alternatively TRUE or FALSE. Map controls always sample independently
  with replacement; `replace = FALSE` is an error for that method.

- seed:

  Seed for assignment generation; restores the session RNG state. NULL
  explicitly uses the session RNG.

- draws:

  Number of independent control draws.

- method:

  Observed donor fixations, uniform coordinates, or grid locations
  sampled with probability proportional to a supplied map's nonnegative
  mass.

- xbounds, ybounds:

  Explicit coordinate bounds for uniform controls.

- control_map:

  A density map for `method = "map"`.

- item:

  Item column for observed donor exclusion. Defaults to `"item"` when
  that column exists, otherwise NULL. Pass another column name for
  custom item keys. NULL explicitly retains trial-only exclusion.
  Ignored for uniform/map controls.

## Value

An `eyesim_strata` data frame. Metadata belongs to the case; x/y are the
case or control coordinates. `case_id` and `control_id` index input rows
(control_id is NA for generated locations); `draw`, `stratum`, and
`is_case` identify assignments. Observed controls always exclude the
case's trial and, when `item` is supplied, all trials of the case's
item. The `control_sources` attribute retains input atomic metadata,
indexed by `control_id`, so donor items remain identifiable after
subsetting cases.

## Details

There is no control-only jitter: perturbing controls but not cases
changes their sampling distribution and can induce spurious preference.
Uniform/map controls sample a different baseline from observed controls.
Excluding the case's item avoids attenuating matched contrasts when
items are presented repeatedly. Without item metadata only trial
exclusion is possible. The resolved item column is recorded in the
`sampling` attribute.
