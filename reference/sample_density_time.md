# Sample density maps at fixation locations over time

This function samples template density maps at source fixation locations
across specified time points. It is useful for analyzing how fixation
patterns relate to a reference density map over the course of a trial
(e.g., "reinstatement across time" analysis).

## Usage

``` r
sample_density_time(
  template_tab,
  source_tab,
  match_on,
  times = seq(0, 3000, by = 50),
  time_bins = NULL,
  template_var = "density",
  source_var = "fixgroup",
  permutations = 0,
  permute_on = NULL,
  aggregate_fun = mean
)
```

## Arguments

- template_tab:

  A data frame or tibble containing template density maps (e.g., from an
  encoding phase).

- source_tab:

  A data frame or tibble containing fixation groups (e.g., from a
  retrieval phase).

- match_on:

  A character string specifying the column name used to match template
  density maps to source fixation groups.

- times:

  A numeric vector of time points (in ms) at which to sample the
  density. Default is `seq(0, 3000, by = 50)`.

- time_bins:

  An optional numeric vector specifying bin boundaries for aggregating
  samples. For example, `c(0, 1000, 2000, 3000)` creates 3 bins:
  \[0-1000), \[1000-2000), \[2000-3000). Default is NULL (no binning).

- template_var:

  A character string specifying the name of the density column in
  `template_tab`. Default is "density".

- source_var:

  A character string specifying the name of the fixation group column in
  `source_tab`. Default is "fixgroup".

- permutations:

  An integer specifying the number of permutation samples for computing
  a baseline. Default is 0 (no permutations).

- permute_on:

  An optional character string specifying the column used to stratify
  permutations (e.g., "subject" to permute within subjects).

- aggregate_fun:

  A function used to aggregate density values within time bins. Default
  is `mean`.

## Value

A tibble containing:

- All columns from `source_tab`

- `sampled`: A list column with data frames containing `z` (density
  values) and `time` columns for each row

- If `time_bins` is provided: columns `bin_1`, `bin_2`, etc. with
  aggregated values for each bin

- If `permutations > 0`: `perm_sampled` (list column with mean permuted
  trajectory) and `perm_bin_1`, `perm_bin_2`, etc.

## Details

The function matches each row in `source_tab` to a corresponding density
map in `template_tab` based on the `match_on` column. For each matched
pair, it samples the template density at the fixation coordinates
interpolated at each time point in `times`.

This approach avoids the problem of having too few fixations to compute
density maps within short time windows. Instead, a single density map is
created from all fixations (e.g., during encoding), and the sampling is
done over time (e.g., during retrieval).

If `time_bins` is provided, the sampled values are aggregated within
each bin using `aggregate_fun`. The result includes columns named
`bin_1`, `bin_2`, etc.

If `permutations > 0`, a baseline is computed by sampling from
non-matching density maps. The result includes `perm_sampled` (mean
permuted trajectory) and bin-specific permutation columns if binning is
used.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create encoding density maps
encoding_dens <- eyetab %>%
  filter(condition == "encoding") %>%
  density_by(groups = c("subject", "image"), sigma = 200)

# Get retrieval fixations (as eye_table, not density)
retrieval_fix <- eyetab %>%
  filter(condition == "retrieval")

# Sample encoding density at retrieval fixation locations over time
results <- sample_density_time(
  template_tab = encoding_dens,
  source_tab = retrieval_fix,
  match_on = "subject_image",
  times = seq(0, 5000, by = 50),
  time_bins = c(0, 1000, 2000, 3000, 4000, 5000),
  permutations = 50,
  permute_on = "subject"
)
} # }
```
