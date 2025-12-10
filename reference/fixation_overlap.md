# Fixation Overlap Measure

Calculate the number of overlapping fixations between two fixation
groups, based on a minimum distance threshold, and the percentage of
overlapping fixations.

## Usage

``` r
fixation_overlap(
  x,
  y,
  dthresh = 60,
  time_samples = seq(0, max(x$onset), by = 20),
  dist_method = c("euclidean", "manhattan")
)
```

## Arguments

- x:

  A \`fixation_group\` object representing the first fixation group.

- y:

  A \`fixation_group\` object representing the second fixation group.

- dthresh:

  A numeric value specifying the distance threshold to determine when
  two fixations overlap (default is 60).

- time_samples:

  A numeric vector of points in time at which to evaluate the
  overlapping fixations (default is \`seq(0, max(x\$onset), by = 20)\`).

- dist_method:

  A character string specifying the distance metric to use for measuring
  the distance between fixations. Options are "euclidean" and
  "manhattan" (default is "euclidean").

## Value

A list containing the following elements:

- overlap The number of overlapping fixations between the two fixation
  groups.

- perc The percentage of overlapping fixations.

## Details

This function computes the number of overlapping fixations within a
minimum distance for two fixation groups at specified time points. An
overlap occurs when the distance between two fixations is less than the
specified distance threshold.

## Examples

``` r
# Create two fixation groups
fg1 <- fixation_group(x = runif(50, 0, 100), y = runif(50, 0, 100),
duration = rep(1, 50), onset = seq(1, 50))
fg2 <- fixation_group(x = runif(50, 0, 100), y = runif(50, 0, 100),
duration = rep(1, 50), onset = seq(1, 50))

# Calculate the number of overlapping fixations and the percentage of overlapping fixations
result <- fixation_overlap(fg1, fg2)
overlap <- result$overlap
#> Error in result$overlap: $ operator is invalid for atomic vectors
perc <- result$perc
#> Error in result$perc: $ operator is invalid for atomic vectors
```
