# Concatenate Fixation Groups

Combines multiple fixation_group objects into one by row-binding. Onset
times are shifted so that each subsequent group continues from where the
previous one ended.

## Usage

``` r
# S3 method for class 'fixation_group'
c(..., recursive = FALSE)
```

## Arguments

- ...:

  `fixation_group` objects to concatenate.

- recursive:

  Ignored (present for S3 method compatibility).

## Value

A single `fixation_group` object.

## Examples

``` r
fg1 <- fixation_group(x = c(1, 2), y = c(3, 4),
                       onset = c(0, 100), duration = c(100, 100))
fg2 <- fixation_group(x = c(5, 6), y = c(7, 8),
                       onset = c(0, 100), duration = c(100, 100))
combined <- c(fg1, fg2)
```
