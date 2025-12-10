# Estimate Scaling Parameters for Fixation Data

This function estimates the scaling parameters for two sets of fixation
data using the Hausdorff distance as an optimization objective.

## Usage

``` r
estimate_scale(x, y, lower = c(0.1, 0.1), upper = c(10, 10), window)
```

## Arguments

- x:

  A data frame or matrix containing the x coordinates of the first set
  of fixations.

- y:

  A data frame or matrix containing the y coordinates of the second set
  of fixations.

- lower:

  A numeric vector of length 2 specifying the lower bounds for the
  scaling parameters (default is c(0.1, 0.1)).

- upper:

  A numeric vector of length 2 specifying the upper bounds for the
  scaling parameters (default is c(10, 10)).

- window:

  A numeric vector of length 2 specifying the time window to restrict
  the fixations in \`y\` (default is NULL, which considers all
  fixations).

## Value

A list containing the estimated scaling parameters.
