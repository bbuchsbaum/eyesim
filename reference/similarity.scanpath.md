# Compute Similarity Between Scanpaths

This function computes the similarity between two scanpaths using a
specified method.

## Usage

``` r
# S3 method for class 'scanpath'
similarity(
  x,
  y,
  method = c("multimatch"),
  window = NULL,
  screensize = NULL,
  ...
)
```

## Arguments

- x:

  A scanpath object containing the first scanpath.

- y:

  A scanpath object containing the second scanpath.

- method:

  A character string specifying the method to compute the similarity
  (default is "multimatch").

- window:

  A numeric vector of length 2 specifying the time window to restrict
  the fixations in the input scanpaths (default is NULL, which considers
  all fixations).

- screensize:

  A numeric vector of length 2 specifying the dimensions of the screen
  (e.g., c(1000, 1000)). Required for the "multimatch" method.

- ...:

  Additional arguments passed to the similarity computation method.

## Value

A numeric value representing the similarity between the two input
scanpaths.

## See also

Other similarity:
[`similarity()`](https://bbuchsbaum.github.io/eyesim/reference/similarity.md)

## Examples

``` r
# Example usage of the similarity.scanpath function
scanpath1 <- # first scanpath data
scanpath2 <- # second scanpath data
similarity_value <- similarity.scanpath(scanpath1, scanpath2, method = "multimatch", screensize = c(1000, 1000))
#> Error in similarity.scanpath(scanpath1, scanpath2, method = "multimatch",     screensize = c(1000, 1000)): could not find function "similarity.scanpath"
```
