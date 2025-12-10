# Compute Similarity Between Two Objects

This function computes the similarity between two objects using a
specified similarity metric.

## Usage

``` r
similarity(x, y, method, ...)
```

## Arguments

- x:

  The first object to compare.

- y:

  The second object to compare.

- method:

  A character string specifying the similarity metric to be used.

- ...:

  Additional arguments passed to the similarity computation method.

## Value

A numeric value representing the similarity between the two input
objects.

## See also

Other similarity:
[`similarity.scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/similarity.scanpath.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage of the similarity function
similarity_value <- similarity(object1, object2, method = "pearson")
} # }
```
