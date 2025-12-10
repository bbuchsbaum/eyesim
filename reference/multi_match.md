# Compute MultiMatch Metrics for Scanpath Similarity

This function computes multiple similarity metrics between two
scanpaths, including vector, direction, length, position, duration, and
EMD-based position similarity.

## Usage

``` r
multi_match(x, y, screensize)
```

## Arguments

- x:

  A data frame representing the first scanpath. Must contain at least
  three columns: `x`, `y`, and `onset`, and at least three rows.

- y:

  A data frame representing the second scanpath. Must contain at least
  three columns: `x`, `y`, and `onset`, and at least three rows.

- screensize:

  A numeric vector of length 2 indicating the width and height of the
  screen in pixels.

## Value

A named numeric vector with the following elements:

- mm_vector:

  Similarity based on the 2D vectors between fixations.

- mm_direction:

  Similarity based on the direction (angle) of saccades between
  fixations.

- mm_length:

  Similarity based on the length of saccades between fixations.

- mm_position:

  Similarity based on the spatial position of fixations.

- mm_duration:

  Similarity based on the duration of fixations.

- mm_position_emd:

  Order-insensitive similarity based on the Earth Mover's Distance (EMD)
  between the spatial positions of fixations.

## Details

The function computes six different similarity metrics between the
scanpaths `x` and `y`:

- `mm_vector`: Similarity based on the 2D vectors between fixations.

- `mm_direction`: Similarity based on the direction (angle) of saccades
  between fixations.

- `mm_length`: Similarity based on the length of saccades between
  fixations.

- `mm_position`: Similarity based on the spatial position of fixations.

- `mm_duration`: Similarity based on the duration of fixations.

- `mm_position_emd`: Order-insensitive similarity based on the Earth
  Mover's Distance (EMD) between the spatial positions of fixations.

The function ensures that both scanpaths have strictly increasing onset
times and contain at least three fixations. It also normalizes the
similarity scores to lie between 0 and 1, with higher values indicating
greater similarity.

## References

Dewhurst, R., Nyström, M., Jarodzka, H., Foulsham, T., Johansson, R., &
Holmqvist, K. (2012). It depends on how you look at it: Scanpath
comparison in multiple dimensions with MultiMatch, a vector-based
approach. Behavior research methods, 44, 1079-1100.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
scanpath1 <- data.frame(x = runif(10, 0, 500), y = runif(10, 0, 500), onset = cumsum(runif(10, 1, 5)))
scanpath2 <- data.frame(x = runif(10, 0, 500), y = runif(10, 0, 500), onset = cumsum(runif(10, 1, 5)))
screensize <- c(500, 500)
similarity_scores <- multi_match(scanpath1, scanpath2, screensize)
print(similarity_scores)
} # }
```
