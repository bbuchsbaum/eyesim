# Declare screen geometry for GazeWeave

Declare screen geometry for GazeWeave

## Usage

``` r
gaze_screen(width, height, unit = c("px", "deg", "normalized"), center = NULL)
```

## Arguments

- width, height:

  Positive screen dimensions.

- unit:

  Coordinate unit. One of \`"px"\`, \`"deg"\`, or \`"normalized"\`.

- center:

  Optional two-element screen centre. Defaults to half the screen
  dimensions.

## Value

A \`gaze_screen\` specification.
