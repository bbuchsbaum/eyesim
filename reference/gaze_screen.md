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

  Optional two-element screen centre, in data coordinates. Defaults to
  half the screen dimensions, so the screen spans \`\[0, width\] x \[0,
  height\]\`; with a centre such as \`c(0, 0)\` (for example degrees of
  visual angle measured from fixation) it spans \`center +/- dimensions
  / 2\`.

## Value

A \`gaze_screen\` specification.

## Examples

``` r
gaze_screen(1024, 768)                                # pixels, origin at a corner
#> $width
#> [1] 1024
#> 
#> $height
#> [1] 768
#> 
#> $unit
#> [1] "px"
#> 
#> $center
#> [1] 512 384
#> 
#> $xlim
#> [1]    0 1024
#> 
#> $ylim
#> [1]   0 768
#> 
#> attr(,"class")
#> [1] "gaze_screen" "list"       
gaze_screen(30, 20, unit = "deg", center = c(0, 0))   # degrees from centre
#> $width
#> [1] 30
#> 
#> $height
#> [1] 20
#> 
#> $unit
#> [1] "deg"
#> 
#> $center
#> [1] 0 0
#> 
#> $xlim
#> [1] -15  15
#> 
#> $ylim
#> [1] -10  10
#> 
#> attr(,"class")
#> [1] "gaze_screen" "list"       
```
