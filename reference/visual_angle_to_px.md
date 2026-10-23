# Convert between visual angle and screen pixels

Converts a visual angle to the on-screen extent in pixels, and back, for
a display viewed at a given distance. The extent is centred on the line
of sight: an angle \\\theta\\ spans \\2 d \tan(\theta / 2)\\ on the
screen, where \\d\\ is the viewing distance. Use it to set distance
thresholds such as the `radius` of
[`fixation_revisits`](https://bbuchsbaum.github.io/eyesim/reference/fixation_revisits.md)
in degrees.

## Usage

``` r
visual_angle_to_px(deg, distance, screen_size, screen_res)

px_to_visual_angle(px, distance, screen_size, screen_res)
```

## Arguments

- deg:

  Visual angle in degrees (vectorised).

- distance:

  Eye-to-screen distance.

- screen_size:

  Physical width (or height) of the screen, in the same units as
  `distance`.

- screen_res:

  Screen resolution in pixels along the same dimension as `screen_size`.

- px:

  Extent in pixels (vectorised).

## Value

`visual_angle_to_px()` returns pixels; `px_to_visual_angle()` returns
degrees.

## Examples

``` r
# 2 degrees on a 52 cm wide, 1920 px display viewed from 60 cm
visual_angle_to_px(2, distance = 60, screen_size = 52, screen_res = 1920)
#> [1] 77.33936
px_to_visual_angle(77, distance = 60, screen_size = 52, screen_res = 1920)
#> [1] 1.991226
```
