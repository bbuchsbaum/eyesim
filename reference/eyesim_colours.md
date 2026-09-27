# eyesim plot colours

Named colour roles used by every eyesim plot. Reference (encoding) and
source (recall) gaze keep the same colours in every panel; time, density
and correspondence each have one ramp or colour.

## Usage

``` r
eyesim_colours(role = NULL)
```

## Arguments

- role:

  Optional character vector of role names. \`NULL\` returns all.

## Value

A named character vector of hex colours.

## See also

Other visualization:
[`anim_scanpath()`](https://bbuchsbaum.github.io/eyesim/reference/anim_scanpath.md),
[`element_text_wrap()`](https://bbuchsbaum.github.io/eyesim/reference/element_text_wrap.md),
[`plot.eye_density()`](https://bbuchsbaum.github.io/eyesim/reference/plot.eye_density.md),
[`plot.fixation_group()`](https://bbuchsbaum.github.io/eyesim/reference/plot.fixation_group.md),
[`scale_fill_eyesim_density()`](https://bbuchsbaum.github.io/eyesim/reference/scale_fill_eyesim_density.md),
[`theme_eyesim()`](https://bbuchsbaum.github.io/eyesim/reference/theme_eyesim.md)

## Examples

``` r
eyesim_colours()
#>      reference         source            raw correspondence     diagnostic 
#>      "#1B2A41"      "#C23B4B"      "#8A8F98"      "#5B4FCF"      "#3F4A5A" 
#>            ink          muted           rule           grid 
#>      "#1F2328"      "#5A6270"      "#D5D9DE"      "#ECEEF1" 
eyesim_colours(c("reference", "source"))
#> reference    source 
#> "#1B2A41" "#C23B4B" 
```
