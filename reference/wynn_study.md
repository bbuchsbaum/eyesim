# Study-phase fixations from Wynn et al.

Fixations recorded while participants studied photographs in a
recognition memory experiment, one row per participant and image.

## Usage

``` r
wynn_study
```

## Format

A grouped tibble (grouped by `ImageNumber` and `Subject`) with 12,240
rows and 8 columns:

- Subject:

  Participant number.

- fixgroup:

  A
  [`fixation_group`](https://bbuchsbaum.github.io/eyesim/reference/fixation_group.md)
  with the pooled study fixations (`x`, `y` in pixels; `duration` and
  `onset` in ms).

- ImageVersion:

  Image and version shown, e.g. `"1_A"`.

- Version:

  Image version (`"A"` or `"B"`).

- Trial:

  Trial index of the first viewing.

- ImageNumber:

  Image identifier.

- blockstudy:

  Study block of the first viewing.

- ImageNumberS:

  Participant-image key, `"<Subject>_<ImageNumber>"`, matching
  `wynn_test$ImageNumberS`.

## Source

Wynn et al. eye-tracking study; built from the raw fixation report by
`data-raw/process_dat.R` in the source repository.

## Details

Each image was studied once in each of four study blocks, and a row
pools the fixations of all four viewings: within a row's `fixgroup`,
onsets restart at each viewing. The pooled fixations suit density maps
(for example study templates for
[`template_similarity`](https://bbuchsbaum.github.io/eyesim/reference/template_similarity.md));
for scanpath measures, which join consecutive rows, split the viewings
first. Fixations outside the 800 x 600 pixel image area were removed and
coordinates are relative to its top-left corner.

## See also

[`wynn_test`](https://bbuchsbaum.github.io/eyesim/reference/wynn_test.md)

## Examples

``` r
data(wynn_study)
head(wynn_study)
#> # A tibble: 6 × 8
#> # Groups:   ImageNumber, Subject [6]
#>   Subject fixgroup            ImageVersion Version Trial ImageNumber blockstudy
#>     <dbl> <list>              <chr>        <chr>   <int> <fct>            <int>
#> 1       5 <fxtn_grp [39 × 6]> 1_A          A          18 1                    1
#> 2       6 <fxtn_grp [42 × 6]> 1_A          A          13 1                    1
#> 3       7 <fxtn_grp [26 × 6]> 1_A          A          25 1                    1
#> 4       8 <fxtn_grp [30 × 6]> 1_A          A          11 1                    1
#> 5       9 <fxtn_grp [32 × 6]> 1_A          A          15 1                    1
#> 6      10 <fxtn_grp [41 × 6]> 1_A          A          12 1                    1
#> # ℹ 1 more variable: ImageNumberS <chr>
```
