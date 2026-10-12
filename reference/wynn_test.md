# Test-phase fixations from Wynn et al.

Fixations recorded during the recognition test of the same experiment as
[`wynn_study`](https://bbuchsbaum.github.io/eyesim/reference/wynn_study.md),
one row per participant and test image (a single viewing). Each probe is
either the studied version of an image (`"old"`) or the other version
(`"lure"`).

## Usage

``` r
wynn_test
```

## Format

A grouped tibble (grouped by `ImageNumber` and `Subject`) with 12,240
rows and 11 columns:

- Subject:

  Participant number.

- fixgroup:

  A
  [`fixation_group`](https://bbuchsbaum.github.io/eyesim/reference/fixation_group.md)
  with the test fixations (`x`, `y` in pixels; `duration` and `onset` in
  ms).

- ImageVersion:

  Image and version shown at test.

- Saliency:

  Saliency condition of the test image.

- Accuracy:

  Response accuracy (1 correct, 0 incorrect).

- Version:

  Image version (`"A"` or `"B"`).

- Trial:

  Test trial index.

- Duration:

  Duration condition of the test trial in ms (250, 500 or 750).

- ImageNumber:

  Image identifier.

- ImageNumberS:

  Participant-image key, matching `wynn_study$ImageNumberS`.

- probe_type:

  `"old"` if the test image is the studied version, `"lure"` otherwise.

## Source

Wynn et al. eye-tracking study; built from the raw fixation report by
`data-raw/process_dat.R` in the source repository.

## See also

[`wynn_study`](https://bbuchsbaum.github.io/eyesim/reference/wynn_study.md)

## Examples

``` r
data(wynn_test)
table(wynn_test$probe_type)
#> 
#> lure  old 
#> 6120 6120 
```
