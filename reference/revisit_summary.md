# Summarise revisits for one episode

Summarise revisits for one episode

## Usage

``` r
revisit_summary(x)
```

## Arguments

- x:

  A table returned by
  [`fixation_revisits()`](https://bbuchsbaum.github.io/eyesim/reference/fixation_revisits.md)
  for a `fixation_group`.

## Value

A one-row tibble: `n_fix`; `opportunities`; `revisits`; `revisit_rate`
(revisits per opportunity); `expected_revisits` and `expected_rate`
under the null (`NA` without a reference, or when any opportunity has no
expectation because no reference location satisfies its departure
condition); `excess_revisits` (observed minus expected) and
`excess_rate`.
