# Inspect one reference map

Inspect one reference map

## Usage

``` r
reference_map(refs, row = 1L, ref = refs$refs[1])
```

## Arguments

- refs:

  A reference set from
  [`reference_maps()`](https://bbuchsbaum.github.io/eyesim/reference/reference_maps.md).

- row:

  Row of `refs$keys` identifying the target.

- ref:

  A requested reference name.

## Value

A `eyesim_reference_map` with x/y/z, domain and contributor count, or
NULL when insufficient contributors exist. Log maps contain scores, not
probabilities, and must not be transformed again.
