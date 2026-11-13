# Construct leave-out reference maps for two crossed factors

Stores source maps and grouped sums once. References are evaluated on
demand, without storing four full maps for every subject/item pair.

## Usage

``` r
reference_maps(
  dens_tab,
  subject = "subject",
  item = "item",
  refs = c("self_item", "self_other_items", "others_item", "others_other_items"),
  space = c("log", "density"),
  floor = NULL,
  min_n = 1,
  density_var = "density",
  normalize = c("sum", "none"),
  transform = c("none", "log")
)
```

## Arguments

- dens_tab:

  Table with one map per subject/item pair. NULL maps are missing.

- subject, item:

  Names of the two crossed-factor columns.

- refs:

  Requested reference names (see Details).

- space:

  Average log values or density values.

- floor:

  Explicit positive additive offset required when `space = "log"` (the
  default) or `transform = "log"`. For unit-sum maps, a starting choice
  is `1 / (10 * ncells)`, where `ncells` is the number of grid cells.
  There is no automatic offset because this choice changes the reference
  scores.

- min_n:

  Minimum contributing units per reference.

- density_var:

  Column containing density objects.

- normalize:

  Normalize each source map to unit sum, or retain its values.

- transform:

  Optional log transformation after density-space aggregation. Log-space
  aggregation already returns log scores and requires `"none"`.

## Value

An `eyesim_references` object with keys, counts, source values, grouped
sums and a specification. Use
[`reference_map()`](https://bbuchsbaum.github.io/eyesim/reference/reference_map.md)
to inspect one map, or
[`sample_density()`](https://bbuchsbaum.github.io/eyesim/reference/sample_density.md)
to score a table of locations carrying subject/item keys.

## Details

`self_item` uses the target map; `self_other_items` averages that
subject's other maps; `others_item` averages other subjects' maps of
that item. `others_other_items` averages other items within each other
subject, then averages those subject means equally, including in
incomplete designs. `min_n` counts maps for self references and subjects
for other references. Insufficient references retain their keys but
return missing scores. Normalization precedes log transformation, which
precedes aggregation.
