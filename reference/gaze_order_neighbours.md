# Declare fixed-neighbour ordinal chronology for Transport

Declare fixed-neighbour ordinal chronology for Transport

## Usage

``` r
gaze_order_neighbours(
  neighbours = 2L,
  coalesce_distance = 0,
  unit = c("native", "px", "deg", "normalized")
)
```

## Arguments

- neighbours:

  Number of forward fixation neighbours represented in the directed
  relation graph.

- coalesce_distance:

  Maximum distance for merging immediately adjacent fixation atoms. The
  default zero merges only exactly identical locations.

- unit:

  Coordinate unit for \`coalesce_distance\`.

## Value

A \`gaze_chronology_spec\` object.

## Examples

``` r
gaze_order_neighbours(neighbours = 2)
#> $type
#> [1] "order_neighbours"
#> 
#> $clock
#> [1] "order"
#> 
#> $neighbours
#> [1] 2
#> 
#> $coalesce_distance
#> [1] 0
#> 
#> $unit
#> [1] "native"
#> 
#> attr(,"class")
#> [1] "gaze_chronology_spec" "list"                
```
