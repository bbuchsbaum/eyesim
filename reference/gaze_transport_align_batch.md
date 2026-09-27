# Batch candidate alignments with Transport

Batch candidate alignments with Transport

## Usage

``` r
gaze_transport_align_batch(
  references,
  source,
  spec,
  batch_size = length(references)
)
```

## Arguments

- references:

  Named list of candidate fixation paths or prepared paths.

- source:

  One source fixation path or prepared path.

- spec:

  A \[gaze_transport_spec()\].

- batch_size:

  Positive number of candidates scheduled per batch.

## Value

Named candidate \`gaze_engine_result\` objects in input order.
