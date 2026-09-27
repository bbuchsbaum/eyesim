# Align gaze paths with edge-normalized Transport

Align gaze paths with edge-normalized Transport

## Usage

``` r
gaze_transport_align(
  reference,
  source,
  spec,
  warp_model = NULL,
  candidate_key = "candidate"
)
```

## Arguments

- reference, source:

  Fixation paths or prepared gaze measures.

- spec:

  A \[gaze_transport_spec()\].

- warp_model:

  Candidate-invariant fitted warp. A non-identity warp must be learned
  outside the scored pair.

- candidate_key:

  Candidate identifier.

## Value

A symmetric \`gaze_engine_result\` with a coverage-integrated Transport
alignment.
