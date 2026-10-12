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

## Examples

``` r
# Example data: two participants study and later recall six images; each
# recall roughly repeats the study path of the same image.
set.seed(1)
make_path <- function(anchor, jitter) {
  fixation_group(
    x = anchor[1] + c(0, 2, 4, 3) + rnorm(4, sd = jitter),
    y = anchor[2] + c(0, 1, -1, 2) + rnorm(4, sd = jitter),
    onset = c(0, 250, 500, 750), duration = c(200, 220, 180, 240)
  )
}
study <- expand.grid(participant = c("p1", "p2"), image = 1:6,
                     stringsAsFactors = FALSE)
anchors <- lapply(study$image, function(i) c(3 * i, 2 * (i %% 3)))
study$fixgroup <- lapply(anchors, make_path, jitter = 0.2)
recall <- study[c("participant", "image")]
recall$fixgroup <- lapply(anchors, make_path, jitter = 0.6)
# Align one recall to the study path of the same image.
alignment <- gaze_transport_align(study$fixgroup[[1]], recall$fixgroup[[1]],
                                  gaze_transport_spec(coverage_nodes = 4))
alignment$log_score
#> [1] -0.07317749
```
