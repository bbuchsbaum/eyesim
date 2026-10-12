# Align recall gaze to an encoding path with a fitted Replay model

Align recall gaze to an encoding path with a fitted Replay model

## Usage

``` r
gaze_replay_align(
  reference,
  source,
  model,
  candidate_key = "candidate",
  warp_group = NULL,
  background_key = NULL
)
```

## Arguments

- reference:

  Encoding fixation path.

- source:

  Recall fixation path.

- model:

  A training-only \[fit_gaze_replay_model()\] result.

- candidate_key:

  Candidate identifier retained in the result.

- warp_group:

  Required fitted group key for grouped models.

- background_key:

  Optional values of the \`background_by\` columns (in their declared
  order) for the recall path; models fitted with \`background_by\` only.
  \`NULL\` uses the pooled training background density.

## Value

A \`gaze_engine_result\` and a \`gaze_replay_alignment\` posterior
object. The score is the total trial log likelihood over recall
fixations.

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
# \donttest{
model <- fit_gaze_replay_model(study, recall,
                               match_on = c("participant", "image"),
                               contrast_on = "participant")
alignment <- gaze_replay_align(study$fixgroup[[1]], recall$fixgroup[[1]],
                               model)
alignment$log_score
#> [1] -21.66163
# }
```
