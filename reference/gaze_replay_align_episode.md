# Align recall gaze to a mixture of encoding presentations

Each encoding path is aligned independently. Their scores (total trial
log likelihoods) are combined as an equal-weight mixture, so no
transition is created between study presentations. The returned
alignment retains all component posteriors and uses the
highest-responsibility component for the ordinary overlay and braid
fields.

## Usage

``` r
gaze_replay_align_episode(
  references,
  source,
  model,
  candidate_key = "candidate",
  template_key = NULL,
  warp_group = NULL,
  background_key = NULL
)
```

## Arguments

- references:

  A list of encoding fixation paths for one candidate item.

- source:

  One recall fixation path.

- model:

  A training-only \[fit_gaze_replay_model()\] result.

- candidate_key:

  Candidate identifier retained in the result.

- template_key:

  Optional unique identifiers for \`references\`.

- warp_group:

  Required fitted group key for grouped models.

- background_key:

  Optional values of the \`background_by\` columns (in their declared
  order) for the recall path; models fitted with \`background_by\` only.
  \`NULL\` uses the pooled training background density.

## Value

A \`gaze_engine_result\` whose alignment inherits from
\`gaze_replay_episode_alignment\` and \`gaze_replay_alignment\`.

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
# Two study presentations of the same image as separate templates.
result <- gaze_replay_align_episode(
  list(first = study$fixgroup[[1]], second = study$fixgroup[[3]]),
  recall$fixgroup[[1]], model
)
result$log_score
#> [1] -22.35478
# }
```
