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
