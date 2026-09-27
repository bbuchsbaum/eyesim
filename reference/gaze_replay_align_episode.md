# Align recall gaze to a mixture of encoding presentations

Each encoding path is aligned independently. Their scores (total log
likelihoods under revision \`"2026.10"\`, mean log likelihoods per
duration bin under \`"2026.08"\`) are combined as an equal-weight
mixture, so no transition is created between study presentations. The
returned alignment retains all component posteriors and uses the
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
  order) for the recall path; revision \`"2026.10"\` models fitted with
  \`background_by\` only. \`NULL\` uses the pooled training background
  density.

## Value

A \`gaze_engine_result\` whose alignment inherits from
\`gaze_replay_episode_alignment\` and \`gaze_replay_alignment\`.
