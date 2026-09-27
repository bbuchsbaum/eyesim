# Elastic consensus matching for encoding and recall fixations

This is an inspectable reference implementation of the order-free
consensus and moving-least-squares relocation described by Wang,
Holmqvist, and Alexa (2021). It is a comparator, not a GazeWeave engine.
The iterative relocation follows their Gaussian weighting equations.
Their method returns relocated recall and matched encoding fixations
rather than a pairwise scalar; \`score\` is therefore an
eyesim-specific, duration-weighted soft match used only as a frozen
compatibility comparator.

## Usage

``` r
elastic_consensus_align(reference, source, spec)
```

## Arguments

- reference:

  Encoding fixation path.

- source:

  Recall fixation path.

- spec:

  A \[gaze_baseline_spec()\].

## Value

An \`elastic_consensus_alignment\` with relocated points, soft spatial
similarity, matching coverage, deformation, and convergence diagnostics.

## References

Wang, X., Holmqvist, K., & Alexa, M. (2021). The recorded trajectories
of visually guided eye movements can be used to recover the spatial
location of previously attended objects. \*Attention, Perception, &
Psychophysics\*, 83, 1700-1716.
