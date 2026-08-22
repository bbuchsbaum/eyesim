# Transport v2 estimator-lock fixture

These public deterministic paths and Transport-v2 outputs were frozen
before Transport-v3 algorithm edits. Later backends regenerate the tables
and compare scientific values within 1e-8. The coupling is an optimized
correspondence, not a posterior probability.

Regenerate from the repository root with:

```sh
Rscript inst/validation/gaze-weave-transport-v3-golden.R
```
