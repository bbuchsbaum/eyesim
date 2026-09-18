# eyesim 0.1.0.9000

* Canonicalized edge-normalized Transport as the sole public Transport method.
  The API is now `gaze_transport_spec()`, `gaze_transport_align()`,
  `gaze_transport_cv()`, and related unversioned helpers; `gaze_weave_cv()`
  accepts `engine = "transport"`. The unreleased original and v2
  implementations, exports, S3 classes, tests, and documentation were removed.
  Their labels remain only in internal and immutable development-court
  provenance.

* Added edge-normalized GazeWeave Transport as the symmetric explanatory
  alignment engine. It integrates bounded coverage, scores separate study
  presentations as fixed equal-prior likelihoods, uses nested response-blind
  calibration, and retains a pure-R oracle alongside the optimized native
  backend. Default tidy output contains only identifiers and `gaze_info_bits`;
  coverage, coordinate-unit spatial RMSE, local order, contraction, rank,
  episode weighting, and solver stability are nested diagnostics.

* Added Transport registered overlays, optimized-correspondence gaze braids,
  raw-versus-registered paths, and a one-score evidence ledger. Labels state
  that Transport correspondence is an optimized alignment, not posterior
  replay probability.

* Froze and ran the Transport simulation, full-cohort repeated-viewing, and
  exact 0--3000 ms retrieval courts. Invariants, runtime, and the study positive
  control passed. Retrieval information and held-out old-item response
  prediction did not pass their uncertainty gates, so Transport is retained
  as an explanatory engine without a retrieval or superiority claim.

* Added GazeWeave Replay, a directional probability model of locally ordered,
  restartable encoding-to-recall replay mixed with candidate-independent
  background gaze. Candidate probability is calibrated on inner-held-out items
  within each outer-training set from mean log predictive density per normalized
  duration bin, with log-domain evidence and a weak fixed residual-temperature
  prior. Raw total HMM likelihood remains available as a diagnostic. The primary endpoint is
  `gaze_info_bits = log2(p_true / prior_true)`.

* Added engine-specific registered overlays, gaze braids, and diagnostics.
  Replay braids show posterior correspondence probabilities. Transport braids
  are labelled as regularized optimized correspondences, not posterior
  uncertainty.

* Added frozen synthetic and real-data validation courts under
  `inst/validation`. The real court compared Replay and the then-current
  Transport development estimator with
  nested-ridge MultiMatch, registered density, and elastic-matching baselines
  under shared participant-by-item folds. The synthetic court passed all gates.
  The real court failed its repeated-viewing gate after calibration exposed
  unstable Replay probabilities, so `gaze_weave_cv()` requires an explicit
  engine choice and no GazeWeave superiority is claimed. Registered density was
  strongest overall and blank-screen imagery evidence was weak for every
  method.

* Added a prospectively frozen participant- and item-disjoint replication of
  the resolution-invariant Replay calibration. The correction restored
  positive held-out Replay information in repeated viewing, but neither
  GazeWeave engine had positive mean information in the blank-screen control.
  The package therefore retains explicit engine selection and makes no power or
  superiority claim.

* `template_similarity()`, `fixation_similarity()`, and `scanpath_similarity()`
  now return an `n_perm` column when `permutations > 0`. It records the number of
  permuted non-matching comparisons that actually contributed to `perm_sim` for
  each row, which varies with small `permute_on` strata or when fewer candidates
  than requested are available, and is `0` when no baseline could be computed
  (`perm_sim = NA`). Use it to exclude rows with too thin a permutation baseline,
  e.g. `dplyr::filter(res, n_perm >= k)`. The column is additive: output for
  `permutations = 0` is unchanged.

## Correctness fixes from the eyes4s parity audit

* `eye_density()` now honours an explicit `weights` vector. It was previously
  ignored, so the map was always unweighted unless `duration_weighted = TRUE`.
  Weights must be finite, non-negative, and supplied one per fixation; they are
  subset along with `window`, and they take precedence over
  `duration_weighted`.

* `eye_density()` now forwards extra named arguments in `...` to `ks::kde()`,
  as documented. Previously any extra argument failed with "unused argument".
  Arguments that `eye_density()` sets itself, and unknown names, are rejected
  with a clear error, as are extra arguments under `kde_pkg = "MASS"`.

* Weighted densities under `kde_pkg = "MASS"` (explicit `weights` or
  `duration_weighted = TRUE`) now work. The internal weighted kernel failed
  with "invalid 'times' argument" and `eye_density()` returned `NULL`.

* `eye_density()` and the geometric warp behind `affine_transform()` and
  `contract_transform()` now round maps with a fixed `zapsmall(digits = 7)`.
  They previously used `getOption("digits")`, so the same call returned maps
  differing by up to about 5e-9 when a session changed `options(digits)`.
  Results under the default `digits = 7` are unchanged.

* The permutation baseline in `template_similarity()`,
  `template_similarity_cv()`, `fixation_similarity()`, and
  `scanpath_similarity()` now removes the true match before sampling
  candidates. It previously sampled first and then dropped the match, so with
  three candidates and `permutations = 2` most rows received one control
  instead of the documented two. `n_perm` is now `min(permutations, available
  non-matching candidates)`. Whenever `permutations` is smaller than the
  candidate pool, the sampled controls, and hence `perm_sim` and
  `eye_sim_diff`, differ from earlier versions for the same seed.

* When several source rows share a `match_on` key, every copy of that key is
  now removed from each such row's permutation candidates. Previously one copy
  remained, so a row could be compared with its own template in the baseline.
  This matches `sample_density_time()`.

* Documentation: `?template_similarity` no longer claims that permutation
  sampling uses a "fixed future seed". Sampling uses the session RNG, so
  `set.seed()` before the call makes the baseline reproducible. No behaviour
  change.

* `template_similarity_cv()` no longer resets the caller's random number
  stream. Fold assignment and permutation draws still run under `seed`, so its
  results are unchanged, but the session RNG state is restored on exit. The
  documentation now states that permutation controls come only from the
  held-out fold.

* `similarity(method = "fisherz")` now returns `atanh(1 - .Machine$double.eps)`
  (about 18.37) for identical constant maps, the same value it already
  returned for identical non-constant maps. It previously returned 1. The
  clamp is documented in `?template_similarity`.

* `similarity()` on two density maps now refuses maps whose lattices (x and y
  grid coordinates) differ, instead of comparing the `z` matrices cell by cell
  as if they were aligned. A numeric `y` must have one value per grid cell.
  `template_similarity()` and related wrappers raise the same error on both
  the fast cosine and the general path.

* `similarity()` for fixation groups with `method = "overlap"`, and therefore
  `fixation_similarity(method = "overlap")`, now defaults to `dthresh = 60`,
  the value documented for `fixation_overlap()`. It previously defaulted to 40,
  so the same pair scored differently through the two entry points.

* `fixation_overlap()` now builds its default `time_samples` grid from the
  onsets of both fixation groups, `seq(0, max(c(x$onset, y$onset)), by = 20)`.
  It previously used `x` only, so swapping the arguments could change the
  result substantially (0.098 versus 0.833 in one case).
