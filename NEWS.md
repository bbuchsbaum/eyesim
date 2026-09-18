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
  Weights must be finite, non-negative, not all zero, and supplied one per
  fixation; they are subset along with `window`, and they take precedence over
  `duration_weighted`.

* `eye_density()` now forwards extra named arguments in `...` to `ks::kde()`,
  as documented. Previously any extra argument failed with "unused argument".
  Arguments that `eye_density()` sets itself (including `eval.points`), and
  unknown names, are rejected with a clear error, as are extra arguments under
  `kde_pkg = "MASS"`.

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

* Permutation candidates are now distinct templates. A template matched by
  several source rows previously appeared once per matching row, so it could
  be drawn twice and was counted twice in `perm_sim` and `n_perm`. This
  applies to `template_similarity()`, `template_similarity_cv()`,
  `fixation_similarity()`, `scanpath_similarity()`, and
  `sample_density_time()`, and changes their permutation columns whenever
  keys repeat within a stratum. Reference rows that no source row matches
  remain outside the candidate set, as before; the documentation now says so.

* Documentation: `?template_similarity` no longer claims that permutation
  sampling uses a "fixed future seed". Sampling uses the session RNG, so
  `set.seed()` before the call makes the baseline reproducible. No behaviour
  change.

* `template_similarity_cv()` no longer resets the caller's random number
  stream. Fold assignment and permutation draws still run under `seed`, so its
  results are unchanged, but the session RNG state is restored on exit. The
  documentation now states that permutation controls come only from the
  held-out fold.

* `similarity(method = "fisherz")` now returns exactly
  `atanh(1 - .Machine$double.eps)` (about 18.37) for every pair of identical
  maps, and for any correlation within `64 * .Machine$double.eps` of 1, such as
  a rescaled copy of a map; correlations that close to -1 give its negative.
  Identical constant maps previously returned 1, and identical non-constant
  maps returned values between about 17.3 and 18.37, depending on rounding
  error in `cor()`. The clamp is documented in `?template_similarity`.

* `similarity()` on two density maps now refuses maps whose lattices (x and y
  grid coordinates) differ, instead of comparing the `z` matrices cell by cell
  as if they were aligned. A numeric `y` must have one value per grid cell.
  `template_similarity()` and related wrappers raise the same error on both
  the fast cosine and the general path, and for multiscale maps, where the
  per-scale error handling previously would have turned it into `NA`.

* `similarity()` for fixation groups with `method = "overlap"`, and therefore
  `fixation_similarity(method = "overlap")`, now defaults to `dthresh = 60`,
  the value documented for `fixation_overlap()`. It previously defaulted to 40,
  so the same pair scored differently through the two entry points.

* `fixation_overlap()` now builds its default `time_samples` grid from the
  onsets of both fixation groups, `seq(0, max(c(x$onset, y$onset)), by = 20)`.
  It previously used `x` only, so swapping the arguments could change the
  result substantially (0.098 versus 0.833 in one case).

* `sample_fixations()` with the default `fast = TRUE` now holds the last
  fixation for time points after its onset, as `fast = FALSE` always did. It
  previously returned `NA` there. It also no longer fails for a group with a
  single fixation. Downstream, `sample_density()` with `times`,
  `sample_density_time()`, `template_sample()` with `time`, and
  `fixation_overlap()` now score time points after the last onset instead of
  treating them as missing. Two smaller changes on the fast path: fixations
  with tied onsets previously had their coordinates averaged
  (`approx(ties = mean)`) and now yield the last tied fixation, as
  `fast = FALSE` does; and a fixation with `NA` coordinates was previously
  skipped, carrying the preceding fixation forward, and now yields `NA`. The
  `fast = FALSE` path now orders fixations by onset first, as the fast path
  does; for groups whose onsets were not in increasing order it previously
  returned coordinates by row position rather than by onset.

* `rep_fixations()` now counts `floor(duration * resolution)` replicates with a
  floating-point tolerance. A duration of 0.29 at resolution 100 previously
  gave 28 copies because `0.29 / 0.01` evaluates to 28.999999999999996; it now
  gives 29.

* `sample_density_time()` now treats every `time_bins` interval as half-open,
  `[lower, upper)`, including the last, as documented. The final bin was
  previously closed, so a time point equal to the last break (for example
  t = 3000 with the default `times` and bins ending at 3000) was averaged
  into the last bin. It now falls outside all bins, which changes the last
  `bin_*`, `perm_bin_*`, and `diff_bin_*` values whenever a sample lies on the
  final break.

* `multi_match()` now computes `mm_position_emd` over all fixations. It
  previously used the saccade table, which omits the final fixation, so two
  three-fixation paths differing only in their last fixation scored 1. Every
  `mm_position_emd` value changes, including the `_perm` and `_diff` columns
  from `scanpath_similarity()`.

* `template_multireg()` now defaults to `method = "lm"`, as documented.
  Omitting `method` previously failed with "the condition has length > 1".
  Unknown methods are rejected by `match.arg()`.

* `template_regression()` now stops with a clear message when a
  `baseline_key` value used by `source_tab` appears in more than one row of
  `baseline_tab`. It previously failed with "$ operator is invalid for atomic
  vectors".

* `template_sample()` now drops rows whose template or fixation group is
  `NULL`, with a warning, as its code comment intended. The `is.null()` filter
  on list columns was a no-op, so a single `NULL` template made the whole call
  fail.

* `fixation_entropy()` on a fixation group with a single fixation now returns
  `NA` for `method = "density"` when `sigma` is not supplied, because
  `suggest_sigma()` cannot estimate a bandwidth from one point. It previously
  failed with "missing values present in assertion". `eye_density()` now
  reports an `NA` or non-finite `sigma` with a clear message.

* `fixation_entropy()` now rejects maps with negative values, such as the
  difference of two densities, with a clear error. Previously an exact
  difference map (total zero) gave `NA` while a signed map with a positive
  total gave a meaningless number. Maps with zero total mass still give `NA`.

* Documentation only, no behaviour change:
  - `?eye_density.fixation_group` now states that `sigma` is the kernel
    standard deviation under `kde_pkg = "ks"` but the `MASS::kde2d()`
    bandwidth under `"MASS"`, whose kernel standard deviation is `sigma / 4`.
  - `?sample_density` now states that points are looked up at the nearest
    lattice point with `round()`, so exact midpoints resolve half to even on
    the index scale, and that off-lattice points are clamped to the edge.
  - `?fixation_entropy.default` now states that the grid method counts
    fixations outside `xbounds`/`ybounds` in the nearest edge cell.
