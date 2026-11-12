test_that("scoring, fitted scaling and contrasts match independent calculations", {
  tab <- reference_fixture()
  refs <- reference_maps(tab, floor = .01)
  strata <- sample_controls(reference_fixations(), n = 2, draws = 2)
  set.seed(99)
  before <- .Random.seed
  raw <- score_strata(strata, refs, standardize = "none")
  expect_identical(.Random.seed, before)
  for (nm in refs$refs) {
    expected <- vapply(seq_len(nrow(strata)), function(i) {
      row <- which(tab$subject == strata$subject[i] & tab$item == strata$item[i])
      map <- reference_oracle(tab, row, nm)
      map[strata$x[i]+1, strata$y[i]+1]
    }, numeric(1))
    expect_equal(raw[[nm]], expected, tolerance = 1e-12)
  }
  scaled <- score_strata(strata, refs)
  fit <- fit_score_scaling(raw, refs$refs, c("subject", "draw"))
  expect_identical(scaled, apply_score_scaling(raw, fit))
  expect_identical(score_strata(strata, refs, scaling = fit), scaled)
  expect_error(score_strata(strata, refs, standardize = "none", scaling = fit), "conflicts")
  for (nm in refs$refs) for (s in unique(strata$subject)) for (d in 1:2) {
    rows <- which(strata$subject == s & strata$draw == d)
    expect_equal(scaled[[nm]][rows], as.vector(scale(raw[[nm]][rows])), tolerance = 1e-11)
  }
  contrasts <- case_control_contrast(scaled)
  for (i in seq_len(nrow(contrasts))) for (nm in refs$refs) {
    rows <- which(scaled$stratum == contrasts$stratum[i])
    expected <- scaled[[nm]][rows[1]] - sum(scaled[[nm]][rows[-1]]) / (length(rows)-1)
    expect_equal(contrasts[[nm]][i], expected, tolerance = 1e-12)
  }
  one <- score_strata(sample_controls(reference_fixations(), n=2, draws=1), refs)
  for (nm in refs$refs) expect_equal(one[[nm]], scaled[[nm]][scaled$draw == 1])
  expect_error(apply_score_scaling(scaled, fit), "already standardized")
  new <- raw; new$subject <- "unseen"
  expect_error(apply_score_scaling(new, fit), "no fitted")
  shuffled <- strata[rev(seq_len(nrow(strata))), ]
  replay <- case_control_contrast(score_strata(shuffled, refs))
  for (nm in refs$refs) expect_equal(replay[[nm]], contrasts[[nm]], tolerance = 1e-12)
})

test_that("constant and missing scores have explicit scaling and contrast behavior", {
  tab <- reference_fixture()
  tab$density <- rep(list(gen_density(0:2,0:1,matrix(1,3,2))),9)
  refs <- reference_maps(tab, floor=.01)
  strata <- sample_controls(reference_fixations(), n=2)
  scaled <- score_strata(strata,refs)
  expect_equal(scaled$self_item, rep(0,nrow(scaled)))
  raw <- score_strata(strata,refs,standardize="none")
  raw$self_item <- NA_real_
  fit <- fit_score_scaling(raw)
  expect_true(all(is.na(apply_score_scaling(raw,fit)$self_item)))
  expect_true(all(is.na(case_control_contrast(raw)$self_item)))
  broken <- strata; broken$is_case[2] <- TRUE
  expect_error(score_strata(broken,refs), "exactly one")
  broken <- strata; broken$item[2] <- "different"
  expect_error(score_strata(broken,refs), "case's reference keys")
})

test_that("log-reference defaults explain the required offset", {
  expect_error(reference_maps(reference_fixture()), "explicit positive floor.*1 / \\(10 \\* ncells\\)")
})

test_that("mismatch evaluation reuses assignments and matches explicit rescoring", {
  refs <- reference_maps(reference_fixture(), floor=.01)
  strata <- sample_controls(reference_fixations(), n=2, draws=2)
  assignments <- mismatch_assignments(strata,refs,n=3)
  set.seed(103)
  before <- .Random.seed
  result <- mismatch_draws(strata,refs,assignments=assignments,replace_ref="all",standardize="none")
  expect_identical(.Random.seed,before)
  expect_equal(result, mismatch_draws(strata,refs,assignments=assignments,replace_ref="all",standardize="none",seed=999))
  for (d in 1:3) {
    alt <- strata
    for (i in seq_len(nrow(alt))) {
      a <- assignments[assignments$subject==alt$subject[i] & assignments$trial==alt$trial[i] & assignments$null_draw==d,]
      alt$item[i] <- a$replacement_item
    }
    keep <- alt$is_case | reference_fixations()$item[alt$control_id] != alt$item
    alt <- alt[keep, ]
    manual <- case_control_contrast(score_strata(alt,refs,standardize="none"))
    for (nm in refs$refs) expect_equal(result$contrasts[[nm]][result$contrasts$null_draw==d],manual[[nm]])
    matched <- case_control_contrast(score_strata(strata[keep, ],refs,standardize="none"))
    for (nm in refs$refs) expect_equal(result$paired_observed[[nm]][result$paired_observed$null_draw==d],matched[[nm]])
  }
  self <- mismatch_draws(strata,refs,assignments=assignments,standardize="none")
  for (nm in setdiff(refs$refs,"self_item")) expect_equal(self$contrasts[[nm]],self$paired_observed[[nm]])
  scaled <- mismatch_draws(strata,refs,assignments=assignments,restandardize=FALSE)
  unscaled <- mismatch_draws(strata,refs,assignments=assignments,standardize="none")
  for (s in unique(strata$subject)) for (draw in 1:2) {
    rows <- scaled$contrasts$subject==s & scaled$contrasts$draw==draw
    g <- which(scaled$scaling$groups$subject==s & scaled$scaling$groups$draw==draw)
    expect_equal(scaled$contrasts$self_item[rows],unscaled$contrasts$self_item[rows]/scaled$scaling$scale[g,"self_item"],tolerance=1e-12)
  }
  bad <- assignments; bad$replacement_item <- bad$item
  expect_error(mismatch_draws(strata,refs,assignments=bad), "different items")
})

test_that("known-answer null, planted preference, generic effects and control jitter", {
  tab <- reference_fixture()
  # Nine distinct peaks: subject/item-specific preferences.
  tab$density <- lapply(1:9,function(i) {
    z <- matrix(.01,3,3); z[i] <- 1
    gen_density(0:2,0:2,z)
  })
  refs <- reference_maps(tab,floor=.001)
  # Every trial visits every cell once: exhaustive controls give exactly zero
  # mean preference, including after any coherent mismatch assignment.
  null <- tab[rep(1:9,each=9),c("subject","item")]
  null$trial <- rep(1:9,each=9)
  null$x <- rep(rep(0:2,3),9); null$y <- rep(rep(0:2,each=3),9)
  strata <- sample_controls(null,n=18)
  result <- mismatch_draws(strata,refs,n=3,by="trial",standardize="none")
  expect_equal(result$observed$self_item,rep(0,9),tolerance=1e-12)
  expect_equal(result$contrasts$self_item,rep(0,27),tolerance=1e-12)
  planted <- tab[c("subject","item")]
  planted$trial <- 1:9; planted$x <- rep(0:2,3); planted$y <- rep(0:2,each=3)
  signal <- mismatch_draws(sample_controls(planted,n=2),refs,n=3,standardize="none")
  expect_true(all(signal$observed$self_item > 4))
  expect_equal(signal$contrasts$self_item, rep(0,27), tolerance=1e-12)
  expect_equal(signal$paired_observed$self_item, rep(signal$observed$self_item,3))
  fixed <- mismatch_draws(sample_controls(planted,n=2),refs,assignments=signal$assignments,
                          standardize="none",control_policy="fixed")
  expect_equal(fixed$contrasts$self_item, -rep(fixed$observed$self_item,3)/2)
  # Pure content preference: all subjects have the same map for each item.
  tab$density <- rep(tab$density[c(1,4,7)],each=3)
  generic <- reference_maps(tab,floor=.001)
  score <- score_strata(sample_controls(planted,n=2),generic,standardize="none")
  expect_equal(score$self_item,score$others_item,tolerance=1e-12)
  # Known failing procedure: jittering only controls creates a preference even
  # when cases and observed controls have exactly the same location.
  tab$density <- rep(tab$density[1],9)
  centered <- reference_maps(tab,floor=.001)
  planted$x <- 0; planted$y <- 0
  same <- sample_controls(planted,n=2)
  expect_equal(case_control_contrast(score_strata(same,centered,standardize="none"))$self_item,rep(0,9))
  jittered <- same; jittered$x[!jittered$is_case] <- 1
  expect_true(all(case_control_contrast(score_strata(jittered,centered,standardize="none"))$self_item > 4))
})
