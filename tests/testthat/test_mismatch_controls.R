test_that("mismatch exclusions use donor items across repeated trials and case subsets", {
  refs <- reference_maps(reference_fixture(), floor = .01)
  fix <- reference_fixations()
  fix <- rbind(fix, fix)
  fix$trial <- paste0("episode-", seq_len(nrow(fix)))
  # All five eligible donor fixations are retained in the assignment table.
  strata <- sample_controls(fix, n = 5, draws = 2)
  selected <- strata[strata$item == "a", ]
  assignments <- mismatch_assignments(selected, refs, n = 3)
  result <- mismatch_draws(selected, refs, assignments = assignments, standardize = "none")
  expect_true(all(result$contrasts$n_controls == 3))
  expect_equal(result$paired_observed$n_controls, result$contrasts$n_controls)
  for (d in 1:3) {
    a <- assignments[assignments$null_draw == d, ]
    # Independent explicit subject/trial matching; trials are not item labels.
    replacement <- vapply(seq_len(nrow(selected)), function(i) {
      a$replacement_item[a$subject == selected$subject[i] & a$trial == selected$trial[i]]
    }, character(1))
    keep <- selected$is_case | fix$item[selected$control_id] != replacement
    alt <- selected[keep, ]; alt$item <- replacement[keep]
    expected <- case_control_contrast(score_strata(alt, refs, standardize = "none"))
    expect_equal(result$contrasts$self_item[result$contrasts$null_draw == d], expected$self_item)
  }
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(selected, path)
  expect_identical(mismatch_draws(readRDS(path), refs, assignments = assignments, standardize = "none"), result)
  # Complete old assignment tables can recover provenance from their case rows.
  old <- strata; attr(old, "control_sources") <- NULL
  expect_equal(mismatch_draws(old, refs, n = 2), mismatch_draws(strata, refs, n = 2))
  incomplete <- selected; attr(incomplete, "control_sources") <- NULL
  expect_error(mismatch_draws(incomplete, refs, assignments = assignments), "Missing donor provenance")
  expect_silent(mismatch_draws(incomplete, refs, assignments = assignments, control_policy = "fixed"))
})

test_that("paired comparisons share retained rows and a single scaling fit", {
  refs <- reference_maps(reference_fixture(), floor = .01)
  fix <- reference_fixations()
  strata <- sample_controls(fix, n = 2, draws = 2)
  assignments <- mismatch_assignments(strata, refs, n = 2)
  for (replace_ref in c("self_item", "all")) {
    result <- mismatch_draws(strata, refs, assignments = assignments, replace_ref = replace_ref)
    for (d in 1:2) {
      a <- assignments[assignments$null_draw == d, ]
      replacement <- vapply(seq_len(nrow(strata)), function(i) {
        a$replacement_item[a$subject == strata$subject[i] & a$trial == strata$trial[i]]
      }, character(1))
      keep <- strata$is_case | fix$item[strata$control_id] != replacement
      matched <- score_strata(strata[keep, ], refs, standardize = "none")
      alt <- strata[keep, ]; alt$item <- replacement[keep]
      mismatched <- score_strata(alt, refs, standardize = "none")
      if (replace_ref == "self_item") {
        rest <- setdiff(refs$refs, "self_item")
        mismatched[rest] <- matched[rest]
      }
      fit <- fit_score_scaling(mismatched)
      for (resolution in c("fixation", "trial")) {
        actual <- if (resolution == "fixation") result else
          mismatch_draws(strata, refs, assignments = assignments, replace_ref = replace_ref, by = "trial")
        expected <- case_control_contrast(apply_score_scaling(matched, fit), by = resolution)
        other <- case_control_contrast(apply_score_scaling(mismatched, fit), by = resolution)
        for (nm in refs$refs) {
          expect_equal(actual$paired_observed[[nm]][actual$paired_observed$null_draw == d], expected[[nm]], tolerance = 1e-12)
          expect_equal(actual$contrasts[[nm]][actual$contrasts$null_draw == d], other[[nm]], tolerance = 1e-12)
        }
      }
    }
  }
})

test_that("exclusion never redraws empty strata or removes generated controls", {
  refs <- reference_maps(reference_fixture(), floor = .01)
  fix <- reference_fixations()
  fix <- fix[fix$item != "c", ]
  strata <- sample_controls(fix, n = 1)
  a <- mismatch_assignments(strata, refs, n = 1)
  a$replacement_item <- ifelse(a$item == "a", "b", "a")
  set.seed(123); before <- .Random.seed
  expect_error(mismatch_draws(strata, refs, assignments = a), "without controls")
  expect_identical(.Random.seed, before)
  for (method in c("uniform", "map")) {
    generated <- sample_controls(fix, n = 2, method = method,
                                 xbounds = c(0,2), ybounds = c(0,1),
                                 control_map = reference_fixture()$density[[1]])
    excluded <- mismatch_draws(generated, refs, assignments = a)
    fixed <- mismatch_draws(generated, refs, assignments = a, control_policy = "fixed")
    expect_identical(excluded$contrasts, fixed$contrasts)
    expect_identical(excluded$paired_observed, fixed$paired_observed)
  }
  broken <- strata; broken$control_id <- NULL
  expect_error(mismatch_draws(broken, refs, assignments = a), "donor provenance")
  broken <- strata; broken$control_id[2] <- NA_integer_
  expect_error(mismatch_draws(broken, refs, assignments = a), "donor provenance")
  broken <- strata; broken$control_id[2] <- .5
  expect_error(mismatch_draws(broken, refs, assignments = a), "index saved donor rows")
})

test_that("paired mismatch reductions preserve keys after shuffling", {
  refs <- reference_maps(reference_fixture(), floor = .01)
  strata <- sample_controls(reference_fixations(), n = 2, draws = 2)
  assignments <- mismatch_assignments(strata, refs, n = 2)
  for (by in c("fixation", "trial")) {
    result <- mismatch_draws(strata, refs, assignments = assignments, by = by)
    shuffled <- mismatch_draws(strata[rev(seq_len(nrow(strata))), ], refs,
                               assignments = assignments, by = by)
    for (part in c("observed", "paired_observed", "contrasts")) {
      expect_equal(shuffled[[part]], result[[part]], tolerance = 1e-12)
    }
  }
})

test_that("reduction preserves standardized contrasts of nearly constant maps", {
  tab <- reference_fixture()
  for (i in seq_along(tab$density)) {
    tab$density[[i]]$z <- matrix(1 + seq_len(6) * 1e-14, 3, 2)
  }
  refs <- reference_maps(tab, floor = .01, refs = "self_item")
  fix <- reference_fixations()
  strata <- sample_controls(fix, n = 2)
  a <- mismatch_assignments(strata, refs, n = 1)
  for (policy in c("fixed", "exclude_replacement")) {
    result <- mismatch_draws(strata, refs, assignments = a, control_policy = policy)
    replacement <- vapply(seq_len(nrow(strata)), function(i) {
      a$replacement_item[a$subject == strata$subject[i] & a$trial == strata$trial[i]]
    }, character(1))
    keep <- policy == "fixed" | strata$is_case | fix$item[strata$control_id] != replacement
    alt <- strata[keep, ]; alt$item <- replacement[keep]
    # Row-level centering/scaling before reduction is the independent order.
    expected <- case_control_contrast(score_strata(alt, refs))
    expect_equal(result$contrasts$self_item, expected$self_item, tolerance = 1e-12)
    expect_equal(result$observed$self_item,
                 case_control_contrast(score_strata(strata, refs))$self_item, tolerance = 1e-12)
  }
})

test_that("paired mismatch scoring preserves unavailable references", {
  tab <- reference_fixture()
  refs <- expect_warning(reference_maps(tab[tab$subject == "s1", ], floor = .01), "min_n")
  fix <- reference_fixations()
  strata <- sample_controls(fix[fix$subject == "s1", ], n = 2)
  result <- mismatch_draws(strata, refs, n = 2, replace_ref = "all")
  for (part in c("observed", "paired_observed", "contrasts")) {
    expect_true(all(is.finite(result[[part]]$self_item)))
    expect_true(all(is.na(result[[part]]$others_item)))
    expect_true(all(is.na(result[[part]]$others_other_items)))
  }
})
