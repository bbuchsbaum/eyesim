test_that("mismatch exclusions use donor items across repeated trials and case subsets", {
  refs <- reference_maps(reference_fixture(), floor = .01)
  fix <- reference_fixations()
  fix <- rbind(fix, fix)
  fix$trial <- paste0("episode-", seq_len(nrow(fix)))
  # All five eligible donor fixations are retained in the assignment table.
  strata <- sample_controls(fix, n = 5, draws = 2, item = NULL)
  selected <- strata[strata$item == "a", ]
  assignments <- mismatch_assignments(selected, refs, n = 3)
  result <- mismatch_draws(selected, refs, assignments = assignments, standardize = "none")
  expect_true(all(result$contrasts$n_controls == 2))
  expect_equal(result$paired_observed$n_controls, result$contrasts$n_controls)
  for (d in 1:3) {
    a <- assignments[assignments$null_draw == d, ]
    # Independent explicit subject/trial matching; trials are not item labels.
    replacement <- vapply(seq_len(nrow(selected)), function(i) {
      a$replacement_item[a$subject == selected$subject[i] & a$trial == selected$trial[i]]
    }, character(1))
    donor_item <- fix$item[selected$control_id]
    keep <- selected$is_case | (donor_item != replacement & donor_item != selected$item)
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
  fix <- rbind(fix, fix)
  fix$trial <- seq_len(nrow(fix))
  strata <- sample_controls(fix, n = 5, draws = 2, item = NULL)
  assignments <- mismatch_assignments(strata, refs, n = 2)
  for (replace_ref in c("self_item", "all")) {
    result <- mismatch_draws(strata, refs, assignments = assignments, replace_ref = replace_ref)
    for (d in 1:2) {
      a <- assignments[assignments$null_draw == d, ]
      replacement <- vapply(seq_len(nrow(strata)), function(i) {
        a$replacement_item[a$subject == strata$subject[i] & a$trial == strata$trial[i]]
      }, character(1))
      donor_item <- fix$item[strata$control_id]
      keep <- strata$is_case | (donor_item != replacement & donor_item != strata$item)
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

test_that("repeated-item planted contrasts are symmetric on other-item donors", {
  tab <- reference_fixture()
  tab$density <- lapply(1:9, function(i) {
    z <- matrix(.01, 3, 3); z[i] <- 1
    gen_density(0:2, 0:2, z)
  })
  refs <- reference_maps(tab, floor = .001)
  fix <- tab[rep(1:9, each = 3), c("subject", "item")]
  fix$trial <- seq_len(nrow(fix))
  fix$x <- rep(rep(0:2, 3), each = 3)
  fix$y <- rep(rep(0:2, each = 3), each = 3)
  # Each case is on its own peak; every other-item donor is on a low cell.
  signal <- log(1 / 1.08 + .001) - log(.01 / 1.08 + .001)
  old <- sample_controls(fix, n = 8, item = NULL, replace = FALSE)
  result <- mismatch_draws(old, refs, n = 2, standardize = "none")
  expect_equal(result$observed$self_item, rep(signal * 6 / 8, 27), tolerance = 1e-12)
  expect_equal(result$paired_observed$self_item, rep(signal, 54), tolerance = 1e-12)
  expect_equal(result$contrasts$self_item, rep(0, 54), tolerance = 1e-12)
  expect_equal(result$paired_observed$n_controls, rep(3L, 54))
  expect_equal(result$paired_observed$n_controls, result$contrasts$n_controls)
  current <- sample_controls(fix, n = 6, replace = FALSE)
  clean <- mismatch_draws(current, refs, assignments = result$assignments, standardize = "none")
  expect_equal(clean$observed$self_item, rep(signal, 27), tolerance = 1e-12)
  expect_equal(clean$paired_observed$self_item, result$paired_observed$self_item, tolerance = 1e-12)
  expect_equal(clean$contrasts$self_item, result$contrasts$self_item, tolerance = 1e-12)

  # No planted preference: every repeated trial visits all cells equally.
  null <- fix[rep(seq_len(nrow(fix)), each = 9), ]
  null$x <- rep(rep(0:2, 3), nrow(fix))
  null$y <- rep(rep(0:2, each = 3), nrow(fix))
  balanced <- mismatch_draws(sample_controls(null, n = 54, replace = FALSE),
                             refs, n = 2, standardize = "none", by = "trial")
  for (part in c("observed", "paired_observed", "contrasts")) {
    expect_equal(balanced[[part]]$self_item, rep(0, nrow(balanced[[part]])), tolerance = 1e-12)
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
    # Exact oracle: scores differ by ~1e-14 around -1.7, where stats::sd() is
    # only accurate to ~1e-5 because its mean cannot be represented exactly.
    # Shifting by one score per group is exact here (Sterbenz: all scores lie
    # within a factor of two), and the SD is shift invariant, so row-level
    # scaling of the exact differences gives contrasts to full precision.
    raw <- score_strata(alt, refs, standardize = "none")
    group <- interaction(raw$subject, raw$draw, drop = TRUE)
    first <- raw$self_item[match(group, group)]
    expect_true(all(raw$self_item / first >= 0.5 & raw$self_item / first <= 2))
    shifted <- raw$self_item - first
    raw$self_item <- shifted / stats::ave(shifted, group, FUN = stats::sd)
    expected <- case_control_contrast(raw)
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

test_that("refitted scale is accurate when a group has no observed scores", {
  # Large offset, tiny spread: the per-draw variance needs a shift near the
  # group mean. s1's own maps are missing, so its observed centre is NA.
  set.seed(3)
  tab <- expand.grid(subject = c("s1", "s2"), item = paste0("i", 1:4), stringsAsFactors = FALSE)
  tab$density <- lapply(seq_len(nrow(tab)), function(i) {
    gen_density(0:3, 0:2, matrix(1e4 + stats::runif(12) * 1e-4, 4))
  })
  tab$density[tab$subject == "s1" & tab$item %in% c("i1", "i2")] <- list(NULL)
  cases <- tab[tab$item %in% c("i1", "i2"), c("subject", "item")]
  fix <- cases[rep(seq_len(nrow(cases)), each = 6), ]
  fix$trial <- paste(fix$subject, fix$item, rep(1:2, each = 3))
  fix$x <- sample(0:3, nrow(fix), TRUE)
  fix$y <- sample(0:2, nrow(fix), TRUE)
  refs <- expect_warning(reference_maps(tab, refs = "self_item", space = "density",
                                        normalize = "none"), "min_n")
  strata <- sample_controls(fix, n = 6, seed = 3)
  a <- mismatch_assignments(strata, refs, n = 1)
  result <- mismatch_draws(strata, refs, assignments = a, control_policy = "fixed")
  alt <- strata
  alt$item <- a$replacement_item[match(paste(strata$subject, strata$trial),
                                       paste(a$subject, a$trial))]
  raw <- score_strata(alt, refs, standardize = "none")
  expected <- case_control_contrast(apply_score_scaling(raw, fit_score_scaling(raw)))
  s1 <- expected$subject == "s1"
  expect_true(all(is.finite(expected$self_item[s1])))
  expect_equal(result$contrasts$self_item, expected$self_item, tolerance = 1e-10)
})


test_that("retained near-constant references are refitted without cancellation", {
  # Dropping the score 100 leaves two almost identical scores. The original
  # centre is far from both, so subtracting shifted moments loses their variance.
  tab <- expand.grid(subject = c("tiny", "small", "constant"), item = c("a", "b", "c"),
                     stringsAsFactors = FALSE)
  gaps <- c(tiny = 1e-8, small = 1e-4, constant = 0)
  tab$density <- lapply(tab$subject, function(s) {
    gen_density(0:2, 0:1, matrix(rep(c(1, 100, 1 + gaps[[s]]), 2), 3))
  })
  refs <- reference_maps(tab, refs = c("self_item", "self_other_items"),
                         space = "density", normalize = "none")
  fix <- tab[c("subject", "item")]
  fix$trial <- match(fix$item, c("a", "b", "c"))
  fix$x <- fix$trial - 1
  fix$y <- 0
  strata <- sample_controls(fix, n = 2, draws = 2, replace = FALSE)
  strata <- strata[strata$item == "a", ]
  assignments <- data.frame(subject = names(gaps), item = "a", trial = 1L,
                            replacement_item = "b", null_draw = 1L)
  for (resolution in c("fixation", "trial")) {
    result <- mismatch_draws(strata, refs, assignments = assignments, by = resolution)
    for (part in c("paired_observed", "contrasts")) {
      expected <- ifelse(result[[part]]$subject == "constant", 0, -sqrt(2))
      for (ref in refs$refs) {
        expect_equal(result[[part]][[ref]], expected, tolerance = 1e-12)
      }
    }
  }
  # The matched numerator must also be stable without refitting or scaling.
  raw <- score_strata(strata, refs, standardize = "none")
  result <- mismatch_draws(strata, refs, assignments = assignments, standardize = "none")
  for (ref in refs$refs) {
    cases <- raw[raw$is_case, ]
    retained <- raw[!raw$is_case & raw$x == 2, ]
    expected <- cases[[ref]] - retained[[ref]]
    expect_equal(result$paired_observed[[ref]], expected, tolerance = 1e-12)
  }
})

test_that("stable moment fallback preserves missing and singleton groups", {
  values <- cbind(a = c(1, 1 + 1e-8, NA, 5), b = c(2, 2, NA, NA))
  ids <- c(1L, 1L, 2L, 3L)
  shift <- matrix(34, 3, 2)
  moments <- shifted_moments(values, ids, 3, shift)
  scale <- refit_moment_scale(moments, c("a", "b"), values, ids)$scale
  expect_equal(scale[1, "a"], (values[2, "a"] - values[1, "a"]) / sqrt(2),
               tolerance = 1e-15)
  expect_equal(unname(scale[1, "b"]), 1)
  expect_true(all(is.na(scale[2, ])))
  expect_equal(unname(scale[3, "a"]), 1)
  expect_true(is.na(scale[3, "b"]))
  well_conditioned <- matrix(c(3, 0, 2), 1)
  expect_equal(refit_moment_scale(well_conditioned, "a", stop("eager values"),
                                  stop("eager ids"))$scale, matrix(1, 1, dimnames = list(NULL, "a")))
})
