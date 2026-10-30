fg_abca <- function() {
  fixation_group(x = c(100, 400, 400, 105), y = c(100, 100, 400, 98),
                 onset = c(0, 250, 500, 750), duration = rep(200, 4))
}

test_that("A-B-C-A scores one revisit of the first fixation", {
  rv <- fixation_revisits(fg_abca(), radius = 50)
  expect_s3_class(rv, "revisit_table")
  expect_equal(rv$opportunity, c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(rv$revisit, c(FALSE, FALSE, FALSE, TRUE))
  expect_equal(rv$match_index[4], 1L)
  expect_equal(rv$lag[4], 3L)
  expect_equal(rv$distance[4], sqrt(5^2 + 2^2))
  s <- revisit_summary(rv)
  expect_equal(s$n_fix, 4L)
  expect_equal(s$opportunities, 2)
  expect_equal(s$revisits, 1)
  expect_equal(s$revisit_rate, 0.5)
  expect_true(is.na(s$expected_rate))
})

test_that("a dwell cluster is a revisit only without the departure rule", {
  fg <- fixation_group(x = c(100, 105, 110), y = c(100, 100, 100), onset = c(0, 1, 2), duration = c(1, 1, 1))
  strict <- fixation_revisits(fg, radius = 50)
  expect_equal(sum(strict$opportunity), 0)
  expect_equal(sum(strict$revisit), 0)
  loose <- fixation_revisits(fg, radius = 50, require_departure = FALSE)
  expect_equal(loose$revisit, c(FALSE, FALSE, TRUE))
  expect_equal(loose$lag[3], 2L)
  # min_lag = 1 also lets a fixation revisit the one immediately before it
  lag1 <- fixation_revisits(fg, radius = 50, min_lag = 1, require_departure = FALSE)
  expect_equal(lag1$revisit, c(FALSE, TRUE, TRUE))
  expect_equal(lag1$lag[3], 1L)
})

test_that("the most recent matching fixation is reported and max_lag limits the look-back", {
  # A B A' C A'': A'' matches A (lag 4) and A' (lag 2); the most recent is A'
  fg <- fixation_group(x = c(0, 300, 5, 300, 10), y = c(0, 0, 0, 300, 0), onset = 1:5, duration = rep(1, 5))
  rv <- fixation_revisits(fg, radius = 50)
  expect_equal(rv$revisit, c(FALSE, FALSE, TRUE, FALSE, TRUE))
  expect_equal(rv$match_index[5], 3L)
  expect_equal(rv$lag[5], 2L)
  # a look-back of exactly two fixations still finds both returns
  rv2 <- fixation_revisits(fg, radius = 50, min_lag = 2, max_lag = 2)
  expect_equal(rv2$revisit, c(FALSE, FALSE, TRUE, FALSE, TRUE))
  # exactly three back: C is compared with A and A'' with B, neither within the radius
  rv3 <- fixation_revisits(fg, radius = 50, min_lag = 3, max_lag = 3)
  expect_equal(sum(rv3$revisit), 0)
})

test_that("the spatial null gives the conditional proportion of reference locations", {
  fg <- fg_abca()
  # reference: one point near A, one near B, one near C and one far away. A reference point only counts
  # for fixation i if it is also >= radius from fixation i - 1 (the departure condition).
  ref <- rbind(c(102, 101), c(398, 102), c(401, 399), c(800, 800))
  rv <- fixation_revisits(fg, radius = 50, reference = ref)
  # fixation 3: eligible {A}; departure from B excludes the near-B point -> 3 valid, 1 near A
  expect_equal(rv$expected[3], 1 / 3)
  # fixation 4: eligible {A, B}; departure from C excludes the near-C point -> 3 valid, near A or B: 2
  expect_equal(rv$expected[4], 2 / 3)
  s <- revisit_summary(rv)
  expect_equal(s$expected_revisits, 1)
  expect_equal(s$expected_rate, 0.5)
  expect_equal(s$excess_revisits, 0)
  # data-frame references with x/y columns are accepted
  rv_df <- fixation_revisits(fg, radius = 50, reference = data.frame(x = ref[, 1], y = ref[, 2]))
  expect_equal(rv_df$expected, rv$expected)
})

test_that("results do not depend on row order, only on onset order", {
  fg <- fg_abca()
  shuffled <- fg[c(3, 1, 4, 2), ]
  class(shuffled) <- class(fg)
  expect_equal(fixation_revisits(shuffled, radius = 50)$revisit, fixation_revisits(fg, radius = 50)$revisit)
})

test_that("manhattan distance changes the location criterion", {
  fg <- fixation_group(x = c(0, 300, 30), y = c(0, 0, 30), onset = 1:3, duration = rep(1, 3))
  expect_true(fixation_revisits(fg, radius = 50)$revisit[3])                           # euclidean 42.4
  expect_false(fixation_revisits(fg, radius = 50, dist_method = "manhattan")$revisit[3]) # manhattan 60
})

test_that("empty and single-fixation episodes have no opportunities", {
  empty <- fixation_group(numeric(0), numeric(0), numeric(0), numeric(0))
  rv <- fixation_revisits(empty, radius = 50)
  expect_equal(nrow(rv), 0)
  s <- revisit_summary(rv)
  expect_equal(s$opportunities, 0)
  expect_true(is.na(s$revisit_rate))
  one <- fixation_group(100, 100, onset = 0, duration = 200)
  expect_equal(revisit_summary(fixation_revisits(one, radius = 50))$opportunities, 0)
})

test_that("invalid arguments are rejected", {
  fg <- fg_abca()
  expect_error(fixation_revisits(fg), "radius")
  expect_error(fixation_revisits(fg, radius = -1), "radius")
  expect_error(fixation_revisits(fg, radius = 50, min_lag = 0), "min_lag")
  expect_error(fixation_revisits(fg, radius = 50, min_lag = 3, max_lag = 2), "max_lag")
  expect_error(fixation_revisits(fg, radius = 50, reference = matrix(1:3, 1)), "two columns")
  expect_error(revisit_summary(data.frame()), "fixation_revisits")
})

test_that("the eye_table method summarises every episode with a leave-episode-out subject null", {
  set.seed(3)
  et <- simulate_eye_table(n_fixations = 80, n_groups = 4)
  et$subject <- c("s1", "s1", "s2", "s2")
  res <- fixation_revisits(et, radius = 200, null = "subject", subject_var = "subject", keep_fixations = TRUE)
  expect_s3_class(res, "eye_table")
  expect_equal(nrow(res), 4)
  expect_true(all(c("n_fix", "opportunities", "revisits", "revisit_rate", "expected_revisits",
                    "excess_revisits") %in% names(res)))
  # row 1's reference is exactly row 2's fixations (same subject, other episode)
  ref <- cbind(et$fixgroup[[2]]$x, et$fixgroup[[2]]$y)
  direct <- revisit_summary(fixation_revisits(et$fixgroup[[1]], radius = 200, reference = ref))
  expect_equal(res$expected_revisits[1], direct$expected_revisits)
  expect_equal(res$revisits[1], direct$revisits)
  expect_gt(sum(res$opportunities), 0)
  expect_true(all(is.finite(res$expected_rate[res$opportunities > 0])))
  expect_s3_class(res$revisit_fixations[[1]], "revisit_table")
  # without a null, expected columns are NA
  plain <- fixation_revisits(et, radius = 200)
  expect_true(all(is.na(plain$expected_rate)))
  expect_error(fixation_revisits(et, radius = 200, null = "subject"), "subject_var")
})

test_that("the item null pools other subjects viewing the same item", {
  set.seed(4)
  et <- simulate_eye_table(n_fixations = 120, n_groups = 6)
  et$subject <- rep(c("s1", "s2", "s3"), each = 2)
  et$image <- rep(c("a", "b"), times = 3)
  res <- fixation_revisits(et, radius = 200, null = "item", item_var = "image", subject_var = "subject")
  # row 1 (s1, image a): reference = rows 3 and 5 (s2, s3 viewing a)
  ref <- do.call(rbind, lapply(et$fixgroup[c(3, 5)], function(f) cbind(f$x, f$y)))
  direct <- revisit_summary(fixation_revisits(et$fixgroup[[1]], radius = 200, reference = ref))
  expect_equal(res$expected_revisits[1], direct$expected_revisits)
  # without subject_var only the episode itself is left out (rows 3 and 5 for row 1, same here; check row 3)
  et2 <- et; et2$subject[1] <- "s2"   # row 1 now shares subject s2 with row 3
  with_subj <- fixation_revisits(et2, radius = 200, null = "item", item_var = "image", subject_var = "subject")
  no_subj <- fixation_revisits(et2, radius = 200, null = "item", item_var = "image")
  ref3_subj <- do.call(rbind, lapply(et2$fixgroup[5], function(f) cbind(f$x, f$y)))           # excludes row 1 (same subject)
  ref3_all <- do.call(rbind, lapply(et2$fixgroup[c(1, 5)], function(f) cbind(f$x, f$y)))     # leave-episode-out only
  expect_equal(with_subj$expected_revisits[3], revisit_summary(fixation_revisits(et2$fixgroup[[3]], radius = 200, reference = ref3_subj))$expected_revisits)
  expect_equal(no_subj$expected_revisits[3], revisit_summary(fixation_revisits(et2$fixgroup[[3]], radius = 200, reference = ref3_all))$expected_revisits)
  expect_error(fixation_revisits(et, radius = 200, null = "item"), "item_var")
  expect_error(fixation_revisits(et, radius = 200, null = "item", item_var = "image", subject_var = "nope"), "subject_var")
})

test_that("re-scoring a scored eye_table replaces the summary columns", {
  set.seed(5)
  et <- simulate_eye_table(n_fixations = 60, n_groups = 3)
  first <- fixation_revisits(et, radius = 100, keep_fixations = TRUE)
  again <- fixation_revisits(first, radius = 300)
  expect_equal(names(again), setdiff(names(first), "revisit_fixations"))
  expect_equal(again$revisits, fixation_revisits(et, radius = 300)$revisits)
  expect_equal(attr(again, "origin"), attr(et, "origin"))
  expect_equal(attr(again, "revisit_columns"), attr(fixation_revisits(et, radius = 300), "revisit_columns"))
})

test_that("a caller's own column with a summary name is never overwritten", {
  set.seed(6)
  et <- simulate_eye_table(n_fixations = 60, n_groups = 3)
  et$n_fix <- c(99L, 98L, 97L)
  expect_error(fixation_revisits(et, radius = 100), "n_fix")
  expect_equal(et$n_fix, c(99L, 98L, 97L))
  # renamed, it is kept alongside the summary
  names(et)[names(et) == "n_fix"] <- "tracker_n_fix"
  res <- fixation_revisits(et, radius = 100)
  expect_equal(res$tracker_n_fix, c(99L, 98L, 97L))
})

test_that("an empty reference gives NA expectations, not zero", {
  rv <- fixation_revisits(fg_abca(), radius = 50, reference = matrix(numeric(0), 0, 2))
  s <- revisit_summary(rv)
  expect_true(is.na(s$expected_revisits))
  expect_true(is.na(s$excess_rate))
  # an item seen by nobody else has no null either
  et <- simulate_eye_table(n_fixations = 40, n_groups = 2)
  et$subject <- c("s1", "s2"); et$image <- c("a", "b")
  res <- fixation_revisits(et, radius = 200, null = "item", item_var = "image", subject_var = "subject")
  expect_true(all(is.na(res$expected_revisits)))
})

test_that("invalid fixation values are rejected without closing gaps", {
  fg <- fixation_group(c(0, 100, 200, 0), rep(0, 4), rep(1, 4), 1:4)
  for (column in c("x", "y", "onset")) {
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
      bad <- fg
      bad[[column]][3] <- value
      expect_error(fixation_revisits(bad, radius = 10), paste0("column '", column, "'.*finite"))
    }
    bad <- fg
    bad[[column]] <- as.character(bad[[column]])
    expect_error(fixation_revisits(bad, radius = 10), "finite numeric")
  }
  bad <- fg
  bad$x[3] <- NA_real_
  et <- as_eye_table(tibble::tibble(fixgroup = list(fg, bad)))
  expect_error(fixation_revisits(et, radius = 10), "column 'x'.*finite")
})

test_that("lag and logical validation is consistent for empty and nonempty inputs", {
  empty <- fixation_group(numeric(), numeric(), numeric(), numeric())
  et <- as_eye_table(tibble::tibble(fixgroup = list()))
  for (input in list(fg_abca(), empty, et)) {
    for (value in list(2.9, 0, -Inf, NA_real_, NaN, "3", TRUE, 2 + 1i, numeric(), c(2, 3))) {
      expect_error(fixation_revisits(input, radius = 50, max_lag = value), "max_lag")
    }
    for (value in list(1.5, Inf, NA_real_, "2", TRUE, numeric())) {
      expect_error(fixation_revisits(input, radius = 50, min_lag = value), "min_lag")
    }
    for (value in list(Inf, NA_real_, "50", TRUE, 1 + 1i, numeric())) {
      expect_error(fixation_revisits(input, radius = value), "radius")
    }
    for (value in list(NA, 1, "TRUE", logical(), c(TRUE, FALSE))) {
      expect_error(fixation_revisits(input, radius = 50, require_departure = value), "require_departure")
    }
  }
  expect_error(fixation_revisits(et, radius = 50, keep_fixations = NA), "keep_fixations")
  unlimited <- fixation_revisits(fg_abca(), radius = 50)
  large <- fixation_revisits(fg_abca(), radius = 50, max_lag = 3e9)
  expect_equal(large$revisit, unlimited$revisit)
  expect_equal(large$lag, unlimited$lag)
  expect_false(any(fixation_revisits(fg_abca(), radius = 50, min_lag = 3e9)$opportunity))
})

test_that("invalid references are rejected even for an empty episode", {
  empty <- fixation_group(numeric(), numeric(), numeric(), numeric())
  for (fg in list(fg_abca(), empty)) {
    for (value in c(NA_real_, NaN, Inf, -Inf)) {
      expect_error(fixation_revisits(fg, radius = 50, reference = rbind(c(0, 0), c(value, 0))),
                   "reference.*finite")
    }
    expect_error(fixation_revisits(fg, radius = 50, reference = matrix("x", 1, 2)), "finite numeric")
    expect_error(fixation_revisits(fg, radius = 50, reference = matrix(1:3, 1)), "two columns")
  }
})

test_that("empty eye tables retain a typed summary and scoring specification", {
  empty <- as_eye_table(tibble::tibble(subject = character(), image = character(), fixgroup = list()))
  attr(empty, "origin") <- c(12, 34)
  for (null in c("none", "subject", "item")) {
    res <- fixation_revisits(empty, radius = 50, null = null, subject_var = "subject",
                             item_var = "image", keep_fixations = TRUE, dist_method = "manhattan")
    expect_s3_class(res, "eye_table")
    expect_equal(nrow(res), 0L)
    expect_identical(res$n_fix, integer())
    expect_identical(res$opportunities, integer())
    expect_identical(res$expected_rate, double())
    expect_identical(res$revisit_fixations, list())
    expect_equal(attr(res, "origin"), c(12, 34))
    expect_equal(attr(res, "revisit_spec")$null, null)
    again <- fixation_revisits(res, radius = 20)
    expect_false("revisit_fixations" %in% names(again))
    expect_equal(attr(again, "revisit_spec")$radius, 20)
  }
  singleton <- as_eye_table(tibble::tibble(fixgroup = list(fg_abca())))
  prototype <- fixation_revisits(singleton, radius = 50)[0, ]
  result <- fixation_revisits(singleton[0, ], radius = 50)
  expect_equal(result, prototype)
})

test_that("scoring provenance is retained and updated on summary-only tables", {
  rv <- fixation_revisits(fg_abca(), radius = 50, dist_method = "manhattan")
  expect_equal(attr(rv, "dist_method"), "manhattan")
  empty <- fixation_group(numeric(), numeric(), numeric(), numeric())
  expect_equal(attr(fixation_revisits(empty, radius = 50), "dist_method"), "euclidean")
  et <- as_eye_table(tibble::tibble(subject = c("s1", "s2"), image = c("a", "a"),
                                   scans = list(fg_abca(), fg_abca())))
  res <- fixation_revisits(et, radius = 40, min_lag = 1, max_lag = 3,
                           require_departure = FALSE, dist_method = "manhattan",
                           null = "item", subject_var = "subject", item_var = "image", var = "scans")
  expect_false("revisit_fixations" %in% names(res))
  expect_equal(attr(res, "revisit_spec"),
               list(radius = 40, min_lag = 1, max_lag = 3, require_departure = FALSE,
                    dist_method = "manhattan", null = "item", subject_var = "subject",
                    item_var = "image", var = "scans"))
  again <- fixation_revisits(res, radius = 50, var = "scans")
  expect_equal(attr(again, "revisit_spec"), attr(fixation_revisits(et, radius = 50, var = "scans"), "revisit_spec"))
})

test_that("nearest distance and most recent match distance are distinct", {
  fg <- fixation_group(c(0, 100, 8, 100, 1), rep(0, 5), rep(1, 5), 1:5)
  rv <- fixation_revisits(fg, radius = 10)
  expect_equal(rv$match_index[5], 3L)
  expect_equal(rv$nearest_distance[5], 1)
  expect_equal(rv$match_distance[5], 7)
  expect_identical(rv$distance, rv$nearest_distance)
  expect_true(all(is.na(rv$match_distance[!rv$revisit])))
})

test_that("incoming steps and strict radius boundaries define returns", {
  fg <- fixation_group(c(0, -9, 9), rep(0, 3), rep(1, 3), 1:3)
  expect_true(fixation_revisits(fg, radius = 10)$revisit[3])
  # The incoming step equals radius: eligible. A target exactly radius away is not a match.
  edge <- fixation_group(c(0, 20, 10), rep(0, 3), rep(1, 3), 1:3)
  rv <- fixation_revisits(edge, radius = 10, reference = rbind(c(10, 0), c(0, 0)))
  expect_true(rv$opportunity[3])
  expect_false(rv$revisit[3])
  expect_equal(rv$expected[3], 0.5)
})

test_that("window expiry preserves overlap counts and updates during non-opportunities", {
  # At i=5, fixation 1 expires but fixation 2 still covers the reference. At i=6,
  # fixation 2 also expires and no eligible fixation covers it.
  fg <- fixation_group(c(0, 0, 100, 100, 100, 0), rep(0, 6), rep(1, 6), 1:6)
  rv <- fixation_revisits(fg, radius = 10, min_lag = 2, max_lag = 3,
                          require_departure = FALSE, reference = matrix(c(0, 0), 1))
  expect_equal(rv$expected[3:6], c(1, 1, 1, 0))
  strict <- fixation_revisits(fg, radius = 10, min_lag = 2, max_lag = 3,
                              reference = matrix(c(0, 0), 1))
  expect_false(any(strict$opportunity[4:5]))
  expect_true(strict$opportunity[6])
  expect_equal(strict$expected[6], 0)
})

test_that("one unsupported opportunity keeps all aggregate expectations unavailable", {
  fg <- fixation_group(c(0, 100, 200, 0), rep(0, 4), rep(1, 4), 1:4)
  rv <- fixation_revisits(fg, radius = 10, reference = matrix(c(200, 0), 1))
  expect_equal(rv$expected[3:4], c(0, NA_real_))
  s <- revisit_summary(rv)
  expect_true(all(is.na(unlist(s[c("expected_revisits", "expected_rate", "excess_revisits", "excess_rate")]))))
})

test_that("Euclidean distances avoid avoidable overflow and underflow", {
  for (scale in c(1e200, 1e-200)) {
    fg <- fixation_group(c(0, 20, 3) * scale, c(0, 0, 4) * scale, rep(1, 3), 1:3)
    rv <- fixation_revisits(fg, radius = 6 * scale, reference = matrix(c(3, 4) * scale, 1))
    expect_true(rv$revisit[3])
    expect_equal(rv$nearest_distance[3] / scale, 5, tolerance = 1e-14)
    expect_equal(rv$expected[3], 1)
    outside <- fixation_revisits(fg, radius = 4 * scale)
    expect_false(outside$revisit[3])
  }
  # A true distance beyond double range is Inf, not NaN; coincident points stay zero.
  huge <- fixation_group(c(-1e308, 1e308, -1e308), c(0, 0, 0), rep(1, 3), 1:3)
  expect_equal(fixation_revisits(huge, radius = 1e307)$match_distance[3], 0)
  integer_coords <- fixation_group(c(-2000000000L, 2000000000L, -2000000000L),
                                   c(0L, 0L, 0L), rep(1, 3), 1:3)
  expect_true(fixation_revisits(integer_coords, radius = 10)$revisit[3])
})

# Deliberately scalar: reconstruct history and test every candidate/reference
# location independently, without production distance or window helpers.
scalar_revisits <- function(fg, radius, min_lag, max_lag, require_departure, reference, method) {
  fg <- fg[order(fg$onset), ]
  n <- nrow(fg)
  dist <- function(a, b) {
    delta <- abs(a - b)
    if (method == "euclidean") sqrt(sum(delta^2)) else sum(delta)
  }
  xy <- cbind(fg$x, fg$y)
  result <- data.frame(index = fg$index, onset = fg$onset,
                       opportunity = rep(FALSE, n), revisit = rep(FALSE, n),
                       lag = rep(NA_integer_, n), match_index = rep(NA_integer_, n),
                       distance = rep(NA_real_, n), nearest_distance = rep(NA_real_, n),
                       match_distance = rep(NA_real_, n), expected = rep(NA_real_, n))
  for (i in seq_len(n)) {
    history <- Filter(function(j) i - j >= min_lag && i - j <= max_lag, seq_len(i - 1L))
    if (!length(history)) next
    if (require_departure && dist(xy[i, ], xy[i - 1L, ]) < radius) next
    result$opportunity[i] <- TRUE
    distances <- vapply(history, function(j) dist(xy[i, ], xy[j, ]), numeric(1))
    result$distance[i] <- result$nearest_distance[i] <- min(distances)
    for (j in history) {
      if (dist(xy[i, ], xy[j, ]) < radius) {
        result$revisit[i] <- TRUE
        result$match_index[i] <- fg$index[j]
        result$lag[i] <- i - j
        result$match_distance[i] <- dist(xy[i, ], xy[j, ])
      }
    }
    if (is.null(reference)) next
    numerator <- denominator <- 0L
    for (r in seq_len(nrow(reference))) {
      if (require_departure && dist(reference[r, ], xy[i - 1L, ]) < radius) next
      denominator <- denominator + 1L
      hit <- FALSE
      for (j in history) {
        if (dist(reference[r, ], xy[j, ]) < radius) hit <- TRUE
      }
      numerator <- numerator + hit
    }
    if (denominator > 0L) result$expected[i] <- numerator / denominator
  }
  result
}

test_that("incremental scoring agrees with a scalar oracle across 500 seeded cases", {
  set.seed(20261029)
  for (case in seq_len(500)) {
    n <- sample(0:20, 1)
    fg <- fixation_group(sample(-8:8, n, replace = TRUE), sample(-8:8, n, replace = TRUE),
                          rep(1, n), sample(0:20, n, replace = TRUE))
    m <- sample(0:20, 1)
    ref <- if (case %% 7L == 0L) NULL else matrix(sample(-8:8, 2 * m, replace = TRUE), ncol = 2)
    min_lag <- sample(1:6, 1)
    max_lag <- sample(c(min_lag:(min_lag + 8), Inf, 3e9), 1)
    departure <- sample(c(TRUE, FALSE), 1)
    method <- sample(c("euclidean", "manhattan"), 1)
    radius <- sample(c(0.5, 1, 3, 6.25, 12, 50), 1)
    rv <- fixation_revisits(fg, radius = radius, min_lag = min_lag, max_lag = max_lag,
                            require_departure = departure, reference = ref, dist_method = method)
    expected <- scalar_revisits(fg, radius, min_lag, max_lag, departure, ref, method)
    # Relative 1e-14 tolerance covers rounding differences in scaled hypot;
    # logical decisions, indices and NA support must agree exactly.
    discrete <- c("index", "onset", "opportunity", "revisit", "lag", "match_index")
    expect_identical(lapply(rv[discrete], identity), lapply(expected[discrete], identity),
                      info = paste("case", case))
    expect_equal(lapply(rv, identity), lapply(expected, identity), tolerance = 1e-14,
                  info = paste("case", case))
  }
})
