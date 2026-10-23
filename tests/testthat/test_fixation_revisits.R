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
