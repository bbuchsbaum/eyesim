test_that("observed controls match an independent seeded sampler and retain case metadata", {
  fix <- reference_fixations()
  fix$epoch <- "test"
  set.seed(19)
  expected <- list()
  for (d in 1:2) for (i in seq_len(nrow(fix))) {
    pool <- which(fix$subject == fix$subject[i] & fix$trial != fix$trial[i])
    expected[[length(expected) + 1L]] <- pool[sample.int(length(pool), 4, replace = TRUE)]
  }
  set.seed(888)
  before <- .Random.seed
  strata <- sample_controls(fix, n = 4, draws = 2, stratify_by = "epoch", seed = 19)
  expect_identical(.Random.seed, before)
  expect_identical(strata, sample_controls(fix, n = 4, draws = 2, stratify_by = "epoch", seed = 19))
  expect_equal(strata$control_id[!strata$is_case], unlist(expected))
  controls <- strata[!strata$is_case, ]
  expect_equal(controls$x, fix$x[controls$control_id])
  expect_equal(controls$y, fix$y[controls$control_id])
  expect_equal(controls$item, fix$item[controls$case_id])
  expect_true(all(fix$subject[controls$control_id] == controls$subject))
  expect_true(all(fix$trial[controls$control_id] != controls$trial))
  one <- sample_controls(fix, n = 2, seed = 19)
  expect_true(all(vapply(split(one$control_id[!one$is_case], one$stratum[!one$is_case]),
                         function(x) !anyDuplicated(x), logical(1))))
  expect_error(sample_controls(fix, n = 3, replace = FALSE), "exceed pool")
  expect_error(sample_controls(fix[fix$trial == 1, ]), "No eligible")
})

test_that("uniform and map controls have explicit sampling distributions", {
  fix <- reference_fixations()
  uniform <- sample_controls(fix, n = 5, method = "uniform", xbounds = c(-1,1), ybounds = c(10,20))
  controls <- uniform[!uniform$is_case, ]
  expect_true(all(controls$x >= -1 & controls$x <= 1 & controls$y >= 10 & controls$y <= 20))
  expect_true(all(is.na(controls$control_id)))
  map <- gen_density(0:2, 0:1, matrix(c(0,0,1,0,0,0),3))
  sampled <- sample_controls(fix, n = 5, method = "map", control_map = map)
  expect_true(all(sampled$x[!sampled$is_case] == 2))
  expect_true(all(sampled$y[!sampled$is_case] == 0))
  expect_error(sample_controls(fix, method = "uniform"), "bounds")
  expect_error(sample_controls(fix, method = "map", control_map = map, replace = FALSE), "require replacement")
})

test_that("map controls reproduce concentrated probability mass independently", {
  z <- matrix(.001, 10, 10); z[5,5] <- 1; z[2,2] <- .5
  map <- gen_density(1:10, 1:10, z)
  fix <- data.frame(subject = "s", trial = 1:200, x = 5, y = 5)
  strata <- sample_controls(fix, n = 20, method = "map", control_map = map, seed = 13)
  controls <- strata[!strata$is_case, ]
  expect_equal(mean(controls$x == 5 & controls$y == 5), z[5,5] / sum(z), tolerance = .02)
  explicit <- sample_controls(fix, n = 20, method = "map", control_map = map, seed = 13, replace = TRUE)
  expect_equal(strata$x, explicit$x)
  expect_equal(strata$y, explicit$y)
  expect_error(sample_controls(fix, n = 1, method = "map", control_map = map, replace = FALSE), "require replacement")
})

test_that("stratification separates donor pools even with repeated trial labels", {
  fix <- reference_fixations()
  a <- fix; a$epoch <- "first"
  b <- fix; b$epoch <- "second"; b$x <- b$x + 10
  input <- rbind(a, b)
  out <- sample_controls(input, n = 2, stratify_by = "epoch")
  donor <- out$control_id[!out$is_case]
  expect_equal(input$epoch[donor], out$epoch[!out$is_case])
  expect_true(all(input$trial[donor] != out$trial[!out$is_case]))
  expect_true(all(input$subject[donor] == out$subject[!out$is_case]))
})

test_that("mismatch assignments are trial coherent, replayable and preserve RNG", {
  refs <- reference_maps(reference_fixture(), floor = .01)
  strata <- sample_controls(reference_fixations(), n = 2, draws = 2)
  set.seed(44)
  before <- .Random.seed
  assignments <- mismatch_assignments(strata, refs, n = 4, seed = 91)
  expect_identical(.Random.seed, before)
  expect_equal(assignments, mismatch_assignments(strata, refs, n = 4, seed = 91))
  expect_equal(nrow(assignments), 36)
  expect_true(all(assignments$item != assignments$replacement_item))
  expect_false(anyDuplicated(assignments[c("subject","trial","null_draw")]) > 0)
  sparse <- expect_warning(reference_maps(reference_fixture()[1:3, ], space="density"), "min_n")
  expect_error(mismatch_assignments(strata[strata$item == "a", ], sparse),
               "No other")
})

test_that("seeded constructors restore an initially absent RNG state even on error", {
  saved <- snapshot_session_rng()
  on.exit(saved())
  if (exists(".Random.seed", globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  sample_controls(reference_fixations(), n = 1)
  expect_false(exists(".Random.seed", globalenv(), inherits = FALSE))
  expect_error(sample_controls(reference_fixations(), n = 4, replace = FALSE), "exceed pool")
  expect_false(exists(".Random.seed", globalenv(), inherits = FALSE))
})
