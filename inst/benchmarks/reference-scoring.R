# Synthetic scale benchmark for issue #53. Run against an installed package,
# or devtools::load_all() and then source this file. No research data are used.
# Default: 98,880 test fixations, 20 observed controls, and 100 null draws.
library(eyesim)

benchmark_reference_scoring <- function(subjects = 103L, items = 120L,
                                        grid = 20L, null_draws = 100L, seed = 53L,
                                        fixations = 8L, controls = 20L) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else
      rm(".Random.seed", envir = globalenv())
  })
  set.seed(seed)
  maps <- expand.grid(subject = sprintf("s%03d", seq_len(subjects)),
                       item = sprintf("i%03d", seq_len(items)), stringsAsFactors = FALSE)
  axis <- seq(0, 100, length.out = grid)
  maps$density <- lapply(seq_len(nrow(maps)), function(i) {
    gen_density(axis, axis, matrix(rexp(grid^2), grid))
  })
  fix <- maps[rep(seq_len(nrow(maps)), each = fixations), c("subject", "item")]
  fix$trial <- fix$item
  fix$x <- runif(nrow(fix), 0, 100)
  fix$y <- runif(nrow(fix), 0, 100)
  timings <- list()
  timings$references <- system.time(refs <- reference_maps(maps, floor = 1 / (10 * grid^2)))
  timings$controls <- system.time(strata <- sample_controls(fix, n = controls))
  timings$scoring <- system.time(scored <- score_strata(strata, refs))
  timings$mismatch <- system.time(null <- mismatch_draws(strata, refs, n = null_draws, by = "trial"))
  stopifnot(nrow(null$contrasts) == nrow(maps) * null_draws,
            nrow(null$paired_observed) == nrow(null$contrasts),
            all(is.finite(null$contrasts$self_item)))
  list(timings = timings, reference_bytes = as.numeric(object.size(refs)),
       cases = nrow(fix), strata_rows = nrow(strata), grid = c(grid, grid),
       strata_bytes = as.numeric(object.size(strata)), controls = controls,
       fixations_per_trial = fixations,
       null_draws = null_draws, session = sessionInfo())
}

if (sys.nframe() == 0L) print(benchmark_reference_scoring())
