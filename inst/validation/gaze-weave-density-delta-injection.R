# Post-review (2026-09-24) semi-synthetic power on the REAL template structure
# for the density Delta court.
#
# This file sources the frozen court script unchanged; it lives in a separate
# file so that the 1.1.0 freeze hashes of gaze-weave-density-delta.R remain
# valid. It is not part of the frozen protocol.
#
# Retrieval positions are redrawn from the real-data null (1) generator (the
# fold-specific fitted M0 mixture on the real group/background templates) and
# a share s of fixations is replaced by draws from the participant's own
# template O_si (h = 40 px). Real fixation counts and durations are kept. The
# own and pseudo-own designs are then fully cross-fitted and the paired
# own - pseudo decision rule is applied. Participant-linked replicate output
# stays in the ignored results directory; only the aggregate is written to
# the committed simulation directory.
#
#   devtools::load_all()
#   source("inst/validation/gaze-weave-density-delta-injection.R")
#   run_density_delta_injection()

density_delta_injection_file <- system.file(
  "validation", "gaze-weave-density-delta.R", package = "eyesim"
)
if (!nzchar(density_delta_injection_file)) {
  density_delta_injection_file <- file.path(
    "inst", "validation", "gaze-weave-density-delta.R"
  )
}
source(density_delta_injection_file, local = TRUE)

density_delta_inject_own <- function(y, own_templates, share, h, screen, seed) {
  set.seed(seed)
  for (t in seq_along(y)) {
    inject <- stats::runif(nrow(y[[t]])) < share
    if (any(inject)) {
      xy <- density_delta_sample_template(own_templates[[t]], sum(inject), h,
                                          screen)
      y[[t]]$x[inject] <- xy$x
      y[[t]]$y[inject] <- xy$y
    }
  }
  y
}

run_density_delta_injection <- function(
    shares = c(0, 0.05, 0.10, 0.15, 0.20), replicates = 40L, cores = 6L,
    injection_bandwidth = 40, config = density_delta_config(),
    output_dir = density_delta_result_dir,
    summary_dir = density_delta_simulation_dir) {
  density_delta_verify_freeze(config = config)
  context <- density_delta_context()
  folds <- density_delta_real_folds(context)
  own_design <- density_delta_real_design(context, "old_lure", "primary",
                                          "combined", config)
  pseudo_design <- density_delta_real_design(context, "old_lure", "pseudo_own",
                                             "combined", config)
  if (!identical(own_design$trials[c("participant", "item")],
                 pseudo_design$trials[c("participant", "item")])) {
    stop("Own and pseudo-own designs do not align.")
  }
  gen <- density_delta_generating_fits(own_design, folds)
  own_design$fold_of <- gen$fold_of
  jobs <- expand.grid(replicate = seq_len(replicates), share = shares)
  started <- proc.time()[["elapsed"]]
  reps <- parallel::mclapply(seq_len(nrow(jobs)), function(k) {
    share <- jobs$share[[k]]
    seed <- config$seed + 700000L + k
    y <- density_delta_null1_draw(own_design, gen$fits, seed)
    y <- density_delta_inject_own(y, own_design$templates$own, share,
                                  injection_bandwidth, config$screen, seed + 1L)
    own <- density_delta_crossfit(density_delta_redesign(own_design, y), folds)
    pseudo <- density_delta_crossfit(density_delta_redesign(pseudo_design, y),
                                     folds)
    if (!identical(own$scores$trial, pseudo$scores$trial)) {
      stop("Own and pseudo-own scores do not align.")
    }
    draws <- config$simulation_bootstrap_draws
    alpha <- config$alpha_one_sided
    raw <- density_delta_crossed_test(own$scores, own$scores$delta, "raw",
                                      draws, seed, alpha)
    omp <- density_delta_crossed_test(
      own$scores, own$scores$delta - pseudo$scores$delta, "own_minus_pseudo",
      draws, seed, alpha
    )
    data.frame(
      share = share, replicate = jobs$replicate[[k]],
      raw_estimate = raw$estimate, raw_reject = raw$reject_one_sided,
      contrast_estimate = omp$estimate, contrast_lower = omp$lower_95,
      contrast_upper = omp$upper_95, contrast_reject = omp$reject_one_sided,
      own_weight = mean(own$fits$w1_own),
      pseudo_weight = mean(pseudo$fits$w1_own)
    )
  }, mc.cores = cores, mc.preschedule = FALSE)
  if (any(vapply(reps, inherits, logical(1), "try-error"))) {
    stop("Injection replicates failed.")
  }
  reps <- do.call(rbind, reps)
  summary <- do.call(rbind, lapply(split(reps, reps$share), function(r) {
    ci <- density_delta_rate_interval(sum(r$contrast_reject), nrow(r))
    data.frame(
      share = r$share[[1L]], replicates = nrow(r),
      contrast_power = mean(r$contrast_reject),
      power_lower_95 = ci[["lower"]], power_upper_95 = ci[["upper"]],
      mean_contrast = mean(r$contrast_estimate),
      sd_contrast = stats::sd(r$contrast_estimate),
      mean_bootstrap_se = mean((r$contrast_upper - r$contrast_lower) /
                                 (2 * stats::qnorm(0.975))),
      raw_rejection_rate = mean(r$raw_reject),
      mean_raw_delta = mean(r$raw_estimate),
      mean_own_weight = mean(r$own_weight),
      mean_pseudo_weight = mean(r$pseudo_weight)
    )
  }))
  summary$injection_bandwidth <- injection_bandwidth
  summary$elapsed_seconds_total <- proc.time()[["elapsed"]] - started
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(reps, file.path(output_dir, "injection-replicates.csv"),
                   row.names = FALSE)
  utils::write.csv(summary, file.path(summary_dir, "injection-summary.csv"),
                   row.names = FALSE)
  invisible(summary)
}
