# GazeWeave A3 shared calibration study (revision 2026.10). Synthetic data only.
#
# Produces the before/after numbers reported in NEWS.md for the A3 gates:
#   Rscript inst/validation/gaze-weave-calibration-a3.R <job>
# with <job> one of tnull, rnull, tsparse, rsparse, ttyp, rtyp, tdisc, rdisc,
# run from the package root (or set EYESIM_ROOT). "old" is the pre-A3
# 2026.10 calibration (method = "global"), "new" the evidence-scaled
# calibration without a typicality offset, "new+typ-mean" and "new+typ-std"
# add the mean and standardized offsets. The simulators and summaries are
# those of tests/testthat/test_gaze_weave_calibration_revision.R, whose slow
# block (EYESIM_SLOW_TESTS=true) asserts the gates.
job <- commandArgs(TRUE)[[1]]
WT <- Sys.getenv("EYESIM_ROOT", ".")
suppressMessages(devtools::load_all(WT, quiet = TRUE))
library(testthat)
env <- environment()
exprs <- parse(file.path(WT, "tests/testthat/test_gaze_weave_calibration_revision.R"))
for (e in exprs) if (!(is.call(e) && identical(e[[1]], as.name("test_that")))) eval(e, env)
OLD_T <- list(reliability = "effective_fixations", calibration_control = gaze_calibration_control(method = "global"))
OLD_R <- list(calibration_control = gaze_calibration_control(method = "global"))
NEW <- list(calibration_control = gaze_calibration_control(typicality = "none"))
TYP <- list(calibration_control = gaze_calibration_control(typicality = "mean"))
TYPS <- list(calibration_control = gaze_calibration_control(typicality = "standardized"))
CFGS_T <- list(list("old", OLD_T), list("new", NEW), list("new+typ-mean", TYP), list("new+typ-std", TYPS))
CFGS_R <- list(list("old", OLD_R), list("new", NEW), list("new+typ-mean", TYP), list("new+typ-std", TYPS))
rt <- function(data, args, seed) do.call(run_calibration_transport, c(list(data), args, list(seed = seed)))
rr <- function(data, args, seed) do.call(run_calibration_replay, c(list(data), args, list(seed = seed)))
summ <- function(fits, label) {
  r <- do.call(rbind, lapply(fits, function(f) f$results[, c("top1_credit", "gaze_info_bits", "log_loss", "template_rank", "candidate_count")]))
  chance <- 1 / r$candidate_count
  auc <- (r$candidate_count - r$template_rank) / (r$candidate_count - 1)
  cat(sprintf("%-28s n=%4d top1=%.3f (chance %.3f, mc %.3f) bits=%.4f (se %.4f) logloss=%.4f auc=%.3f\n", label, nrow(r),
    mean(r$top1_credit), mean(chance), sqrt(sum(chance * (1 - chance))) / nrow(r), mean(r$gaze_info_bits),
    sd(r$gaze_info_bits) / sqrt(nrow(r)), mean(r$log_loss), mean(auc)))
  invisible(r)
}
temps <- function(fits, label) {
  t <- do.call(rbind, lapply(fits, function(f) if (!is.null(f$calibration$fold_calibration)) f$calibration$fold_calibration[, c("temperature", "gamma")] else
    data.frame(temperature = sapply(f$folds, function(x) x$calibration$temperature), gamma = NA)))
  cat(sprintf("   %s fold T: %s | gamma: %s\n", label, paste(signif(t$temperature, 3), collapse = " "), paste(round(t$gamma, 2), collapse = " ")))
}
if (job == "tnull") {
  data <- lapply(1:4, function(s) simulate_calibration_transport(8, 6, recall = "centre", seed = 20 + s))
  for (cfg in list(list("old(global+kappa)", OLD_T), list("new", NEW), list("new+typ-mean", TYP), list("new+typ-std", TYPS))) {
    fits <- lapply(1:4, function(s) rt(data[[s]], cfg[[2]], s)); summ(fits, cfg[[1]]); temps(fits, cfg[[1]]) }
}
if (job == "rnull") {
  data <- lapply(1:4, function(s) simulate_calibration_replay(8, 8, recall = "centre", seed = 30 + s))
  for (cfg in list(list("old(global)", OLD_R), list("new", NEW), list("new+typ-mean", TYP), list("new+typ-std", TYPS))) {
    fits <- lapply(1:4, function(s) rr(data[[s]], cfg[[2]], s)); summ(fits, cfg[[1]]); temps(fits, cfg[[1]]) }
}
strata <- function(fit, label) {
  r <- fit$results
  n <- r$raw_fixation_count
  for (k in sort(unique(n))) {
    rows <- n == k
    maxp <- sapply(r$candidates[rows], function(c) max(c$posterior))
    acc <- sapply(r$candidates[rows], function(c) c$is_true[which.max(c$posterior)])
    cat(sprintf("   %s n_fix=%d rows=%d logloss=%.4f mean_maxp=%.3f acc=%.3f gap=%.3f\n", label, k, sum(rows), mean(r$log_loss[rows]), mean(maxp), mean(acc), mean(maxp) - mean(acc)))
  }
}
if (job == "tsparse") {
  fits <- list()
  for (s in 1:2) {
    data <- simulate_calibration_transport(10, 6, keep = c(2, 8), seed = 40 + s)
    for (cfg in list(list("old", OLD_T), list("gamma0", list(calibration_control = gaze_calibration_control(gamma_bounds = c(0, 0), typicality = "none"))), list("new", NEW))) {
      f <- rt(data, cfg[[2]], s); fits[[cfg[[1]]]] <- c(fits[[cfg[[1]]]], list(f)); strata(f, paste(cfg[[1]], "seed", s))
    }
  }
  for (k in names(fits)) { summ(fits[[k]], k); temps(fits[[k]], k) }
}
if (job == "rsparse") {
  fits <- list()
  for (s in 1:2) {
    data <- simulate_calibration_replay(8, 8, n_fixations = c(6, 30), rho = 0.3, seed = 45 + s)
    for (cfg in list(list("old", OLD_R), list("gamma0", list(calibration_control = gaze_calibration_control(gamma_bounds = c(0, 0), typicality = "none"))), list("new", NEW))) {
      f <- rr(data, cfg[[2]], s); fits[[cfg[[1]]]] <- c(fits[[cfg[[1]]]], list(f)); strata(f, paste(cfg[[1]], "seed", s))
    }
  }
  for (k in names(fits)) { summ(fits[[k]], k); temps(fits[[k]], k) }
}
disp <- function(fits, label) {
  d <- do.call(rbind, lapply(fits, function(f) f$calibration$heldout_score_dispersion))
  if (is.null(d)) return(invisible())
  cat(sprintf("   %s held-out: between-candidate SD of mean raw=%.4f ranking=%.4f | within-candidate SD raw median=%.4f (max/min %.2f) ranking median=%.4f (max/min %.2f)\n", label,
    sd(d$raw_mean), sd(d$ranking_mean), median(d$raw_sd, na.rm = TRUE), max(d$raw_sd, na.rm = TRUE) / min(d$raw_sd, na.rm = TRUE),
    median(d$ranking_sd, na.rm = TRUE), max(d$ranking_sd, na.rm = TRUE) / min(d$ranking_sd, na.rm = TRUE)))
  typ <- do.call(rbind, lapply(fits, function(f) f$calibration$typicality))
  if (!is.null(typ)) cat(sprintf("   %s typicality sources: offset SD=%.4f, source SD median=%.4f max/min=%.2f, status=%s\n", label,
    sd(typ$offset), median(typ$source_sd), max(typ$source_sd) / min(typ$source_sd), paste(unique(typ$status), collapse = ",")))
}
typ_group <- function(fits, central, label) {
  s <- central_share(fits, central)
  cat(sprintf("%-20s central share per candidate=%.3f (chance %.3f, mc %.3f); central wins=%.3f expected %.3f\n", label, s$share_per_candidate, s$chance_per_candidate, s$mc_error, s$wins, s$expected))
}
if (job == "ttyp") {
  data <- lapply(1:3, function(s) simulate_calibration_transport(8, 6, recall = "centre", layout = "mixed", seed = 50 + s))
  for (cfg in CFGS_T) {
    fits <- lapply(1:3, function(s) rt(data[[s]], cfg[[2]], s))
    typ_group(fits, function(c) as.integer(sub(".*:", "", c$candidate_key)) <= 2, cfg[[1]]); summ(fits, cfg[[1]]); disp(fits, cfg[[1]])
  }
}
if (job == "rtyp") {
  data <- lapply(1:3, function(s) simulate_calibration_replay(6, 12, recall = "centre", layout = "mixed", seed = 60 + s))
  for (cfg in CFGS_R) {
    fits <- lapply(1:3, function(s) rr(data[[s]], cfg[[2]], s))
    typ_group(fits, function(c) c$image_id %% 3 == 1, cfg[[1]]); summ(fits, cfg[[1]]); disp(fits, cfg[[1]])
  }
}
if (job == "tdisc") {
  data <- lapply(1:2, function(s) simulate_calibration_transport(8, 6, keep = c(2, 8), seed = 70 + s))
  for (cfg in CFGS_T) {
    fits <- lapply(1:2, function(s) rt(data[[s]], cfg[[2]], s)); summ(fits, cfg[[1]]); disp(fits, cfg[[1]]) }
}
if (job == "rdisc") {
  data <- lapply(1:2, function(s) simulate_calibration_replay(6, 12, seed = 72 + s))
  for (cfg in CFGS_R) {
    fits <- lapply(1:2, function(s) rr(data[[s]], cfg[[2]], s)); summ(fits, cfg[[1]]); disp(fits, cfg[[1]]) }
}
