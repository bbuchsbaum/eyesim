# Reconstruct the public pcmri catalog and compare its estimand with the
# completed common-support GazeWeave recognition court.
#
# Protocol: inst/validation/GAZEWEAVE-PCMRI-CATALOG-REPLICATION.md
# Run from the eyesim source root after devtools::load_all(). Participant-level
# outputs belong only in the Git-ignored output directory.

pcmri_catalog_protocol_version <- 1L
pcmri_catalog_seed <- 20260820L
pcmri_catalog_sigma <- 80
pcmri_catalog_grid <- c(80L, 60L)
pcmri_catalog_output_dir <- file.path(
  "inst", "validation", "gaze-weave-pcmri-catalog-replication-results"
)

pcmri_catalog_full_file <- system.file(
  "validation", "gaze-weave-recognition-full-cohort.R", package = "eyesim"
)
pcmri_catalog_this_file <- tryCatch(
  normalizePath(sys.frame(1L)$ofile, mustWork = TRUE),
  error = function(error) ""
)
if (!nzchar(pcmri_catalog_full_file) && nzchar(pcmri_catalog_this_file)) {
  pcmri_catalog_full_file <- file.path(
    dirname(pcmri_catalog_this_file),
    "gaze-weave-recognition-full-cohort.R"
  )
}
if (!nzchar(pcmri_catalog_full_file)) {
  pcmri_catalog_full_file <- file.path(
    "inst", "validation", "gaze-weave-recognition-full-cohort.R"
  )
}
if (!file.exists(pcmri_catalog_full_file)) {
  stop("Run the pcmri catalog reconstruction from the eyesim source root.")
}
source(pcmri_catalog_full_file, local = TRUE)

pcmri_catalog_response <- function(response) {
  response <- suppressWarnings(as.integer(response))
  ifelse(response == 0L | is.na(response), NA_integer_,
         as.integer(response %in% c(1L, 2L)))
}

pcmri_catalog_pool_study <- function(study) {
  key <- interaction(
    study$participant, study$item, drop = TRUE, lex.order = TRUE
  )
  pooled <- lapply(split(seq_len(nrow(study)), key), function(index) {
    part <- study[index, , drop = FALSE]
    presentations <- sort(unique(part$presentation))
    if (!identical(presentations, 1:4)) return(NULL)
    fixation_rows <- do.call(rbind, part$fixgroup)
    fixation_rows$onset <- seq_len(nrow(fixation_rows)) - 1
    data.frame(
      participant = part$participant[[1L]],
      item = part$item[[1L]],
      presentations = length(presentations),
      nfix = nrow(fixation_rows),
      fixgroup = I(list(fixation_rows)),
      stringsAsFactors = FALSE
    )
  })
  tibble::as_tibble(dplyr::bind_rows(pooled))
}

pcmri_catalog_density_vectors <- function(paths, sigma = pcmri_catalog_sigma,
                                          grid = pcmri_catalog_grid) {
  vectors <- lapply(paths$fixgroup, function(path) {
    density <- eye_density(
      path, sigma = sigma,
      xbounds = c(0, 800), ybounds = c(0, 600), outdim = grid,
      duration_weighted = TRUE
    )
    if (is.null(density)) return(rep(NA_real_, prod(grid)))
    as.numeric(density$z)
  })
  matrix <- do.call(rbind, vectors)
  means <- rowMeans(matrix)
  centered <- matrix - means
  norms <- sqrt(rowSums(centered^2))
  valid <- is.finite(norms) & norms > 0
  centered[valid, ] <- centered[valid, , drop = FALSE] / norms[valid]
  centered[!valid, ] <- NA_real_
  list(raw = matrix, pearson = centered, valid = valid)
}

pcmri_catalog_pearson <- function(reference, source, reference_vectors,
                                  source_vectors) {
  if (nrow(reference) != nrow(reference_vectors) ||
      nrow(source) != nrow(source_vectors)) {
    stop("Path tables and density matrices have inconsistent row counts.")
  }
  rows <- lapply(sort(unique(source$participant)), function(participant) {
    ri <- which(reference$participant == participant)
    si <- which(source$participant == participant)
    if (!length(ri) || !length(si)) return(NULL)
    similarity <- source_vectors[si, , drop = FALSE] %*%
      t(reference_vectors[ri, , drop = FALSE])
    match_index <- match(source$item[si], reference$item[ri])
    keep <- which(!is.na(match_index))
    if (!length(keep)) return(NULL)
    matched <- similarity[cbind(keep, match_index[keep])]
    permuted <- vapply(keep, function(row_index) {
      alternatives <- seq_along(ri) != match_index[[row_index]]
      mean(similarity[row_index, alternatives], na.rm = TRUE)
    }, numeric(1))
    part <- source[si[keep], , drop = FALSE]
    part$sim_self <- as.numeric(matched)
    part$sim_self_perm <- permuted
    part$sim_self_diff <- part$sim_self - part$sim_self_perm
    part$n_perm <- length(ri) - 1L
    part
  })
  tibble::as_tibble(dplyr::bind_rows(rows))
}

pcmri_catalog_density_reconstruction <- function(verify = TRUE) {
  raw <- full_recognition_read_inputs(verify = verify)
  reference <- pcmri_catalog_pool_study(raw$study)
  source <- raw$retrieval[
    raw$retrieval$probe_type %in% c("old", "lure"), , drop = FALSE
  ]
  retrieval_meta <- raw$trials$retrieval[c(
    "Subject", "ImageNumber", "Response"
  )]
  names(retrieval_meta) <- c("participant", "item", "response")
  retrieval_meta$participant <- as.character(retrieval_meta$participant)
  retrieval_meta$item <- as.integer(retrieval_meta$item)
  retrieval_meta <- unique(retrieval_meta)
  source <- merge(
    source, retrieval_meta, by = c("participant", "item"),
    all.x = TRUE, sort = FALSE
  )
  source <- source[order(source$participant, source$item), ]
  reference <- reference[order(reference$participant, reference$item), ]
  reference_density <- pcmri_catalog_density_vectors(reference)
  source_density <- pcmri_catalog_density_vectors(source)
  scores <- pcmri_catalog_pearson(
    reference, source, reference_density$pearson, source_density$pearson
  )
  scores$condition <- factor(
    scores$probe_type, levels = c("old", "lure")
  )
  scores$sal10 <- (as.numeric(scores$degradation) - 60) / 10
  scores$said_old <- pcmri_catalog_response(scores$response)
  scores
}

pcmri_catalog_fixed_effect <- function(model, term, analysis, method) {
  coefficients <- as.data.frame(summary(model)$coefficients)
  if (!term %in% rownames(coefficients)) {
    stop("Term was not present in fitted model: ", term)
  }
  row <- coefficients[term, , drop = FALSE]
  p_column <- grep("^Pr\\(", names(row), value = TRUE)
  data.frame(
    analysis = analysis,
    method = method,
    term = term,
    estimate = row[["Estimate"]],
    std_error = row[["Std. Error"]],
    statistic = row[[grep("value$", names(row), value = TRUE)[[1L]]]],
    p_value = if (length(p_column)) row[[p_column[[1L]]]] else NA_real_,
    n = stats::nobs(model),
    singular = lme4::isSingular(model),
    stringsAsFactors = FALSE
  )
}

pcmri_catalog_models <- function(scores) {
  scores$participant <- factor(scores$participant)
  scores$item <- factor(scores$item)
  old <- scores[
    scores$condition == "old" & !is.na(scores$said_old) &
      is.finite(scores$sim_self_diff), , drop = FALSE
  ]
  old$z_metric <- as.numeric(scale(old$sim_self_diff))
  above <- lmerTest::lmer(
    sim_self_diff ~ 1 + (1 | participant) + (1 | item), data = scores
  )
  condition <- lmerTest::lmer(
    sim_self_diff ~ condition + (1 | participant) + (1 | item),
    data = scores
  )
  saliency <- lmerTest::lmer(
    sim_self_diff ~ condition * sal10 +
      (1 | participant) + (1 | item), data = scores
  )
  behavior <- lme4::glmer(
    said_old ~ z_metric + sal10 +
      (1 | participant) + (1 | item),
    family = stats::binomial(), data = old,
    control = lme4::glmerControl(optimizer = "bobyqa")
  )
  dplyr::bind_rows(
    pcmri_catalog_fixed_effect(
      above, "(Intercept)", "mean_reinstatement", "density_sigma_80_all4"
    ),
    pcmri_catalog_fixed_effect(
      condition, "conditionlure", "lure_minus_old", "density_sigma_80_all4"
    ),
    pcmri_catalog_fixed_effect(
      saliency, "sal10", "old_saliency_per10", "density_sigma_80_all4"
    ),
    pcmri_catalog_fixed_effect(
      saliency, "conditionlure:sal10", "lure_saliency_interaction",
      "density_sigma_80_all4"
    ),
    pcmri_catalog_fixed_effect(
      behavior, "z_metric", "old_said_old", "density_sigma_80_all4"
    )
  )
}

pcmri_catalog_common_support <- function(
    result_dir = file.path(
      "inst", "validation", "gaze-weave-recognition-full-cohort-results"
    )) {
  trial_path <- file.path(result_dir, "trial-scores.csv")
  candidate_path <- file.path(result_dir, "candidate-scores.csv")
  if (!file.exists(trial_path) || !file.exists(candidate_path)) {
    stop("Completed GW-13 trial and candidate scores were not found.")
  }
  trials <- utils::read.csv(trial_path, stringsAsFactors = FALSE)
  candidates <- utils::read.csv(candidate_path, stringsAsFactors = FALSE)
  response <- utils::read.csv(
    probe_delay_data_files()[["retrieval"]], stringsAsFactors = FALSE,
    check.names = FALSE
  )[c("Subject", "ImageNumber", "Response")]
  names(response) <- c("participant", "item", "response")
  response <- unique(response)
  trials <- merge(
    trials, response, by = c("participant", "item"),
    all.x = TRUE, sort = FALSE
  )
  trials$said_old <- pcmri_catalog_response(trials$response)
  trials$sal10 <- (trials$degradation - 60) / 10
  catalog <- candidates |>
    dplyr::group_by(participant, target_item, method) |>
    dplyr::summarise(
      catalog_diff = mean(log_score[is_true]) - mean(log_score[!is_true]),
      .groups = "drop"
    )
  names(catalog)[names(catalog) == "target_item"] <- "item"
  trials <- merge(
    trials, catalog, by = c("participant", "item", "method"),
    all.x = TRUE, sort = FALSE
  )
  trials
}

pcmri_catalog_compare_methods <- function(trials) {
  methods <- c(
    "replay", "transport_v2", "density_sigma_80_raw_calibrated",
    "density_ridge_registered", "multimatch_ridge_registered",
    "multimatch_mm_position_raw_calibrated"
  )
  rows <- lapply(methods, function(method) {
    part <- trials[
      trials$method == method & trials$probe_type == "old" &
        !is.na(trials$said_old) & is.finite(trials$gaze_info_bits),
      , drop = FALSE
    ]
    part$participant <- factor(part$participant)
    part$item <- factor(part$item)
    part$z_metric <- as.numeric(scale(part$gaze_info_bits))
    model <- lme4::glmer(
      said_old ~ z_metric + sal10 +
        (1 | participant) + (1 | item),
      family = stats::binomial(), data = part,
      control = lme4::glmerControl(optimizer = "bobyqa")
    )
    pcmri_catalog_fixed_effect(
      model, "z_metric", "old_said_old_information", method
    )
  })
  dplyr::bind_rows(rows)
}

pcmri_catalog_run <- function(
    output_dir = pcmri_catalog_output_dir, verify = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  density_scores <- pcmri_catalog_density_reconstruction(verify = verify)
  density_models <- pcmri_catalog_models(density_scores)
  common_scores <- pcmri_catalog_common_support()
  common_models <- pcmri_catalog_compare_methods(common_scores)
  aggregate <- common_scores |>
    dplyr::filter(is.finite(gaze_info_bits)) |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      n = dplyr::n(),
      participants = dplyr::n_distinct(participant),
      mean_gaze_info_bits = mean(gaze_info_bits),
      mean_rank = mean(template_rank),
      top1_credit = mean(top1_credit),
      mean_catalog_diff = mean(catalog_diff),
      .groups = "drop"
    )
  density_export <- density_scores[
    !vapply(density_scores, is.list, logical(1))
  ]
  utils::write.csv(
    density_export, file.path(output_dir, "density-trial-scores.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    density_models, file.path(output_dir, "density-models.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    common_models, file.path(output_dir, "common-support-behavior-models.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    aggregate, file.path(output_dir, "common-support-summary.csv"),
    row.names = FALSE
  )
  invisible(list(
    density_scores = density_scores,
    density_models = density_models,
    common_scores = common_scores,
    common_models = common_models,
    aggregate = aggregate
  ))
}
