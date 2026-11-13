# Matching with column vectors avoids delimiter collisions in crossed keys.
# The one-column legacy path deliberately retains base match() semantics.
match_map_rows <- function(source, target, by, strict = FALSE) {
  if (strict) {
    check_key_columns(source, by)
    check_key_columns(target, by)
    if (anyDuplicated(target[by])) stop("Reference keys must be unique.")
  }
  if (length(by) == 1L) return(match(source[[by]], target[[by]]))
  ids <- key_group_ids(rbind(as.data.frame(target[by]), as.data.frame(source[by])), by)
  match(ids[nrow(target) + seq_len(nrow(source))], ids[seq_len(nrow(target))])
}

check_key_columns <- function(tab, keys) {
  if (!is.data.frame(tab) || !is.character(keys) || !length(keys) ||
      anyDuplicated(keys) || !all(keys %in% names(tab))) {
    stop("keys must name distinct columns in a data frame.")
  }
  if (anyNA(tab[keys])) stop("Key columns must not contain missing values.")
  if (any(vapply(tab[keys], is.list, logical(1)))) stop("Key columns must be atomic.")
}

# Factor levels are storage metadata, not part of a key's identity.
equal_key_values <- function(x, y) {
  if (is.factor(x)) x <- as.character(x)
  if (is.factor(y)) y <- as.character(y)
  x == y
}

key_group_ids <- function(tab, keys) {
  as.integer(dplyr::group_indices(dplyr::group_by(as.data.frame(tab),
                                               dplyr::across(dplyr::all_of(keys)))))
}

# Exclusion, deduplication and bounded sampling are shared by legacy baselines
# and new assignment constructors. Preserve candidate order and session RNG.
reference_candidates <- function(pool, exclude) unique(pool[!pool %in% exclude])

draw_candidates <- function(pool, n, replace = FALSE, cap = TRUE) {
  if (cap && !replace && n >= length(pool)) return(pool)
  if (!length(pool)) stop("No eligible controls or replacement templates.")
  pool[sample.int(length(pool), n, replace = replace)]
}

with_sampling_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  check_count(seed, "seed", zero = TRUE)
  restore <- snapshot_session_rng()
  on.exit(restore())
  set.seed(seed)
  force(code)
}

check_count <- function(x, name, zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x > .Machine$integer.max || x != trunc(x) || x < if (zero) 0 else 1) {
    stop(name, " must be a ", if (zero) "non-negative" else "positive", " integer.")
  }
}

#' Draw replayable case-control locations
#'
#' Creates explicit assignments. Scoring this table never samples again.
#' @param fix_tab A flat fixation table with finite x/y and case metadata.
#' @param n Number of controls per case.
#' @param subject,trial Columns identifying subjects and trials (trial is nested
#'   within subject).
#' @param stratify_by Additional columns restricting observed donor pools.
#' @param replace `"auto"` uses replacement only if an observed pool is too small;
#'   alternatively TRUE or FALSE. Map controls always sample independently with
#'   replacement; `replace = FALSE` is an error for that method.
#' @param seed Seed for assignment generation; restores the session RNG state.
#'   NULL explicitly uses the session RNG.
#' @param draws Number of independent control draws.
#' @param method Observed donor fixations, uniform coordinates, or grid locations
#'   sampled with probability proportional to a supplied map's nonnegative mass.
#' @param xbounds,ybounds Explicit coordinate bounds for uniform controls.
#' @param control_map A density map for `method = "map"`.
#' @param item Item column for observed donor exclusion. Defaults to `"item"`
#'   when that column exists, otherwise NULL. Pass another column name for
#'   custom item keys. NULL explicitly retains trial-only exclusion. Ignored
#'   for uniform/map controls.
#' @return An `eyesim_strata` data frame. Metadata belongs to the case; x/y are
#'   the case or control coordinates. `case_id` and `control_id` index input rows
#'   (control_id is NA for generated locations); `draw`, `stratum`, and `is_case`
#'   identify assignments. Observed controls always exclude the case's trial
#'   and, when `item` is supplied, all trials of the case's item.
#'   The `control_sources` attribute retains input atomic metadata, indexed by
#'   `control_id`, so donor items remain identifiable after subsetting cases.
#' @details There is no control-only jitter: perturbing controls but not cases
#'   changes their sampling distribution and can induce spurious preference.
#'   Uniform/map controls sample a different baseline from observed controls.
#'   Excluding the case's item avoids attenuating matched contrasts when items
#'   are presented repeatedly. Without item metadata only trial exclusion is
#'   possible. The resolved item column is recorded in the `sampling` attribute.
#' @md
#' @export
sample_controls <- function(fix_tab, n = 20, subject = "subject", trial = "trial",
                            stratify_by = NULL, replace = "auto", seed = 1, draws = 1,
                            method = c("observed", "uniform", "map"),
                            xbounds = NULL, ybounds = NULL, control_map = NULL,
                            item = if ("item" %in% names(fix_tab)) "item" else NULL) {
  method <- match.arg(method)
  check_count(n, "n")
  check_count(draws, "draws")
  check_key_columns(fix_tab, unique(c(subject, trial, stratify_by)))
  check_locations(fix_tab)
  reserved <- c("case_id", "control_id", "draw", "stratum", "is_case")
  if (any(reserved %in% names(fix_tab))) stop("Input contains reserved assignment column names.")
  if (!nrow(fix_tab)) stop("fix_tab must contain at least one case.")
  if (!identical(replace, "auto") && !identical(replace, TRUE) && !identical(replace, FALSE)) {
    stop("replace must be 'auto', TRUE, or FALSE.")
  }
  pool_ids <- key_group_ids(fix_tab, unique(c(subject, stratify_by)))
  trial_ids <- key_group_ids(fix_tab, unique(c(subject, trial)))
  if (method != "observed") item <- NULL
  item_ids <- integer(nrow(fix_tab))
  if (!is.null(item)) {
    if (!is.character(item) || length(item) != 1L || is.na(item)) {
      stop("item must name one column or be NULL for trial-only exclusion.")
    }
    check_key_columns(fix_tab, item)
    item_ids <- key_group_ids(fix_tab, unique(c(subject, item)))
  }
  pools <- split(seq_len(nrow(fix_tab)), pool_ids)
  # Cache eligible donors per pool/trial/item, not once per fixation.
  group <- key_group_ids(data.frame(pool = pool_ids, trial = trial_ids, item = item_ids),
                         c("pool", "trial", "item"))
  examples <- match(seq_len(max(group)), group)
  candidates <- lapply(examples, function(i) {
    pool <- pools[[as.character(pool_ids[i])]]
    exclude <- trial_ids[pool] == trial_ids[i] |
      (item_ids[pool] != 0L & item_ids[pool] == item_ids[i])
    reference_candidates(pool, pool[exclude])
  })
  if (method == "observed" && any(lengths(candidates) == 0L)) {
    stop("No eligible controls outside the case's trial",
         if (!is.null(item)) " and item", " in one or more pools.")
  }
  if (method == "uniform") {
    valid_bounds <- function(b) is.numeric(b) && length(b) == 2L &&
      all(is.finite(b)) && b[2] > b[1]
    if (!valid_bounds(xbounds) || !valid_bounds(ybounds)) stop("Uniform controls require increasing finite xbounds and ybounds.")
  }
  if (method == "map") {
    if (identical(replace, FALSE)) stop("Map controls require replacement; use replace = TRUE or 'auto'.")
    validate_density_grid(control_map)
    if (!inherits(control_map, "density") || any(!is.finite(control_map$z)) ||
        any(control_map$z < 0) || sum(control_map$z) <= 0) {
      stop("control_map must be a density with finite non-negative mass.")
    }
    cells <- which(control_map$z > 0)
  }
  use_replacement <- function(size) {
    yes <- if (identical(replace, "auto")) size < n else replace
    if (!yes && size < n) stop("Requested controls exceed pool size with replace = FALSE.")
    yes
  }
  out <- with_sampling_seed(seed, {
    assignments <- vector("list", nrow(fix_tab) * draws)
    for (d in seq_len(draws)) for (i in seq_len(nrow(fix_tab))) {
      donor <- rep(NA_integer_, n)
      if (method == "observed") {
        pool <- candidates[[group[i]]]
        donor <- draw_candidates(pool, n, use_replacement(length(pool)), cap = FALSE)
        cx <- fix_tab$x[donor]
        cy <- fix_tab$y[donor]
      } else if (method == "uniform") {
        cx <- stats::runif(n, xbounds[1], xbounds[2])
        cy <- stats::runif(n, ybounds[1], ybounds[2])
      } else {
        selected <- cells[sample.int(length(cells), n, replace = TRUE,
                                      prob = control_map$z[cells])]
        cx <- control_map$x[(selected - 1L) %% length(control_map$x) + 1L]
        cy <- control_map$y[(selected - 1L) %/% length(control_map$x) + 1L]
      }
      sid <- (d - 1L) * nrow(fix_tab) + i
      row <- as.data.frame(fix_tab[rep(i, n + 1L), , drop = FALSE])
      row$x <- c(fix_tab$x[i], cx)
      row$y <- c(fix_tab$y[i], cy)
      row$case_id <- i
      row$control_id <- c(i, donor)
      row$draw <- d
      row$stratum <- sid
      row$is_case <- c(TRUE, rep(FALSE, n))
      assignments[[sid]] <- row
    }
    dplyr::bind_rows(assignments)
  })
  attr(out, "sampling") <- list(subject = subject, trial = trial, stratify_by = stratify_by,
                                 method = method, replace = replace, seed = seed, n = n, draws = draws,
                                 item = item)
  if (method == "observed") {
    metadata <- setdiff(names(fix_tab)[vapply(fix_tab, is.atomic, logical(1))], c("x", "y"))
    attr(out, "control_sources") <- as.data.frame(fix_tab[metadata])
  }
  class(out) <- c("eyesim_strata", "data.frame")
  out
}

#' Draw replayable mismatched-template assignments
#' @param strata A case-control table from [sample_controls()].
#' @param refs A reference set.
#' @param n Number of null draws.
#' @param trial Trial column, nested within the reference subject column.
#' @param seed Seed, or NULL to use the session RNG.
#' @return A data frame containing original subject/item/trial keys, `null_draw`,
#'   and `replacement_item`. Each trial gets one other valid item from the same
#'   subject per draw. Assignments can be saved and passed to [mismatch_draws()].
#' @md
#' @export
mismatch_assignments <- function(strata, refs, n = 100, trial = "trial", seed = 1) {
  check_count(n, "n")
  if (!inherits(refs, "eyesim_references")) stop("refs must be a reference set.")
  subject <- refs$spec$subject
  item <- refs$spec$item
  keys <- unique(c(subject, trial, item))
  if (any(keys %in% c("replacement_item", "null_draw"))) stop("Key names collide with assignment columns.")
  check_key_columns(strata, keys)
  trials <- unique(as.data.frame(strata[keys]))
  if (anyDuplicated(trials[unique(c(subject, trial))])) stop("Each subject/trial must identify exactly one item.")
  original <- match_map_rows(trials, refs$keys, c(subject, item), strict = TRUE)
  if (anyNA(original)) {
    stop("Some trials have no matching reference keys.")
  }
  pools <- lapply(seq_len(nrow(trials)), function(i) {
    reference_candidates(which(refs$si == refs$si[original[i]] & refs$valid),
                         which(refs$ii == refs$ii[original[i]]))
  })
  if (any(lengths(pools) == 0L)) stop("No other valid item for one or more trials.")
  out <- with_sampling_seed(seed, {
    dplyr::bind_rows(lapply(seq_len(n), function(d) {
      rows <- vapply(pools, function(pool) draw_candidates(pool, 1L, cap = FALSE), integer(1))
      ret <- trials
      ret$null_draw <- d
      ret$replacement_item <- refs$keys[[item]][rows]
      ret
    }))
  })
  attr(out, "seed") <- seed
  out
}
