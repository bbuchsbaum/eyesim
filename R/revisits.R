#' Revisitation of previously fixated locations
#'
#' Scores, fixation by fixation, whether gaze returns to a location already
#' fixated earlier in the same viewing episode, and summarises each episode as
#' counts, rates and, optionally, rates expected under a spatial reference
#' distribution. This is a spatial-return measure related to Kragel et al.
#' (2021), not the backtracking-sequence classification of Voss et al. (2011),
#' which classified A-B-C-A as "other" rather than revisitation.
#'
#' Fixation \eqn{i} is a \emph{revisit} when it lies within \code{radius} of at
#' least one earlier fixation \eqn{j} that is at least \code{min_lag} fixations
#' back (and at most \code{max_lag}). The default \code{min_lag = 2} means a
#' fixation never revisits the fixation immediately before it. With
#' \code{require_departure = TRUE} (default), the incoming step must be at least
#' \code{radius}: the preceding fixation must lie at least \code{radius} away
#' from fixation \eqn{i}. This does not require leaving and re-entering a region
#' centred on the earlier fixation. For example, with radius 10, positions
#' 0, -9, 9 on one axis count as a return, although both later positions are
#' within radius of the original position.
#' A fixation is an \emph{opportunity} when it could have been a revisit:
#' it has an eligible earlier fixation and, under \code{require_departure},
#' arrives from at least \code{radius} away.
#'
#' Revisit counts rise mechanically when an observer makes many fixations within
#' a small area, so raw counts partly restate fixation number and spatial
#' dispersion. Supplying \code{reference} locations gives a spatial null: for each
#' opportunity, the expected probability of a revisit is the proportion of
#' reference locations that would count as a revisit of the same eligible earlier
#' fixations, among those reference locations that also satisfy the departure
#' condition. A natural reference is the observer's fixations from other trials,
#' which preserves the observer's spatial bias but not the order within the
#' episode (see the \code{eye_table} method, \code{null = "subject"}).
#' Observed minus expected is a descriptive correction under this reference
#' distribution, not an inferential test or a generative model of scanpaths.
#'
#' Coordinates and onsets must be finite numeric values. Invalid values cause
#' an error; dropping fixations would change adjacency, lags and denominators.
#' Fixations are scored in onset order, preserving input order for tied onsets.
#' Lag bounds are whole numbers; large finite bounds are supported without
#' conversion to 32-bit integers.
#'
#' Reference coverage is updated as fixations enter and leave the eligible lag
#' window. For \eqn{n} fixations, \eqn{m} reference points and maximum eligible
#' history size \eqn{L}, scoring uses \eqn{O(n + nL + nm)} time and \eqn{O(n + m)}
#' storage, including result rows. No pairwise distance matrix is allocated.
#'
#' \code{radius} is in the units of the fixation coordinates. Kragel et al.
#' (2021) used about 2 degrees of visual angle; convert degrees to pixels for
#' your display geometry before calling.
#'
#' @param x A \code{fixation_group} or an \code{eye_table}.
#' @param radius Distance within which two fixations count as the same location.
#' @param min_lag Minimum number of fixations back for an eligible earlier
#'   fixation (default 2).
#' @param max_lag Maximum number of fixations back, a whole number at least
#'   \code{min_lag}, or \code{Inf} (default) for no limit.
#' @param require_departure Logical; require the preceding fixation to lie at
#'   least \code{radius} from the current fixation (default \code{TRUE}).
#' @param reference Optional two-column matrix or data frame of \code{x}, \code{y}
#'   finite numeric locations used for the spatial reference distribution.
#' @param dist_method \code{"euclidean"} (default) or \code{"manhattan"}.
#' @param ... Further arguments passed to methods.
#'
#' @return For a \code{fixation_group}, a tibble with one row per fixation in
#'   onset order: \code{index}, \code{onset}, \code{opportunity},
#'   \code{revisit}, \code{lag} (fixations back to the most recent matched
#'   fixation), \code{match_index} (that fixation's original index),
#'   \code{nearest_distance} (to the nearest eligible earlier fixation),
#'   \code{match_distance} (to \code{match_index}, which need not be nearest),
#'   \code{distance} (a compatibility alias for \code{nearest_distance}) and
#'   \code{expected} (null probability of a revisit, or \code{NA} without
#'   reference support). Distances are \code{NA} outside opportunities;
#'   \code{match_distance} is also \code{NA} without a match. Attributes retain
#'   \code{radius}, \code{min_lag}, \code{max_lag}, \code{require_departure},
#'   \code{dist_method} and \code{null} (whether a reference was supplied).
#'   Use \code{revisit_summary()} for one row per episode. For an
#'   \code{eye_table}, see \code{fixation_revisits.eye_table}.
#'
#' @references
#' Voss, J. L., Warren, D. E., Gonsalves, B. D., Federmeier, K. D., Tranel, D.,
#' & Cohen, N. J. (2011). Spontaneous revisitation during visual exploration as a
#' link among strategic behavior, learning, and the hippocampus. \emph{PNAS},
#' 108(31), E402-E409.
#'
#' Kragel, J. E., Schuele, S., VanHaerents, S., Rosenow, J. M., & Voss, J. L.
#' (2021). Rapid coordination of effective learning by the human hippocampus.
#' \emph{Science Advances}, 7(25), eabf7144.
#'
#' @examples
#' # A-B-C-A: the fourth fixation returns to the first location
#' fg <- fixation_group(x = c(100, 400, 400, 105), y = c(100, 100, 400, 98),
#'                      onset = c(0, 250, 500, 750), duration = rep(200, 4))
#' rv <- fixation_revisits(fg, radius = 50)
#' rv
#' revisit_summary(rv)
#' @export
fixation_revisits <- function(x, ...) {
  UseMethod("fixation_revisits", x)
}

#' @rdname fixation_revisits
#' @export
fixation_revisits.fixation_group <- function(x, radius, min_lag = 2, max_lag = Inf,
                                             require_departure = TRUE, reference = NULL,
                                             dist_method = c("euclidean", "manhattan"), ...) {
  dist_method <- match.arg(dist_method)
  check_revisit_args(radius, min_lag, max_lag, require_departure)

  for (col in c("x", "y", "onset")) {
    if (!is.numeric(x[[col]]) || is.complex(x[[col]]) || any(!is.finite(x[[col]]))) {
      stop("Fixation column '", col, "' must contain only finite numeric values; ",
           "invalid fixations cannot be dropped without changing adjacency and lags.")
    }
  }
  ref <- if (is.null(reference)) NULL else as_reference_matrix(reference)
  x <- x[order(x$onset, method = "radix"), , drop = FALSE]
  n <- nrow(x)
  opportunity <- rep(FALSE, n)
  revisit <- rep(FALSE, n)
  lag <- rep(NA_integer_, n)
  match_index <- rep(NA_integer_, n)
  nearest_distance <- rep(NA_real_, n)
  match_distance <- rep(NA_real_, n)
  expected <- rep(NA_real_, n)

  xy <- cbind(x$x, x$y)
  storage.mode(xy) <- "double"
  # Each count records how many currently eligible fixations cover a reference
  # point. Update even when the current fixation is not an opportunity.
  covered <- integer(if (is.null(ref)) 0L else nrow(ref))

  for (i in seq_len(n)) {
    if (length(covered)) {
      entering <- i - min_lag
      leaving <- i - max_lag - 1
      if (entering >= 1) {
        covered <- covered + (revisit_dist_to(ref, xy[entering, ], dist_method) < radius)
      }
      if (leaving >= 1) {
        covered <- covered - (revisit_dist_to(ref, xy[leaving, ], dist_method) < radius)
      }
    }
    elig <- revisit_eligible(i, min_lag, max_lag)
    if (!length(elig)) next
    if (require_departure &&
        revisit_dist_to(xy[i - 1L, , drop = FALSE], xy[i, ], dist_method) < radius) next
    opportunity[i] <- TRUE
    distances <- revisit_dist_to(xy[elig, , drop = FALSE], xy[i, ], dist_method)
    nearest_distance[i] <- min(distances)
    near <- which(distances < radius)
    if (length(near)) {
      recent <- near[length(near)]
      j <- elig[recent]
      revisit[i] <- TRUE
      match_index[i] <- x$index[j]
      lag[i] <- as.integer(i - j)
      match_distance[i] <- distances[recent]
    }
    if (length(covered)) {
      ok <- if (require_departure) {
        revisit_dist_to(ref, xy[i - 1L, ], dist_method) >= radius
      } else {
        rep(TRUE, length(covered))
      }
      if (any(ok)) {
        expected[i] <- mean(covered[ok] > 0L)
      }
    }
  }

  out <- tibble::tibble(index = x$index, onset = x$onset, opportunity = opportunity,
                        revisit = revisit, lag = lag, match_index = match_index,
                        distance = nearest_distance, nearest_distance = nearest_distance,
                        match_distance = match_distance, expected = expected)
  new_revisit_table(out, radius, min_lag, max_lag, require_departure,
                    !is.null(reference), dist_method)
}

#' Revisits for every episode of an eye_table
#'
#' Applies \code{\link{fixation_revisits}} to each row (episode) of an
#' \code{eye_table} and appends per-episode summary columns.
#'
#' @param x An \code{eye_table}.
#' @param radius,min_lag,max_lag,require_departure,dist_method See
#'   \code{\link{fixation_revisits}}.
#' @param null The spatial null model:
#'   \describe{
#'     \item{\code{"none"}}{(default) no null; expected columns are \code{NA}.}
#'     \item{\code{"subject"}}{every fixation of the same subject from the
#'       \emph{other} rows of the table (leave-episode-out). Preserves the observer's
#'       spatial bias but not the content of the viewed image.}
#'     \item{\code{"item"}}{every fixation from other rows showing the same item
#'       (\code{item_var}). Preserves where the image draws gaze, but not the
#'       observer's own bias. If \code{subject_var} is given, all of the same subject's
#'       rows are excluded, which matters when an observer sees an item more than once.}
#'   }
#' @param subject_var Column identifying the subject; required when
#'   \code{null = "subject"}, optional for \code{null = "item"}.
#' @param item_var Column identifying the viewed item (image); required when
#'   \code{null = "item"}.
#' @param var Name of the fixation-group column (default \code{"fixgroup"}).
#' @param keep_fixations Logical; keep the per-fixation tables in a list column
#'   \code{revisit_fixations} (default \code{FALSE}).
#' @param ... Ignored.
#'
#' @details Re-scoring a table that \code{fixation_revisits()} already scored (for example at another radius)
#'   replaces the columns the earlier call wrote. A column of the same name that the caller added is never
#'   overwritten: it is an error, so rename it first.
#'
#' @return \code{x} with columns from \code{\link{revisit_summary}}
#'   (\code{n_fix}, \code{opportunities}, \code{revisits}, \code{revisit_rate},
#'   \code{expected_revisits}, \code{expected_rate}, \code{excess_revisits},
#'   \code{excess_rate}). The \code{revisit_spec} attribute retains
#'   \code{radius}, \code{min_lag}, \code{max_lag}, \code{require_departure},
#'   \code{dist_method}, \code{null}, \code{subject_var}, \code{item_var} and
#'   \code{var}, including when \code{keep_fixations = FALSE} or the table has
#'   zero rows. Re-scoring replaces this specification with the new arguments.
#' @examples
#' set.seed(1)
#' et <- simulate_eye_table(n_fixations = 120, n_groups = 6)
#' et$subject <- rep(c("s1", "s2"), each = 3)
#' et$image <- rep(c("a", "b", "c"), times = 2)
#' res <- fixation_revisits(et, radius = 150, null = "subject", subject_var = "subject")
#' tibble::as_tibble(res)[, c("subject", "n_fix", "revisits", "expected_revisits", "excess_revisits")]
#' # same-image null from the other subject's viewing
#' res_item <- fixation_revisits(et, radius = 150, null = "item", item_var = "image", subject_var = "subject")
#' tibble::as_tibble(res_item)[, c("subject", "image", "revisits", "expected_revisits")]
#' @export
fixation_revisits.eye_table <- function(x, radius, min_lag = 2, max_lag = Inf,
                                        require_departure = TRUE,
                                        null = c("none", "subject", "item"), subject_var = NULL,
                                        item_var = NULL, var = "fixgroup", keep_fixations = FALSE,
                                        dist_method = c("euclidean", "manhattan"), ...) {
  null <- match.arg(null)
  dist_method <- match.arg(dist_method)
  check_revisit_args(radius, min_lag, max_lag, require_departure)
  if (!is.logical(keep_fixations) || length(keep_fixations) != 1L || is.na(keep_fixations)) {
    stop("keep_fixations must be TRUE or FALSE.")
  }
  if (!var %in% names(x)) stop("Column '", var, "' not found in the eye_table.")
  has_col <- function(v) !is.null(v) && length(v) == 1L && v %in% names(x)
  if (null == "subject" && !has_col(subject_var)) {
    stop("null = \"subject\" needs subject_var naming a column of the eye_table.")
  }
  if (null == "item" && !has_col(item_var)) {
    stop("null = \"item\" needs item_var naming a column of the eye_table.")
  }
  if (null == "item" && !is.null(subject_var) && !has_col(subject_var)) {
    stop("subject_var '", subject_var, "' is not a column of the eye_table.")
  }
  # Columns written by an earlier call are recorded in an attribute and replaced; any other column with one of
  # these names belongs to the caller and is never overwritten.
  ours <- c(revisit_summary_cols, "revisit_fixations")
  prev <- attr(x, "revisit_columns") %||% character(0)
  clash <- setdiff(intersect(names(x), ours), prev)
  if (length(clash)) {
    stop("The eye_table already has column(s) ", paste(clash, collapse = ", "), " that fixation_revisits() writes. ",
         "Rename or drop them first.")
  }
  fgs <- x[[var]]
  xy_of <- function(fg) if (is.null(fg) || !nrow(fg)) matrix(numeric(0), 0, 2) else cbind(fg$x, fg$y)

  per <- lapply(seq_along(fgs), function(k) {
    fg <- fgs[[k]]
    if (is.null(fg)) fg <- fixation_group(numeric(0), numeric(0), numeric(0), numeric(0))
    ref <- NULL
    if (null != "none") {
      pool <- if (null == "subject") {
        setdiff(which(x[[subject_var]] == x[[subject_var]][k]), k)
      } else {
        same_item <- which(x[[item_var]] == x[[item_var]][k])
        if (is.null(subject_var)) setdiff(same_item, k) else same_item[x[[subject_var]][same_item] != x[[subject_var]][k]]
      }
      ref <- do.call(rbind, lapply(fgs[pool], xy_of))
      if (is.null(ref)) ref <- matrix(numeric(0), 0, 2)
    }
    fixation_revisits(fg, radius = radius, min_lag = min_lag, max_lag = max_lag,
                      require_departure = require_departure, reference = ref,
                      dist_method = dist_method)
  })
  summ <- if (length(per)) {
    dplyr::bind_rows(lapply(per, revisit_summary))
  } else {
    tibble::tibble(n_fix = integer(), opportunities = integer(), revisits = integer(),
                   revisit_rate = double(), expected_revisits = double(),
                   expected_rate = double(), excess_revisits = double(), excess_rate = double())
  }
  stopifnot(identical(names(summ), revisit_summary_cols))
  # Re-scoring a scored table (e.g. at another radius) replaces the columns the earlier call wrote.
  extra <- attributes(x)[setdiff(names(attributes(x)), c("names", "row.names", "class", "revisit_columns"))]
  old <- intersect(names(x), prev)
  if (length(old)) x <- x[, setdiff(names(x), old), drop = FALSE]
  out <- dplyr::bind_cols(x, summ)
  if (keep_fixations) out$revisit_fixations <- per
  for (a in names(extra)) attr(out, a) <- extra[[a]]
  attr(out, "revisit_columns") <- c(revisit_summary_cols, if (keep_fixations) "revisit_fixations")
  attr(out, "revisit_spec") <- list(radius = radius, min_lag = min_lag, max_lag = max_lag,
                                    require_departure = require_departure, dist_method = dist_method,
                                    null = null, subject_var = subject_var, item_var = item_var, var = var)
  class(out) <- class(x)
  out
}

revisit_summary_cols <- c("n_fix", "opportunities", "revisits", "revisit_rate", "expected_revisits",
                          "expected_rate", "excess_revisits", "excess_rate")

#' Summarise revisits for one episode
#'
#' @param x A table returned by \code{fixation_revisits()} for a
#'   \code{fixation_group}.
#'
#' @return A one-row tibble: \code{n_fix}; \code{opportunities}; \code{revisits};
#'   \code{revisit_rate} (revisits per opportunity); \code{expected_revisits} and
#'   \code{expected_rate} under the null (\code{NA} without a reference, or when
#'   any opportunity has no expectation because no reference location satisfies
#'   its departure condition);
#'   \code{excess_revisits} (observed minus expected) and \code{excess_rate}.
#' @export
revisit_summary <- function(x) {
  if (!inherits(x, "revisit_table")) stop("revisit_summary() needs the output of fixation_revisits() for a fixation_group.")
  opp <- sum(x$opportunity)
  rev <- sum(x$revisit)
  has_null <- isTRUE(attr(x, "null"))
  e <- x$expected[x$opportunity]
  exp_n <- if (has_null && opp > 0 && all(is.finite(e))) sum(e) else NA_real_
  rate <- if (opp > 0) rev / opp else NA_real_
  exp_rate <- if (has_null && opp > 0) exp_n / opp else NA_real_
  tibble::tibble(n_fix = nrow(x), opportunities = opp, revisits = rev, revisit_rate = rate,
                 expected_revisits = exp_n, expected_rate = exp_rate,
                 excess_revisits = rev - exp_n, excess_rate = rate - exp_rate)
}

#' @export
print.revisit_table <- function(x, ...) {
  cat(sprintf("Revisits: %d of %d opportunities (%d fixations; radius %g, lag %s-%s%s)\n",
              sum(x$revisit), sum(x$opportunity), nrow(x), attr(x, "radius"), attr(x, "min_lag"),
              format(attr(x, "max_lag")), if (isTRUE(attr(x, "require_departure"))) ", departure required" else ""))
  print(tibble::as_tibble(unclass_revisit(x)), ...)
  invisible(x)
}

new_revisit_table <- function(out, radius, min_lag, max_lag, require_departure, has_null, dist_method) {
  attr(out, "radius") <- radius
  attr(out, "min_lag") <- min_lag
  attr(out, "max_lag") <- max_lag
  attr(out, "require_departure") <- require_departure
  attr(out, "null") <- has_null
  attr(out, "dist_method") <- dist_method
  class(out) <- c("revisit_table", class(out))
  out
}

unclass_revisit <- function(x) {
  class(x) <- setdiff(class(x), "revisit_table")
  x
}

check_revisit_args <- function(radius, min_lag, max_lag, require_departure) {
  scalar_number <- function(x) is.numeric(x) && !is.complex(x) && length(x) == 1L && !is.na(x)
  if (missing(radius) || !scalar_number(radius) || !is.finite(radius) || radius <= 0) {
    stop("radius must be a single positive number (in fixation coordinate units).")
  }
  if (!scalar_number(min_lag) || !is.finite(min_lag) || min_lag < 1 || min_lag != floor(min_lag)) {
    stop("min_lag must be a positive integer.")
  }
  if (!scalar_number(max_lag) || max_lag < min_lag || max_lag != floor(max_lag)) {
    stop("max_lag must be a whole number >= min_lag (Inf for no limit).")
  }
  if (!is.logical(require_departure) || length(require_departure) != 1L || is.na(require_departure)) {
    stop("require_departure must be TRUE or FALSE.")
  }
}

revisit_eligible <- function(i, min_lag, max_lag) {
  hi <- i - min_lag
  lo <- max(1, i - max_lag)
  if (hi < lo) integer(0) else seq.int(lo, hi)
}

revisit_dist_to <- function(xy, point, method) {
  dx <- abs(xy[, 1] - point[1])
  dy <- abs(xy[, 2] - point[2])
  if (method == "euclidean") {
    # Scale before squaring to avoid overflow/underflow of finite distances.
    scale <- pmax(dx, dy)
    nonzero <- is.finite(scale) & scale > 0
    scale[nonzero] <- scale[nonzero] * sqrt((dx[nonzero] / scale[nonzero])^2 +
                                           (dy[nonzero] / scale[nonzero])^2)
    scale
  } else {
    dx + dy
  }
}

as_reference_matrix <- function(reference) {
  m <- if (is.data.frame(reference) && all(c("x", "y") %in% names(reference))) {
    cbind(reference$x, reference$y)
  } else {
    as.matrix(reference)
  }
  if (ncol(m) != 2L) stop("reference must have two columns (x, y).")
  if (!is.numeric(m) || is.complex(m) || any(!is.finite(m))) {
    stop("reference coordinates must contain only finite numeric values.")
  }
  storage.mode(m) <- "double"
  m
}
