#' Revisitation of previously fixated locations
#'
#' Scores, fixation by fixation, whether gaze returns to a location already
#' fixated earlier in the same viewing episode ("revisitation"; Voss et al.,
#' 2011; Kragel et al., 2021), and summarises each episode as counts, rates and,
#' optionally, rates expected under a spatial null model.
#'
#' Fixation \eqn{i} is a \emph{revisit} when it lies within \code{radius} of at
#' least one earlier fixation \eqn{j} that is at least \code{min_lag} fixations
#' back (and at most \code{max_lag}). The default \code{min_lag = 2} means a
#' fixation never revisits the fixation immediately before it. With
#' \code{require_departure = TRUE} (default) gaze must also have left: the
#' preceding fixation must lie at least \code{radius} away from fixation
#' \eqn{i}, so a run of small refixations on one spot is not counted as a
#' return. A fixation is an \emph{opportunity} when it could have been a revisit:
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
#'
#' \code{radius} is in the units of the fixation coordinates. Kragel et al.
#' (2021) used about 2 degrees of visual angle; convert degrees to pixels for
#' your display geometry before calling.
#'
#' @param x A \code{fixation_group} or an \code{eye_table}.
#' @param radius Distance within which two fixations count as the same location.
#' @param min_lag Minimum number of fixations back for an eligible earlier
#'   fixation (default 2).
#' @param max_lag Maximum number of fixations back (default \code{Inf}).
#' @param require_departure Logical; require the preceding fixation to lie at
#'   least \code{radius} away (default \code{TRUE}).
#' @param reference Optional two-column matrix or data frame of \code{x}, \code{y}
#'   locations used for the spatial null model.
#' @param dist_method \code{"euclidean"} (default) or \code{"manhattan"}.
#' @param ... Further arguments passed to methods.
#'
#' @return For a \code{fixation_group}, a tibble with one row per fixation in
#'   onset order: \code{index}, \code{onset}, \code{opportunity},
#'   \code{revisit}, \code{lag} (fixations back to the most recent matched
#'   fixation), \code{match_index}, \code{distance} (to the nearest eligible
#'   earlier fixation) and \code{expected} (null probability of a revisit, or
#'   \code{NA} without \code{reference}). Use \code{revisit_summary()} for one row
#'   per episode. For an \code{eye_table}, see \code{fixation_revisits.eye_table}.
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
  check_revisit_args(radius, min_lag, max_lag)

  keep <- is.finite(x$x) & is.finite(x$y) & !is.na(x$onset)
  x <- x[keep, , drop = FALSE]
  x <- x[order(x$onset), , drop = FALSE]
  n <- nrow(x)
  out <- tibble::tibble(index = x$index, onset = x$onset,
                        opportunity = rep(FALSE, n), revisit = rep(FALSE, n),
                        lag = rep(NA_integer_, n), match_index = rep(NA_integer_, n),
                        distance = rep(NA_real_, n), expected = rep(NA_real_, n))
  if (n == 0L) {
    return(new_revisit_table(out, radius, min_lag, max_lag, require_departure, !is.null(reference)))
  }

  xy <- cbind(x$x, x$y)
  d <- revisit_dist(xy, xy, dist_method)
  ref <- if (is.null(reference)) NULL else as_reference_matrix(reference)
  dref <- if (is.null(ref) || nrow(ref) == 0L) NULL else revisit_dist(ref, xy, dist_method)

  for (i in seq_len(n)) {
    elig <- revisit_eligible(i, min_lag, max_lag)
    if (!length(elig)) next
    departed <- i == 1L || d[i, i - 1L] >= radius
    if (require_departure && !departed) next
    out$opportunity[i] <- TRUE
    out$distance[i] <- min(d[i, elig])
    near <- elig[d[i, elig] < radius]
    if (length(near)) {
      j <- max(near)
      out$revisit[i] <- TRUE
      out$match_index[i] <- out$index[j]
      out$lag[i] <- as.integer(i - j)
    }
    if (!is.null(dref)) {
      ok <- if (require_departure && i > 1L) dref[, i - 1L] >= radius else rep(TRUE, nrow(dref))
      if (any(ok)) {
        hit <- apply(dref[ok, elig, drop = FALSE] < radius, 1, any)
        out$expected[i] <- mean(hit)
      }
    }
  }

  new_revisit_table(out, radius, min_lag, max_lag, require_departure, !is.null(reference))
}

#' Revisits for every episode of an eye_table
#'
#' Applies \code{\link{fixation_revisits}} to each row (episode) of an
#' \code{eye_table} and appends per-episode summary columns.
#'
#' @param x An \code{eye_table}.
#' @param radius,min_lag,max_lag,require_departure,dist_method See
#'   \code{\link{fixation_revisits}}.
#' @param null \code{"none"} (default) or \code{"subject"}: the spatial null uses
#'   every fixation of the same subject from the \emph{other} rows of the table
#'   (leave-episode-out), preserving the observer's spatial bias.
#' @param subject_var Column identifying the subject; required when
#'   \code{null = "subject"}.
#' @param var Name of the fixation-group column (default \code{"fixgroup"}).
#' @param keep_fixations Logical; keep the per-fixation tables in a list column
#'   \code{revisit_fixations} (default \code{FALSE}).
#' @param ... Ignored.
#'
#' @return \code{x} with columns from \code{\link{revisit_summary}}
#'   (\code{n_fix}, \code{opportunities}, \code{revisits}, \code{revisit_rate},
#'   \code{expected_revisits}, \code{expected_rate}, \code{excess_revisits},
#'   \code{excess_rate}).
#' @examples
#' set.seed(1)
#' et <- simulate_eye_table(n_fixations = 120, n_groups = 6)
#' et$subject <- rep(c("s1", "s2"), each = 3)
#' res <- fixation_revisits(et, radius = 150, null = "subject", subject_var = "subject")
#' tibble::as_tibble(res)[, c("subject", "n_fix", "revisits", "expected_revisits", "excess_revisits")]
#' @export
fixation_revisits.eye_table <- function(x, radius, min_lag = 2, max_lag = Inf,
                                        require_departure = TRUE,
                                        null = c("none", "subject"), subject_var = NULL,
                                        var = "fixgroup", keep_fixations = FALSE,
                                        dist_method = c("euclidean", "manhattan"), ...) {
  null <- match.arg(null)
  dist_method <- match.arg(dist_method)
  if (!var %in% names(x)) stop("Column '", var, "' not found in the eye_table.")
  if (null == "subject" && (is.null(subject_var) || !subject_var %in% names(x))) {
    stop("null = \"subject\" needs subject_var naming a column of the eye_table.")
  }
  fgs <- x[[var]]
  xy_of <- function(fg) if (is.null(fg) || !nrow(fg)) matrix(numeric(0), 0, 2) else cbind(fg$x, fg$y)

  per <- lapply(seq_along(fgs), function(k) {
    fg <- fgs[[k]]
    if (is.null(fg)) fg <- fixation_group(numeric(0), numeric(0), numeric(0), numeric(0))
    ref <- NULL
    if (null == "subject") {
      same <- which(x[[subject_var]] == x[[subject_var]][k])
      ref <- do.call(rbind, lapply(fgs[setdiff(same, k)], xy_of))
      if (is.null(ref)) ref <- matrix(numeric(0), 0, 2)
    }
    fixation_revisits(fg, radius = radius, min_lag = min_lag, max_lag = max_lag,
                      require_departure = require_departure, reference = ref,
                      dist_method = dist_method)
  })
  summ <- dplyr::bind_rows(lapply(per, revisit_summary))
  out <- dplyr::bind_cols(x, summ)
  if (keep_fixations) out$revisit_fixations <- per
  class(out) <- class(x)
  out
}

#' Summarise revisits for one episode
#'
#' @param x A table returned by \code{fixation_revisits()} for a
#'   \code{fixation_group}.
#'
#' @return A one-row tibble: \code{n_fix}; \code{opportunities}; \code{revisits};
#'   \code{revisit_rate} (revisits per opportunity); \code{expected_revisits} and
#'   \code{expected_rate} under the null (\code{NA} without a reference);
#'   \code{excess_revisits} (observed minus expected) and \code{excess_rate}.
#' @export
revisit_summary <- function(x) {
  if (!inherits(x, "revisit_table")) stop("revisit_summary() needs the output of fixation_revisits() for a fixation_group.")
  opp <- sum(x$opportunity)
  rev <- sum(x$revisit)
  has_null <- isTRUE(attr(x, "null"))
  exp_n <- if (has_null && opp > 0) sum(x$expected[x$opportunity], na.rm = TRUE) else NA_real_
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

new_revisit_table <- function(out, radius, min_lag, max_lag, require_departure, has_null) {
  attr(out, "radius") <- radius
  attr(out, "min_lag") <- min_lag
  attr(out, "max_lag") <- max_lag
  attr(out, "require_departure") <- require_departure
  attr(out, "null") <- has_null
  class(out) <- c("revisit_table", class(out))
  out
}

unclass_revisit <- function(x) {
  class(x) <- setdiff(class(x), "revisit_table")
  x
}

check_revisit_args <- function(radius, min_lag, max_lag) {
  if (missing(radius) || length(radius) != 1L || !is.finite(radius) || radius <= 0) {
    stop("radius must be a single positive number (in fixation coordinate units).")
  }
  if (length(min_lag) != 1L || !is.finite(min_lag) || min_lag < 1 || min_lag != round(min_lag)) {
    stop("min_lag must be a positive integer.")
  }
  if (length(max_lag) != 1L || is.na(max_lag) || max_lag < min_lag) {
    stop("max_lag must be >= min_lag (Inf for no limit).")
  }
}

revisit_eligible <- function(i, min_lag, max_lag) {
  hi <- i - min_lag
  lo <- if (is.finite(max_lag)) max(1L, i - as.integer(max_lag)) else 1L
  if (hi < lo) integer(0) else seq.int(lo, hi)
}

revisit_dist <- function(a, b, method) {
  if (method == "euclidean") {
    sqrt(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2)
  } else {
    abs(outer(a[, 1], b[, 1], "-")) + abs(outer(a[, 2], b[, 2], "-"))
  }
}

as_reference_matrix <- function(reference) {
  m <- if (is.data.frame(reference) && all(c("x", "y") %in% names(reference))) {
    cbind(reference$x, reference$y)
  } else {
    as.matrix(reference)
  }
  if (ncol(m) != 2L) stop("reference must have two columns (x, y).")
  m <- m[is.finite(m[, 1]) & is.finite(m[, 2]), , drop = FALSE]
  storage.mode(m) <- "double"
  m
}
