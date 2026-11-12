#' Construct leave-out reference maps for two crossed factors
#'
#' Stores source maps and grouped sums once. References are evaluated on demand,
#' without storing four full maps for every subject/item pair.
#' @param dens_tab Table with one map per subject/item pair. NULL maps are missing.
#' @param subject,item Names of the two crossed-factor columns.
#' @param refs Requested reference names (see Details).
#' @param space Average log values or density values.
#' @param floor Explicit positive additive offset required when `space = "log"`
#'   (the default) or `transform = "log"`. For unit-sum maps, a starting choice
#'   is `1 / (10 * ncells)`, where `ncells` is the number of grid cells. There is
#'   no automatic offset because this choice changes the reference scores.
#' @param min_n Minimum contributing units per reference.
#' @param density_var Column containing density objects.
#' @param normalize Normalize each source map to unit sum, or retain its values.
#' @param transform Optional log transformation after density-space aggregation.
#'   Log-space aggregation already returns log scores and requires `"none"`.
#' @details `self_item` uses the target map; `self_other_items` averages that
#'   subject's other maps; `others_item` averages other subjects' maps of that
#'   item. `others_other_items` averages other items within each other subject,
#'   then averages those subject means equally, including in incomplete designs.
#'   `min_n` counts maps for self references and subjects for other references.
#'   Insufficient references retain their keys but return missing scores.
#'   Normalization precedes log transformation, which precedes aggregation.
#' @return An `eyesim_references` object with keys, counts, source values, grouped
#'   sums and a specification. Use [reference_map()] to inspect one map, or
#'   [sample_density()] to score a table of locations carrying subject/item keys.
#' @md
#' @export
reference_maps <- function(dens_tab, subject = "subject", item = "item",
                           refs = c("self_item", "self_other_items", "others_item",
                                    "others_other_items"),
                           space = c("log", "density"), floor = NULL, min_n = 1,
                           density_var = "density", normalize = c("sum", "none"),
                           transform = c("none", "log")) {
  space <- match.arg(space)
  normalize <- match.arg(normalize)
  transform <- match.arg(transform)
  if (space == "log" && transform != "none") stop("Log-space references are already log scores.")
  check_count(min_n, "min_n")
  keys <- c(subject, item)
  check_key_columns(dens_tab, keys)
  if (!nrow(dens_tab) || anyDuplicated(dens_tab[keys])) stop("Reference keys must be nonempty and unique.")
  allowed <- c("self_item", "self_other_items", "others_item", "others_other_items")
  if (!length(refs) || anyNA(refs) || anyDuplicated(refs) || !all(refs %in% allowed)) {
    stop("refs must contain distinct supported reference names.")
  }
  if (!density_var %in% names(dens_tab) || !is.list(dens_tab[[density_var]])) {
    stop("density_var must identify a list column of density maps.")
  }
  if (space == "log" || transform == "log") {
    if (is.null(floor)) {
      stop("Log references require an explicit positive floor. For unit-sum maps, consider floor = 1 / (10 * ncells), or use space = 'density', transform = 'none'.")
    }
    check_log_floor(floor)
  }
  maps <- dens_tab[[density_var]]
  valid <- !vapply(maps, is.null, logical(1))
  if (!any(valid)) stop("At least one non-NULL density map is required.")
  first <- maps[[which(valid)[1]]]
  validate_density_grid(first)
  ncells <- length(first$z)
  values <- matrix(0, ncells, length(maps))
  for (j in which(valid)) {
    map <- maps[[j]]
    validate_density_grid(map)
    check_same_density_lattice(first, map)
    if (!inherits(map, "density") || any(!is.finite(map$z)) || any(map$z < 0) || sum(map$z) <= 0) {
      stop("Reference inputs must be finite non-negative density maps with positive mass.")
    }
    values[, j] <- as.vector(prepare_map_values(map$z, normalize,
                                               if (space == "log") "log" else "none", floor))
  }
  subjects <- unique(dens_tab[[subject]])
  items <- unique(dens_tab[[item]])
  si <- match(dens_tab[[subject]], subjects)
  ii <- match(dens_tab[[item]], items)
  ns <- tabulate(si[valid], length(subjects))
  ni <- tabulate(ii[valid], length(items))
  sum_columns <- function(ids, n) {
    out <- matrix(0, ncells, n)
    grouped <- rowsum(t(values), ids, reorder = FALSE)
    out[, as.integer(rownames(grouped))] <- t(grouped)
    out
  }
  subject_sums <- sum_columns(si, length(subjects))
  item_sums <- sum_columns(ii, length(items))
  # One subject-balanced other-item sum per item, not per subject/item pair.
  cross_sums <- matrix(0, ncells, length(items))
  cross_n <- integer(length(items))
  for (j in seq_along(items)) {
    rows <- which(ii == j & valid)
    counts <- ns
    counts[si[rows]] <- counts[si[rows]] - 1L
    sums <- subject_sums
    sums[, si[rows]] <- sums[, si[rows], drop = FALSE] - values[, rows, drop = FALSE]
    keep <- counts > 0L
    cross_n[j] <- sum(keep)
    if (any(keep)) cross_sums[, j] <- rowSums(sweep(sums[, keep, drop = FALSE], 2, counts[keep], "/"))
  }
  counts <- data.frame(self_item = as.integer(valid),
                       self_other_items = ns[si] - valid,
                       others_item = ni[ii] - valid,
                       others_other_items = cross_n[ii] - (ns[si] - valid > 0))
  out <- structure(list(keys = as.data.frame(dens_tab[keys]), refs = refs,
                        counts = counts[refs], x = first$x, y = first$y,
                        values = values, valid = valid, si = si, ii = ii,
                        subject_sums = subject_sums, item_sums = item_sums,
                        cross_sums = cross_sums, ns = ns, ni = ni, cross_n = cross_n,
                        spec = list(subject = subject, item = item, space = space,
                                    floor = floor, normalize = normalize, min_n = min_n,
                                    transform = transform,
                                    domain = if (space == "log" || transform == "log") "log" else "density")),
                   class = "eyesim_references")
  if (any(as.matrix(out$counts) < min_n)) {
    warning("Some references have fewer than min_n contributors; their scores will be NA.")
  }
  out
}

# Evaluate sufficient statistics only at requested cells. No full reference
# maps or fixation-by-candidate Cartesian products are materialized.
reference_values <- function(refs, target, cell, reference) {
  si <- refs$si[target]
  ii <- refs$ii[target]
  own <- refs$values[cbind(cell, target)]
  own_n <- refs$ns[si] - refs$valid[target]
  own_other <- (refs$subject_sums[cbind(cell, si)] - own) / pmax(own_n, 1)
  value <- switch(reference,
    self_item = own,
    self_other_items = own_other,
    others_item = (refs$item_sums[cbind(cell, ii)] - own) /
      pmax(refs$ni[ii] - refs$valid[target], 1),
    others_other_items = (refs$cross_sums[cbind(cell, ii)] - own_other) /
      pmax(refs$cross_n[ii] - (own_n > 0), 1))
  count <- refs$counts[[reference]][target]
  value[is.na(count) | count < refs$spec$min_n] <- NA_real_
  # Density sums can lose a few ulps when excluded maps exhaust a group.
  if (refs$spec$space == "density") value <- pmax(value, 0)
  if (refs$spec$transform == "log") value <- log(value + refs$spec$floor)
  value
}

#' Inspect one reference map
#' @param refs A reference set from [reference_maps()].
#' @param row Row of `refs$keys` identifying the target.
#' @param ref A requested reference name.
#' @return A `eyesim_reference_map` with x/y/z, domain and contributor count,
#'   or NULL when insufficient contributors exist. Log maps contain scores,
#'   not probabilities, and must not be transformed again.
#' @md
#' @export
reference_map <- function(refs, row = 1L, ref = refs$refs[1]) {
  if (!inherits(refs, "eyesim_references")) stop("refs must be a reference set.")
  check_count(row, "row")
  if (row > nrow(refs$keys) || length(ref) != 1L || !ref %in% refs$refs) stop("Unknown reference row or name.")
  if (refs$counts[[ref]][row] < refs$spec$min_n) return(NULL)
  z <- reference_values(refs, rep(row, length(refs$x) * length(refs$y)),
                        seq_len(length(refs$x) * length(refs$y)), ref)
  structure(list(x = refs$x, y = refs$y, z = matrix(z, length(refs$x)),
                 domain = refs$spec$domain, n = refs$counts[[ref]][row]),
            class = "eyesim_reference_map")
}

#' @rdname sample_density
#' @md
#' @export
sample_density.eyesim_reference_map <- function(x, fix, times = NULL, ...,
                                               normalize = "none", transform = "none", floor = NULL) {
  if (normalize != "none" || transform != "none" || !is.null(floor)) {
    stop("Reference maps are already prepared; normalization and transformation are fixed.")
  }
  sample_density.density(x, fix, times = times)
}

#' @rdname sample_density
#' @md
#' @export
sample_density.eyesim_references <- function(x, fix, ..., format = c("wide", "long")) {
  format <- match.arg(format)
  if (length(list(...))) stop("Reference sampling accepts locations and format only; map preparation is fixed.")
  target <- match_map_rows(fix, x$keys, names(x$keys), strict = TRUE)
  if (anyNA(target)) stop("Some locations have no matching reference keys.")
  check_locations(fix)
  idx <- density_grid_index(fix, x$x, x$y)
  cell <- idx[, 1] + (idx[, 2] - 1L) * length(x$x)
  values <- evaluate_references(x, target, cell, x$refs)
  if (format == "wide") return(values)
  dplyr::bind_rows(lapply(x$refs, function(ref) {
    data.frame(fixation = seq_len(nrow(fix)), reference = rep(ref, nrow(fix)), z = values[[ref]])
  }))
}

evaluate_references <- function(refs, target, cell, columns) {
  values <- lapply(columns, function(ref) reference_values(refs, target, cell, ref))
  as.data.frame(stats::setNames(values, columns))
}

check_locations <- function(x) {
  if (!all(c("x", "y") %in% names(x)) || !is.numeric(x$x) || !is.numeric(x$y) ||
      any(!is.finite(x$x)) || any(!is.finite(x$y))) stop("Locations require finite numeric x and y columns.")
}

#' @export
#' @noRd
print.eyesim_references <- function(x, ...) {
  cat("Reference set:", nrow(x$keys), "target keys;", length(x$x), "x", length(x$y), "grid\n")
  cat("Aggregation:", x$spec$space, "| output:", x$spec$domain, "\n")
  print(data.frame(reference = x$refs,
                   supported_targets = vapply(x$counts, function(n) sum(n >= x$spec$min_n), integer(1))),
        row.names = FALSE)
  invisible(x)
}
