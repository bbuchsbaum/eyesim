reference_fixture <- function() {
  tab <- expand.grid(subject = c("s1", "s2", "s3"), item = c("a", "b", "c"),
                     stringsAsFactors = FALSE)
  tab$density <- lapply(seq_len(nrow(tab)), function(i) {
    gen_density(0:2, 0:1, matrix((1:6 + i)^(1 + i %% 3), 3))
  })
  tab
}

reference_fixations <- function() {
  tab <- reference_fixture()[c("subject", "item")]
  tab$trial <- rep(1:3, each = 3)
  tab$x <- rep(0:2, 3)
  tab$y <- rep(c(0, 1, 0), each = 3)
  tab$onset <- seq_len(nrow(tab))
  tab
}

# Independent, deliberately slow membership-based oracle. No production sums,
# matching, coordinate lookup, normalization or aggregation helpers are used.
reference_oracle <- function(tab, row, ref, space = "log", offset = .01) {
  valid <- !vapply(tab$density, is.null, logical(1))
  subject <- tab$subject[row]
  item <- tab$item[row]
  prepare <- function(m) {
    z <- m$z / sum(m$z)
    if (space == "log") log(z + offset) else z
  }
  average <- function(rows) {
    if (!length(rows)) return(NULL)
    Reduce(`+`, lapply(tab$density[rows], prepare)) / length(rows)
  }
  if (ref == "self_item") return(average(which(valid & tab$subject == subject & tab$item == item)))
  if (ref == "self_other_items") return(average(which(valid & tab$subject == subject & tab$item != item)))
  if (ref == "others_item") return(average(which(valid & tab$subject != subject & tab$item == item)))
  means <- lapply(setdiff(unique(tab$subject), subject), function(s) {
    average(which(valid & tab$subject == s & tab$item != item))
  })
  means <- Filter(Negate(is.null), means)
  if (!length(means)) NULL else Reduce(`+`, means) / length(means)
}
