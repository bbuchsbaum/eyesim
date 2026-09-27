# Text elements that wrap to the drawn width ---------------------------------
#
# ggplot2 text elements do not wrap. This subclass of element_text measures
# the width available when the plot is drawn and wraps its label to it, so
# titles, subtitles and captions fit any figure size. theme_eyesim() uses it
# for the plot title, subtitle and caption.

element_text_wrap_class <- S7::new_class(
  "eyesim_element_text_wrap",
  parent = ggplot2::element_text,
  package = NULL
)

#' Wrapping text element
#'
#' Like [ggplot2::element_text()], but the text is wrapped to the width
#' available when the plot is drawn. [theme_eyesim()] uses it for the plot
#' title, subtitle and caption.
#'
#' @param ... Arguments passed to [ggplot2::element_text()].
#' @return A theme element.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() +
#'   labs(caption = paste(rep("A long caption that wraps.", 8), collapse = " ")) +
#'   theme(plot.caption = element_text_wrap(hjust = 0))
#' @family visualization
#' @export
element_text_wrap <- function(...) {
  do.call(element_text_wrap_class, text_element_args(ggplot2::element_text(...)))
}

# Arguments that rebuild a text element through element_text()'s constructor.
text_element_args <- function(element) {
  props <- S7::props(element)
  keep <- intersect(names(props), names(formals(ggplot2::element_text)))
  props[keep]
}

#' @exportS3Method ggplot2::element_grob
element_grob.eyesim_element_text_wrap <- function(element, label = "", ...) {
  if (is.null(label) || length(label) == 0L || !nzchar(paste(label, collapse = ""))) {
    return(ggplot2::zeroGrob())
  }
  plain <- do.call(ggplot2::element_text, text_element_args(element))
  grid::gTree(
    label = paste(label, collapse = " "), plain = plain, args = list(...),
    cache = new.env(parent = emptyenv()),
    cl = "eyesim_wrapped_text"
  )
}

# Build the ggplot2 title grob for the label wrapped to the current width.
#
# Layout asks for the height before the final drawing width is known (for
# example, patchwork measures at the full figure width), so the height is
# computed for the width then in view less `slack` and remembered. At draw
# time the label is wrapped to the real width; if it would then need more
# height than was reserved, the text is set smaller until it fits, so wrapped
# text never overprints what lies below it.
wrapped_title_grob <- function(x, for_height = FALSE, slack = 0.7) {
  width_in <- grid::convertWidth(grid::unit(1, "npc"), "in", valueOnly = TRUE)
  if (for_height) {
    width_in <- max(width_in - slack, 0.5)
  }
  el <- x$plain
  size <- if (is.null(el@size)) 11 else el@size
  build <- function(scale) {
    gp <- grid::gpar(fontsize = size * scale, fontfamily = el@family %||% "",
                     fontface = el@face %||% "plain")
    margin <- el@margin
    side <- if (is.null(margin)) 0 else
      grid::convertWidth(margin[2] + margin[4], "in", valueOnly = TRUE)
    text <- wrap_to_width(x$label, max(width_in - side, 0.5), gp)
    args <- text_element_args(el)
    args$size <- size * scale
    scaled <- do.call(ggplot2::element_text, args)
    do.call(ggplot2::element_grob, c(list(scaled, label = text), x$args))
  }
  grob <- build(1)
  height <- function(g) grid::convertHeight(grid::grobHeight(g), "in", valueOnly = TRUE)
  if (for_height) {
    x$cache$height <- height(grob)
    return(grob)
  }
  reserved <- x$cache$height
  if (!is.null(reserved)) {
    for (scale in seq(0.95, 0.5, by = -0.05)) {
      if (height(grob) <= reserved + 1e-3) break
      grob <- build(scale)
    }
    x$cache$drawn <- height(grob)
  }
  grob
}

# Greedy wrapping to `width_in`, then balanced: the narrowest width that
# keeps the same number of lines, so the last line is not an orphan.
# Newlines in the label are kept as line breaks, as in ggplot2.
wrap_to_width <- function(text, width_in, gp) {
  paragraphs <- strsplit(text, "\n", fixed = TRUE)[[1]]
  measure <- function(s) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = gp)), "in",
                       valueOnly = TRUE)
  }
  greedy <- function(words, w) {
    out <- character(0)
    current <- words[[1]]
    for (word in words[-1]) {
      candidate <- paste(current, word)
      if (measure(candidate) <= w) {
        current <- candidate
      } else {
        out <- c(out, current)
        current <- word
      }
    }
    c(out, current)
  }
  lines <- unlist(lapply(paragraphs, function(p) {
    words <- strsplit(p, " ", fixed = TRUE)[[1]]
    words <- words[nzchar(words)]
    if (length(words) == 0L) return("")
    first <- greedy(words, width_in)
    k <- length(first)
    if (k < 2L) return(first)
    lo <- width_in / 2
    hi <- width_in
    best <- first
    for (iter in seq_len(10L)) {
      mid <- (lo + hi) / 2
      trial <- greedy(words, mid)
      if (length(trial) == k && all(vapply(trial, measure, 1) <= mid + 1e-6)) {
        best <- trial
        hi <- mid
      } else {
        lo <- mid
      }
    }
    best
  }))
  paste(lines, collapse = "\n")
}

#' @export
#' @importFrom grid makeContent
makeContent.eyesim_wrapped_text <- function(x) {
  grid::setChildren(x, grid::gList(wrapped_title_grob(x)))
}

#' @export
#' @importFrom grid heightDetails
heightDetails.eyesim_wrapped_text <- function(x) {
  grid::grobHeight(wrapped_title_grob(x, for_height = TRUE))
}

#' @export
#' @importFrom grid widthDetails
widthDetails.eyesim_wrapped_text <- function(x) {
  grid::unit(1, "null")
}

# Merging a plain element_text into a wrapping element (for example
# `theme_eyesim() + theme(plot.title = element_text(size = 14))`) keeps the
# wrapping and takes the new element's non-NULL properties.
merge_into_wrap <- function(new, old, ...) {
  new_args <- text_element_args(new)
  old_args <- text_element_args(old)
  for (nm in names(new_args)) {
    if (!is.null(new_args[[nm]])) old_args[[nm]] <- new_args[[nm]]
  }
  do.call(element_text_wrap_class, old_args)
}

ggplot2_merge_element <- S7::new_external_generic("ggplot2", "merge_element",
                                                  c("new", "old"))
S7::method(ggplot2_merge_element, list(ggplot2::element_text, element_text_wrap_class)) <-
  merge_into_wrap
