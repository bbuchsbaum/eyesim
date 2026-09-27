make_viz_fixations <- function() {
  set.seed(11)
  centres <- rbind(c(256, 476), c(696, 422), c(890, 225))
  who <- sample(1:3, 24, TRUE)
  fixation_group(
    x = centres[who, 1] + stats::rnorm(24, sd = 35),
    y = centres[who, 2] + stats::rnorm(24, sd = 35),
    duration = round(stats::runif(24, 100, 450)),
    onset = cumsum(stats::runif(24, 150, 400))
  )
}

layer_with <- function(plot, column) {
  built <- ggplot2::ggplot_build(plot)
  Filter(function(d) column %in% names(d), built$data)[[1]]
}

opacity <- function(colours) {
  grDevices::col2rgb(colours, alpha = TRUE)[4, ] / 255
}

test_that("every plot.fixation_group type builds on the eyesim theme", {
  fg <- make_viz_fixations()
  for (type in c("points", "contour", "filled_contour", "density", "raster")) {
    p <- plot(fg, type = type)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
    expect_true(isTRUE(p$coordinates$ratio == 1))
  }
})

boxes_overlap <- function(placed, w, h) {
  idx <- which(placed$shown)
  if (length(idx) < 2L) return(FALSE)
  pairs <- utils::combn(idx, 2)
  any(abs(placed$lx[pairs[1, ]] - placed$lx[pairs[2, ]]) < (w[pairs[1, ]] + w[pairs[2, ]]) / 2 &
        abs(placed$ly[pairs[1, ]] - placed$ly[pairs[2, ]]) < (h[pairs[1, ]] + h[pairs[2, ]]) / 2)
}

test_that("fixation numbers are unambiguous and never overlap at any panel size", {
  fg <- make_viz_fixations()
  for (panel_w in c(150, 47)) {
    s <- panel_w / diff(range(fg$x))
    px <- (fg$x - min(fg$x)) * s
    py <- (fg$y - min(fg$y)) * s
    w <- 1 + 1.3 * nchar(as.character(fg$index))
    h <- rep(3.1, nrow(fg))
    r <- 2.5 * sqrt(fg$duration / max(fg$duration)) + 0.3
    placed <- place_labels_mm(px, py, w, h, r, panel_w, max(py),
                              start = which.min(fg$onset), start_r = 3.4)
    expect_false(boxes_overlap(placed, w, h))
    for (i in which(placed$shown & !placed$leader)) {
      d <- (px - placed$lx[i])^2 + (py - placed$ly[i])^2
      expect_identical(which.min(d), i)
    }
    # a leader never passes over another fixation's marker
    for (i in which(placed$shown & placed$leader)) {
      others <- setdiff(seq_along(px), i)
      clearance <- seg_point_dist(placed$x0[i], placed$y0[i], placed$lx[i], placed$ly[i],
                                  px[others], py[others])
      expect_true(all(clearance >= r[others]))
    }
    # At regular size the start of the scanpath is labelled (at very small
    # sizes it may be omitted; the start ring still marks it).
    if (panel_w == 150) expect_true(placed$shown[[1]])
  }
})

test_that("unlabelled fixation numbers are grouped or noted inside the panel", {
  set.seed(3)
  crowd <- fixation_group(x = stats::rnorm(40, 100, 2), y = stats::rnorm(40, 100, 2),
                          duration = rep(200, 40), onset = seq(0, 3900, by = 100))
  p <- plot(crowd, xlim = c(0, 200), ylim = c(0, 200)) + ggplot2::labs(caption = "mine")
  grDevices::pdf(NULL, width = 3, height = 3)
  on.exit(grDevices::dev.off())
  print(p)
  forced <- grid::grid.force()
  texts <- unlist(lapply(grid::grid.ls(print = FALSE, recursive = TRUE)$name, function(n) {
    g <- tryCatch(grid::grid.get(n), error = function(e) NULL)
    if (inherits(g, "text")) g$label else NULL
  }))
  # every number is accounted for: in a group label or in the note
  expect_true(any(grepl("not labelled", texts)) || any(grepl(", ", texts)))
  expect_true("mine" %in% texts)
})

test_that("density bands above zero are visible and the zero band is empty", {
  fg <- make_viz_fixations()
  for (type in c("density", "filled_contour")) {
    bands <- layer_with(plot(fg, type = type), "level_low")
    expect_true(all(opacity(bands$fill[bands$level_low <= 0]) == 0))
    expect_gt(min(opacity(bands$fill[bands$level_low > 0])), 0.1)
  }
})

test_that("density above shared limits takes the top colour", {
  fg <- make_viz_fixations()
  peak <- max(layer_with(plot(fg, type = "raster"), "density")$density)
  for (type in c("density", "filled_contour")) {
    bands <- layer_with(plot(fg, type = type, limits = c(0, peak / 2)), "level_high")
    top <- bands$level_high == max(bands$level_high)
    expect_true(is.infinite(max(bands$level_high)))
    expect_gt(min(opacity(bands$fill[top])), 0.5)
  }
  ed <- eye_density(fg, sigma = 50, xbounds = c(0, 1024), ybounds = c(0, 768))
  cells <- layer_with(plot(ed, limits = c(0, max(ed$z) / 2)), "fill")
  expect_false(any(is.na(cells$fill)))
  expect_false(any(toupper(cells$fill) %in% c("GREY50", "#7F7F7F")))
})

test_that("rank transform keeps the zero band empty and prints no ticks", {
  fg <- make_viz_fixations()
  for (type in c("density", "filled_contour")) {
    bands <- suppressWarnings(layer_with(plot(fg, type = type, transform = "rank"), "level_low"))
    expect_true(all(opacity(bands$fill[bands$level_low <= 0]) == 0))
    expect_gt(min(opacity(bands$fill[bands$level_low > 0])), 0.05)
  }
  guide <- ggplot2::get_guide_data(plot(fg, type = "raster", transform = "rank"), "fill")
  expect_true(is.null(guide) || nrow(guide) == 0L)
  expect_warning(plot(fg, type = "density", transform = "rank", limits = c(0, 1e-4)),
                 "ignored")
})

test_that("transform changes the density fill", {
  fg <- make_viz_fixations()
  fill_at_quarter <- function(transform) {
    cells <- layer_with(plot(fg, type = "raster", transform = transform), "density")
    cells$fill[which.min(abs(cells$density - 0.25 * max(cells$density)))]
  }
  expect_false(identical(fill_at_quarter("identity"), fill_at_quarter("sqroot")))
})

test_that("colours and legend titles can be overridden", {
  fg <- make_viz_fixations()
  default <- layer_with(plot(fg), "fill")$fill
  custom <- layer_with(plot(fg, colours = c("red", "blue")), "fill")$fill
  expect_false(identical(default, custom))
  ed <- eye_density(fg, sigma = 50, xbounds = c(0, 1024), ybounds = c(0, 768))
  expect_identical((plot(ed) + ggplot2::labs(fill = "X"))$labels$fill, "X")
})

test_that("background images keep their stored colours", {
  skip_if_not_installed("imager")
  skip_if_not_installed("png")
  file <- tempfile(fileext = ".png")
  img <- array(0, c(4, 4, 3))
  img[, , 1] <- 0xD9 / 255
  img[, , 2] <- 0xD4 / 255
  img[, , 3] <- 0xC7 / 255
  png::writePNG(img, file)
  layer <- eyesim_bg_layer(file, c(0, 4), c(0, 4))
  raster <- layer$geom_params$raster
  expect_length(raster, 16L)
  expect_true(all(toupper(as.vector(raster)) == "#D9D4C7"))
})

test_that("theme, colours and wrapping helpers behave", {
  expect_s3_class(theme_eyesim(), "theme")
  expect_s3_class(theme_eyesim_spatial(), "theme")
  expect_named(eyesim_colours(c("reference", "source")), c("reference", "source"))
  expect_error(eyesim_colours("nope"), "Unknown")
})

test_that("titles, subtitles and captions wrap to the drawn width", {
  long <- paste(rep("a caption that keeps going", 12), collapse = " ")
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + ggplot2::labs(caption = long) + theme_eyesim()
  for (w in c(3, 8)) {
    grDevices::pdf(NULL, width = w, height = 4)
    print(p)
    grid::grid.force()
    texts <- list()
    for (n in grid::grid.ls(print = FALSE, recursive = TRUE)$name) {
      g <- tryCatch(grid::grid.get(n), error = function(e) NULL)
      if (inherits(g, "text") && any(grepl("caption that keeps", g$label))) texts <- c(texts, list(g))
    }
    expect_true(length(texts) >= 1L)
    lines <- strsplit(texts[[1]]$label, "\n")[[1]]
    widest <- max(vapply(lines, function(l) {
      grid::convertWidth(grid::grobWidth(grid::textGrob(l, gp = texts[[1]]$gp)), "in",
                         valueOnly = TRUE)
    }, numeric(1)))
    grDevices::dev.off()
    expect_gt(length(lines), 1L)
    expect_lte(widest, w)
  }
})

test_that("the stepped density bar shows the zero band as transparent", {
  fg <- make_viz_fixations()
  steps <- ggplot2::get_guide_data(plot(fg, type = "filled_contour"), "fill")
  steps <- steps[!is.na(steps$fill), , drop = FALSE]
  expect_equal(unname(opacity(steps$fill[[1]])), 0)
  expect_true(all(opacity(steps$fill[-1]) > 0))
})

test_that("both braids size points by mass and fade by matched share", {
  ribbons <- data.frame(x = 0.2, y = 1, xend = 0.3, yend = 0, relative_mass = 1)
  reference <- data.frame(time = c(0.2, 0.6), rail = 1, mass = c(0.5, 0.5))
  source <- data.frame(time = c(0.3, 0.7), rail = 0, mass = c(0.3, 0.7),
                       share = c(0.1, 1))
  p <- gaze_braid_rails(ribbons, reference, source, c("source", "reference"),
                        "t", "s", "x", share_note = "share")
  pts <- Filter(function(d) "fill" %in% names(d) && !anyNA(d$fill),
                ggplot2::ggplot_build(p)$data)[[1]]
  expect_lt(opacity(pts$fill[[1]]), opacity(pts$fill[[2]]))
  expect_lt(pts$size[[1]], pts$size[[2]])
})

test_that("user theme tweaks merge into the wrapping text elements", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point() + ggplot2::labs(title = "t") + theme_eyesim() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, colour = "red"))
  el <- ggplot2::calc_element("plot.title", ggplot2::complete_theme(p$theme))
  expect_s3_class(el, "eyesim_element_text_wrap")
  expect_equal(el@size, 16)
  expect_equal(el@colour, "red")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(print(p))
})

test_that("fixations are numbered and joined in onset order", {
  fg <- fixation_group(x = c(10, 50, 90), y = c(10, 60, 20),
                       duration = c(200, 200, 200), onset = c(400, 0, 200))
  p <- plot(fg)
  path <- Filter(function(d) "group" %in% names(d) && nrow(d) == 3L && !"size" %in% names(d),
                 ggplot2::ggplot_build(p)$data)[[1]]
  expect_equal(path$x, c(50, 90, 10))
  lab <- Filter(function(d) "label" %in% names(d), ggplot2::ggplot_build(p)$data)[[1]]
  expect_equal(as.integer(lab$label[match(c(50, 90, 10), lab$x)]), 1:3)
})

test_that("overlapping fixations are labelled as a group with a clear leader", {
  set.seed(7)
  px <- c(stats::rnorm(6, 40, 1.2), 90, 120)
  py <- c(stats::rnorm(6, 40, 1.2), 70, 20)
  n <- length(px)
  r <- rep(2, n)
  labels <- as.character(seq_len(n))
  placed <- place_labels_mm(px, py, rep(4, n), rep(3.1, n), r, 150, 90,
                            labels = labels)
  groups <- attr(placed, "clusters")
  grouped <- as.integer(unlist(strsplit(groups$members, ",")))
  # every number is labelled or in a group label
  expect_setequal(c(which(placed$shown), grouped), seq_len(n))
  for (k in seq_len(nrow(groups))) {
    members <- as.integer(strsplit(groups$members[k], ",")[[1]])
    expect_identical(groups$text[k], paste(sort(members), collapse = ", "))
    outside <- setdiff(seq_len(n), members)
    clearance <- seg_point_dist(groups$x0[k], groups$y0[k], groups$lx[k], groups$ly[k],
                                px[outside], py[outside])
    expect_true(all(clearance >= r[outside]))
  }
})

test_that("newlines in wrapped labels stay line breaks", {
  gp <- grid::gpar(fontsize = 11)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_identical(wrap_to_width("Condition A\n(n = 20)", 5, gp), "Condition A\n(n = 20)")
})
