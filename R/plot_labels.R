# Fixation-number labels placed at draw time ------------------------------
#
# Placement needs the physical panel size (label size is fixed in mm while
# the data scale is not), so labels are computed in makeContent(), which
# runs whenever the plot is drawn at a given size.

# Place labels for points at (px, py) in mm. `w`, `h` are label box sizes in
# mm and `r` the marker radii in mm. Returns a data frame with the label
# centre (lx, ly), whether a leader line is drawn, and whether the label is
# shown. Labels are unambiguous: a label next to its point must be clearly
# nearer its own fixation than any other; a label further out gets a leader
# line to its fixation. Boxes never overlap each other or other markers.
place_labels_mm <- function(px, py, w, h, r, panel_w, panel_h,
                            gap = 0.6, ownership = 0.85, reserve = NULL,
                            start = NULL, start_r = 0, labels = NULL,
                            measure = function(s) 1 + 1.3 * nchar(s),
                            label_h = NULL) {
  n <- length(px)
  out <- data.frame(lx = rep(NA_real_, n), ly = rep(NA_real_, n),
                    leader = rep(FALSE, n), shown = rep(FALSE, n),
                    x0 = px, y0 = py)
  if (n == 0L) {
    return(out)
  }
  dirs <- rbind(c(1, 1), c(-1, 1), c(1, -1), c(-1, -1),
                c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  dirs <- dirs / sqrt(rowSums(dirs^2))
  angles <- seq(0, 2 * pi, length.out = 17L)[-17L] + pi / 4
  dirs16 <- cbind(cos(angles), sin(angles))
  # occupied boxes (xmin, xmax, ymin, ymax); `reserve` blocks a region
  boxes <- if (is.null(reserve)) matrix(numeric(0), ncol = 4) else matrix(reserve, ncol = 4)
  leaders <- matrix(numeric(0), ncol = 4) # x0, y0, x1, y1

  box_of <- function(i, cx, cy) {
    c(cx - w[[i]] / 2, cx + w[[i]] / 2, cy - h[[i]] / 2, cy + h[[i]] / 2)
  }
  fits <- function(i, b) {
    if (b[[1]] < 0 || b[[2]] > panel_w || b[[3]] < 0 || b[[4]] > panel_h) {
      return(FALSE)
    }
    if (nrow(boxes) > 0L && any(
      boxes[, 1] < b[[2]] + gap & boxes[, 2] > b[[1]] - gap &
        boxes[, 3] < b[[4]] + gap & boxes[, 4] > b[[3]] - gap
    )) {
      return(FALSE)
    }
    # the box must not cover any marker (closest box point vs radius)
    qx <- pmin(pmax(px, b[[1]]), b[[2]])
    qy <- pmin(pmax(py, b[[3]]), b[[4]])
    clear <- sqrt((qx - px)^2 + (qy - py)^2) >= r + gap / 2
    clear[[i]] <- TRUE
    all(clear)
  }
  # A leader must stay clear of every other marker, cross no other leader
  # and pass through no placed label, so it can only be read one way.
  # The start fixation's leader begins at the edge of its ring, which marks
  # it uniquely even when other markers overlap it.
  leader_origin <- function(i, x1, y1) {
    if (!is.null(start) && i == start && start_r > 0) {
      len <- sqrt((x1 - px[[i]])^2 + (y1 - py[[i]])^2)
      if (len > start_r) {
        return(c(px[[i]] + (x1 - px[[i]]) * start_r / len,
                 py[[i]] + (y1 - py[[i]]) * start_r / len))
      }
    }
    c(px[[i]], py[[i]])
  }
  leader_clear <- function(i, x1, y1) {
    origin <- leader_origin(i, x1, y1)
    x0 <- origin[[1]]
    y0 <- origin[[2]]
    others <- setdiff(seq_len(n), i)
    if (length(others) > 0L &&
        any(seg_point_dist(x0, y0, x1, y1, px[others], py[others]) <
            r[others] + gap)) {
      return(FALSE)
    }
    if (nrow(leaders) > 0L && any(vapply(seq_len(nrow(leaders)), function(k) {
      segments_cross(x0, y0, x1, y1, leaders[k, 1], leaders[k, 2],
                     leaders[k, 3], leaders[k, 4])
    }, logical(1)))) {
      return(FALSE)
    }
    t <- seq(0.05, 0.95, length.out = 12)
    sx <- x0 + t * (x1 - x0)
    sy <- y0 + t * (y1 - y0)
    if (nrow(boxes) > 0L) {
      for (k in seq_len(nrow(boxes))) {
        if (any(sx > boxes[k, 1] & sx < boxes[k, 2] &
                sy > boxes[k, 3] & sy < boxes[k, 4])) {
          return(FALSE)
        }
      }
    }
    TRUE
  }
  take <- function(i, cx, cy, leader) {
    boxes <<- rbind(boxes, box_of(i, cx, cy))
    if (leader) {
      origin <- leader_origin(i, cx, cy)
      leaders <<- rbind(leaders, c(origin, cx, cy))
      out$x0[[i]] <<- origin[[1]]
      out$y0[[i]] <<- origin[[2]]
    }
    out$lx[[i]] <<- cx
    out$ly[[i]] <<- cy
    out$leader[[i]] <<- leader
    out$shown[[i]] <<- TRUE
  }

  # The start and end of the scanpath get first choice of position.
  first <- if (is.null(start)) 1L else start
  priority <- unique(c(first, n, seq_len(n)))
  for (i in priority) {
    reach <- r[[i]] + gap + sqrt((w[[i]] / 2)^2 + (h[[i]] / 2)^2)
    placed <- FALSE
    for (d in seq_len(nrow(dirs))) {
      cx <- px[[i]] + dirs[d, 1] * reach
      cy <- py[[i]] + dirs[d, 2] * reach
      dist <- sqrt((px - cx)^2 + (py - cy)^2)
      other <- if (n > 1L) min(dist[-i]) else Inf
      if (dist[[i]] <= ownership * other && fits(i, box_of(i, cx, cy))) {
        take(i, cx, cy, leader = FALSE)
        placed <- TRUE
        break
      }
    }
    if (placed) next
    # Further out, a leader line makes ownership explicit; search 16
    # directions at increasing distance for a line that stays clear.
    for (scale in c(1.6, 2.2, 3, 4, 5.5, 7, 9)) {
      for (d in seq_len(nrow(dirs16))) {
        cx <- px[[i]] + dirs16[d, 1] * reach * scale
        cy <- py[[i]] + dirs16[d, 2] * reach * scale
        if (fits(i, box_of(i, cx, cy)) && leader_clear(i, cx, cy)) {
          take(i, cx, cy, leader = TRUE)
          placed <- TRUE
          break
        }
      }
      if (placed) break
    }
  }
  attr(out, "clusters") <- place_cluster_labels(
    out, px, py, r, labels, measure, label_h %||% (if (n > 0L) h[[1]] else 3),
    boxes, leaders, panel_w, panel_h, gap
  )
  out
}

# Numbers that could not be placed individually sit on markers overlapping
# other markers. Such a group is outlined and labelled once ("4, 9, 17"),
# with a leader from the outline that stays clear of every marker outside
# the group. Returns one row per placed group label.
place_cluster_labels <- function(out, px, py, r, labels, measure, label_h,
                                 boxes, leaders, panel_w, panel_h, gap) {
  n <- length(px)
  empty <- data.frame(cx = numeric(0), cy = numeric(0), cr = numeric(0),
                      lx = numeric(0), ly = numeric(0), x0 = numeric(0),
                      y0 = numeric(0), w = numeric(0), text = character(0),
                      members = character(0))
  missing <- which(!out$shown)
  if (length(missing) == 0L || is.null(labels)) {
    return(empty)
  }
  # connected components of overlapping markers
  touch <- outer(seq_len(n), seq_len(n), function(i, j) {
    sqrt((px[i] - px[j])^2 + (py[i] - py[j])^2) < r[i] + r[j] + gap
  })
  group <- seq_len(n)
  repeat {
    new_group <- vapply(seq_len(n), function(i) min(group[touch[i, ]]), integer(1))
    if (identical(new_group, group)) break
    group <- new_group
  }
  angles <- seq(0, 2 * pi, length.out = 17L)[-17L] + pi / 4
  dirs <- cbind(cos(angles), sin(angles))
  rows <- list()
  handled <- integer(0)
  for (g in unique(group[missing])) {
    members <- which(group == g)
    todo <- setdiff(intersect(members, missing), handled)
    if (length(members) < 2L || length(todo) == 0L) next
    # The outline encloses only the unlabelled numbers it lists (a single
    # number becomes a ring on its own marker). Any other unlabelled marker
    # the outline would enclose joins the list, so none sits inside unnamed.
    repeat {
      cx <- mean(px[todo])
      cy <- mean(py[todo])
      cr <- max(sqrt((px[todo] - cx)^2 + (py[todo] - cy)^2) + r[todo]) + gap / 2
      extra <- setdiff(missing, c(todo, handled))
      extra <- extra[sqrt((px[extra] - cx)^2 + (py[extra] - cy)^2) < cr]
      if (length(extra) == 0L) break
      todo <- c(todo, extra)
    }
    handled <- c(handled, todo)
    text <- paste(sort(as.integer(labels[todo])), collapse = ", ")
    bw <- measure(text)
    bh <- label_h
    outside <- setdiff(seq_len(n), todo)
    done <- FALSE
    for (scale in c(1.3, 1.8, 2.5, 3.5, 5)) {
      reach <- cr + gap + sqrt((bw / 2)^2 + (bh / 2)^2) * scale
      for (d in seq_len(nrow(dirs))) {
        lx <- cx + dirs[d, 1] * reach
        ly <- cy + dirs[d, 2] * reach
        b <- c(lx - bw / 2, lx + bw / 2, ly - bh / 2, ly + bh / 2)
        if (b[[1]] < 0 || b[[2]] > panel_w || b[[3]] < 0 || b[[4]] > panel_h) next
        if (nrow(boxes) > 0L && any(
          boxes[, 1] < b[[2]] + gap & boxes[, 2] > b[[1]] - gap &
            boxes[, 3] < b[[4]] + gap & boxes[, 4] > b[[3]] - gap
        )) next
        qx <- pmin(pmax(px, b[[1]]), b[[2]])
        qy <- pmin(pmax(py, b[[3]]), b[[4]])
        if (any(sqrt((qx - px)^2 + (qy - py)^2) < r + gap / 2)) next
        x0 <- cx + dirs[d, 1] * cr
        y0 <- cy + dirs[d, 2] * cr
        if (length(outside) > 0L &&
            any(seg_point_dist(x0, y0, lx, ly, px[outside], py[outside]) <
                r[outside] + gap)) next
        if (nrow(leaders) > 0L && any(vapply(seq_len(nrow(leaders)), function(k) {
          segments_cross(x0, y0, lx, ly, leaders[k, 1], leaders[k, 2],
                         leaders[k, 3], leaders[k, 4])
        }, logical(1)))) next
        boxes <- rbind(boxes, b)
        leaders <- rbind(leaders, c(x0, y0, lx, ly))
        rows[[length(rows) + 1L]] <- data.frame(
          cx = cx, cy = cy, cr = cr, lx = lx, ly = ly, x0 = x0, y0 = y0,
          w = bw, text = text, members = paste(todo, collapse = ",")
        )
        done <- TRUE
        break
      }
      if (done) break
    }
  }
  if (length(rows) == 0L) empty else do.call(rbind, rows)
}

# Distance from points (qx, qy) to the segment (x0, y0)-(x1, y1).
seg_point_dist <- function(x0, y0, x1, y1, qx, qy) {
  dx <- x1 - x0
  dy <- y1 - y0
  len2 <- dx^2 + dy^2
  t <- if (len2 > 0) pmin(pmax(((qx - x0) * dx + (qy - y0) * dy) / len2, 0), 1) else 0
  sqrt((x0 + t * dx - qx)^2 + (y0 + t * dy - qy)^2)
}

segments_cross <- function(ax, ay, bx, by, cx, cy, dx, dy) {
  orient <- function(px, py, qx, qy, rx, ry) sign((qx - px) * (ry - py) - (qy - py) * (rx - px))
  o1 <- orient(ax, ay, bx, by, cx, cy)
  o2 <- orient(ax, ay, bx, by, dx, dy)
  o3 <- orient(cx, cy, dx, dy, ax, ay)
  o4 <- orient(cx, cy, dx, dy, bx, by)
  o1 != o2 && o3 != o4 && o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0
}

GeomFixationLabel <- ggplot2::ggproto(
  "GeomFixationLabel", ggplot2::Geom,
  required_aes = c("x", "y", "label"),
  default_aes = ggplot2::aes(size = 2.3, radius = 1.5, start = FALSE),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(data, panel_params, coord, colour = "#1F2328",
                        start_r = 0) {
    coords <- coord$transform(data, panel_params)
    grid::gTree(
      pos = data.frame(x = coords$x, y = coords$y,
                       label = as.character(coords$label),
                       size = coords$size, radius = coords$radius,
                       start = as.logical(coords$start)),
      colour = colour, start_r = start_r,
      cl = "eyesim_fixation_labels"
    )
  }
)

#' @export
#' @importFrom grid makeContent
makeContent.eyesim_fixation_labels <- function(x) {
  pos <- x$pos
  pw <- grid::convertWidth(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
  ph <- grid::convertHeight(grid::unit(1, "npc"), "mm", valueOnly = TRUE)
  fontsize <- pos$size[[1]] * ggplot2::.pt
  gp <- grid::gpar(fontsize = fontsize)
  pad <- 0.5
  measure <- function(s) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = gp)), "mm",
                       valueOnly = TRUE) + 2 * pad
  }
  w <- vapply(pos$label, measure, numeric(1))
  label_h <- fontsize / ggplot2::.pt * 1.15 + 2 * pad
  h <- rep(label_h, nrow(pos))
  px <- pos$x * pw
  py <- pos$y * ph
  start <- which(pos$start)[1]
  if (is.na(start)) start <- NULL
  place <- function(reserve = NULL) {
    place_labels_mm(px, py, w, h, pos$radius, pw, ph, reserve = reserve,
                    start = start, start_r = x$start_r, labels = pos$label,
                    measure = measure, label_h = label_h)
  }
  unlabelled <- function(placed) {
    grouped <- as.integer(unlist(strsplit(attr(placed, "clusters")$members, ",")))
    setdiff(which(!placed$shown), grouped)
  }
  note_gp <- grid::gpar(fontsize = fontsize * 0.9, col = eyesim_colours("muted"))
  note_text <- function(missing) {
    if (length(missing) <= 6L) {
      paste0("not labelled: ", paste(sort(as.integer(pos$label[missing])), collapse = ", "))
    } else {
      paste0(length(missing), " of ", nrow(pos), " numbers not labelled")
    }
  }
  placed <- place()
  note_w <- note_h <- 0
  note_x <- note_y <- 0.5
  if (length(unlabelled(placed)) > 0L) {
    # Put the note in the corner with the fewest markers under it, reserve
    # that corner so it never covers a number, then re-place the labels.
    probe <- grid::textGrob(note_text(unlabelled(placed)), gp = note_gp)
    note_w <- grid::convertWidth(grid::grobWidth(probe), "mm", valueOnly = TRUE) + 2
    note_h <- grid::convertHeight(grid::grobHeight(probe), "mm", valueOnly = TRUE) + 1.5
    corners <- rbind(c(0.5, 0.5), c(pw - note_w - 0.5, 0.5),
                     c(0.5, ph - note_h - 0.5), c(pw - note_w - 0.5, ph - note_h - 0.5))
    covered <- apply(corners, 1, function(cn) {
      qx <- pmin(pmax(px, cn[[1]]), cn[[1]] + note_w)
      qy <- pmin(pmax(py, cn[[2]]), cn[[2]] + note_h)
      sum(sqrt((qx - px)^2 + (qy - py)^2) < pos$radius)
    })
    best <- which.min(covered)
    note_x <- corners[best, 1]
    note_y <- corners[best, 2]
    placed <- place(reserve = c(note_x - 0.5, note_x + note_w + 0.5,
                                note_y - 0.5, note_y + note_h + 0.5))
  }
  shown <- which(placed$shown)
  groups <- attr(placed, "clusters")
  missing <- unlabelled(placed)
  kids <- list()
  if (nrow(groups) > 0L) {
    kids$groups <- grid::circleGrob(
      x = grid::unit(groups$cx, "mm"), y = grid::unit(groups$cy, "mm"),
      r = grid::unit(groups$cr, "mm"),
      gp = grid::gpar(col = x$colour, fill = NA, lwd = 0.6, lty = "22", alpha = 0.8)
    )
  }
  lead <- shown[placed$leader[shown]]
  lx0 <- c(placed$x0[lead], groups$x0)
  if (length(lx0) > 0L) {
    kids$leaders <- grid::segmentsGrob(
      x0 = grid::unit(lx0, "mm"),
      y0 = grid::unit(c(placed$y0[lead], groups$y0), "mm"),
      x1 = grid::unit(c(placed$lx[lead], groups$lx), "mm"),
      y1 = grid::unit(c(placed$ly[lead], groups$ly), "mm"),
      gp = grid::gpar(col = x$colour, lwd = 0.5, alpha = 0.7)
    )
  }
  box_x <- c(placed$lx[shown], groups$lx)
  if (length(box_x) > 0L) {
    kids$boxes <- grid::rectGrob(
      x = grid::unit(box_x, "mm"), y = grid::unit(c(placed$ly[shown], groups$ly), "mm"),
      width = grid::unit(c(w[shown], groups$w), "mm"),
      height = grid::unit(label_h, "mm"),
      gp = grid::gpar(fill = grDevices::adjustcolor("white", 0.85), col = NA)
    )
    kids$text <- grid::textGrob(
      c(pos$label[shown], groups$text), name = "eyesim_fixation_numbers",
      x = grid::unit(box_x, "mm"), y = grid::unit(c(placed$ly[shown], groups$ly), "mm"),
      gp = grid::gpar(fontsize = fontsize, col = x$colour)
    )
  }
  if (length(missing) > 0L) {
    # The note lives in the panel, so a user caption cannot remove it.
    kids$note_bg <- grid::rectGrob(
      x = grid::unit(note_x, "mm"), y = grid::unit(note_y, "mm"),
      width = grid::unit(note_w, "mm"), height = grid::unit(note_h, "mm"),
      just = c(0, 0), gp = grid::gpar(fill = grDevices::adjustcolor("white", 0.9), col = NA)
    )
    kids$note <- grid::textGrob(
      note_text(missing),
      x = grid::unit(note_x + 1, "mm"), y = grid::unit(note_y + 0.7, "mm"), just = c(0, 0),
      gp = note_gp
    )
  }
  grid::setChildren(x, do.call(grid::gList, unname(kids)))
}

geom_fixation_label <- function(mapping = NULL, data = NULL, colour = "#1F2328",
                                ..., inherit.aes = TRUE) {
  ggplot2::layer(
    geom = GeomFixationLabel, stat = "identity", position = "identity",
    data = data, mapping = mapping, inherit.aes = inherit.aes,
    show.legend = FALSE, params = list(colour = colour, ...)
  )
}
