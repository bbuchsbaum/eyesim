test_that("create_graph path matches Python dijkstra path", {
  testthat::skip_on_cran()
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    testthat::skip("reticulate not available")
  }
  if (!reticulate::py_module_available("multimatch_gaze")) {
    testthat::skip("Python module 'multimatch_gaze' not available")
  }

  mmgaze <- reticulate::import("multimatch_gaze")
  np <- reticulate::import("numpy")

  if (!reticulate::py_has_attr(mmgaze, "createdirectedgraph") ||
      !reticulate::py_has_attr(mmgaze, "dijkstra")) {
    testthat::skip("multimatch_gaze build does not expose graph helpers")
  }

  fg1 <- fixation_group(
    x = c(100, 200, 260, 300),
    y = c(100, 140, 120, 180),
    duration = c(0.2, 0.25, 0.2, 0.3),
    onset = c(0.1, 0.3, 0.6, 1.0)
  )
  fg2 <- fixation_group(
    x = c(120, 210, 255, 310),
    y = c(110, 150, 130, 175),
    duration = c(0.22, 0.24, 0.21, 0.28),
    onset = c(0.05, 0.35, 0.7, 1.2)
  )

  sp1 <- scanpath(fg1)
  sp2 <- scanpath(fg2)
  gout <- eyesim:::create_graph(sp1[1:(nrow(sp1)-1),], sp2[1:(nrow(sp2)-1),])
  r_path <- as.integer(igraph::as_ids(gout$vpath))

  # Also compare M distance matrices across R and Python for sanity
  M_R <- as.matrix(proxy::dist(cbind(sp1$lenx[1:(nrow(sp1)-1)], sp1$leny[1:(nrow(sp1)-1)]),
                               cbind(sp2$lenx[1:(nrow(sp2)-1)], sp2$leny[1:(nrow(sp2)-1)])))

  # Build path via Python reference
  # Build Python-style path dicts using our sacx/sacy fields
  fix1 <- data.frame(start_x = fg1$x, start_y = fg1$y, duration = fg1$duration)
  fix2 <- data.frame(start_x = fg2$x, start_y = fg2$y, duration = fg2$duration)
  sacx <- sp1[1:(nrow(sp1)-1),]
  sacy <- sp2[1:(nrow(sp2)-1),]
  p1 <- reticulate::dict(
    fix = reticulate::dict(dur = as.numeric(fix1$duration)),
    sac = reticulate::dict(
      x = as.numeric(sacx$x), y = as.numeric(sacx$y),
      lenx = as.numeric(sacx$lenx), leny = as.numeric(sacx$leny),
      theta = as.numeric(sacx$theta), rho = as.numeric(sacx$rho)
    )
  )
  p2 <- reticulate::dict(
    fix = reticulate::dict(dur = as.numeric(fix2$duration)),
    sac = reticulate::dict(
      x = as.numeric(sacy$x), y = as.numeric(sacy$y),
      lenx = as.numeric(sacy$lenx), leny = as.numeric(sacy$leny),
      theta = as.numeric(sacy$theta), rho = as.numeric(sacy$rho)
    )
  )
  nr <- nrow(M_R); nc <- ncol(M_R)
  # Recreate M_assignment with row-major ordering as in Python
  M_assignment <- matrix(0:(nr*nc - 1), nrow = nr, ncol = nc, byrow = TRUE)
  g <- mmgaze$createdirectedgraph(list(nr, nc), reticulate::r_to_py(M_R), M_assignment)
  numVert <- reticulate::py_to_r(g[[1]])
  rows <- reticulate::py_to_r(g[[2]])
  cols <- reticulate::py_to_r(g[[3]])
  weight <- reticulate::py_to_r(g[[4]])
  # Compare edge sets between R and Python
  # R edges from igraph
  r_edges <- igraph::as_edgelist(gout$g, names = TRUE)
  r_edges <- cbind(as.integer(r_edges[,1]), as.integer(r_edges[,2]))
  r_w <- igraph::E(gout$g)$weight
  r_df <- data.frame(from = r_edges[,1], to = r_edges[,2], w = as.numeric(r_w))
  r_df <- r_df[order(r_df$from, r_df$to), ]

  # Python edges (0-based -> 1-based)
  py_df <- data.frame(from = as.integer(rows) + 1L,
                      to = as.integer(cols) + 1L,
                      w = as.numeric(weight))
  py_df <- py_df[order(py_df$from, py_df$to), ]

  testthat::expect_equal(r_df$from, py_df$from)
  testthat::expect_equal(r_df$to, py_df$to)
  testthat::expect_equal(r_df$w, py_df$w, tolerance = 1e-12)
  end_idx <- nr * nc - 1
  dres <- mmgaze$dijkstra(as.integer(numVert), rows, cols, weight, 0L, as.integer(end_idx))
  py_path0 <- reticulate::py_to_r(dres[[1]])
  py_path <- as.integer(py_path0) + 1L

  testthat::expect_equal(r_path, py_path)

  # Also compare per-path elementwise differences across implementations
  # Compute cds (row, col) indices from node ids (1-based)
  r_i <- (r_path - 1L) %/% nc + 1L
  r_j <- (r_path - 1L) %% nc + 1L
  cds <- cbind(r_i, r_j)

  # R differences
  sacx <- sp1[1:(nrow(sp1)-1),]
  sacy <- sp2[1:(nrow(sp2)-1),]
  vector_d_R <- eyesim:::vector_diff_2d(sacx, sacy, "lenx", "leny", cds)
  direction_d_R <- eyesim:::angle_diff_1d(sacx$theta, sacy$theta, cds)
  duration_d_R <- eyesim:::duration_diff_1d(sacx$duration, sacy$duration, cds)
  length_d_R <- abs(eyesim:::vector_diff_1d(sacx, sacy, "rho", metric = "l1", cds))
  position_d_R <- eyesim:::vector_diff_2d(sacx, sacy, "x", "y", cds)

  # Skip per-path elementwise Python helpers (not exported in all builds);
  # we validate edges, path indices, and total path cost instead.

  # Check total path cost equals Python's dijkstra distance
  dist_py <- as.numeric(reticulate::py_to_r(dres[[2]]))
  # Sum weights along the path using move-dependent weights
  total_R <- 0
  for (k in seq_len(length(r_path) - 1L)) {
    i1 <- r_i[k]; j1 <- r_j[k]
    i2 <- r_i[k+1]; j2 <- r_j[k+1]
    if (i2 == i1 && j2 == j1 + 1L) {
      total_R <- total_R + M_R[i1, j1 + 1L]
    } else if (i2 == i1 + 1L && j2 == j1) {
      total_R <- total_R + M_R[i1 + 1L, j1]
    } else if (i2 == i1 + 1L && j2 == j1 + 1L) {
      total_R <- total_R + M_R[i1 + 1L, j1 + 1L]
    } else {
      stop("Unexpected move in path; not right/down/diag")
    }
  }
  testthat::expect_equal(total_R, dist_py, tolerance = 1e-8)
})
