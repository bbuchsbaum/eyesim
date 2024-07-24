
#' @noRd
merge_fix <- function(fg) {
  v_x = sum(fg$lenx)
  v_y = sum(fg$leny)
  polar <- cart2pol(v_x, v_y)
  tibble(index=fg$index[1], x=fg$x[1], y=fg$y[1], duration=sum(fg$duration), onset=fg$onset[1],
         group_index=fg$group_index[1], lenx=v_x, leny=v_y, rho=polar[,1], theta=polar[,2])
}

#' @noRd
simplify_dir <- function(x, TDir, TDur) {
  G <- do.call(rbind, purrr::map(1:(nrow(x)-1), function(i) {
    ang <- calcangle(c(x$lenx[i], x$leny[i]), c(x$lenx[i+1], x$leny[i+1]))
    if (ang < TDir & x$duration[i+1] < TDur) {
      cbind(i, i+1)
    } else {
      NULL
    }
  }) )

  if (is.null(G) || nrow(G) == 0) {
    return(x)
  }


  gs <- igraph::groups(igraph::components(igraph::graph_from_data_frame(G)))
  m <- as.integer(unlist(gs))

  merged <- do.call(rbind, lapply(gs, function(fg) {
    merge_fix(x[as.integer(fg),])
  }))

  if (length(m) != nrow(x)) {
    singletons <- (1:nrow(x))[!(1:nrow(x) %in% m)]
    rem <- x[singletons,]
    out <- rbind(merged, rem) %>% arrange(onset) %>% mutate(index=1:n())
  } else {
    out <- merged %>% arrange(onset) %>% mutate(index=1:n())
  }

}


#' @noRd
duration_diff_1d <- function(dur1, dur2, path, cds) {
  dur1 <- dur1[cds[,1]]
  dur2 <- dur2[cds[,2]]

  adiff <- abs(dur1-dur2)
  denom <- sapply(1:length(dur1), function(i) max(dur1[i], dur2[i]))
  adiff/denom
}


#' @noRd
angle_diff_1d <- function(theta1, theta2, path, cds) {
  theta1 <- theta1[cds[,1]]
  theta2 <- theta2[cds[,2]]

  theta1 <- ifelse(theta1 < 0, pi + (pi + theta1), theta1)
  theta2 <- ifelse(theta2 < 0, pi + (pi + theta2), theta2)

  adiff <- abs(theta1-theta2)
  adiff <- ifelse(adiff > pi, 2 * pi - adiff, adiff)
  adiff
}

#' @noRd
vector_diff_1d <- function(x,y,path, v1, metric=c("l1", "l2"), cds) {
  metric <- match.arg(metric)
  x1 = x[[v1]]
  x2 = y[[v1]]

  #cds <- arrayInd(as.integer(path), c(nrow(x),nrow(y)))
  #cds <- cbind(cds[,2], cds[,1])

  if (metric == "l2") {
    sqrt((x1[cds[,1]] - x2[cds[,2]])^2)
  } else {
    x1[cds[,1]] - x2[cds[,2]]
  }

}

#' @noRd
vector_diff_2d <- function(x,y,path, v1,v2, cds) {
  x1 = x[[v1]]
  x2 = y[[v1]]
  y1 = x[[v2]]
  y2 = y[[v2]]

  #cds <- arrayInd(as.integer(path), c(nrow(x),nrow(y)))
  #cds <- cbind(cds[,2], cds[,1])

  sqrt((x1[cds[,1]] - x2[cds[,2]])^2 + (y1[cds[,1]] - y2[cds[,2]])^2)

}


#' @noRd
create_graph <- function(x,y) {
  M <- proxy::dist(cbind(x$lenx, x$leny), cbind(y$lenx, y$leny))
  M_assignment <- matrix(seq(nrow(M) * ncol(M)), nrow(M), ncol(M), byrow=TRUE)

  # loop through every node rowwise
  out <- do.call(rbind, lapply(seq(1, nrow(M)), function(i) {
     # loop through every node columnwise
     do.call(rbind, lapply(seq(1, ncol(M)), function(j) {
       currentNode = (i-1) * ncol(M) + j
        # if in the last (bottom) row, only go right
        if (i == nrow(M) && (j < ncol(M))) {
          cbind(currentNode, currentNode+1,(M[i, j + 1]))
        } else if (i < nrow(M) && j == ncol(M)) {
          cbind(currentNode, currentNode+ncol(M),M[i + 1, j])
        } else if ((i == nrow(M)) && (j == ncol(M))) {
          cbind(currentNode, currentNode,0)
        } else {
            cbind(rep(currentNode,3),
                  c(currentNode+1,currentNode+ncol(M), currentNode+ncol(M)+1),
                  c(M[i, j + 1],M[i + 1, j], M[i + 1, j + 1]))
        }
     }))
  }))

  g <- igraph::graph_from_data_frame(out[,1:2])
  igraph::E(g)$weight <- out[,3]
  spath <- igraph::shortest_paths(g, 1, max(igraph::V(g)), weights=igraph::E(g)$weight, output="both", predecessors = TRUE)
  list(g=g, vpath=spath$vpath[[1]], epath=spath$epath[[1]], M=M, M_assignment=M_assignment, pred=spath$predecessors)

}


#' @export
install_multimatch <- function() {
  reticulate::py_install("multimatch_gaze", pip=TRUE)
}


#' multi_match
#'
#' @description
#' \code{multi_match} compares two scanpaths based on vector, direction, length, position, and duration.
#'
#' @param x the first \code{scanpath}
#' @param y the second \code{scanpath}
#' @param screensize a two element vector indicating screen size
#' @export
multi_match <- function(x,y, screensize) {

  if (any(diff(x$onset) <= 0)) {
    stop("multi_match: x `onset` vector must be strictly increasing")
  }

  if (any(diff(x$onset) <= 0)) {
    stop("multi_match: y `onset` vector must be strictly increasing")
  }

  if (nrow(x) < 3 || nrow(y) < 3) {
    warning("multi_match requires 3 or more coordinates in each scanpath, returning NAs")
    c(mm_vector=NA, mm_direction=NA,
      mm_length=NA, mm_position=NA,
      mm_duration=NA)

  }

  sacx <- x[1:(nrow(x)-1),]
  sacy <- y[1:(nrow(y)-1),]

  gout <- create_graph(sacx,sacy)
  p <- as.integer(gout$vpath)
  rnum <- ceiling(p / ncol(gout$M))
  cnum <-  p %% ncol(gout$M)
  cnum[cnum==0] <- ncol(gout$M)

  cds <- cbind(rnum, cnum)

  vector_d <- vector_diff_2d(sacx, sacy, as.integer(gout$vpath), "lenx", "leny", cds)
  vector_sim <-  1- median(vector_d) / (2 * sqrt(screensize[1]^2 + (screensize[2]^2)))

  direction_d <- angle_diff_1d(sacx$theta,sacy$theta, as.integer(gout$vpath), cds)
  direction_sim <- 1 - median(direction_d) / pi

  duration_d <- duration_diff_1d(sacx$duration,sacy$duration, as.integer(gout$vpath), cds)
  duration_sim <- 1 - median(duration_d)

  length_d <- abs(vector_diff_1d(sacx, sacy, as.integer(gout$vpath), "rho", metric="l1", cds))
  length_sim <- 1 - (median(length_d)) / (sqrt(screensize[1]^2 + (screensize[2]^2)))

  position_d <- vector_diff_2d(sacx, sacy, as.integer(gout$vpath), "x", "y", cds)
  position_sim <- 1 - (median(position_d)) / (sqrt(screensize[1]^2 + (screensize[2]^2)))

  c(mm_vector=vector_sim, mm_direction=direction_sim,
    mm_length=length_sim, mm_position=position_sim,
    mm_duration=duration_sim)
}


#' @noRd
py_multi_match <- function(fg1, fg2,
                        screensize,
                        grouping=FALSE,
                        tdir=25,
                        tdur=.05,
                        tamp=100) {

  if (!requireNamespace("reticulate")) {
    stop("multi_match requires access to python library `multimatch_gaze` via `reticulate`, please install")
  }

  if (!exists("mmgaze")) {
    mmgaze <<- try(reticulate::import("multimatch_gaze"))
    if (inherits(mmgaze, "try-error")) {
      stop("cannot load python module `multimatch_gaze`")
    }
  }

  fg1 <- fg1 %>% arrange(onset)
  fg2 <- fg2 %>% arrange(onset)

  fix1 <- fg1[, c("x", "y", "duration")]
  fix2 <- fg2[, c("x", "y", "duration")]

  colnames(fix1) <- c("start_x", "start_y", "duration")
  colnames(fix2) <- c("start_x", "start_y", "duration")

  ret <- mmgaze$docomparison(fix1, fix2, as.integer(screensize), grouping=grouping, TDir=tdir,TDur=tdur,TAmp=tamp)
  names(ret) <- c("mm_vector", "mm_direction", "mm_length", "mm_position", "mm_duration")
  ret
}

