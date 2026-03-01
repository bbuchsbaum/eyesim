


#' Construct an Eye-Movement Data Frame
#'
#' The `eye_table` function creates a `data.frame` to store eye-movement data
#' with columns for x and y coordinates, fixation duration, onset time, and
#' an optional grouping variable. Additional variables can also be retained.
#' The function filters the data based on the clip bounds and can compute
#' relative coordinates.
#'
#' @importFrom dplyr group_by select filter as_tibble rename_with all_of across summarise left_join
#' @importFrom magrittr %>%
#' @importFrom assertthat assert_that
#'
#' @param x A character or symbol representing the column for the x coordinates in the source data.
#' @param y A character or symbol representing the column for the y coordinates in the source data.
#' @param duration A character or symbol representing the column for fixation durations in the source data.
#' @param onset A character or symbol representing the column for fixation onset times in the source data.
#' @param groupvar A character or symbol representing the column for the grouping variable in the source data.
#' @param vars A character vector of additional variable names to retain, or NULL (default) if no additional variables are needed.
#' @param data The source `data.frame` containing the eye-movement data.
#' @param clip_bounds A numeric vector of length 4 representing the clip bounds for the field of view in the form c(xmin, xmax, ymin, ymax). Default is c(0, 1280, 0, 1280).
#' @param relative_coords A logical value indicating whether to compute relative coordinates (TRUE by default). If TRUE, x and y coordinates will be transformed based on the clip_bounds.
#'
#' @return A `data.frame` of class "eye_table" with columns for x and y coordinates, fixation duration, onset time, the grouping variable, and any additional specified variables. The data frame will also have an "origin" attribute containing the center coordinates of the field of view.
#' @export
#'
#' @details
#' The `eye_table` function first checks that the input `data` is a `data.frame` and then renames the columns specified by x, y, duration, and onset to their standard names. The function then filters the data based on the specified clip_bounds, ensuring that all x and y coordinates fall within the bounds. If relative_coords is TRUE, the x and y coordinates will be transformed to be relative to the clip_bounds.
#'
#' The function groups the data by the specified grouping variable and constructs a fixation_group object for each group, which is added to the output data frame as a new "fixgroup" column. The output data frame retains the specified additional variables and is assigned a class of "eye_table". The "origin" attribute of the output data frame contains the center coordinates of the field of view, which are computed based on the clip_bounds and whether relative_coords is TRUE or FALSE.
eye_table <- function(x, y, duration, onset, groupvar, vars=NULL, data,
                      clip_bounds=c(0,1280, 0,1280), relative_coords=TRUE) {

  assertthat::assert_that(inherits(data, "data.frame"))

  # Validate that required columns exist
  required <- c(x, y, duration, onset, groupvar)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }

  # Validate that coordinate/timing columns are numeric
  for (col in c(x, y, duration, onset)) {
    if (!is.numeric(data[[col]])) {
      stop("Column '", col, "' must be numeric, but is ", class(data[[col]])[1])
    }
  }

  # Validate non-negative duration and onset
  if (any(data[[duration]] < 0, na.rm = TRUE)) {
    stop("Column '", duration, "' contains negative values. Durations must be non-negative.")
  }
  if (any(data[[onset]] < 0, na.rm = TRUE)) {
    stop("Column '", onset, "' contains negative values. Onset times must be non-negative.")
  }

  colmapping <- c("x","y","duration", "onset")
  names(colmapping) <- c(x,y,duration,onset)


  data <- data %>% dplyr::rename_with(.cols=all_of(c(x,y,duration,onset)),
                               .fn = function(x){colmapping[x]}) %>% as_tibble()

  data <- if (is.null(vars)) {
    data %>% select(all_of(c("x","y","duration", "onset", groupvar)))
  } else {
    data %>% select(all_of(c("x","y","duration", "onset", vars, groupvar)))
  }


  xdir <- sign(clip_bounds[2] - clip_bounds[1])
  ydir <- sign(clip_bounds[4] - clip_bounds[3])
  xr <- sort(c(clip_bounds[1], clip_bounds[2]))
  yr <- sort(c(clip_bounds[3], clip_bounds[4]))

  data <- data %>% filter(x >= xr[1] & x <= xr[2]
                        & y >= yr[1] & y <= yr[2])


  if (relative_coords) {
    data <- data %>% mutate(x = (x - clip_bounds[1]) * xdir, y = (y - clip_bounds[3]) * ydir)
  }


  res <- data %>%
    group_by(across(all_of(groupvar))) %>%
    summarise(fixgroup = list(fixation_group(.data[["x"]], .data[["y"]], .data[["duration"]], .data[["onset"]])),
              .groups = "drop")

  if (!is.null(vars)) {
    # Re-attach the first value of each kept variable per group
    var_data <- data %>%
      group_by(across(all_of(groupvar))) %>%
      summarise(across(all_of(vars), ~ .x[1]), .groups = "drop")
    res <- dplyr::left_join(res, var_data, by = groupvar)
  }

  class(res) <- c("eye_table", class(res))

  if (relative_coords) {
    xr1 <- (xr[1] - clip_bounds[1]) * xdir
    xr2 <- (xr[2] - clip_bounds[1]) * xdir
    yr1 <- (yr[1] - clip_bounds[3]) * ydir
    yr2 <- (yr[2] - clip_bounds[3]) * ydir
    attr(res, "origin") <- c((xr1+xr2)/2, (yr1+yr2)/2)
  } else {
    attr(res, "origin") <- c((xr[1]+xr[2])/2, (yr[1]+yr[2])/2)
  }
  res

}


#' @export
print.eye_table <- function(x, ...) {
  ngroups <- nrow(x)
  nfix <- sum(vapply(x$fixgroup, nrow, integer(1)))
  origin <- attr(x, "origin")
  cat("Eye table:", ngroups, "groups,", nfix, "total fixations\n")
  if (!is.null(origin)) {
    cat("  origin: (", round(origin[1], 1), ", ", round(origin[2], 1), ")\n", sep = "")
  }
  NextMethod()
}

#' Reapply the 'eye_table' Class to an Object
#'
#' The `as_eye_table` function is a simple utility function that reapplies
#' the 'eye_table' class to a given object if it is not already part of
#' its class attribute. This function is not intended for serious use.
#'
#' @param x The input object to which the 'eye_table' class should be (re)applied.
#'
#' @return The input object with the 'eye_table' class (re)applied to its class attribute.
#' @export
#'
#' @examples
#' # Create a simple data.frame
#' df <- data.frame(x = 1:5, y = 6:10)
#'
#' # Apply the 'eye_table' class to df
#' eye_table_df <- as_eye_table(df)
#' class(eye_table_df)
as_eye_table <- function(x) {
  cls <- class(x)
  if (!"eye_table" %in% cls) {
    class(x) <- c("eye_table", cls)
  }
  x
}

#' Subset an 'eye_table' Object
#'
#' The `[.eye_table` function is an S3 method for subsetting 'eye_table' objects.
#' It uses `NextMethod()` to perform the subsetting operation and then reapplies
#' the 'eye_table' class to the result using the `as_eye_table` function. This
#' ensures that the returned object still has the 'eye_table' class after subsetting.
#'
#' @param x The 'eye_table' object to be subsetted.
#' @param i Row indices for subsetting.
#' @param j Column indices for subsetting.
#' @param drop A logical value indicating whether to drop dimensions that have
#'   only one level after subsetting. Defaults to `FALSE`.
#'
#' @return An 'eye_table' object after subsetting.
#' @export
#'
#' @examples
#' # Create an 'eye_table' object
#' df <- data.frame(x = 1:5, y = 6:10)
#' eye_table_df <- as_eye_table(df)
#'
#' # Subset the 'eye_table' object
#' subset_eye_table <- eye_table_df[1:3,]
#' class(subset_eye_table)
`[.eye_table` <- function(x, i, j, drop = FALSE) {
  as_eye_table(NextMethod())
}

#' Generate a Simulated Eye-Movement Data Frame
#'
#' The `simulate_eye_table` function generates a simulated `eye_table` object
#' containing eye-movement data with columns for x and y coordinates, fixation duration,
#' onset time, and an optional grouping variable.
#'
#' @param n_fixations The number of fixations to simulate.
#' @param n_groups The number of groups to simulate.
#' @param clip_bounds A numeric vector of length 4 representing the clip bounds for the field of view in the form c(xmin, xmax, ymin, ymax). Default is c(0, 1280, 0, 1280).
#' @param relative_coords A logical value indicating whether to compute relative coordinates (TRUE by default). If TRUE, x and y coordinates will be transformed based on the clip_bounds.
#'
#' @return A `data.frame` of class "eye_table" with simulated data.
#' @export
#'
#' @examples
#' sim_eye_table <- simulate_eye_table(n_fixations = 100, n_groups = 10)
#' @importFrom stats runif rnorm
simulate_eye_table <- function(n_fixations, n_groups, clip_bounds=c(0,1280, 0,1280), relative_coords=TRUE) {

  # Simulate eye-movement data with per-group onset times
  fix_per_group <- as.integer(n_fixations / n_groups)
  data <- do.call(rbind, lapply(seq_len(n_groups), function(g) {
    data.frame(
      x = runif(fix_per_group, min = clip_bounds[1], max = clip_bounds[2]),
      y = runif(fix_per_group, min = clip_bounds[3], max = clip_bounds[4]),
      duration = abs(rnorm(fix_per_group, mean = 300, sd = 50)),
      onset = cumsum(abs(rnorm(fix_per_group, mean = 400, sd = 100))),
      groupvar = factor(g, levels = seq_len(n_groups))
    )
  }))

  # Create an eye_table object from the simulated data
  sim_eye_table <- eye_table(
    x = "x", y = "y", duration = "duration", onset = "onset",
    groupvar = "groupvar", data = data,
    clip_bounds = clip_bounds, relative_coords = relative_coords
  )

  return(sim_eye_table)
}




