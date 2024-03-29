% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{anim_scanpath}
\alias{anim_scanpath}
\title{Animate a Fixation Scanpath with gganimate}
\usage{
anim_scanpath(
  x,
  bg_image = NULL,
  xlim = range(x$x),
  ylim = range(x$y),
  alpha = 1,
  anim_over = c("index", "onset"),
  type = c("points", "raster"),
  time_bin = 1
)
}
\arguments{
\item{x}{A `fixation_group` object.}

\item{bg_image}{An optional image file name to use as the background.}

\item{xlim}{The range in x coordinates (default: range of x values in the fixation group).}

\item{ylim}{The range in y coordinates (default: range of y values in the fixation group).}

\item{alpha}{The opacity of each dot (default: 1).}

\item{anim_over}{Animate over index (ordered) or onset (real time) (default: c("index", "onset")).}

\item{type}{Display as points or a raster (default: c("points", "raster")).}

\item{time_bin}{The size of the time bins (default: 1).}
}
\value{
A gganimate object representing the animated scanpath.
}
\description{
This function creates an animated visualization of a fixation scanpath using gganimate.
}
\examples{
# Create a fixation group
fg <- fixation_group(x=c(.1,.5,1), y=c(1,.5,1), onset=1:3, duration=rep(1,3))
# Animate the scanpath for the fixation group
anim_sp <- anim_scanpath(fg)
}
\seealso{
Other visualization: 
\code{\link{plot.eye_density}()},
\code{\link{plot.fixation_group}()}
}
\concept{visualization}
