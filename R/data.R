#' Study-phase fixations from Wynn et al.
#'
#' Fixations recorded while participants studied photographs in a recognition
#' memory experiment, one row per participant and image.
#'
#' Each image was studied once in each of four study blocks, and a row pools
#' the fixations of all four viewings: within a row's \code{fixgroup}, onsets
#' restart at each viewing. The pooled fixations suit density maps (for
#' example study templates for \code{\link{template_similarity}}); for
#' scanpath measures, which join consecutive rows, split the viewings first.
#' Fixations outside the 800 x 600 pixel image area were removed and
#' coordinates are relative to its top-left corner.
#'
#' @format A grouped tibble (grouped by \code{ImageNumber} and \code{Subject}) with 12,240 rows and 8 columns:
#' \describe{
#'   \item{Subject}{Participant number.}
#'   \item{fixgroup}{A \code{\link{fixation_group}} with the pooled study
#'     fixations (\code{x}, \code{y} in pixels; \code{duration} and
#'     \code{onset} in ms).}
#'   \item{ImageVersion}{Image and version shown, e.g. \code{"1_A"}.}
#'   \item{Version}{Image version (\code{"A"} or \code{"B"}).}
#'   \item{Trial}{Trial index of the first viewing.}
#'   \item{ImageNumber}{Image identifier.}
#'   \item{blockstudy}{Study block of the first viewing.}
#'   \item{ImageNumberS}{Participant-image key, \code{"<Subject>_<ImageNumber>"},
#'     matching \code{wynn_test$ImageNumberS}.}
#' }
#'
#' @source Wynn et al. eye-tracking study; built from the raw fixation report
#'   by \code{data-raw/process_dat.R} in the source repository.
#' @seealso \code{\link{wynn_test}}
#'
#' @examples
#' data(wynn_study)
#' head(wynn_study)
#'
"wynn_study"

#' Test-phase fixations from Wynn et al.
#'
#' Fixations recorded during the recognition test of the same experiment as
#' \code{\link{wynn_study}}, one row per participant and test image (a single
#' viewing). Each probe is either the studied version of an image
#' (\code{"old"}) or the other version (\code{"lure"}).
#'
#' @format A grouped tibble (grouped by \code{ImageNumber} and \code{Subject}) with 12,240 rows and 11 columns:
#' \describe{
#'   \item{Subject}{Participant number.}
#'   \item{fixgroup}{A \code{\link{fixation_group}} with the test fixations
#'     (\code{x}, \code{y} in pixels; \code{duration} and \code{onset} in ms).}
#'   \item{ImageVersion}{Image and version shown at test.}
#'   \item{Saliency}{Saliency condition of the test image.}
#'   \item{Accuracy}{Response accuracy (1 correct, 0 incorrect).}
#'   \item{Version}{Image version (\code{"A"} or \code{"B"}).}
#'   \item{Trial}{Test trial index.}
#'   \item{Duration}{Duration condition of the test trial in ms (250, 500 or
#'     750).}
#'   \item{ImageNumber}{Image identifier.}
#'   \item{ImageNumberS}{Participant-image key, matching
#'     \code{wynn_study$ImageNumberS}.}
#'   \item{probe_type}{\code{"old"} if the test image is the studied version,
#'     \code{"lure"} otherwise.}
#' }
#'
#' @source Wynn et al. eye-tracking study; built from the raw fixation report
#'   by \code{data-raw/process_dat.R} in the source repository.
#' @seealso \code{\link{wynn_study}}
#'
#' @examples
#' data(wynn_test)
#' table(wynn_test$probe_type)
#'
"wynn_test"
