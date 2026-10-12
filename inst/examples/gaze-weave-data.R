# Example data: two participants study and later recall six images; each
# recall roughly repeats the study path of the same image.
set.seed(1)
make_path <- function(anchor, jitter) {
  fixation_group(
    x = anchor[1] + c(0, 2, 4, 3) + rnorm(4, sd = jitter),
    y = anchor[2] + c(0, 1, -1, 2) + rnorm(4, sd = jitter),
    onset = c(0, 250, 500, 750), duration = c(200, 220, 180, 240)
  )
}
study <- expand.grid(participant = c("p1", "p2"), image = 1:6,
                     stringsAsFactors = FALSE)
anchors <- lapply(study$image, function(i) c(3 * i, 2 * (i %% 3)))
study$fixgroup <- lapply(anchors, make_path, jitter = 0.2)
recall <- study[c("participant", "image")]
recall$fixgroup <- lapply(anchors, make_path, jitter = 0.6)
