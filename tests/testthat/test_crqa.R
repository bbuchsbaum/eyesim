context("crqa wrapper")

skip_if_not_installed("crqa")

library(testthat)

# Create simple fixation groups
fg1 <- data.frame(x = 1:5, y = 1:5)
fg2 <- data.frame(x = 2:6, y = 3:7)

# Use package function directly for reference
nr <- min(nrow(fg1), nrow(fg2))
ref <- crqa::crqa(as.matrix(fg1[1:nr, 1:2]),
                  as.matrix(fg2[1:nr, 1:2]),
                  method = "mdcrqa", radius = 60)

# Call wrapper
res <- crqa(fg1, fg2, radius = 60)

expect_equal(res, ref)

