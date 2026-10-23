test_that("visual angle and pixels convert exactly and invert", {
  # 2 * 60 * tan(1 deg) cm on a 52 cm / 1920 px screen
  expect_equal(visual_angle_to_px(2, 60, 52, 1920), 2 * 60 * tan(pi / 180) * 1920 / 52)
  px <- visual_angle_to_px(c(0, 1, 2, 10), distance = 57, screen_size = 40, screen_res = 1024)
  expect_equal(px[1], 0)
  expect_equal(px_to_visual_angle(px, distance = 57, screen_size = 40, screen_res = 1024), c(0, 1, 2, 10))
  # at 57.3 cm, 1 cm on screen is close to 1 degree
  expect_equal(px_to_visual_angle(1024 / 40, 57.2958, 40, 1024), 1, tolerance = 1e-4)
})

test_that("invalid display geometry is rejected", {
  expect_error(visual_angle_to_px(2, distance = 0, screen_size = 40, screen_res = 1024), "positive")
  expect_error(visual_angle_to_px(2, distance = 60, screen_size = c(40, 30), screen_res = 1024), "single")
  expect_error(visual_angle_to_px(-1, 60, 40, 1024), "deg")
  expect_error(visual_angle_to_px(180, 60, 40, 1024), "deg")
  expect_error(px_to_visual_angle(-5, 60, 40, 1024), "px")
})
