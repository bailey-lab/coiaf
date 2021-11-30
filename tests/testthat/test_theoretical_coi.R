test_that("variant method works", {
  # Define coi and interval
  coi <- 2
  interval <- seq(0, 0.5, 0.25)

  curve <- single_theoretical_coi(coi, interval)
  expect_equal(curve, c(0.00, 0.375, 0.5))
})

test_that("frequency method works", {
  # Define coi and interval
  coi <- 4
  interval <- seq(0.01, 0.50, 0.24)

  curve <- single_theoretical_coi(coi, interval, coi_method = "frequency")
  expect_equal(round(curve, 3), c(0.254, 0.362, 0.494))
})
