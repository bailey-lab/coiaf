context("theoretical coi")

#------------------------------------------------
test_that("single_theoretical_coi works", {
  # define coi and interval
  coi = 2
  interval = seq(0, 0.5, 0.25)

  curve <- single_theoretical_coi(coi, interval)
  expect_equal(curve, c(0.00, 0.375, 0.5))
})
