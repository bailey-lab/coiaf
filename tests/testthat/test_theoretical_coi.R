context("theoretical coi")

#------------------------------------------------
test_that("theoretical_coi works", {
  # define coi and interval
  coi = 2
  interval = seq(0, 0.5, 0.25)

  curve <- theoretical_coi(coi, interval)
  expect_equal(curve, c(0.00, 0.375, 0.5))
})
