test_that("helper returns a vector", {
  expect_vector(
    single_theoretical_coi(
      coi = 2,
      plmaf = seq(0, 0.5, 0.01),
      coi_method = "variant"
    ),
    ptype = double(),
    size = 51
  )
})

test_that("variant method works", {
  curve <- single_theoretical_coi(2, seq(0, 0.5, 0.25), "variant")
  expect_equal(curve, c(0.00, 0.375, 0.5))
})

test_that("frequency method works", {
  curve <- single_theoretical_coi(4, seq(0.01, 0.50, 0.24), "frequency")
  expect_equal(round(curve, 3), c(0.254, 0.362, 0.494))
})

test_that("coi_method argument is restricted", {
  expect_snapshot(theoretical_coi(1:5, seq(0, 0.5, 0.1)))
  expect_snapshot(
    theoretical_coi(1:5, seq(0, 0.5, 0.1), coi_method = "frequency")
  )
  expect_snapshot_error(theoretical_coi(1:5, coi_method = "wrong method"))
})

test_that("theoretical_coi returns a data frane", {
  expect_s3_class(theoretical_coi(1:2), class = "data.frame")
})

test_that("theoretical_coi column names are correct", {
  expect_equal(
    colnames(theoretical_coi(1:10)),
    c(paste0("coi_", 1:10), "plmaf")
  )
})
