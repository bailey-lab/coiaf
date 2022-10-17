set.seed(375)
data <- sim_biallelic(7)

test_that("compute_coi() args deprecated", {
  expect_snapshot(compute_coi(data, "sim", comparison = "end", use_bins = TRUE))
  expect_snapshot(
    compute_coi(data, "sim", distance = "abs_sum", use_bins = TRUE)
  )
  expect_snapshot_warning(compute_coi(data, "sim", use_bins = TRUE))
  expect_snapshot(compute_coi(data, "sim", bin_size = 100, use_bins = TRUE))
})

test_that("optimize_coi() args deprecated", {
  expect_snapshot(
    optimize_coi(data, "sim", distance = "abs_sum", use_bins = TRUE)
  )
  expect_snapshot_warning(optimize_coi(data, "sim", use_bins = TRUE))
  expect_snapshot(optimize_coi(data, "sim", bin_size = 100, use_bins = TRUE))
})

test_that("compute_coi_regression() args deprecated", {
  expect_snapshot_warning(
    compute_coi_regression(data, "sim", distance = "abs_sum")
  )
})

test_that("optimize_coi_regression() args deprecated", {
  expect_snapshot_warning(
    optimize_coi_regression(data, "sim", distance = "abs_sum")
  )
})

test_that("bootstrap_ci() args deprecated", {
  expect_snapshot_warning(bootstrap_ci(data, use_bins = TRUE, replicates = 50))
  expect_snapshot(
    bootstrap_ci(data, "sim", bin_size = 100, use_bins = TRUE, replicates = 10)
  )
})
