test_that("check_freq_method catches COI = 1", {
  set.seed(101)
  data <- sim_biallelic(1, runif(1000, 0, 0.5), coverage = 200)

  expect_false(check_freq_method(data$data$wsmaf, data$data$plmaf, 0))
})

test_that("when no loci return NaN with note and estimated coi of 1", {
  set.seed(202)
  data <- sim_biallelic(1, runif(1000, 0, 0.5), coverage = 200)

  compute <- compute_coi(data, "sim", coi_method = "frequency", use_bins = TRUE)
  optim <- optimize_coi(data, "sim", coi_method = "frequency", use_bins = TRUE)

  # Discrete
  expect_true(is.nan(compute$coi))
  expect_named(compute, c("coi", "probability", "notes", "estimated_coi"))
  expect_equal(compute$estimated_coi, 1)

  # Continuous
  expect_true(is.nan(optim))
  expect_named(attributes(optim), c("notes", "estimated_coi"))
  expect_equal(attributes(optim)$estimated_coi, 1)
})

test_that("when too few loci returns NaN with note and estimated coi", {
  set.seed(303)
  data <- sim_biallelic(1, runif(1000, 0, 0.5), coverage = 200, epsilon = 0.01)

  compute <- compute_coi(
    data,
    "sim",
    seq_error = 0.015,
    coi_method = "frequency",
    use_bins = TRUE
  )

  optim <- optimize_coi(
    data,
    "sim",
    seq_error = 0.015,
    coi_method = "frequency",
    use_bins = TRUE
  )

  # Discrete
  expect_true(is.nan(compute$coi))
  expect_named(compute, c("coi", "probability", "notes", "estimated_coi"))
  expect_gt(compute$estimated_coi, 1)

  # Continuous
  expect_true(is.nan(optim))
  expect_named(attributes(optim), c("notes", "estimated_coi"))
  expect_gt(attributes(optim)$estimated_coi, 1)
})
