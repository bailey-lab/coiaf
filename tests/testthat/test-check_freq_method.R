test_that("check_freq_method catches COI = 1", {
  withr::local_seed(101)
  data <- sim_biallelic(1, runif(1000, 0, 0.5), coverage = 200)

  check <- check_freq_method(data$data$wsmaf, data$data$plmaf, 0)
  expect_lt(check["variant"], check["lower_ci"])
})

test_that("when no loci return NaN with note and estimated coi of 1", {
  withr::local_seed(202)
  data <- sim_biallelic(1, runif(1000, 0, 0.5), coverage = 200)

  # Buckets
  compute <- compute_coi(data, "sim", coi_method = "frequency", use_bins = TRUE)
  optim <- optimize_coi(data, "sim", coi_method = "frequency", use_bins = TRUE)

  # Regression
  compute_reg <- compute_coi(
    data,
    "sim",
    coi_method = "frequency",
    use_bins = FALSE
  )
  optim_reg <- optimize_coi(
    data,
    "sim",
    coi_method = "frequency",
    use_bins = FALSE
  )

  # Discrete buckets
  expect_true(is.nan(compute$coi))
  expect_named(
    compute,
    c("coi", "probability", "notes", "estimated_coi", "num_variant_loci")
  )
  expect_equal(compute$estimated_coi, 1)

  # Continuous buckets
  expect_true(is.nan(optim))
  expect_named(
    attributes(optim),
    c("notes", "estimated_coi", "num_variant_loci")
  )
  expect_equal(attributes(optim)$estimated_coi, 1)

  # Discrete regression
  expect_true(is.nan(compute_reg$coi))
  expect_named(
    compute_reg,
    c("coi", "probability", "notes", "estimated_coi", "num_variant_loci")
  )
  expect_equal(compute_reg$estimated_coi, 1)

  # Continuous regression
  expect_true(is.nan(optim_reg))
  expect_named(
    attributes(optim_reg),
    c("notes", "estimated_coi", "num_variant_loci")
  )
  expect_equal(attributes(optim_reg)$estimated_coi, 1)
})

test_that("when too few loci returns NaN with note and estimated coi", {
  withr::local_seed(303)
  data <- sim_biallelic(1, runif(1000, 0, 0.5), coverage = 200, epsilon = 0.01)

  # Buckets
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

  # Regression
  compute_reg <- compute_coi(
    data,
    "sim",
    seq_error = 0.015,
    coi_method = "frequency",
    use_bins = FALSE
  )
  optim_reg <- optimize_coi(
    data,
    "sim",
    seq_error = 0.015,
    coi_method = "frequency",
    use_bins = FALSE
  )

  # Discrete buckets
  expect_true(is.nan(compute$coi))
  expect_named(
    compute,
    c("coi", "probability", "notes", "estimated_coi", "num_variant_loci", "expected_num_loci")
  )
  expect_gt(compute$estimated_coi, 1)

  # Continuous buckets
  expect_true(is.nan(optim))
  expect_named(
    attributes(optim),
    c("notes", "estimated_coi", "num_variant_loci", "expected_num_loci")
  )
  expect_gt(attributes(optim)$estimated_coi, 1)

  # Discrete regression
  expect_true(is.nan(compute_reg$coi))
  expect_named(
    compute_reg,
    c("coi", "probability", "notes", "estimated_coi", "num_variant_loci", "expected_num_loci")
  )
  expect_gt(compute_reg$estimated_coi, 1)

  # Continuous regression
  expect_true(is.nan(optim_reg))
  expect_named(
    attributes(optim_reg),
    c("notes", "estimated_coi", "num_variant_loci", "expected_num_loci")
  )
  expect_gt(attributes(optim_reg)$estimated_coi, 1)
})
