test_that("compute_coi for sims", {
  withr::local_options(lifecycle_verbosity = "quiet")
  withr::local_seed(11L)
  data <- sim_biallelic(3, runif(1000, 0, 0.5), coverage = 200)

  out1 <- compute_coi(data, "sim", 25, 0.01, coi_method = "frequency", use_bins = TRUE)
  out2 <- compute_coi(data, "sim", 25, 0.01, coi_method = "variant", use_bins = TRUE)
  out3 <- compute_coi(data, "sim", 25, 0.01, coi_method = "frequency", use_bins = FALSE)
  out4 <- compute_coi(data, "sim", 25, 0.01, coi_method = "variant", use_bins = FALSE)

  # test that they return
  expect_equal(names(out1), c("coi", "probability"))
  expect_equal(names(out2), c("coi", "probability"))
  expect_equal(names(out3), c("coi", "probability"))
  expect_equal(names(out4), c("coi", "probability"))
})

test_that("optimize_coi for sims", {
  withr::local_options(lifecycle_verbosity = "quiet")
  withr::local_seed(11L)
  data <- sim_biallelic(3, runif(1000, 0, 0.5), coverage = 200)

  out1 <- optimize_coi(data, "sim", 25, 0.01, coi_method = "frequency", use_bins = TRUE)
  out2 <- optimize_coi(data, "sim", 25, 0.01, coi_method = "variant", use_bins = TRUE)
  out3 <- optimize_coi(data, "sim", 25, 0.01, coi_method = "frequency", use_bins = FALSE)
  out4 <- optimize_coi(data, "sim", 25, 0.01, coi_method = "variant", use_bins = FALSE)

  # test that they return
  expect_true(length(out1) == 1)
  expect_true(class(out1) == "numeric")
  expect_true(length(out2) == 1)
  expect_true(class(out2) == "numeric")
  expect_true(length(out3) == 1)
  expect_true(class(out3) == "numeric")
  expect_true(length(out4) == 1)
  expect_true(class(out4) == "numeric")
})

test_that("compute_coi for real", {
  withr::local_options(lifecycle_verbosity = "quiet")
  withr::local_seed(11L)
  data <- readRDS(system.file("testdata", "test_real_data.rds", package = "coiaf"))

  out1 <- compute_coi(data, "real", 25, 0.01, coi_method = "frequency", use_bins = TRUE)
  out2 <- compute_coi(data, "real", 25, 0.01, coi_method = "variant", use_bins = TRUE)
  out3 <- compute_coi(data, "real", 25, 0.01, coi_method = "frequency", use_bins = FALSE)
  out4 <- compute_coi(data, "real", 25, 0.01, coi_method = "variant", use_bins = FALSE)

  # test that they return
  expect_equal(names(out1), c("coi", "probability"))
  expect_equal(names(out2), c("coi", "probability"))
  expect_equal(names(out3), c("coi", "probability"))
  expect_equal(names(out4), c("coi", "probability"))
})

test_that("optimize_coi for real", {
  withr::local_options(lifecycle_verbosity = "quiet")
  withr::local_seed(11L)
  data <- readRDS(system.file("testdata", "test_real_data.rds", package = "coiaf"))

  out1 <- optimize_coi(data, "real", 25, 0.01, coi_method = "frequency", use_bins = TRUE)
  out2 <- optimize_coi(data, "real", 25, 0.01, coi_method = "variant", use_bins = TRUE)
  out3 <- optimize_coi(data, "real", 25, 0.01, coi_method = "frequency", use_bins = FALSE)
  out4 <- optimize_coi(data, "real", 25, 0.01, coi_method = "variant", use_bins = FALSE)

  # test that they return
  expect_true(length(out1) == 1)
  expect_true(class(out1) == "numeric")
  expect_true(length(out2) == 1)
  expect_true(class(out2) == "numeric")
  expect_true(length(out3) == 1)
  expect_true(class(out3) == "numeric")
  expect_true(length(out4) == 1)
  expect_true(class(out4) == "numeric")
})
