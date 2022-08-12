# A simple simulation we will leverage to test our function
simple_sim <- withr::with_preserve_seed(sim_biallelic(3))

test_that("if coverage is a single number, expand across all loci", {
  sim <- sim_biallelic(coi = 5, coverage = 75)
  coverage <- sim$data$coverage
  expect_true(unique(coverage) == 75)
})

test_that("simulator returns object of correct class", {
  expect_type(simple_sim, "list")
  expect_s3_class(simple_sim, "sim", exact = TRUE)
})

test_that("simulator returns a list of tibbles", {
  expect_s3_class(simple_sim$parameters, c("tbl_df", "tbl", "dta.frame"))
  expect_s3_class(
    simple_sim$strain_proportions,
    c("tbl_df", "tbl", "dta.frame")
  )
  expect_s3_class(simple_sim$phased_haplotypes, c("tbl_df", "tbl", "dta.frame"))
  expect_s3_class(simple_sim$data, c("tbl_df", "tbl", "dta.frame"))
})

test_that("parameter tibble reflects inputs", {
  expect_equal(
    sim_biallelic(
      coi = 3,
      alpha = 0.5,
      overdispersion = 0.7,
      relatedness = 0.2,
      epsilon = 0.1
    )$parameters,
    tibble::tribble(
      ~parameter, ~value,
      "coi", 3,
      "alpha", 0.5,
      "overdispersion", 0.7,
      "relatedness", 0.2,
      "epsilon", 0.1
    )
  )
})

test_that("length of strain proportions is the COI", {
  coi <- 10
  sim <- sim_biallelic(coi)
  expect_length(sim$strain_proportions$proportion, coi)
})

test_that("strain proporitons sum to 1", {
  expect_true(sum(sim_biallelic(5)$strain_proportions$proportion) == 1)
  expect_true(sum(sim_biallelic(15)$strain_proportions$proportion) == 1)
})

test_that("phased haplotypes dimensions are correct", {
  sim <- sim_biallelic(coi = 5, plmaf = runif(20, 0, 0.5))
  expect_equal(dim(sim$phased_haplotypes), c(5, 20))
})

test_that("data columns are correct", {
  expect_equal(
    colnames(simple_sim$data),
    c("plmaf", "coverage", "counts", "wsmaf")
  )

  plmaf <- runif(20, 0, 0.5)
  expect_identical(sim_biallelic(coi = 2, plmaf = plmaf)$data$plmaf, plmaf)
})

test_that("relatedness works as expected", {
  withr::local_seed(1)

  # Define number of loci, and PLMAF
  L <- 1e3
  p <- rbeta(L, 1, 5)
  p[p > 0.5] <- 1 - p[p > 0.5]
  k <- 3

  # Two simulations: one with and one without relatedness
  nonrelated_sim <- sim_biallelic(coi = k, plmaf = p)
  related_sim <- sim_biallelic(coi = k, plmaf = p, relatedness = 0.75)

  # Table up the phased counts
  nonrelated_tbl <- table(colMeans(nonrelated_sim$phased_haplotypes))
  related_tbl <- table(colMeans(related_sim$phased_haplotypes))

  # We would expect in related samples there will be more homozyogous calls
  expect_gt(related_tbl[1], nonrelated_tbl[1])
  expect_gt(related_tbl[4], nonrelated_tbl[4])
})

# Plotting test cases ----------------------------------------------------------
plot_sim <- withr::with_seed(500, sim_biallelic(3, runif(100, 0, 0.5)))

test_that("plot and autoplot methods work", {
  vdiffr::expect_doppelganger("plot method works", plot(plot_sim))
  vdiffr::expect_doppelganger("autoplot method works", autoplot(plot_sim))
})

test_that("can pass along ggplot2 parameters", {
  vdiffr::expect_doppelganger(
    "can pass additional params",
    plot(plot_sim, color = "red")
  )
})
