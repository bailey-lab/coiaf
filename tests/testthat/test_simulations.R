context("simulation")

#------------------------------------------------
test_that("sim_biallelic works", {

  set.seed(1)

  # define number of loci, and distribution of minor allele frequencies
  L <- 1e5
  p <- rbeta(L, 1, 5)
  p[p > 0.5] <- 1 - p[p > 0.5]
  k <- 2

  sim1 <- sim_biallelic(COI = k, PLAF = p, overdispersion = 0.01)
  expect_identical(names(sim1), c("COI", "strain_proportions", "phased", "data"))
  expect_equal(round(sim1$strain_proportions, 2), c(0.62, 0.38))

})
