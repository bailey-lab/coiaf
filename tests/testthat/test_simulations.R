test_that("sim_biallelic works", {
  set.seed(1)

  # define number of loci, and distribution of minor allele frequencies
  L <- 1e5
  p <- rbeta(L, 1, 5)
  p[p > 0.5] <- 1 - p[p > 0.5]
  k <- 2

  sim1 <- sim_biallelic(coi = k, plaf = p, overdispersion = 0.01)
  expect_identical(names(sim1), c("coi", "strain_proportions", "phased", "data", "inputs"))
  expect_equal(round(sim1$strain_proportions, 2), c(0.62, 0.38))
})

test_that("sim_biallelic relatedness works", {
  set.seed(1)

  # define number of loci, and distribution of minor allele frequencies
  L <- 1e3
  p <- rbeta(L, 1, 5)
  p[p > 0.5] <- 1 - p[p > 0.5]
  k <- 3

  # two sims one with and one without relatedness
  sim1 <- sim_biallelic(coi = k, plaf = p)
  sim2 <- sim_biallelic(coi = k, plaf = p, relatedness = 0.75)

  # table up the phased counts
  tbl1 <- table(colMeans(sim1$phased))
  tbl2 <- table(colMeans(sim2$phased))

  # we would expect in related samples there will be more homozyogous calls
  expect_gt(tbl2[1], tbl1[1])
  expect_gt(tbl2[4], tbl1[4])
})
