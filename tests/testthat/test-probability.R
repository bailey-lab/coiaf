test_that("rdirichlet works", {
  withr::local_seed(1)
  out <- rdirichlet(rep(1, 3))
  expect_equal(round(out, 3), c(0.04, 0.49, 0.47))
})

test_that("rbetabinomial works", {
  withr::local_seed(1)
  out <- rbetabinom(1, 10, 1, 1)
  expect_equal(out, 7)
})
