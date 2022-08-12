computed_midpoints <- function(data, cuts) {
  find_cut_midpoints(data, cuts) %>%
    dplyr::distinct() %>%
    dplyr::pull(midpoints) %>%
    sort()
}

test_that("can find cut midpoints", {
  withr::local_seed(100)
  plmaf <- runif(100, 0, 0.5)
  data <- data.frame(cuts = Hmisc::cut2(plmaf, m = 25))

  expect_equal(
    computed_midpoints(data, cuts),
    c(0.08455, 0.20625, 0.31635, 0.43475)
  )
})

test_that("can find cut midpoints even if there is a cut with one value", {
  data <- data.frame(cuts = as.factor(c(
    "[0.1861,0.311)", "[0.0219,0.186)", "[0.3110,0.355]", "0.0000", "0.4750"
  )))

  expect_equal(
    computed_midpoints(data, cuts),
    c(0.0000, 0.10395, 0.24855, 0.333, 0.4750),
    tolerance = 0.001
  )
})

test_that("cuts can end with a `(` or a `[`", {
  data <- data.frame(cuts = as.factor(c("[0.1861,0.311)", "[0.0219,0.186]")))
  expect_equal(computed_midpoints(data, cuts), c(0.10395, 0.24855))
})

test_that("cut midpoints work for a large data set (#18)", {
  # The data set we use here used to fail when running our algorithms
  data <- readRDS("issue_18.rds")
  expect_snapshot(
    compute_coi(
      tibble::as_tibble(data),
      "real",
      seq_error = data$seq_error,
      bin_size = data$bin_size,
      coi_method = data$coi_method,
      use_bins = TRUE
    )
  )
})
