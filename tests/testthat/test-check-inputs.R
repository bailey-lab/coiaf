input <- tibble::tibble(plmaf = c(0.26, NA), wsmaf = c(0.35, NA))

# Simulate data
set.seed(101)
sim_input <- sim_biallelic(coi = 3)
sim_input$data <- input

test_that("check_input_data drops NAs", {
  expect_length(nrow(check_input_data(sim_input, "sim")$data), 1)
  expect_length(nrow(check_input_data(input, "real")), 1)
})

test_that("check_input_data adds coverage", {
  expect_true(rlang::has_name(
    check_input_data(sim_input, "sim")$data,
    "coverage"
  ))

  expect_true(rlang::has_name(check_input_data(input, "real"), "coverage"))
})

test_that("check_input_data coerces to minor allele", {
  input <- tibble::tibble(
    plmaf = c(0.5, 0.7),
    wsmaf = c(0.3, 0.2),
    coverage = c(10, 10)
  )

  minor <- tibble::tibble(
    plmaf = c(0.5, 0.3),
    wsmaf = c(0.3, 0.8),
    coverage = c(10, 10)
  )

  expect_equal(check_input_data(input, "real"), minor)
})
