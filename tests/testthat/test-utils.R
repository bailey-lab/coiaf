test_that("quiet() silences output", {
  expect_silent(quiet(cat("a")))
  expect_silent(quiet(print("b")))
})

test_that("quiet() shows messages, warnings, and errors", {
  expect_message(quiet(message("a")))
  expect_warning(quiet(warning("b")))
  expect_error(quiet(error("c")))
})
