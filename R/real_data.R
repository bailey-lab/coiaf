#------------------------------------------------
#' Run real data
#'
#' Run the algorithm on real data.
#'
#' @param data The dataset.
#' @param max_coi The maximum COI that the model will look at. Looks at the
#' theoretical COIs from `1` till `max_coi`.
#' @inheritParams sensitivity
#'
#' @return A list of samples. Each sample contains:
#' * The predicted COI.
#' * The probability distribution for the predictions.
#' @family real data functions
#' @export

run_real <- function(data,
                     max_coi = 25,
                     seq_error = 0.01,
                     cut = seq(0, 0.5, 0.01),
                     comparison = "overall",
                     distance = "squared",
                     weighted = TRUE) {

  # Check inputs
  assert_single_pos_int(max_coi)
  assert_single_bounded(seq_error)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)
  assert_single_string(comparison)
  assert_in(comparison, c("end", "ideal", "overall"))
  assert_single_string(distance)
  assert_in(distance, c("abs_sum", "sum_abs", "squared"))
  assert_logical(weighted)

  # Function to determine if pbapply is installed. If it is installed, it will
  # display a progress bar
  list_apply <- function(x, fun, ...) {
    if (requireNamespace("pbapply", quietly = TRUE)) {
      pbapply::pblapply(x, fun, ...)
    } else {
      lapply(x, fun, ...)
    }
  }

  # PLAF is the same for the entire area so we can compute it once outside the
  # lapply
  plaf <- colMeans(data, na.rm = T)

  # Run each sample in the data input
  coi_pred <- lapply(seq_len(nrow(data)), function(x) {
    # Get wsaf and remove any missing data
    wsaf  <- data[x,]
    input <- data.frame(wsaf = wsaf, plaf = plaf) %>% tidyr::drop_na()

    # Format data in the proper way
    processed_data <- process_real(input$wsaf, input$plaf, seq_error, cut)

    # Compute the coi
    sample_coi <- compute_coi(processed_data, 1:max_coi, cut,
                              comparison, distance, weighted)
  })

  # Add names to the output
  names(coi_pred) <- rownames(data)

  # Return predicted COIs and probabilities for each sample
  return(coi_pred)
}

