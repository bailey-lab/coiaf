#------------------------------------------------
#' @title Generate Real COI Curve
#'
#' @description Generate the COI curve for real data.
#'
#' @param wsaf WSAF.
#' @param plaf PLAF.
#' @param seq_error The sequencing error.
#' @param cut How often the data is summarized.
#'
#' @return Real COI curve.
#'
#' @export

process_real_data <- function(wsaf, plaf, seq_error = 0.01,
                              cut = seq(0, 0.5, 0.01)){
  # Check inputs
  assert_vector(wsaf)
  assert_bounded(wsaf)
  assert_vector(plaf)
  assert_bounded(plaf)
  assert_single_bounded(seq_error)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)

  # Ensure that the PLAF is at most 0.5
  plaf[plaf > 0.5] <- 1 - plaf[plaf > 0.5]
  assert_bounded(plaf, left = 0, right = 0.5)

  # Isolate PLAF, determine the PLAF cuts, and whether a site is a variant
  df <- data.frame(
    plaf_cut = cut(plaf, cut, include.lowest = TRUE),
    variant  = ifelse(wsaf < seq_error | wsaf > (1 - seq_error), 0, 1))

  # Average over intervals of PLAF
  df_grouped <- df %>%
    dplyr::group_by(.data$plaf_cut) %>%
    dplyr::summarise(m_variant   = mean(.data$variant),
                     bucket_size = dplyr::n()) %>%
    as.data.frame()

  # Include midpoints
  df_grouped$midpoints <- cut[-length(cut)] + diff(cut)/2

  return(df_grouped)
}


#------------------------------------------------
#' @title Run Real Data COIs
#'
#' @description Run the algorithm on real data
#'
#' @param data The dataset.
#' @inheritParams coi_test
#'
#'
#' @return A list of the following dataframes:
#' \describe{
#'   \item{\code{predicted_coi}}{A dataframe of the predicted COIs. COIs are
#'   predicted using \link{compute_coi}. Each column represents a separate set
#'   of parameters. Each row represents a predicted COI. Predictions are done
#'   many times, depending on the value of \code{repetitions}.}
#'   \item{\code{param_grid}}{The parameter grid. The parameter grid is all
#'   possible combinations of the parameters inputted. Each row represents a
#'   unique combination.}
#'   \item{\code{error_bias}}{A dataframe containing any parameter that was
#'   varied and the associated mean absolute error and bias (mean error). By
#'   showing only parameters that were varied, the output is easier to interpret
#'   and does not have information about parameters that were held constant.}
#' }
#'
#' @export

run_real_data <- function(data,
                          seq_error = 0.01,
                          cut = seq(0, 0.5, 0.01),
                          max_COI = 25,
                          method = "overall",
                          dist_method = "squared",
                          weighted = TRUE){

  # Check inputs
  assert_single_pos_int(max_COI)
  assert_single_bounded(seq_error)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)
  assert_single_string(method)
  assert_in(method, c("end", "ideal", "overall"))
  assert_single_string(dist_method)
  assert_in(dist_method, c("abs_sum", "sum_abs", "squared", "KL"))
  assert_logical(weighted)

  # Function to determine if pbapply is installed. If it is installed, it will
  # display a progress bar
  list_apply <- function(x, fun, ...){
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
  coi_pred <- list_apply(seq_len(nrow(data)), function(x) {
    # Get wsaf and remove any missing data
    wsaf  <- data[x,]
    input <- data.frame(wsaf = wsaf, plaf = plaf) %>% tidyr::drop_na

    # Format data in the proper way
    processed_data <- process_real_data(input$wsaf, input$plaf, seq_error, cut)

    # Compute the coi
    sample_coi <- compute_coi(processed_data, 1:max_COI, cut,
                              method, dist_method, weighted)

    return(sample_coi)
  })

  # Add names to the output
  names(coi_pred) <- rownames(data)

  # Return predicted COIs and probabilities for each sample
  return (coi_pred)
}

