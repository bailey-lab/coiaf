#------------------------------------------------
#' Theoretical COI
#'
#' @description
#' \loadmathjax
#' Generate the theoretical relationship between the WSMAF (\mjseqn{\bf{w}}),
#' the PLMAF (\mjseqn{\bf{p}}), and the COI (\mjseqn{k}).
#'
#' @param coi_range The COIs for which the relationship will be generated.
#' @param plmaf The population-level minor allele frequency over which the
#'   relationship will be generated.
#' @param coi_method The method we will use to generate the theoretical
#'   relationship. The method is either "variant" or "frequency". The default
#'   value is "variant".
#'
#' @return
#' A [`tibble()`][tibble::tibble-package] containing the generated values. Each
#' column is named with the COI used. The last column of the tibble contains the
#' PLMAF.
#'
#' @export
#' @examples
#' theoretical_coi(1:5)
#' theoretical_coi(1:5, coi_method = "frequency")
theoretical_coi <- function(coi_range,
                            plmaf = seq(0, 0.5, length.out = 101),
                            coi_method = c("variant", "frequency")) {
  # Check inputs
  assert_pos(coi_range, zero_allowed = FALSE)
  assert_vector(plmaf)
  assert_bounded(plmaf, left = 0, right = 0.5)
  # assert_increasing(plmaf)

  # Argument match coi_method
  coi_method <- rlang::arg_match(coi_method)

  # Compute curve for the COIs
  named_coi_range <- rlang::set_names(coi_range, ~ paste0("coi_", coi_range))
  curves <- purrr::map_dfc(named_coi_range, function(x) {
    single_theoretical_coi(x, plmaf, coi_method)
  })

  # Include the PLMAF in the final output and return
  curves$plmaf <- plmaf
  curves
}

single_theoretical_coi <- function(coi, plmaf, coi_method) {
  switch(coi_method,
    variant = 1 - plmaf^coi - (1 - plmaf)^coi,
    frequency = (plmaf - plmaf^coi) / (1 - plmaf^coi - (1 - plmaf)^coi)
  )
}
