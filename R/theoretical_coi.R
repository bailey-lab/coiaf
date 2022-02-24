#------------------------------------------------
#' Theoretical COI
#'
#' Generate theoretical COI curves.
#'
#' \loadmathjax
#' The theoretical curve can be visualized as the WSMAF of an individual at
#' various PLMAFs. At each locus, the theoretical COI is defined as:
#' \mjsdeqn{1-p^{coi}-q^{coi}} where \mjseqn{p} is the PLMAF and \mjseqn{q} is
#' 1 - PLMAF.
#'
#' @param coi_range The COIs for which the curve will be calculated.
#' @param plmaf The PLMAF over which the curve will be calculated.
#' @param coi_method The method we will use to calculate the theoretical COI.
#'   The method is either "variant" or "frequency". The default value is
#'   "variant".
#'
#' @return The theoretical COI curves for the specified COIs and PLMAF.
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
  for (i in coi_range) {
    if (i == coi_range[1]) {
      curves <- data.frame(first = single_theoretical_coi(i, plmaf, coi_method))
      colnames(curves) <- paste0("coi_", i)
    } else {
      curves[[paste0("coi_", i)]] <- single_theoretical_coi(
        i,
        plmaf,
        coi_method
      )
    }
  }

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
