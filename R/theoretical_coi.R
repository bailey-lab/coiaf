#------------------------------------------------
#' Theoretical COI
#'
#' Generate theoretical COI curves.
#'
#' \loadmathjax
#' The theoretical curve can be visualized as the WSAF of an individual at
#' various PLAFs. At each locus, the theoretical COI is defined as:
#' \mjsdeqn{1-p^{coi}-q^{coi}} where \mjseqn{p} is the PLAF and \mjseqn{q} is
#' 1 - PLAF.
#'
#' @param coi_range The COIs for which the curve will be calculated.
#' @param plaf The PLAF over which the curve will be calculated.
#' @param coi_method The method we will use to calculate the theoretical COI.
#' The method is either "1" or "2". The default value is "1".
#'
#' @return The theoretical COI curves for the specified COIs and PLAF.
#'
#' @export

theoretical_coi <- function(coi_range, plaf = seq(0, 0.5, l = 101),
                            coi_method = "1") {
  # Check inputs
  assert_pos(coi_range, zero_allowed = FALSE)
  assert_vector(plaf)
  assert_bounded(plaf, left = 0, right = 0.5)
  assert_increasing(plaf)
  assert_single_string(coi_method)
  assert_in(coi_method, c("1", "2"))

  # Compute curve for the COIs
  for (i in coi_range) {
    if (i == coi_range[1]) {
      curves <- data.frame(first = single_theoretical_coi(i, plaf, coi_method))
      colnames(curves) <- paste0("coi_", i)
    } else {
      curves[[paste0("coi_", i)]] = single_theoretical_coi(i, plaf, coi_method)
    }
  }

  # Include the PLAF in the final output and return
  curves$plaf <- plaf
  return(curves)
}

#------------------------------------------------
#' Single Theoretical COI
#'
#' Generate the theoretical COI curve for a particular COI value.
#'
#' \loadmathjax
#' The theoretical curve can be visualized as the WSAF of an individual at
#' various PLAFs. At each locus, the theoretical COI is defined as:
#' \mjsdeqn{1-p^{coi}-q^{coi}}
#'
#' @param coi The COI for which the curve will be calculated.
#' @inheritParams theoretical_coi
#'
#' @return The theoretical COI curve for the specified PLAF.
#'
#' @keywords internal

single_theoretical_coi <- function(coi,
                                   plaf = seq(0, 0.5, l = 101),
                                   coi_method = "1") {
  # Check inputs
  assert_single_pos(coi)
  assert_vector(plaf)
  assert_bounded(plaf, left = 0, right = 0.5)
  assert_increasing(plaf)
  assert_single_string(coi_method)
  assert_in(coi_method, c("1", "2"))

  # Determine the curve
  if (coi_method == "1") {
    curve <- 1 - plaf^coi - (1 - plaf)^coi
  } else if (coi_method == "2") {
    if (coi == 1) {
      curve <- rep(0, length(plaf))
    } else {
      curve <- (plaf - plaf^coi)/(1 - plaf^coi - (1 - plaf)^coi)
    }
  }
}
