#------------------------------------------------
#' @title Generate Theoretical COI Curves
#'
#' @description Generate theoretical COI curves. At each locus defined as:
#' \eqn{1-p^{coi}-q^{coi}}. This expression is derived in the corresponding
#' manuscript.
#'
#' @param coi_range The COIs for which the curve will be calculated.
#' @param plaf The PLAF over which the curve will be calculated.
#' @param method The method we will use to calculate the theoretical COI.
#' The method is either "1" or "2". The default value is "1".
#'
#' @return The theoretical COI curves for the specified COIs and PLAF.
#'
#' @export

theoretical_coi <- function(coi_range, plaf = seq(0, 0.5, l = 101),
                            method = "1"){
  # Check inputs
  assert_pos_int(coi_range, zero_allowed = FALSE)
  assert_vector(plaf)
  assert_bounded(plaf, left = 0, right = 0.5)
  assert_increasing(plaf)
  assert_single_string(method)
  assert_in(method, c("1", "2"))

  # Compute curve for the COIs
  for (i in coi_range){
    if (i == coi_range[1]){
      curves <- data.frame(first = single_theoretical_coi(i, plaf, method))
      colnames(curves) <- paste0("coi_", i)
    } else {
      curves[[paste0("coi_", i)]] = single_theoretical_coi(i, plaf, method)
    }
  }

  # Include the PLAF in the final output and return
  curves$plaf <- plaf
  return(curves)
}

#------------------------------------------------
#' @title Generate Theoretical COI Curve
#'
#' @description Generate theoretical COI curve. At each locus defined as:
#' \eqn{1-p^coi-q^coi}. This expression is derived in the corresponding
#' manuscript.
#'
#' @param coi The COI for which the curve will be calculated.
#' @inheritParams theoretical_coi
#'
#' @return The theoretical COI curve for the specified PLAF.
#'
#' @keywords internal

single_theoretical_coi <- function(coi, plaf = seq(0, 0.5, l = 101),
                                   method = "1"){
  # Check inputs
  assert_single_pos_int(coi)
  assert_vector(plaf)
  assert_bounded(plaf, left = 0, right = 0.5)
  assert_increasing(plaf)
  assert_single_string(method)
  assert_in(method, c("1", "2"))

  # Determine the curve
  if (method == "1"){
    curve <- 1 - plaf^coi - (1 - plaf)^coi
  } else if (method == "2"){
    curve <- (plaf - plaf^coi)/(1 - plaf^coi - (1 - plaf)^coi)
  }
  return(curve)
}
