#------------------------------------------------
#' @title Generate Theoretical COI Curves
#'
#' @description Generate theoretical COI curves. At each locus defined as:
#' \eqn{1-p^coi-q^coi}. This expression is derived in the corresponding
#' manuscript.
#'
#' @param cois The COIs for which the curve will be calculated
#' @param interval The interval over which the curve will be calculated
#'
#' @return The theoretical COI curves for the specified COIs and interval
#'
#' @export

theoretical_coi <- function(cois, interval){
  # Check inputs
  assert_pos_int(cois, zero_allowed = FALSE)
  assert_bounded(interval, left = 0, right = 0.5)
  assert_increasing(interval)

  # Compute curve for the COIs
  for (i in cois){
    if (i == cois[1]){
      curves <- data.frame(first = single_theoretical_coi(i, interval))
      colnames(curves) <- paste("coi_", i, sep="")
    } else {
      curves[[paste("coi_", i, sep="")]] = single_theoretical_coi(i, interval)
    }
  }

  curves$PLAF <- interval
  return(curves)
}

#------------------------------------------------
#' @title Generate Theoretical COI Curve
#'
#' @description Generate theoretical COI curve. At each locus defined as:
#' \eqn{1-p^coi-q^coi}. This expression is derived in the corresponding
#' manuscript.
#'
#' @param coi The COI
#' @param interval The interval over which the curve will be calculated
#'
#' @return The theoretical COI curve for the specified interval
#'
#' @keywords internal

single_theoretical_coi <- function(coi, interval){
  # Check inputs
  assert_single_pos_int(coi)
  assert_bounded(interval, left = 0, right = 0.5)
  assert_increasing(interval)

  # Determine the curve
  curve <- 1 - interval^coi - (1 - interval)^coi
  return(curve)
}
