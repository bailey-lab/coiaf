#------------------------------------------------
#' @title Generate Theoretical COI Curves
#'
#' @description Generate theoretical COI curves. At each locus defined as:
#' \eqn{1-p^coi-q^coi}. This expression is derived in the corresponding
#' manuscript.
#'
#' @param coi_range The COIs for which the curve will be calculated.
#' @param interval The interval over which the curve will be calculated.
#' @param method The method we will use to calculate the theoretical COI.
#' The method is either "1" or "2". The default value is "1".
#'
#' @return The theoretical COI curves for the specified COIs and interval.
#'
#' @export

theoretical_coi <- function(coi_range, interval = seq(0, 0.5, l = 101),
                            method = "1"){
  # Check inputs
  assert_pos_int(coi_range, zero_allowed = FALSE)
  assert_vector(interval)
  assert_bounded(interval, left = 0, right = 0.5)
  assert_increasing(interval)
  assert_single_string(method)
  assert_in(method, c("1", "2"))

  # Compute curve for the COIs
  for (i in coi_range){
    if (i == coi_range[1]){
      curves <- data.frame(first = single_theoretical_coi(i, interval, method))
      colnames(curves) <- paste0("coi_", i)
    } else {
      curves[[paste0("coi_", i)]] = single_theoretical_coi(i, interval, method)
    }
  }

  # Include the PLAF in the final output and return
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
#' @param coi The COI for which the curve will be calculated.
#' @inheritParams theoretical_coi
#'
#' @return The theoretical COI curve for the specified interval.
#'
#' @keywords internal

single_theoretical_coi <- function(coi, interval = seq(0, 0.5, l = 101),
                                   method = "1"){
  # Check inputs
  assert_single_pos_int(coi)
  assert_vector(interval)
  assert_bounded(interval, left = 0, right = 0.5)
  assert_increasing(interval)
  assert_single_string(method)
  assert_in(method, c("1", "2"))

  # Determine the curve
  if (method == "1"){
    curve <- 1 - interval^coi - (1 - interval)^coi
  } else if (method == "2"){
    curve <- (interval - interval^coi)/(1 - interval^coi - (1 - interval)^coi)
  }
  return(curve)
}
