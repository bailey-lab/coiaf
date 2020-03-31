#------------------------------------------------
#' @title Generate Theoretical COI Curve
#'
#' @description Generate theoretical COI curve. At each locus defined as:
#' \eqn{1-p^coi-q^coi}. This expression is derived in the corresponding
#' manuscript.
#'
#' @param coi The COI
#' @param interval The interval over which the curve will be calculated
#' @return The theoretical COI curve for the specified interval.
#'
#' @export

theoretical_coi <- function(coi, interval){
  curve <- 1 - interval^coi - (1 - interval)^coi
  return(curve)
}
