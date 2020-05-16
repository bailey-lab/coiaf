#------------------------------------------------
#' @title Get dataframe of P. falciparum chromosome lengths
#'
#' @description Get dataframe of P. falciparum chromosome lengths.
#'
#' @export

Pf_chrom_lengths <- function() {
  ret <- data.frame(chrom = 1:14,
                    length = c(643292, 947102, 1060087,
                               1204112, 1343552, 1418244,
                               1501717, 1419563, 1541723,
                               1687655, 2038337, 2271478,
                               2895605, 3291871))
  return(ret)
}
