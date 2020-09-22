#------------------------------------------------
#' Process data
#'
#' Helper function to process data
#'
#' The function computes whether a SNP is a variant site or not, based on the
#' simulated WSAF at that SNP. This process additionally accounts for potential
#' sequencing error.
#'
#' @inheritParams process_real
#'
#' @return Real COI curve.
#' @keywords internal

process <- function(wsaf,
                    plaf,
                    seq_error = 0.01,
                    cut = seq(0, 0.5, 0.01),
                    coi_method = "1") {

  if (coi_method == "1") {
    # Isolate PLAF, determine the PLAF cuts, and whether a site is a variant,
    # accounting for sequence error
    df <- data.frame(
      plaf_cut = cut(plaf, cut, include.lowest = TRUE),
      variant = ifelse(wsaf <= seq_error | wsaf >= (1 - seq_error), 0, 1)
      )

  } else if (coi_method == "2") {
    # Isolate PLAF, and keep WSAF as is
    df <- data.frame(
      plaf_cut = cut(plaf, cut, include.lowest = TRUE),
      variant = wsaf
    )
  }

  # Average over intervals of PLAF
  df_grouped <- df %>%
    dplyr::group_by(.data$plaf_cut, .drop = FALSE) %>%
    dplyr::summarise(m_variant   = mean(.data$variant),
                     bucket_size = dplyr::n()) %>%
    as.data.frame()

  # Include midpoints and remove missing data
  df_grouped$midpoints <- cut[-length(cut)] + diff(cut)/2
  df_grouped <- stats::na.omit(df_grouped)
}

#------------------------------------------------
#' Process simulated data
#'
#' Generate the simulated COI curve.
#'
#' Utilize the output of [sim_biallelic()], which creates simulated
#' data. The PLAF is kept, and the function computes whether a SNP is a
#' variant site or not, based on the simulated WSAF at that SNP. This process
#' additionally accounts for potential sequencing error. To check whether the
#' simulated WSAF correctly indicated a variant site or not, the phased
#' haplotype of the parasites is computed.
#'
#' @param sim Output of [sim_biallelic()].
#' @param seq_error The level of sequencing error that is assumed.
#' @param cut How often the data is summarized.
#' @inheritParams theoretical_coi
#'
#' @return Simulated COI curve.
#' @family simulated data functions
#' @seealso [process_real()] to process real data.
#' @export

process_sim <- function(sim,
                        seq_error = 0.01,
                        cut = seq(0, 0.5, 0.01),
                        coi_method = "1") {
  # Check inputs
  assert_single_bounded(seq_error)
  assert_bounded(cut, left = 0, right = 0.5)
  assert_vector(cut)
  assert_increasing(cut)

  # Subset data to focus only on heterozygous sites
  sim$data <- subset(sim$data, WSAF > 0 & WSAF < 1)

  # Run helper to process
  processed_sim <- process(wsaf       = sim$data$WSAF,
                           plaf       = sim$data$PLAF,
                           seq_error  = seq_error,
                           cut        = cut,
                           coi_method = coi_method)
}

#------------------------------------------------
#' Process real data
#'
#' Generate the COI curve for real data.
#'
#' The function computes whether a SNP is a variant site or not, based on the
#' simulated WSAF at that SNP. This process additionally accounts for potential
#' sequencing error.
#'
#' @param wsaf The within-sample allele frequency.
#' @param plaf The population-level allele frequency.
#' @inheritParams process_sim
#'
#' @return Real COI curve.
#' @family real data functions
#' @seealso [process_sim()] to process simulated data.
#' @export

process_real <- function(wsaf, plaf,
                         seq_error = 0.01,
                         cut = seq(0, 0.5, 0.01),
                         coi_method = "1") {
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

  # Run helper to process
  processed_real <- process(wsaf       = wsaf,
                            plaf       = plaf,
                            seq_error  = seq_error,
                            cut        = cut,
                            coi_method = coi_method)
}
