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
#' @return A list of the following:
#' * `data`: A tibble with
#'  + `plaf_cut`: Breaks of the form `[a, b)`.
#'  + `m_variant`: The average WSAF or proportion of variant sites in each
#'  segment defined by `plaf_cut`.
#'  + `bucket_size`: The number of loci in each bucket.
#'  + `midpoints`: The midpoint of each bucket.
#' * `seq_error`: The sequence error inferred.
#' * `bin_size`: The minimum size of each bin.
#' * `cuts`: The breaks utilized in splitting the data.
#'  of each COI.
#' @keywords internal

process <- function(wsaf,
                    plaf,
                    seq_error = NULL,
                    bin_size = 20,
                    coi_method = "1") {

  # Infer value of seq_error if NULL
  if (is.null(seq_error)) {
    # Cut the data
    bins <- cut(plaf, seq(0, 0.5, 0.05), include.lowest = TRUE)

    # Data points in the lowest bin
    low_plafs <- wsaf[which(bins == levels(bins)[1])]

    # Number of data points with WSAF > 0
    error <- sum(low_plafs > 0, na.rm = T)

    # Expected number of points
    expected <- length(low_plafs) * 0.025

    # Compare expected number and actual number of points, but ensure that
    # seq_error is greater than 1%
    seq_error <- round(max((error - round(expected))/length(low_plafs), 0.01), 4)
  }

  if (coi_method == "1") {
    # Isolate PLAF, determine the PLAF cuts, and whether a site is a variant,
    # accounting for sequence error
    df <- data.frame(
      plaf_cut = Hmisc::cut2(plaf, m = bin_size),
      variant = ifelse(wsaf <= seq_error | wsaf >= (1 - seq_error), 0, 1))

  } else if (coi_method == "2") {
    # Subset to heterozygous sites
    data <- data.frame(wsaf = wsaf, plaf = plaf) %>%
      dplyr::filter(wsaf >= seq_error & wsaf <= (1 - seq_error))
    wsaf <- data$wsaf
    plaf <- data$plaf

    # If remove all data, need to return a pseudo result to not induce errors.
    # Additionally, in order to define a cut, need at least 2 data points
    if (length(plaf) <= 1) {
      vec <- setNames(rep("", 4), c("plaf_cut", "m_variant", "bucket_size", "midpoints"))
      df_grouped <- dplyr::bind_rows(vec)[0, ]
      res <- list(data = df_grouped,
                  seq_error = seq_error,
                  bin_size = bin_size,
                  cuts = NULL)
      return(res)
    }

    # Isolate PLAF, and keep WSAF as is
    df <- data.frame(plaf_cut = Hmisc::cut2(plaf, m = bin_size),
                     variant = wsaf)
  }

  # Average over intervals of PLAF
  df_grouped <- df %>%
    dplyr::group_by(.data$plaf_cut, .drop = FALSE) %>%
    dplyr::summarise(m_variant   = mean(.data$variant),
                     bucket_size = dplyr::n()) %>%
    stats::na.omit()

  cuts <- Hmisc::cut2(plaf, m = bin_size, onlycuts = TRUE)
  df_grouped$midpoints <- cuts[-length(cuts)] + diff(cuts) / 2

  # Return data, seq_error, and cuts
  res <- list(data = df_grouped,
              seq_error = seq_error,
              bin_size = bin_size,
              cuts = cuts)
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
#' @param seq_error The level of sequencing error that is assumed. If no value
#' is inputted, then we infer the level of sequence error.
#' @param bin_size The minimum size of each bin of data.
#' @inheritParams theoretical_coi
#'
#' @return A list of the following:
#' * `data`: A tibble with
#'  + `plaf_cut`: Breaks of the form `[a, b)`.
#'  + `m_variant`: The average WSAF or proportion of variant sites in each
#'  segment defined by `plaf_cut`.
#'  + `bucket_size`: The number of loci in each bucket.
#'  + `midpoints`: The midpoint of each bucket.
#' * `seq_error`: The sequence error inferred.
#' * `bin_size`: The minimum size of each bin.
#' * `cuts`: The breaks utilized in splitting the data.
#'  of each COI.
#' @family simulated data functions
#' @seealso [process_real()] to process real data.
#' @export

process_sim <- function(sim,
                        seq_error = NULL,
                        bin_size = 20,
                        coi_method = "1") {
  # Check inputs
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)

  # Run helper to process
  processed_sim <- process(wsaf       = sim$data$wsaf,
                           plaf       = sim$data$plaf,
                           seq_error  = seq_error,
                           bin_size   = bin_size,
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
#' @return A list of the following:
#' * `data`: A tibble with
#'  + `plaf_cut`: Breaks of the form `[a, b)`.
#'  + `m_variant`: The average WSAF or proportion of variant sites in each
#'  segment defined by `plaf_cut`.
#'  + `bucket_size`: The number of loci in each bucket.
#'  + `midpoints`: The midpoint of each bucket.
#' * `seq_error`: The sequence error inferred.
#' * `bin_size`: The minimum size of each bin.
#' * `cuts`: The breaks utilized in splitting the data.
#'  of each COI.
#' @family real data functions
#' @seealso [process_sim()] to process simulated data.
#' @export

process_real <- function(wsaf, plaf,
                         seq_error = NULL,
                         bin_size = 20,
                         coi_method = "1") {
  # Check inputs
  assert_vector(wsaf)
  assert_bounded(wsaf)
  assert_vector(plaf)
  assert_bounded(plaf)
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)

  # Ensure that the PLAF is at most 0.5
  plaf[plaf > 0.5] <- 1 - plaf[plaf > 0.5]
  assert_bounded(plaf, left = 0, right = 0.5)

  # Run helper to process
  processed_real <- process(wsaf       = wsaf,
                            plaf       = plaf,
                            seq_error  = seq_error,
                            bin_size   = bin_size,
                            coi_method = coi_method)
}
