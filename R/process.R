#------------------------------------------------
#' Process data
#'
#' Helper function to process data
#'
#' The function computes whether a SNP is a variant site or not, based on the
#' simulated WSMAF at that SNP. This process additionally accounts for potential
#' sequencing error.
#'
#' @inheritParams process_real
#'
#' @return A list of the following:
#' * `data`: A tibble with
#'  + `plmaf_cut`: Breaks of the form `[a, b)`.
#'  + `m_variant`: The average WSMAF or proportion of variant sites in each
#'  segment defined by `plmaf_cut`.
#'  + `bucket_size`: The number of loci in each bucket.
#'  + `midpoints`: The midpoint of each bucket.
#' * `seq_error`: The sequence error inferred.
#' * `bin_size`: The minimum size of each bin.
#' * `cuts`: The breaks utilized in splitting the data.
#'  of each COI.
#' @keywords internal

process <- function(wsmaf,
                    plmaf,
                    coverage,
                    seq_error = NULL,
                    bin_size = 20,
                    coi_method = "variant") {

  # Infer value of seq_error if NULL
  if (is.null(seq_error)) {
    seq_error <- estimate_seq_error(wsmaf, plmaf, bin_size)
  }

  if (coi_method == "variant") {
    # Isolate PLMAF, determine the PLMAF cuts, and whether a site is a variant,
    # accounting for sequence error
    df <- data.frame(
      plmaf_cut = suppressWarnings(Hmisc::cut2(plmaf, m = bin_size)),
      variant = ifelse(wsmaf <= seq_error | wsmaf >= (1 - seq_error), 0, 1),
      coverage = coverage
    )
  } else if (coi_method == "frequency") {
    # Subset to heterozygous sites
    data <- data.frame(wsmaf = wsmaf, plmaf = plmaf, coverage = coverage) %>%
      dplyr::filter(wsmaf > seq_error & wsmaf < (1 - seq_error))
    wsmaf <- data$wsmaf
    plmaf <- data$plmaf
    coverage <- data$coverage

    # If remove all data, need to return a pseudo result to not induce errors.
    # Additionally, in order to define a cut, need at least 2 data points
    if (length(plmaf) <= 1) {
      vec <- stats::setNames(
        rep("", 4),
        c("plmaf_cut", "m_variant", "bucket_size", "midpoints")
      )
      df_grouped <- dplyr::bind_rows(vec)[0, ]
      res <- list(
        data = df_grouped,
        seq_error = seq_error,
        bin_size = bin_size,
        cuts = NULL
      )
      return(res)
    }

    # Isolate PLMAF, and keep WSMAF as is
    df <- data.frame(
      plmaf_cut = suppressWarnings(Hmisc::cut2(plmaf, m = bin_size)),
      variant = wsmaf,
      coverage = coverage
    )
  }

  # In some instances, Hmisc::cut2 assigns a cut with only 1 number in it.
  # If this happens and we try to group our data, this can mess up our data.
  # Therefore, to account for this, we find all instances where this occurs
  # and combine these factors with the previous factor.
  one_point <- !stringr::str_starts(levels(df$plmaf_cut), "\\[")

  # We find all places where we only have one point. But, we ignore the case
  # where the one point is the first break (0).
  if (sum(one_point) > 1 | (sum(one_point) == 1 & which(one_point)[1] != 1)) {
    if (which(one_point)[1] == 1) {
      # When 0 is its own break, we ignore it and store all the other locations
      points <- which(one_point)[-1]
    } else {
      # When 0 is not its own break, we store all locations
      points <- which(one_point)
    }

    # We make a list of the factor names of all the points we want to remove. We
    # name the list with the previous factor, and combine the two together. This
    # effectively puts the points in the single factor into the previous one.
    point_list <- c(levels(df$plmaf_cut)[points])
    names(point_list) <- levels(df$plmaf_cut)[points - 1]
    df$plmaf_cut <- forcats::fct_recode(df$plmaf_cut, !!!point_list)
  }

  # Average over intervals of PLMAF
  df_grouped <- df %>%
    dplyr::group_by(.data$plmaf_cut, .drop = FALSE) %>%
    dplyr::summarise(
      m_variant   = stats::weighted.mean(.data$variant, .data$coverage),
      bucket_size = dplyr::n()
    ) %>%
    stats::na.omit()

  # Find the cuts for our data
  cuts <- suppressWarnings(Hmisc::cut2(plmaf, m = bin_size, onlycuts = TRUE))
  if (sum(one_point) > 1 | (sum(one_point) == 1 & which(one_point)[1] != 1)) {
    cuts <- cuts[-points]
  }

  # We then find our midpoints
  df_grouped$midpoints <- cuts[-length(cuts)] + diff(cuts) / 2

  # Return data, seq_error, and cuts
  list(
    data = df_grouped,
    seq_error = seq_error,
    bin_size = bin_size,
    cuts = cuts
  )
}

#------------------------------------------------
#' Process simulated data
#'
#' Generate the simulated COI curve.
#'
#' Utilize the output of [sim_biallelic()], which creates simulated
#' data. The PLMAF is kept, and the function computes whether a SNP is a
#' variant site or not, based on the simulated WSMAF at that SNP. This process
#' additionally accounts for potential sequencing error. To check whether the
#' simulated WSMAF correctly indicated a variant site or not, the phased
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
#'  + `plmaf_cut`: Breaks of the form `[a, b)`.
#'  + `m_variant`: The average WSMAF or proportion of variant sites in each
#'  segment defined by `plmaf_cut`.
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
                        coi_method = "variant") {
  # Check inputs
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)

  # Run helper to process
  process(
    wsmaf = sim$data$wsmaf,
    plmaf = sim$data$plmaf,
    coverage = sim$data$coverage,
    seq_error = seq_error,
    bin_size = bin_size,
    coi_method = coi_method
  )
}

#------------------------------------------------
#' Process real data
#'
#' Generate the COI curve for real data.
#'
#' The function computes whether a SNP is a variant site or not, based on the
#' simulated WSMAF at that SNP. This process additionally accounts for potential
#' sequencing error.
#'
#' @param wsmaf The within-sample allele frequency.
#' @param plmaf The population-level allele frequency.
#' @param coverage The read coverage at each locus.
#' @inheritParams process_sim
#'
#' @return A list of the following:
#' * `data`: A tibble with
#'  + `plmaf_cut`: Breaks of the form `[a, b)`.
#'  + `m_variant`: The average WSMAF or proportion of variant sites in each
#'  segment defined by `plmaf_cut`.
#'  + `bucket_size`: The number of loci in each bucket.
#'  + `midpoints`: The midpoint of each bucket.
#' * `seq_error`: The sequence error inferred.
#' * `bin_size`: The minimum size of each bin.
#' * `cuts`: The breaks utilized in splitting the data.
#'  of each COI.
#' @family real data functions
#' @seealso [process_sim()] to process simulated data.
#' @export

process_real <- function(wsmaf,
                         plmaf,
                         coverage,
                         seq_error = NULL,
                         bin_size = 20,
                         coi_method = "variant") {
  # Check inputs
  assert_vector(wsmaf)
  assert_bounded(wsmaf)
  assert_vector(plmaf)
  assert_bounded(plmaf)
  if (!is.null(seq_error)) assert_single_bounded(seq_error)
  assert_single_pos_int(bin_size)

  input <- tibble::tibble(wsmaf = wsmaf, plmaf = plmaf, coverage = coverage) %>%
    tidyr::drop_na()

  # In some cases we are fed in the major allele so we ensure we only examine
  # the minor allele. If the PLAF is > 0.5, we know it is the major allele so
  # we look at 1 - WSAF and 1 - PLAF.
  minor <- input %>%
    dplyr::mutate(
      wsmaf = ifelse(plmaf > 0.5, 1 - wsmaf, wsmaf),
      plmaf = ifelse(plmaf > 0.5, 1 - plmaf, plmaf)
    )

  assert_bounded(minor$plmaf, left = 0, right = 0.5)

  # Run helper to process
  process(
    wsmaf = minor$wsmaf,
    plmaf = minor$plmaf,
    coverage = minor$coverage,
    seq_error = seq_error,
    bin_size = bin_size,
    coi_method = coi_method
  )
}

#' @noRd
check_real_data <- function(wsmaf, plmaf) {

  # Check inputs
  assert_vector(wsmaf)
  assert_bounded(wsmaf)
  assert_vector(plmaf)
  assert_bounded(plmaf)

  input <- tibble::tibble(wsmaf = wsmaf, plmaf = plmaf) %>%
    tidyr::drop_na()

  # In some cases we are fed in the major allele so we ensure we only examine
  # the minor allele. If the PLAF is > 0.5, we know it is the major allele so
  # we look at 1 - WSAF and 1 - PLAF.
  minor <- input %>%
    dplyr::mutate(
      wsmaf = ifelse(plmaf > 0.5, 1 - wsmaf, wsmaf),
      plmaf = ifelse(plmaf > 0.5, 1 - plmaf, plmaf)
    )

  return(minor)

}

#' @noRd
estimate_seq_error <- function(wsmaf, plmaf, bin_size) {

  # Cut the data. We define error break to allow for flexible bucket sizes.
  error_break <- 0.02
  bins <- cut(plmaf, seq(0, 0.5, error_break), include.lowest = TRUE)

  # We want to ensure that we have at least bin_size points in the first bin
  while (table(bins)[1] < bin_size) {
    error_break <- error_break * 1.25
    bins <- cut(plmaf, seq(0, 0.5, error_break), include.lowest = TRUE)
  }

  # Data points in the lowest bin that are likely sequence error
  low_wsmafs <- wsmaf[which(bins == levels(bins)[1])]
  error <- low_wsmafs[low_wsmafs > 0 & low_wsmafs < 0.5]

  # If wanted to do mixture models, would fit something to error

  # Expected number of points
  expected <- round(length(low_wsmafs) * (error_break / 2), 4)

  # Remove expected number of points from true points
  error_dist <- utils::head(sort(error), -expected)

  # Find 95% error
  seq_error <- as.numeric(stats::quantile(error_dist, 0.95, na.rm = T))

  # Ensure that seq_error is greater than 1%
  seq_error <- round(max(seq_error, 0.01, na.rm = T), 4)

  return(seq_error)

}
