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

  # Average over intervals of PLMAF
  df_grouped <- df %>%
    dplyr::group_by(.data$plmaf_cut, .drop = FALSE) %>%
    dplyr::summarise(
      m_variant   = stats::weighted.mean(.data$variant, .data$coverage),
      bucket_size = dplyr::n()
    ) %>%
    stats::na.omit()

  # Compute midpoints and set coverage to be uniform across each bucket
  df_grouped_mid <- find_cut_midpoints(df_grouped, .data$plmaf_cut) %>%
    tibble::add_column(coverage = rep(100, nrow(.)))

  # Return data, seq_error, and cuts
  list(
    data = df_grouped_mid,
    seq_error = seq_error,
    bin_size = bin_size,
    cuts = suppressWarnings(Hmisc::cut2(plmaf, m = bin_size, onlycuts = TRUE))
  )
}

find_cut_midpoints <- function(data, cuts) {
  # Convert single cuts to the standard format: "[lower,upper)"
  fix_single_cuts <- dplyr::mutate(
    data,
    fixed_cuts = as.character({{ cuts }}),
    fixed_cuts = ifelse(
      !stringr::str_starts(.data$fixed_cuts, "\\["),
      glue::glue("[{.data$fixed_cuts},{.data$fixed_cuts})"),
      .data$fixed_cuts
    )
  )

  # Find lower and upper bounds
  extract_bounds <- tidyr::extract(
    data = fix_single_cuts,
    col = .data$fixed_cuts,
    into = c("lower", "upper"),
    regex = "([[:alnum:]].+),([[:alnum:]].+)[\\]\\)]",
    convert = TRUE
  )

  # Determine cut midpoints
  dplyr::mutate(
    extract_bounds,
    midpoints = (.data$upper + .data$lower) / 2,
    .keep = "unused"
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

check_input_data <- function(data, data_type) {
  if (data_type == "sim") {
    # Remove NA from our data frame
    data$data <- tidyr::drop_na(data$data)

    # Add coverage if it is missing
    if (!"coverage" %in% names(data$data)) {
      data$data$coverage <- rep(100, length(data$data$plmaf))
    }
  } else if (data_type == "real") {
    # Remove NA from our data frame
    data <- tidyr::drop_na(data)

    # Add coverage if it is missing
    if (!"coverage" %in% names(data)) {
      data$coverage <- rep(100, length(data$plmaf))
    }

    # Ensure we are calling the minor allele
    data <- check_minor_allele(data)
  }

  data
}

check_minor_allele <- function(data) {
  # In some cases we are fed in the major allele so we ensure we only examine
  # the minor allele. If the PLAF is > 0.5, we know it is the major allele so
  # we look at 1 - WSAF and 1 - PLAF.
  dplyr::mutate(
    data,
    wsmaf = ifelse(plmaf > 0.5, 1 - wsmaf, wsmaf),
    plmaf = ifelse(plmaf > 0.5, 1 - plmaf, plmaf)
  )
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
