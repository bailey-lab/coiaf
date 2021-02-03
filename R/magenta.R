#------------------------------------------------
#' Simulate biallelic data from a magenta output
#'
#' Simulate biallelic data using a \code{magenta} model output.
#'
#' @param coverage Coverage at each locus. If a single value then the same
#'   coverage is applied over all loci.
#' @param alpha Shape parameter of the symmetric Dirichlet prior on strain
#'   proportions.
#' @param overdispersion The extent to which counts are over-dispersed relative
#'   to the binomial distribution. Counts are Beta-binomially distributed, with
#'   the beta distribution having shape parameters
#'   \mjseqn{\frac{p}{overdispersion}} and
#'   \mjseqn{\frac{1-p}{overdispersion}}.
#' @param epsilon The probability of a single read being miscalled as the other
#'   allele. Applies in both directions.
#'
#' @return A list of:
#' * `coi`: The COI used to simulate the data.
#' * `strain_proportions`: The strain proportion of each strain.
#' * `phased`: The phased haplotype.
#' * `data`: A dataframe of:
#'   + `plaf`: The population-level allele frequency.
#'   + `coverage`: The coverage at each locus.
#'   + `counts`: The count at each locus.
#'   + `wsaf`: The within-sample allele frequency.
#' * `inputs`: A dataframe of function input arguments:
#'   + `alpha`: Shape parameters of Dirichlet controlling strain proportions.
#'   + `overdispersion`: Overdispersion in count data.
#'   + `relatedness`: Within sample relatedness between strains.
#'   + `epsilon`: Probability of a single read being miscalled.
#' * `relatedness`: Within sample relatedness
#'
#' @family simulated data functions
#'
#' @example
#' \dontrun{
#'
#'  devtools::install_github("OJWatson/magenta")
#'  library(magenta)
#'
#'  out <- magenta::pipeline(
#'  EIR=1,N=1000,years=2,num_loci = 1000,full_save=TRUE,
#'  drug_list = drug_list_create(drugs = list(magenta:::drug_create_dhappq()))
#'  )
#'
#'  sims <- magenta_sim_format(out)
#'
#'  sims_with_data <- sim_biallelic_magenta(sims)
#'
#' }
#'
#' @export
sim_biallelic_magenta <- function(out,
                                  coverage = 200,
                                  alpha = 1,
                                  overdispersion = 0,
                                  epsilon = 0) {

  # If a single value was input, then repeat coverage so that coverage is
  # applied over all loci.
  if (length(coverage) == 1) {
    coverage <- rep(coverage, nrow(out$data))
  }

  # Continue to check inputs
  assert_vector(coverage)
  assert_pos_int(coverage)
  assert_single_pos(alpha, zero_allowed = FALSE)
  assert_single_pos(overdispersion, zero_allowed = TRUE)
  assert_single_pos(epsilon, zero_allowed = TRUE)
  assert_bounded(epsilon)

  # Get strain proportions
  w <- out$strain_proportions

  # Get coi
  coi <- out$coi

  # Get plaf
  plaf <- out$data$plaf

  # Get num_loci
  L <- length(plaf)

  # Generate true WSAF levels by summing binomial draws over strain proportions
  m <- mapply(function(x) rbinom(coi, 1, x), x = plaf)

  # Draw with the within sample allele frequencies (p_levels) are
  if (coi == 1) {
    p_levels = m*w
  } else{
    p_levels <- colSums(sweep(m, 1, w, "*"))
  }

  # Rounding errors from multiplying w by m can cause numbers greater than 1
  p_levels[p_levels > 1] <- 1L

  # Add in genotyping error
  p_error <- p_levels * (1 - epsilon) + (1 - p_levels) *epsilon

  # Draw read counts, taking into account overdispersion
  if (overdispersion == 0) {
    counts <- rbinom(L, size = coverage, prob = p_error)
  } else {
    counts <- rbetabinom(L,
                         k = coverage,
                         alpha = p_error/overdispersion,
                         beta = (1 - p_error)/overdispersion)
  }

  # Return list
  ret <- list(coi = coi,
              strain_proportions = w,
              phased = m,
              data = data.frame(plaf     = plaf,
                                coverage = coverage,
                                counts   = counts,
                                wsaf     = counts/coverage),
              inputs = data.frame(alpha          = alpha,
                                  overdispersion = overdispersion,
                                  relatedness    = relatedness,
                                  epsilon        = epsilon),
              relatedness = out$relatedness)

  return(ret)
}

#' @noRd
magenta_sim_format <- function(r) {

  # phased genotypes
  samples <- lapply(r$populations_event_and_strains_List$Strain_barcode_vectors, lapply, as.integer)
  phased_list <- lapply(samples, function(x) { do.call(rbind, x) })

  # cois
  coi_list <- lengths(phased_list)/r$parameters_List$g_num_loci

  # strain proportions
  birth <- r$populations_event_and_strains_List$Strain_day_of_acquisition_vectors
  age <- lapply(birth, function(x) {r$parameters_List$g_current_time - x})
  strain_proportions_list <- lapply(age, function(x) {
    if (length(x)>0) {
      ((x/max(x))/max(x)*1.1)/sum((x/max(x))/max(x)*1.1)
    } else {
      integer(0)
    }
  })

  # subset to those that are parasites
  pos <- which(coi_list > 0)
  phased_list <- phased_list[pos]
  coi_list <- coi_list[pos]
  strain_proportions_list <- strain_proportions_list[pos]
  coiaf_rets <- vector("list", length(pos))

  for(i in seq_along(coi_list)) {

    if(coi_list[i] > 1) {
      relatedness <- mean(dist(phased_list[[i]], method = "manhattan") / r$parameters_List$g_num_loci)
    } else {
      relatedness <- NA
    }

    # Return list
    coiaf_rets[[i]] <- list(coi = coi_list[i],
                            strain_proportions = strain_proportions_list[[i]],
                            phased = phased_list[[i]],
                            data = data.frame(plaf= r$parameters_List$g_plaf),
                            relatedness = relatedness)

  }

  return(coiaf_rets)
}
