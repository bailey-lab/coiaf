
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coiaf

<!-- badges: start -->

[![R build
status](https://github.com/r-lib/usethis/workflows/R-CMD-check/badge.svg)](https://github.com/r-lib/usethis/actions)
<!-- badges: end -->

## Introduction

In malaria infections, individuals can often be infected by multiple
parasites due to repeated mosquito bites or mosquitoes harboring
multiple parasites. Such mixed infections may represent unrelated
strains or parasites that are related. The number of strains found in an
individual is known as the Complexity Of Infection (COI). Previous
methods have utilized statistics to determine the relative COI based on
the assumption that parasites are unrelated. While such assumptions
often do not hold, these methods provide insight into the force of
infection and the parasite population. While multiple methods have been
proposed to measure COI, there have been no direct statistics allowing
for direct calculation of the COI. Probabilistic likelihood and Bayesian
models have been developed but a rapid direct measure such as
heterozygosity or FwS have not been directly tied to the COI. Here we
present two measures that show equivalent power to other methods, but
are directly calculable from allele observed proportions within
infection based on the knowledge of overall population allele
frequencies.

## Installation

To install the package, please follow the code below. In order to
install, devtools must be installed and loaded.

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/OJWatson/coiaf")
```

## Usage

In order to run the two measures, please refer to the analysis folder.
More specifically, the script that is titled:
`02_simple_simulation_comparison_demonstration.Rmd`.

## Development

Please note that this package is still under development and therefore
there may be missing features. This `README` will be updated as the
package is developed accordingly.
