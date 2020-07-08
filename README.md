
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coiaf <a href='https://ojwatson.github.io/coiaf/'><img src='man/figures/logo_v2_light_green_no_background.png' align="right" height="139" /></a>

<!-- badges: start -->

[![R build
status](https://github.com/OJWatson/coiaf/workflows/R-CMD-check/badge.svg)](https://github.com/OJWatson/coiaf/actions)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

# Introduction

In malaria infections, individuals can often be infected by multiple
parasites due to repeated mosquito bites or mosquitoes harboring
multiple parasites. Such mixed infections may represent unrelated
strains or parasites that are related. The number of strains found in an
individual is known as the Complexity Of Infection (COI). Previous
methods have utilized statistics to determine the relative COI based on
the assumption that parasites are unrelated. While such assumptions
often do not hold, these methods provide insight into the force of
infection and the parasite population. Prior methods have been proposed
to measure COI, including probabilistic likelihood models and Bayesian
models. However, a rapid direct measure such as heterozygosity or FwS
has not yet been developed. Here we present two measures that show
equivalent power to other methods, but are directly calculable from
allele observed proportions based on the knowledge of overall population
allele frequencies.

# Installation

To install the package, please follow the code below. In order to
install, devtools must be installed.

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/OJWatson/coiaf")
```

# Usage

In order to run real data, please refer to the analysis folder. More
specifically, the script that is titled, `real_mccoil_data.Rmd` provides
a comprehensive set of instructions on how to run real data. For
questions regarding documentation, please refer to the [online package
website](ojwatson.github.io/coiaf/).

## Development

Please note that this package is still under development and may be
missing features. This `README` will be updated as the package is
developed accordingly.
