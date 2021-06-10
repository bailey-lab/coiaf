
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coiaf <a href='https://ojwatson.github.io/coiaf/'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![R build
status](https://github.com/OJWatson/coiaf/workflows/R-CMD-check/badge.svg)](https://github.com/OJWatson/coiaf/actions)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Introduction

In malaria infections, individuals can often be infected by multiple
parasites due to repeated mosquito bites or mosquitoes harboring
multiple parasites. Such mixed infections may represent unrelated
strains or parasites that are related. The number of strains found in an
individual is known as the Complexity of Infection (COI). These methods
provide insight into the force of infection and the parasite population,
two measures that are becoming increasingly more important as they
provide indirect measures of the effectiveness of malaria control
efforts. Prior methods have been proposed to measure COI, including
probabilistic likelihood models and Bayesian models. However, a rapid
direct measure has not yet been developed that truly represents COI. The
goal of this project is to develop two rapid, direct measures that show
equivalent power to other well-established methods. In addition,
previous methods have utilized statistics to determine the relative COI
based on the assumption that parasites are unrelated, and the second aim
will be to incorporate relatedness to provide a more unbiased measure in
field situations.

## Installation

To install the package, please follow the code below. In order to
install, `devtools` must be installed.

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/OJWatson/coiaf")
```

## Usage

In order to run real data, please refer to the Articles drop down menu.
Several articles are provided which detail how the algorithm works, how
data was simulated to test the algorithm, and importantly how to run
real data. A short example on running real data is included and outlines
the necessary data structure as well as the commands to run.

### Development

Please note that this package is still under development and may be
missing features. This `README` will be updated as the package is
developed accordingly.

Features in development:

  - [ ] Create infrastructure for running real data given a particular
    data format
  - [ ] Account for read depth by weighing coverage
  - [ ] Create a shiny app

To dos:

  - [ ] Clean up package documentation and vignettes for release
  - [ ] Run the full SNP dataset
  - [ ] Statistics for the world map plots (ex: bootstrapping)
  - [ ] Add citation file for the package
  - [ ] Pf7k data will be released summer of 2021. Run coiaf on new data
