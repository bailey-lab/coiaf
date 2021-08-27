
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coiaf <a href='https://bailey-lab.github.io/coiaf/'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/bailey-lab/coiaf/workflows/R-CMD-check/badge.svg)](https://github.com/bailey-lab/coiaf/actions)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

## Introduction

In malaria, individuals are often infected with different parasite
strains. These multi-clonal infections occur either because an
individual received repeated mosquito bites from infectious mosquitoes
or because they received a single bite from an infectious mosquito
harboring multiple parasites. Such mixed infections may be composed of
genetically unrelated strains or parasites that are related. The number
of genetically different parasite strains in an individual, known as the
complexity of infection (COI), provides important insight into the
severity of infection and the parasite population. Probabilistic
likelihood and Bayesian models have been developed to estimate the COI,
but rapid direct measures such as heterozygosity or *FwS* have not been
directly related to the COI. Moreover, current methods have utilized
statistics that rely on the assumption that parasites in mixed
infections are unrelated to determine the relative COI. In this package
we present two new methods that use easily computable measures to
directly estimate the COI from sequence read depth data. Our methods are
computationally efficient and are proven to be comparably accurate to
current methods in literature.

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
  - [ ] Add progress bars using `cli`
  - [ ] Run the full SNP dataset
  - [ ] Statistics for the world map plots (ex: bootstrapping)
  - [ ] Add citation file for the package
  - [ ] Pf7k data will be released summer of 2021. Run coiaf on new data
