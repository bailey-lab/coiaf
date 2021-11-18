<!-- README.md is generated from README.Rmd. Please edit that file -->

# coiaf <a href='https://bailey-lab.github.io/coiaf/'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/bailey-lab/coiaf/workflows/R-CMD-check/badge.svg)](https://github.com/bailey-lab/coiaf/actions)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- badges: end -->

## Introduction

In malaria, individuals are often infected with different parasite
strains, with the complexity of infection (COI) giving the number of
genetically different parasite strains in an individual. Changes in the
mean COI in a population have been shown to be informative of changes in
transmission intensity with a number of probabilistic likelihood and
Bayesian models now developed to estimate COI. However, rapid, direct
measures based on heterozygosity or _FwS_ have not been directly related
to the COI. In this package we present two new methods that use easily
calculated measures to directly estimate the COI from sequence read
depth data.

## Installation

You may install the package from
[Github](https://github.com/bailey-lab/coiaf) using `devtools`. You may
either install the most recent release or the development version of the
software.

```r
# install.packages("devtools")

# Most recent release
devtools::install_github("bailey-lab/coiaf@v0.1.0")

# Development version
devtools::install_github("bailey-lab/coiaf")
```

## Usage

Please consult the [example real data
vignette](https://github.com/bailey-lab/coiaf/blob/main/analysis/vignettes/example_real_data.Rmd)
for a primer on how to run real data using `coiaf`. For a description of
how the estimation methods work, please consult the [example coi
prediction
vignette](https://github.com/bailey-lab/coiaf/blob/main/analysis/vignettes/example_coi_prediction.Rmd).

<!-- In order to run real data, please refer to the Articles drop down menu. Several -->
<!-- articles are provided which detail how the algorithm works, how data was  -->
<!-- simulated to test the algorithm, and importantly how to run real data. -->
<!-- A short example on running real data is included and outlines the necessary -->
<!-- data structure as well as the commands to run. -->

### Development

Please note that this package is still under development and may be
missing features.
