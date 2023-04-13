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
to the COI. This package features two rapid, direct measures for
characterizing polyclonal infections.

## Installation

<div class=".pkgdown-release">

```r
# Install most recent released version
devtools::install_github("bailey-lab/coiaf@v0.1.2")
```

</div>

<div class=".pkgdown-devel">

```r
# Install development version
devtools::install_github("bailey-lab/coiaf")
```

</div>

## Usage

In order to run real data, please refer to the Articles drop-down menu.
Several articles detail how the algorithm works, how data is simulated
to test the algorithm, and, importantly, how to run real data. A short
example of running real data is included and outlines the necessary data
structure and the commands to run.
