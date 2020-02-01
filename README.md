
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SSP development version

## Estimation of sampling effort in community ecology with SSP

Edlin Guerra-Castro, Juan Carlos Cajas, Juan Jose Cruz-Motta, Nuno
Simoes and Maite Mascaro

`SSP` is an R package design to estimate sampling effort in studies of
ecological communities based on the definition of pseudo multivariate
standard error (*MultSE*) (Anderson & Santana-Garcon 2015), simulation
of data and resampling.

`SSP` include seven functions: `assempar` for extrapolation of
assemblage parameters using pilot data; `simdata` for simulation of
several data sets based on extrapolated parameters; `datquality` for
evaluation of plausibility of simulated data; `sampsd` for repeated
estimations of *MultSE* for different sampling designs in simulated data
sets; `summary_sd` for summarizing the behavior of *MultSE* for each
sampling design across all simulated data sets, `ioptimum` for
identification of the optimal sampling effort, and `plot_ssp` to plot
sampling effort and *MultSE*.

## PACKAGE NEEDED IN SSP

  - Required: `vegan`, `sampling`, `stats` `ggplot2`
    [R](https://cran.r-project.org/). These are installed automatically.
  - Suggested: `devtools`, `httr` to build `SSP` from
    [github](https://github.com/edlinguerra/SSP). All these must be
    installed by you.

## HOW TO RUN SSP:

The `SSP` package will be available on
[CRAN](https://cran.r-project.org/) but can be downloaded from github
using the following commands:

``` r
## install the latest version from github
install.packages('devtools')
library(devtools)
install_github('edlinguerra/SSP', build_vignettes = TRUE)
```

For examples about how to use SSP, see the vignette in SSP after
instalation.
