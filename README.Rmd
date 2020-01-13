# Estimation of sampling effort in community ecology with SSP
Edlin Guerra-Castro, Juan Carlos Cajas, Juan Jose Cruz-Motta, Nuno Simoes and Maite Mascaro

`SSP` is an R package design to estimate sampling effort in studies of ecological communities based on the definition of pseudo multivariate standard error (*MultSE*) (Anderson & Santana-Garcon 2015) and simulation of data. This guide will provide you a brief overview in how to use `SSP`. The protocol in `SSP` consists in simulating several extensive data matrices that mimic some of the relevant ecological features of the community of interest using a pilot data set. For each simulated data, several sampling efforts are repeatedly executed and *MultSE* is calculated to each one. The mean value, 0.025 and 0.975 quantiles of *MultSE* for each sampling effort across all simulated data are then estimated, and optionally plotted using `ggplot2`. The mean values are then standardized regarding the lowest sampling effort (consequently, the worst precision), and an optimal sampling effort can be identified as that in which the increase in sample size do not improve the precision beyond a threshold value (e.g. 1 %).

`SSP` include six functions: `assempar` for extrapolation of assemblage parameters using pilot data; `simdata` for simulation of several data sets based on extrapolated parameters; `datquality` for evaluation of plausibility of simulated data; `sampsd` for repeated estimations of *MultSE* for different sampling designs in simulated data sets; `summary_sd` for summarizing the behavior of *MultSE* for each sampling design across all simulated data sets, and `ioptimum` for identification of the optimal sampling effort. The theoretical background are described in a submitted paper by Guerra-Castro et al. (2020).

## PACKAGE NEEDED IN SSP
- Required: `vegan`, `sampling`, `stats` [R](https://cran.r-project.org/)
- Suggested: `ggplot2` [R](https://cran.r-project.org/)
- Also `devtools`, `httr` to build `SSP` from [github](https://github.com/edlinguerra/SSP)

## HOW TO RUN SSP:
The `SSP` package will be available on [CRAN](https://cran.r-project.org/) but can be downloaded from github using the following commands:  

```{r eval=FALSE}
## install the latest version from github
install.packages('devtools')
library(devtools)
install_github('edlinguerra/SSP')
