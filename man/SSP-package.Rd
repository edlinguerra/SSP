\name{SSP-package}
\alias{SSP}

\title{
SSP: Simulated Sampling Procedure for Community Ecology
}
\description{
SSP is an R package designed to estimate sampling effort in studies of ecological communities based on the definition of pseudo multivariate standard error (MultSE) (Anderson & Santana-Garcon 2015) and simulation of data (Guerra-Castro et al., 2020).
}
\details{
The protocol in SSP consists in simulating several extensive data matrices that mimic some of the relevant ecological features of the community of interest using a pilot data set. For each simulated data, several sampling efforts are repeatedly executed and MultSE is calculated to each one. The mean value, 0.025 and 0.975 quantiles of MultSE for each sampling effort across all simulated data are then estimated and plotted. The mean values are standardized in relation to the lowest sampling effort (consequently, the worst precision), and an optimal sampling effort can be identified as that in which the increase in sample size do not improve the precision beyond a threshold value (e.g. 2.5\%).

SSP includes seven functions: \code{\link{assempar}} for extrapolation of assemblage parameters using pilot data; \code{\link{simdata}} for simulation of several data sets based on extrapolated parameters; \code{\link{datquality}} for evaluation of plausibility of simulated data; \code{\link{sampsd}} for repeated estimations of MultSE for different sampling designs in simulated data sets; \code{\link{summary_ssp}} for summarizing the behavior of MultSE for each sampling design across all simulated data sets, \code{\link{ioptimum}} for identification of the optimal sampling effort, and  \code{\link{plot_ssp}} to plot sampling effort vs MultSE of simulated data.

The SSP package is developed at GitHub (\url{https://github.com/edlinguerra/SSP/}).
}

\references{
-Anderson, M. J., & J. Santana-Garcon. (2015). Measures of precision for dissimilarity-based multivariate analysis of ecological communities. Ecology Letters 18:66-73.

-Guerra-Castro, E. J., J. C. Cajas, F. N. Dias Marques Simoes, J. J. Cruz-Motta, and M. Mascaro. (2020). SSP: An R package to estimate sampling effort in studies of ecological communities. bioRxiv:2020.2003.2019.996991.
}
\author{
The SSP development team is Edlin Guerra-Castro, Maite Mascaro, Nuno Simoes, Juan Cruz-Motta and Juan Cajas
}
\examples{
###To speed up the simulation of these examples, the cases, sites and N were set small.

##Single site: micromollusk from Cayo Nuevo (Yucatan, Mexico)
data(micromollusk)

#Estimation of parameters of pilot data
par.mic<-assempar (data = micromollusk,
                    type= "P/A",
                    Sest.method = "average")

#Simulation of 3 data sets, each one with 20 potential sampling units from a single site
sim.mic<-simdata(par.mic, cases= 3, N = 20, sites = 1)

#Sampling and estimation of MultSE for each sample size (few repetitions
#to speed up the example)

sam.mic<-sampsd(dat.sim = sim.mic,
               Par = par.mic,
               transformation = "P/A",
               method = "jaccard",
               n = 10,
               m = 1,
               k = 3)

#Summary of MultSE for each sampling effort
summ.mic<-summary_ssp(results = sam.mic, multi.site = FALSE)

#Cut-off points to identify optimal sampling effort
opt.mic<-ioptimum(xx = summ.mic, multi.site = FALSE)

#Plot
plot_ssp(xx = summ.mic, opt = opt.mic, multi.site = FALSE)

##Multiple sites: Sponges from Alacranes National Park (Yucatan, Mexico).
data(sponges)

#Estimation of parameters of pilot data
par.spo<-assempar(data = sponges,
                  type= "counts",
                  Sest.method = "average")

#Simulation of 3 data sets, each one with 10 potential sampling units in 3 sites.
sim.spo<-simdata(par.spo, cases= 3, N = 10, sites = 3)

#Sampling and estimation of MultSE for each sampling design (few repetitions
#to speed up the example)

sam.spo<-sampsd(dat.sim = sim.spo,
                Par = par.spo,
                transformation = "square root",
                method = "bray",
                n = 10,
                m = 3,
                k = 3)

#Summary of MultSE for each sampling effort
summ.spo<-summary_ssp(results = sam.spo, multi.site = TRUE)

#Cut-off points to identify optimal sampling effort
opt.spo<-ioptimum(xx = summ.spo, multi.site = TRUE)

#Plot
plot_ssp(xx = summ.spo, opt = opt.spo, multi.site = TRUE)
}
