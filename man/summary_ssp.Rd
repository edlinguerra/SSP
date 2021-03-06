\name{summary_ssp}
\alias{summary_ssp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of MultSE for Each Sampling Effort in Simulated Data Sets
}
\description{
For each simulated data set, averages of MultSE are estimated for each sampling size. Then an overall mean, as well as lower and upper intervals of means for each sample size are tabulated. A relativization to the maximum is applied to the average MultSE and a numerical derivative, using a forward finite difference, of the resulting quantity is obtained.

}
\usage{
summary_ssp(results, multi.site)
}

\arguments{
  \item{results}{
A matrix generated by \code{\link{sampsd}}
}
  \item{multi.site}{
Logical argument indicating whether several sites were simulated}
}

\details{
For each set of simulated data, the average of the MultSE in each sampling effort is estimated (Anderson & Santana-Garcon 2015). Then, an overall mean, lower and upper quantiles of means are tabulated for each sampling effort among all simulated data. In order to have a general and comparable criteria to evaluate the rate of change of the average MultSE with respect to the sampling effort, a relativization to the maximum MultSE (obtained with the lower sampling effort) is calculated; then, a standard forward finite derivation is calculated.
}
\value{
\item{mse.results}{
A data frame including the summary of multivariate standard error for each sampling effort.
}
}
\references{
Anderson, M.J. & Santana-Garcon, J. (2015) Measures of precision for dissimilarity-based multivariate analysis of ecological communities. Ecology Letters, 18, 66-73

Guerra-Castro, E. J., J. C. Cajas, F. N. Dias Marques Simoes, J. J. Cruz-Motta, and M. Mascaro. (2020). SSP: An R package to estimate sampling effort in studies of ecological communities. bioRxiv:2020.2003.2019.996991.
}
\author{
Edlin Guerra-Castro (edlinguerra@gmail.com), Juan Carlos Cajas, Juan Jose Cruz-Motta, Nuno Simoes and Maite Mascaro (mmm@ciencias.unam.mx).
}
\note{
This data frame can then be used to plot MultSE with respect to the sampling effort
}

\seealso{
\code{\link[SSP]{sampsd}}, \code{\link[SSP]{ioptimum}}
}
\examples{
###To speed up the simulation of these examples, the cases, sites and n were set small.

##Single site: micromollusk from Cayo Nuevo (Yucatan, Mexico)
data(micromollusk)

#Estimation of parameters of pilot data
par.mic<-assempar (data = micromollusk,
                    type= "P/A",
                    Sest.method = "average")

#Simulation of 3 data sets, each one with 20 potential sampling units from a single site
sim.mic<-simdata(par.mic, cases= 3, N = 10, sites = 1)

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

##Multiple sites: Sponges from Alacranes National Park (Yucatan, Mexico).
data(sponges)

#Estimation of parameters of pilot data
par.spo<-assempar(data = sponges,
                  type = "counts",
                  Sest.method = "average")

#Simulation of 3 data sets, each one with 20 potential sampling units in 3 sites.
sim.spo<-simdata(par.spo, cases= 3, N = 20, sites = 3)

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

}

\keyword{ multivariate	&	Multivariate Techniques }
\keyword{ ecological communities }
\keyword{ numerical ecology }
