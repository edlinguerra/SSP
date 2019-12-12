\name{assempar}
\alias{assempar}

\title{Estimation of Ecological Parameters of the Assemblage
}
\description{
The function extracts the main parameters of the pilot data using basic R functions as well as functions like \code{\link[vegan]{specpool}} and \code{\link[vegan]{dispweight}}
}
\usage{
assempar(data, type, Sest.method)
}

\arguments{
  \item{data}{
Data frame with species names (columns) and samples (rows) information. The first column should indicate the site to which the sample belongs, regardless of whether a single site has been sampled.
}
  \item{type}{
Nature of the data to be processed. It may be presence / absence ("P/A"), counts of individuals ("counts"), or coverage ("cover")
}
 \item{Sest.method}{
Method for estimating species richness. The function \code{\link[vegan]{specpool}} is used for this. Available methods are the incidence-based Chao "chao", first order jackknife "jack1", second order jackknife "jack2" and Bootstrap "boot". By default, the "average" of the four estimates is used.
}
}
\details{
The expected number of species in the assemblage is estimated using non-parametric methods (Gotelli et al. 2011). Due to the variability in the estimates of each approximation (Reese et al. 2014), we recommend using an average of these. The probability detection of each species is estimated among and within sites. The former is calculated as the frequency of occurrences of each of the species against the number of sites sampled, the second as the weighted average frequencies in sites where the species were present. Also, the degree of spatial aggregation of species (only for real counts of individuals), is identified with the index of dispersion D (Clarke et al. 2006). The corresponding properties of unseen species are approximated using information of observed species. Specifically, the probabilities of detection are assumed to be equal to the rarest species of pilot data. The mean and variance of the abundances are defined using random Poisson values with lambda as the overall mean of species abundances.
}
\value{

\item{Par}{
The function returns an object of class List, to be used by \code{\link{simdata}}}
}
\references{
Clarke, K. R., Chapman, M. G., Somerfield, P. J., & Needham, H. R. (2006). Dispersion-based weighting of species counts in assemblage analyses. Journal of Experimental Marine Biology and Ecology, 320, 11-27.

Gotelli, N. J., & Colwell, R. K. (2011). Estimating species richness. Pages 39– 54 in A Magurranand B McGill editors. Biological diversity: frontiers in measurement and assessment. Oxford University Press, Oxford, UK.

Reese, G. C., Wilson, K. R., & Flather, C. H. (2014). Performance of species richness estimators across assemblage types and survey parameters. Global Ecology and Biogeography, 23(5), 585-594.
}
\author{
Edlin Guerra-Castro (edlinguerra@gmail.com), Juan Carlos Cajas, Juan Jose Cruz-Motta, Nuno Simoes and Maite Mascaro (mmm@ciencias.unam.mx).
}
\note{
Be careful about the data file, the first column should indicate the site ID of each sample (as character or numeric), even when a single site was sampled.
}

\seealso{
\code{\link[vegan]{dispweight}}, \code{\link[vegan]{specpool}}, \code{\link{simdata}}}

\examples{
##Single site: micromollusk from Cayo Nuevo (Yucatan, Mexico)
data(micromollusk)
par.mic<-assempar(data = micromollusk,
                  type= "P/A",
                  Sest.method = "average")
par.mic


##Multiple sites: Sponges from Alacranes National Park (Yucatan, Mexico).
data(sponges)
par.spo<-assempar(data = sponges,
                  type= "counts",
                  Sest.method = "average")
par.spo

}

\keyword{ multivariate	&	Multivariate Techniques }
\keyword{ ecological communities }
\keyword{ numerical ecology }
