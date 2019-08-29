\name{assempar}
\alias{assempar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of ecological parameters of the assemblage
}
\description{
The function extracts the main parameters of the pilot data using basic R functions as well as functions like \code{\link[vegan]{specpool}} and \code{\link[ape]{Moran.I}}
}
\usage{
assempar(data, type, spatial.str, coords, moranI.alfa, spatial.unsee.sp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data frame giving P/A or abundance (counts or coverage) of species
}
  \item{type}{
Nature of the data to be processed. It may be presence / absence ("P/A"), counts of individuals ("counts"), or coverage ("cover")
}
  \item{spatial.str}{
Logical argument indicating whether spatial information (coordinates) of samples is available. If set TRUE, a matrix or data frame with X and Y coordinates must be supplied. This file must have the same number of rows than data to be used by the function Moran.I{ape}. The analysis is only valid for individual counts (argument type = "counts"), and will gave an Error if two or more samples share the same X and Y values
}
\item{coords}{
Matrix or data frame given Y and X coordinates of samples
}
\item{moranI.alfa}{
Probability for rejecting random spatial distribution. Default is < 0.05
}
\item{spatial.unsee.sp}{
Spatial distribution to be assumed for unseen species. It might be "random", "positive" or "negative"
}
}
\details{
The observed and expected number of species in the assemblage is estimated with the non-parametric method Chao2. The probability of detection of each species is estimated as the frequency of occurrence of each of the species against the sample size of pilot data. Also, the degree of spatial aggregation of species (only for real counts of individuals), is identified with the index of dispersion D (Clarke et al. 2006). Finally, when the spatial coordinates of sampling units in the pilot data is supplied, the spatial auto-correlation of each species is assessed with the Moran's I test (Legendre & Legendre 1998). The corresponding properties of unseen species are approximated using information of observed species. Specifically, the probabilities of detection are assumed to be equal to the rarest species of pilot data. The mean (and variance) of the abundance is defined as the mean of all species abundances. In addition, the spatial auto-correlation of species was assumed to be random.
}
\value{
%%  If it is a LIST, use
\item{Par}{
The function returns an object of class List, to be used by \code{\link[SSP]{simdat}}}
}
\references{
Clarke, K. R., Chapman, M. G., Somerfield, P. J., & Needham, H. R. (2006). Dispersion-based weighting of species counts in assemblage analyses. Journal of Experimental Marine Biology and Ecology, 320, 11-27.

Legendre, P. & Legendre, L. (1998). Numerical Ecology, 2nd edn. Elsevier, Amsterdam.
}
\author{
Edlin Guerra-Castro (edlinguerra@gmail.com) and Maite Mascaro (mmm@ciencias.unam.mx)
}
\note{
%%  ~~~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## micromollusk from 'Cayo Nuevo', Campeche bank, Gulf of Mexico, Yucatan, Mexico
data(micromollusk)
data(micromolluskGEO)

dat.par<-assempar (data = micromollusk,
                    type= "P/A",
                    spatial.str = TRUE,
                    coords = micromolluskGEO,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "random")
dat.par


##Sponges from Alacranes National Park, Gulf of Mexico
data(sponges)
data(spongesGEO)

dat.par<-assempar (data = sponges,
                    type= "counts",
                    spatial.str = TRUE,
                    coords = spongesGEO,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "random")
dat.par

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate	&	Multivariate Techniques }% use one of  RShowDoc("KEYWORDS")
\keyword{ ecological communities }% __ONLY ONE__ keyword per line
