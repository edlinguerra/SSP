\name{summary.ssp}
\alias{summary.ssp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sumamary of MultSE for each sampling effort in simulated data sets
}
\description{
For each simulated data set, averages of MultSE are estimated at each sampling size. Then an overall mean, as well as lower and upper intervals of means for each sample size are tabulated. DESCRIPCIÓN BÁSICA DE LA ESTANDARIZACIÓN AL MÁXIMO E IDENTIFICACIÓN DE TASA DE CAMBIO.
}
\usage{
mseplot(results, multi.site = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{
A data frame indluding all estimated MultSE for each simulated data, combination of sample replicates and sites for each k repetition
}
  \item{multi.site}{
Logical argument indicating whether several sites were simulated
}
}
\details{
The average of the MultSE generated for each source of variation according to the repetitions of each sampling design are estimated for each simulated data. Values of the overall mean, lower and upper quantiles are tabulated for each sampling effort accross all simulated estimations.  DETALLES SOBRE CÓMO Y PORQUÉ LA ESTÁNDARIZACIÓN. DETALLES SOBRE CÓMO Y PORQUÉ LA DERIVADA FINITA. The values are summarized in a dataframe.
}
\value{
\item{mse.results}{
A data frame indluding the summary of multivariate standard error for each sampling effort.
}
}
\references{
Anderson, M.J. & Santana-Garcon, J. (2015) Measures of precision for dissimilarity-based multivariate analysis of ecological communities. Ecology Letters, 18, 66-73
}
\author{
Edlin Guerra-Castro (edlinguerra@gmail.com), Maite Mascaro (mmm@ciencias.unam.mx), Juan Carlos Cajas (carlos.cajas@enesmerida.unam.mx)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## micromollusk from 'Cayo Nuevo', Campeche bank, Gulf of Mexico, Yucatan, Mexico
data(micromollusk)
data(micromolluskGEO)

#Estimation of parameters of pilot data
dat.par1<-assempar (data = micromollusk,
                    type= "P/A",
                    spatial.str = TRUE,
                    coords = micromolluskGEO,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "random")

#Simulation of 10 data sets, each one with 100 potential sampling units from a single site
dat.sim1<-simdata(dat.par1, cases= 10, n= 100, sites = 1)

#Sampling and estimation of MSE for each sample size
samples1<-sampsd(dat.sim1,
               Sest = 98,
               transformation = "P/A",
               method = "jaccard",
               multi.site = FALSE,
               n=100,
               p.n = 50,
               sites = 1,
               p.s = 1,
               k=10)
#Summary of MSE for each sampling effort
mse.means<-mseplot(results = samples1, multi.site = FALSE)

#plotting MSE and sampling effort
library(ggplot2)
mse.means50<-mse.means[1:48,]
plot1<-ggplot(mse.means50, aes(x=levels, y=mean))+
        geom_point(size = 0.5)+
        geom_errorbar(aes(ymin=lower, ymax=upper), size=0.1, width=.2)+
        theme_bw(base_size=16) +
        ylab ("Multivariate pseudo SE")+
        xlab("Sampling effort (n)")+
        scale_y_continuous(breaks=seq(0.0, 0.3, 0.025))+
        scale_x_discrete(breaks=seq(4, 50, 2))+
        theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
              axis.text.y = element_text(colour="black", size=rel(0.7)),
              axis.title.x = element_text(colour="black", size=rel(0.9)),
              axis.title.y = element_text(colour="black", size=rel(0.9)),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size=0.4),
              axis.ticks= element_line(size=0.2))
plot1

##Sponges from Alacranes National Park, Gulf of Mexico
data(sponges)
data(spongesGEO)

#Estimation of parameters of pilot data
dat.par2<-assempar(data = sponges,
                    type= "counts",
                    spatial.str = TRUE,
                    coords = spongesGEO,
                    moranI.alfa=0.05,
                    spatial.unsee.sp= "random")

#Simulation of 10 data sets, each one with 20 potential sampling units (transects) in 20 sites.
dat.sim2<-simdata(dat.par2, cases= 10, n= 20, sites = 20)

#Sampling and estimation of MSE for each combination of sites and samples
samples2<-sampsd(dat.sim2,
                Sest = 50,
                transformation = "square root",
                method = "bray",
                multi.site = TRUE,
                n=20,
                p.n = 10,
                sites = 20,
                p.s = 10,
                k=10)

#Summary of MSE for each sampling effort
mse.means2<-mseplot(results = samples2, multi.site = TRUE)

#plotting MSE and sampling effort
library(ggplot2)
mse.means2[20:37,5]<-"transects"
plot2<-ggplot(mse.means2, aes(x=levels, y=mean))+
        geom_point()+
        geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)+
        facet_grid(.~sv)+
        theme_bw(base_size=16) +
        ylab ("Multivariate pseudo SE")+
        xlab("Sampling effort")+
        scale_y_continuous(breaks=seq(0.0, 0.3, 0.05))+
        theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
        axis.text.y = element_text(colour="black", size=rel(0.7)),
        axis.title.x = element_text(colour="black", size=rel(0.9)),
        axis.title.y = element_text(colour="black", size=rel(0.9)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0.4),
        axis.ticks= element_line(size=0.2))
plot2

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate	&	Multivariate Techniques }% use one of  RShowDoc("KEYWORDS")
\keyword{ ecological communities }
