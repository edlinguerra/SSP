#' Simulation of Ecological Data Sets
#'
#' Simulates multiple ecological data sets using parameters estimated from a pilot study. The output can be used in
#' downstream SSP functions for quality evaluation and sampling effort estimation.
#'
#' @param Par A list of parameters estimated by \code{\link{assempar}}.
#' @param cases Number of data sets to simulate.
#' @param N Number of samples to simulate in each site.
#' @param sites Number of sites to simulate in each data set.
#'
#' @details
#' Presence/absence data are simulated using Bernoulli trials based on empirical frequencies of occurrence among sites (for site-level presence)
#' and within sites (for local occurrence patterns). These matrices are then converted into abundance matrices using values drawn from
#' Poisson or negative binomial distributions (for count data), or from log-normal distributions (for continuous data like coverage or biomass),
#' depending on the aggregation properties estimated in the pilot data.
#'
#' This process is repeated \code{cases} times, producing a list of simulated data sets that reflect the statistical properties of the original
#' assemblage, but without incorporating environmental constraints or species co-occurrence structures.
#'
#' @return A list of simulated community data sets, to be used by \code{\link{datquality}} and \code{\link{sampsd}}.
#'
#' @note
#' This simulation assumes that differences in composition or abundance are due to spatial aggregation, as captured by the pilot data.
#' It does not incorporate environmental gradients or species associations. For more advanced modeling of species associations,
#' copula-based approaches as suggested by Anderson et al. (2019) may be integrated in future versions of SSP.
#'
#' @references
#' Anderson, M. J., & Walsh, D. C. I. (2013). PERMANOVA, ANOSIM, and the Mantel test in the face of heterogeneous dispersions: What null hypothesis are you testing? Ecological Monographs, 83(4), 557–574.
#'
#' Anderson, M. J., de Valpine, P., Punnett, A., & Miller, A. E. (2019). A pathway for multivariate analysis of ecological communities using copulas. Ecology and Evolution, 9, 3276–3294.
#'
#' Guerra-Castro, E.J., Cajas, J.C., Simões, N., Cruz-Motta, J.J., & Mascaró, M. (2021). SSP: an R package to estimate sampling effort in studies of ecological communities. Ecography 44(4), 561-573. doi: \doi{10.1111/ecog.05284}
#'
#' McArdle, B. H., & Anderson, M. J. (2004). Variance heterogeneity, transformations, and models of species abundance: a cautionary tale. Canadian Journal of Fisheries and Aquatic Sciences, 61, 1294–1302.
#'
#' @seealso \code{\link{sampsd}}, \code{\link{datquality}}
#'
#' @examples
#' ## Single site simulation
#' data(micromollusk)
#' par.mic <- assempar(data = micromollusk, type = "P/A", Sest.method = "average")
#' sim.mic <- simdata(par.mic, cases = 3, N = 10, sites = 1)
#'
#' ## Multiple site simulation
#' data(sponges)
#' par.spo <- assempar(data = sponges, type = "counts", Sest.method = "average")
#' sim.spo <- simdata(par.spo, cases = 3, N = 10, sites = 3)
#'
#' @importFrom stats rbinom rnbinom rnorm
#'
#' @export

simdata <- function(Par, cases, N, sites) {
    Par <- Par
    N <- N
    sites <- sites

# function for simulation
    Simul <- function(Par, N, sites) {

        fs <- Par$par$fs #Probability detection between sites
        fw <- Par$par$fw #Probability detection within sites
        Yh <- array(data = NA, dim = c(N, Par$Sest, sites))

    ## Simulation of presence of species using Bernoulli trials
        for (p in seq_len(sites)) {
            for (i in seq_len(Par$Sest)){
                Yh[, i, p] <- rbinom(N, 1, fs[i]) # Presence in sites
                    if (Yh[1,i,p] == 1) {
                        Yh[, i, p] <- rbinom(N, 1, fw[i]) # Presence within sites
                    } else {Yh[, i, p]<-rep(0, N)
                            }
                    }
                }

# Simulation of abundances
        for (p in seq_len(sites)) {
            for (i in seq_len(Par$Sest)) {
                for (j in seq_len(N)) {
                    if (Par$type == "counts" & Yh[j, i, p] == 1) {
                        if (Par$par$d[i] > 0 & Par$par$prob[i] < 0.05) {
                            Yh[j, i, p] <- rnbinom(1, size = Par$par$d[i], mu = Par$par$mean[i])
                        }
                        if (Par$par$d[i] >= 0 & Par$par$prob[i] > 0.05 | Par$par$d[i] >=
                            0 & is.na(Par$par$prob[i])) {
                            Yh[j, i, p] <- rpois(1, lambda = Par$par$mean[i])
                        }
                    }
                    if (Par$type == "cover" & Yh[j, i, p] == 1) {
                        Yh[j, i, p] <- exp(rnorm(1, mean = Par$par$logmean, sd = sqrt(Par$par$logvar)))
                    }
                }
            }
        }
    #bind the sites in a data frame
        Yh<-data.frame(apply(Yh, MARGIN = 2, rbind))
        colnames(Yh)<- Par$par$Species

    #Add columns for the id of samples and sites
        Yh$N<-rep.int(1:N, sites)
        Yh$sites<-gl(n = sites, k = N)
        return(Yh)
    }

# generation of several simulated data sets
    simulated.data <- vector("list", cases)
    for (i in 1:cases) {
        simulated.data[[i]] <- Simul(Par = Par, N = N, sites = sites)
    }
    return(simulated.data)
}

