#' Estimation of Ecological Parameters of the Assemblage
#'
#' This function extracts the main parameters of the pilot data using base R functions,
#' as well as functions like \code{\link[vegan]{specpool}} and \code{\link[vegan]{dispweight}}.
#'
#' @param data Data frame with species names (columns) and samples (rows). The first column should indicate the site to which the sample belongs, regardless of whether a single site has been sampled.
#' @param type Nature of the data to be processed. It may be presence/absence ("P/A"), counts of individuals ("counts"), or coverage ("cover").
#' @param Sest.method Method for estimating species richness. The function \code{\link[vegan]{specpool}} is used. Available methods are "chao", "jack1", "jack2", and "boot". By default, the "average" of the four estimates is used.
#'
#' @details
#' The expected number of species in the assemblage is estimated using non-parametric methods (Gotelli et al. 2011).
#' Due to variability in the estimates of each approximation (Reese et al. 2014), we recommend using the average.
#' The probability of detection of each species is estimated among and within sites. Among-site detection is calculated
#' as the frequency of occurrences of each species across sampled sites; within-site detection is calculated as the
#' weighted average of frequencies in sites where the species are present. Spatial aggregation (only for count data) is
#' evaluated using the index of dispersion D (Clarke et al. 2006). Properties of unseen species are approximated using
#' information from observed species, assuming their detection probabilities match those of the rarest observed species.
#' Abundance distributions are simulated using random Poisson values with lambda as the overall mean of observed abundances.
#'
#' @return A list (class \code{list}) containing the estimated parameters of the assemblage, to be used by \code{\link{simdata}}.
#'
#' @note Important: The first column should indicate the site ID of each sample (as character or numeric), even when only a single site was sampled.
#'
#' @references
#' Clarke, K. R., Chapman, M. G., Somerfield, P. J., & Needham, H. R. (2006). Dispersion-based weighting of species counts in assemblage analyses. Journal of Experimental Marine Biology and Ecology, 320, 11–27.
#'
#' Gotelli, N. J., & Colwell, R. K. (2011). Estimating species richness. In A. E. Magurran & B. J. McGill (Eds.), Biological diversity: frontiers in measurement and assessment (pp. 39–54). Oxford University Press.
#'
#' Guerra-Castro, E.J., Cajas, J.C., Simões, N., Cruz-Motta, J.J., & Mascaró, M. (2021). SSP: an R package to estimate sampling effort in studies of ecological communities. Ecography 44(4), 561-573. doi: \url{https://doi.org/10.1111/ecog.05284}
#'
#' Reese, G. C., Wilson, K. R., & Flather, C. H. (2014). Performance of species richness estimators across assemblage types and survey parameters. Global Ecology and Biogeography, 23(5), 585–594.
#'
#' @seealso \code{\link[vegan]{dispweight}}, \code{\link[vegan]{specpool}}, \code{\link{simdata}}
#'
#' @examples
#' ## Single site: micromollusk from Cayo Nuevo (Yucatan, Mexico)
#' data(micromollusk)
#' par.mic <- assempar(data = micromollusk, type = "P/A", Sest.method = "average")
#' par.mic
#'
#' ## Multiple sites: Sponges from Alacranes National Park (Yucatan, Mexico)
#' data(sponges)
#' par.spo <- assempar(data = sponges, type = "counts", Sest.method = "average")
#' par.spo
#'
#' @importFrom vegan specpool
#' @importFrom vegan dispweight
#' @importFrom stats rpois
#' @importFrom stats weighted.mean
#' @importFrom stats dist
#' @importFrom stats var
#' @importFrom stats variable.names
#' @importFrom stats aggregate

#' @export

assempar <- function(data, type = c("P/A", "counts", "cover"), Sest.method = "average") {

    # Functions to estimate probability detection of species among and within sites
    freq.site <- function(x,y){
        x<-factor(x)
        a<-length(levels(x))
        yy<-aggregate(y, by = list(x), sum)
        ys<-1 * (yy[,2] >0)
        fs<-sum(ys)/a
        return(fs)
    }
    freq.within <- function(x,y){
        yi <- 1 * (y > 0)
        yyi<-aggregate(yi, by = list(x), sum)
        tt<-aggregate(y, by = list(x), length)
        yw<-yyi[yyi$x!=0,2]/tt[yyi$x!=0,2]
        wt<-tt[yyi$x!=0,2]/sum(tt[yyi$x!=0,2])
        fw<-weighted.mean(yw, wt)
        return(fw)
    }

    #Type of data: P/A
    if (type == "P/A") {
        X <- data[,1] # Vector of sites
        Y <- 1 * (data[,2:length(data)] > 0)  # Ensure data is of presence/absence

        #List to store parameters for each species
        Par <- NULL
        Par$type <- "P/A"

        #Species observed and expected
        richness <- specpool(Y)
        Par$Sobs <- richness$Species
        if (Sest.method == "chao") {
            Par$Sest <- round(richness$chao)  # Extrapoled richness using Chao2
        }
        if (Sest.method == "jack1") {
            Par$Sest <- round(richness$jack1)  # Extrapoled richness using First order jackknife
        }
        if (Sest.method == "jack2") {
            Par$Sest <- round(richness$jack2)  # Extrapoled richness using Second order jackknife
        }
        if (Sest.method == "boot") {
            Par$Sest <- round(richness$boot)  # Extrapoled richness using Bootstrap
        }
        if (Sest.method == "average") {
            Par$Sest <-round(mean(richness$chao, richness$jack1, richness$jack2, richness$boot))  # Average extrapolated richness
        }

        # frequency of occurrence of each species among sites
        fs<-apply(Y, 2, freq.site, x=X)
        fs<-fs[fs!=0]

        # ensure exclude unregistered species in data
        sp.names <- names(fs)
        Y <- subset(Y, select = sp.names)

        # frequency of occurrence of each species within sites
        fw<-apply(Y, 2, freq.within, x=X)

        # data frame with parameters estimated for each species in the assemblage
        if (Par$Sest > Par$Sobs) {
            Par$par <- data.frame(Species = c(rep(NA, Par$Sest)),
                                  fs = c(rep(NA, Par$Sest)),
                                  fw = c(rep(NA, Par$Sest)))
        # labels for unseen species
            delta.sp <- Par$Sest - Par$Sobs
            unsee <- c(rep(NA, delta.sp))
            for(i in 1:delta.sp){
              unsee[i]<- paste ("unseen.species",i)
            }
        # filling data frame
            Par$par[, 1] <- c(sp.names, unsee)
            Par$par[, 2] <- c(fs, rep(min(fs, na.rm = TRUE), Par$Sest - Par$Sobs))
            Par$par[, 3] <- c(fw, rep(min(fw, na.rm = TRUE), Par$Sest - Par$Sobs))
        } else {
            Par$par <- data.frame(Species = c(rep(NA, Par$Sobs)), fw = c(rep(NA, Par$Sobs)))
            Par$par[, 1] <- sp.names
            Par$par[, 2] <- fs
            Par$par[, 3] <- fw
        }
    }

    #Type of data: counts
    if (type == "counts") {
        X <- data[,1] # Vector of sites
        Y <- data[,2:length(data)]

        #List to store parameters for each species
        Par <- NULL
        Par$type <- "counts"

        #Species observed and expected
        richness <- specpool(Y)
        Par$Sobs <- richness$Species
        if (Sest.method == "chao") {
            Par$Sest <- round(richness$chao)  # Extrapoled richness using Chao2
        }
        if (Sest.method == "jack1") {
            Par$Sest <- round(richness$jack1)  # Extrapoled richness using First order jackknife
        }
        if (Sest.method == "jack2") {
            Par$Sest <- round(richness$jack2)  # Extrapoled richness using Second order jackknife
        }
        if (Sest.method == "boot") {
            Par$Sest <- round(richness$boot)  # Extrapoled richness using Bootstrap
        }
        if (Sest.method == "average") {
            Par$Sest <-round(mean(richness$chao, richness$jack1, richness$jack2, richness$boot))  # Average extrapolated richness
        }

        # frequency of occurrence of each species among sites
        fs<-apply(Y, 2, freq.site, x=X)
        fs<-fs[fs!=0]

        # ensure exclude unregistered species in data
        sp.names <- names(fs)
        Y <- subset(Y, select = sp.names)

        # frequency of occurrence of each species within sites
        fw<-apply(Y, 2, freq.within, x=X)

        # data frame with parameters estimated for each species in the assemblage
        if (Par$Sest > Par$Sobs) {
            Par$par <- data.frame(Species = c(rep(NA, Par$Sest)),
                                  fs = c(rep(NA, Par$Sest)),
                                  fw = c(rep(NA, Par$Sest)),
                                  mean = c(rep(NA, Par$Sest)),
                                  var = c(rep(NA, Par$Sest)),
                                  D = c(rep(NA, Par$Sest)),
                                  prob = c(rep(NA, Par$Sest)),
                                  d = c(rep(NA, Par$Sest)))

            # Mean and variance of species abundances, zero values are omitted
            M.0rm <- replace(Y, Y == 0, NA)
            ave<-apply(M.0rm, 2, mean, na.rm = TRUE)
            v<-apply(M.0rm, 2, var, na.rm = TRUE)

            # Asignation of mean and variance abundance of unseen species using random poisson values
            pois<-rpois(n = Par$Sest - Par$Sobs, lambda = mean(ave))

            # labels for unseen species
            delta.sp <- Par$Sest - Par$Sobs
            unsee <- c(rep(NA, delta.sp))
            for(i in 1:delta.sp){
              unsee[i]<- paste ("unseen.species",i)
            }
            # filling data frame
            Par$par[, 1] <- c(sp.names, unsee)
            Par$par[, 2] <- c(fs, rep(min(fs, na.rm = TRUE), delta.sp))
            Par$par[, 3] <- c(fw, rep(min(fw, na.rm = TRUE), delta.sp))
            Par$par[, 4] <- c(ave, pois)
            Par$par[, 5] <- c(v, pois)
            Par$par[is.na(Par$par[,5]),5]<-0

            #Estimation of dispersion index (D) for each species
            DW <- dispweight(Y)
            Par$par[, 6] <- c(attr(DW, "D"), rep(1, Par$Sest - Par$Sobs))
            Par$par[, 7] <- c(attr(DW, "p"), rep(1, Par$Sest - Par$Sobs))
            for (i in 1:Par$Sobs) {
                if (Par$par$D[i] > 1 & Par$par$prob[i] < 0.05 & Par$par$var[i] > Par$par$mean[i] & Par$par$mean[i] > 0) {
                    Par$par$d[i] <- (Par$par$var[i] - Par$par$mean[i])/Par$par$mean[i]^2
                } else {
                    Par$par$d[i] <- 0
                }
            }
            Par$par$d[(Par$Sobs + 1):Par$Sest] <- 0

        } else {
            Par$par <- data.frame(Species = c(rep(NA, Par$Sest)),
                                  fs = c(rep(NA, Par$Sest)),
                                  fw = c(rep(NA, Par$Sest)),
                                  mean = c(rep(NA, Par$Sest)),
                                  var = c(rep(NA, Par$Sest)),
                                  D = c(rep(NA, Par$Sest)),
                                  prob = c(rep(NA, Par$Sest)),
                                  d = c(rep(NA, Par$Sest)))

            Ypa <- 1 * (Y > 0)
            M.0rm <- replace(Y, Y == 0, NA)  # remove zero by NA
            Par$par[, 1] <- sp.names
            Par$par[, 2] <- fs
            Par$par[, 3] <- fw
            Par$par[, 4] <- apply(M.0rm, 2, mean, na.rm = TRUE)
            Par$par[, 5] <- apply(M.0rm, 2, var, na.rm = TRUE)
            Par$par[is.na(Par$par[,5]),5]<-0
            DW <- dispweight(Y)
            Par$par[, 6] <- attr(DW, "D")
            Par$par[, 7] <- attr(DW, "p")
            for (i in 1:Par$Sobs) {
                if (Par$par$D[i] > 1 & Par$par$prob[i] < 0.05 & Par$par$var[i] > Par$par$mean[i] & Par$par$mean[i] > 0) {
                    Par$par$d[i] <- (Par$par$var[i] - Par$par$mean[i])/Par$par$mean[i]^2
                } else {
                    Par$par$d[i] <- 0
                }
            }
        }
    }

    #Type of data: cover
    if (type == "cover") {
        X <- data[,1] # Vector of sites
        Y <- data[,2:length(data)]

        #List to store parameters for each species
        Par <- NULL
        Par$type <- "cover"

        #Species observed and expected
        richness <- specpool(Y)
        Par$Sobs <- richness$Species  #Species observed
        if (Sest.method == "chao") {
            Par$Sest <- round(richness$chao)  # Extrapoled richness using Chao2
        }
        if (Sest.method == "jack1") {
            Par$Sest <- round(richness$jack1)  # Extrapoled richness using First order jackknife
        }
        if (Sest.method == "jack2") {
            Par$Sest <- round(richness$jack2)  # Extrapoled richness using Second order jackknife
        }
        if (Sest.method == "boot") {
            Par$Sest <- round(richness$boot)  # Extrapoled richness using Bootstrap
        }
        if (Sest.method == "average") {
            Par$Sest <-round(mean(richness$chao, richness$jack1, richness$jack2, richness$boot))  # Average extrapolated richness
        }

        # frequency of occurrence of each species among sites
        fs<-apply(Y, 2, freq.site, x=X)
        fs<-fs[fs!=0]

        # ensure exclude unregistered species in data
        sp.names <- names(fs)
        Y <- subset(Y, select = sp.names)

        # frequency of occurrence of each species within sites
        fw<-apply(Y, 2, freq.within, x=X)

        # data frame with parameters estimated for each species in the assemblage
        if (Par$Sest > Par$Sobs) {
            Par$par <- data.frame(Species = c(rep(NA, Par$Sest)),
                                  fs = c(rep(NA, Par$Sest)),
                                  fw = c(rep(NA, Par$Sest)),
                                  logmean = c(rep(NA, Par$Sest)),
                                  logvar = c(rep(NA, Par$Sest)))

            # Mean and variance of species abundances, zero values are omitted
            M.0rm <- replace(Y, Y == 0, NA)
            ave<-apply(log(M.0rm), 2, mean, na.rm = TRUE)
            v<-apply(log(M.0rm), 2, var, na.rm = TRUE)

            # Asignation of mean and variance abundance of unseen species using random poisson values
            pois<-rpois(n = Par$Sest - Par$Sobs, lambda = mean(ave))

            # labels for unseen species
            delta.sp <- Par$Sest - Par$Sobs
            unsee <- c(rep(NA, delta.sp))
            for(i in 1:delta.sp){
              unsee[i]<- paste ("unseen.species",i)
            }
            # filling data frame
            Par$par[, 1] <- c(sp.names, unsee)
            Par$par[, 2] <- c(fs, rep(min(fs, na.rm = TRUE), delta.sp))
            Par$par[, 3] <- c(fw, rep(min(fw, na.rm = TRUE), delta.sp))
            Par$par[, 4] <- c(ave, pois)
            Par$par[, 5] <- c(v, pois)
            Par$par[is.na(Par$par[,5]),5]<-0

        } else {
            Par$par <- data.frame(Species = c(rep(NA, Par$Sest)),
                                  fs = c(rep(NA, Par$Sest)),
                                  fw = c(rep(NA, Par$Sest)),
                                  logmean = c(rep(NA, Par$Sest)),
                                  logvar = c(rep(NA, Par$Sest)))
            Ypa <- 1 * (Y > 0)
            M.0rm <- replace(Y, Y == 0, NA)
            ave<-apply(log(M.0rm), 2, mean, na.rm = TRUE)
            v<-apply(log(M.0rm), 2, var, na.rm = TRUE)
            Par$par[, 1] <- sp.names
            Par$par[, 2] <- fs
            Par$par[, 3] <- fw
            Par$par[, 4] <- ave
            Par$par[, 5] <- v
            Par$par[is.na(Par$par[,5]),5]<-0
        }
    }
    return(Par)
}
