## assempar: Script to estimate parameters of the assemblage
#'@export
#'@importFrom vegan specpool
#'@importFrom vegan dispweight
#'@importFrom stats rpois
#'@importFrom stats weighted.mean
#'@importFrom stats dist
#'@importFrom stats var
#'@importFrom stats variable.names

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
            Par$par[, 1] <- c(sp.names, rep("unseen species", Par$Sest - Par$Sobs))
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

            Par$par[, 1] <- c(sp.names, rep("unseen species", Par$Sest - Par$Sobs))
            Par$par[, 2] <- c(fs, rep(min(fs, na.rm = TRUE), Par$Sest - Par$Sobs))
            Par$par[, 3] <- c(fw, rep(min(fw, na.rm = TRUE), Par$Sest - Par$Sobs))
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

            Par$par[, 1] <- c(sp.names, rep("unseen species", Par$Sest - Par$Sobs))
            Par$par[, 2] <- c(fs, rep(min(fs, na.rm = TRUE), Par$Sest - Par$Sobs))
            Par$par[, 3] <- c(fw, rep(min(fw, na.rm = TRUE), Par$Sest - Par$Sobs))
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
