## simdata: Script to simulate assemblages
#'@importFrom stats rbinom
#'@importFrom stats rnbinom
#'@importFrom stats rpois
#'@importFrom stats rnorm
#'@export
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

