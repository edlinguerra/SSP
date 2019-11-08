## simdata: Script to simulate assemblages
#'@importFrom stats rbinom
#'@importFrom stats rnbinom
#'@importFrom stats rpois
#'@importFrom stats rnorm
#'@export
simdata <- function(Par, cases, n, sites) {
    Par <- Par
    n <- n
    sites <- sites

# function for simulation
    Simul <- function(Par, n, sites) {

        fs <- Par$par$fs
        fw <- Par$par$fw
        Yh <- array(data = NA, dim = c(n, Par$Sest, sites))

    ## Presence of species using Bernoulli trials
        for (p in seq_len(sites)) {
            for (i in seq_len(Par$Sest)){
                Yh[, i, p] <- rbinom(n, 1, fs[i]) # Presence in sites
                    if (Yh[1,i,p] == 1) {
                        Yh[, i, p] <- rbinom(n, 1, fw[i]) # Presence within sites
                    } else {Yh[, i, p]<-rep(0, n)
                            }
                    }
                }

# Simulation of abundances
        for (p in seq_len(sites)) {
            for (i in seq_len(Par$Sest)) {
                for (j in seq_len(n)) {
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

    #Add columns for the id of samples and sites
        Yh$n<-rep.int(1:n, sites)
        Yh$sites<-gl(n=sites,k=n)
        return(Yh)
    }

# generation of several simulated data sets
    simulated.data <- vector("list", cases)
    for (i in 1:cases) {
        simulated.data[[i]] <- Simul(Par = Par, n = n, sites = sites)
    }
    return(simulated.data)
}

