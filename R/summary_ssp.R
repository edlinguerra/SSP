## summaryssp: Summary of MultSE for each sampling effort in simulated data sets
#'@importFrom stats quantile
#'@importFrom stats aggregate
#'@export
summary_ssp<- function(results, multi.site) {
  lower <- function(x) {
    quantile(x, 0.025)
  }
  upper <- function(x) {
    quantile(x, 0.975)
  }
  if (multi.site == TRUE) {
    # General average and 95% quartiles, of the multSE on the scales of sites
    sites.mse <- aggregate(MSE.sites ~ dat.sim * m, data = results, mean)
    sites.mean <- aggregate(MSE.sites ~ m, data = sites.mse, mean)
    colnames(sites.mean) <- c("m", "mean")
    sites.lower <- aggregate(MSE.sites ~ m, data = sites.mse, lower)
    colnames(sites.lower) <- c("m", "lower")
    sites.upper <- aggregate(MSE.sites ~ m, data = sites.mse, upper)
    colnames(sites.upper) <- c("m", "upper")
    sites.results <- cbind(sites.mean, sites.upper[, 2], sites.lower[, 2])
    colnames(sites.results) <- c("samples", "mean", "upper", "lower")
    sites.results$sv <- c(rep("sites", nrow(sites.results)))

    # General average and 95% quartiles, of the multSE on the scales of samples
    n.mse <- aggregate(MSE.n ~ dat.sim * n, data = results, mean)
    n.mean <- aggregate(MSE.n ~ n, data = n.mse, mean)
    colnames(n.mean) <- c("n", "mean")
    n.lower <- aggregate(MSE.n ~ n, data = n.mse, lower)
    colnames(n.lower) <- c("n", "lower")
    n.upper <- aggregate(MSE.n ~ n, data = n.mse, upper)
    colnames(n.upper) <- c("n", "upper")
    n.results <- cbind(n.mean, n.upper[, 2], n.lower[, 2])
    colnames(n.results) <- c("samples", "mean", "upper", "lower")
    n.results$sv <- c(rep("samples", nrow(n.results)))
    xx <- rbind(sites.results, n.results)

    #Relativization of the MultSE to the maximum for the minimum sampling effort

    max.mse.samples<- max(xx[xx$sv=="samples", 2])
    max.mse.sites<- max(xx[xx$sv=="sites", 2])
    xx$rel<-c((xx[xx$sv=="sites", 2]/max.mse.sites)*100, (xx[xx$sv=="samples", 2]/max.mse.samples)*100)
    xx$der<-c(rep(NA, nrow(xx)))
    mse.sites<-xx[xx$sv=="sites",c(1,6,7)]
    mse.residual<-xx[xx$sv=="samples",c(1,6,7)]

    for (i in 1:(nrow(mse.sites)-1)){
      mse.sites$der[i+1]<-(mse.sites$rel[i]-mse.sites$rel[i+1])/(mse.sites$samples[i+1]-mse.sites$samples[i])
    }
    mse.sites$cum<-c(NA,cumsum(mse.sites$der[2:nrow(mse.sites)]))

    for (i in 1:(nrow(mse.residual)-1)){
      mse.residual$der[i+1]<-(mse.residual$rel[i]-mse.residual$rel[i+1])/(mse.residual$samples[i+1]-mse.residual$samples[i])
    }
    mse.residual$cum<-c(NA,cumsum(mse.residual$der[2:nrow(mse.residual)]))

    xx$der<-abs(c(mse.sites$der, mse.residual$der))
    xx$der<-round(xx$der, 3)

    xx$cum<-c(mse.sites$cum, mse.residual$cum)
    xx$cum<-round(xx$cum, 0)
    return(xx)
  }
  if (multi.site == FALSE) {
    # General average and 95% quartiles of the multSE
    n.mse <- aggregate(mSE ~ dat.sim * n, data = results, mean)
    n.mean <- aggregate(mSE ~ n, data = n.mse, mean)
    colnames(n.mean) <- c("n", "mean")
    n.lower <- aggregate(mSE ~ n, data = n.mse, lower)
    colnames(n.lower) <- c("n", "lower")
    n.upper <- aggregate(mSE ~ n, data = n.mse, upper)
    colnames(n.upper) <- c("n", "upper")
    xx <- cbind(n.mean, n.upper[, 2], n.lower[, 2])
    colnames(xx) <- c("samples", "mean", "upper", "lower")

    #Relativization of the MultSE to the maximum for the minimum sampling effort
    xx$rel<-(xx$mean/xx$mean[1])*100
    xx$der<-c(rep(NA, nrow(xx)))

    for (i in 1:(nrow(xx)-1)){
     xx$der[i+1]<-(xx$rel[i]-xx$rel[i+1])/(xx$samples[i+1]-xx$samples[i])
    }
    xx$der<-abs(round(xx$der, 3))
    cum<-c(NA,cumsum(xx$der[2:nrow(xx)]))
    xx$cum<-round(cum, 0)
    return(xx)
  }
}


####



