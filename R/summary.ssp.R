## summary.ssp: Script to estimate the tendency of the MultSE for each potential sampling effort
#'@export
summary.ssp<- function(results, multi.site = TRUE) {
  lower <- function(x) {
    quantile(x, 0.025)
  }
  upper <- function(x) {
    quantile(x, 0.975)
  }
  if (multi.site == TRUE) {
    # General average and 95% quartiles, of the multSE on the scales of sites
    sites.mse <- aggregate(MSE.sites ~ dat.sim * p.s, data = results, mean)
    sites.mean <- aggregate(MSE.sites ~ p.s, data = sites.mse, mean)
    colnames(sites.mean) <- c("p.s", "mean")
    sites.lower <- aggregate(MSE.sites ~ p.s, data = sites.mse, lower)
    colnames(sites.lower) <- c("p.s", "lower")
    sites.upper <- aggregate(MSE.sites ~ p.s, data = sites.mse, upper)
    colnames(sites.upper) <- c("p.s", "upper")
    sites.results <- cbind(sites.mean, sites.upper[, 2], sites.lower[, 2])
    colnames(sites.results) <- c("samples", "mean", "upper", "lower")
    sites.results$sv <- c(rep("sites", nrow(sites.results)))

    # General average and 95% quartiles, of the multSE on the scales of samples
    n.mse <- aggregate(MSE.n ~ dat.sim * p.n, data = results, mean)
    n.mean <- aggregate(MSE.n ~ p.n, data = n.mse, mean)
    colnames(n.mean) <- c("p.n", "mean")
    n.lower <- aggregate(MSE.n ~ p.n, data = n.mse, lower)
    colnames(n.lower) <- c("p.n", "lower")
    n.upper <- aggregate(MSE.n ~ p.n, data = n.mse, upper)
    colnames(n.upper) <- c("p.n", "upper")
    n.results <- cbind(n.mean, n.upper[, 2], n.lower[, 2])
    colnames(n.results) <- c("samples", "mean", "upper", "lower")
    n.results$sv <- c(rep("samples", nrow(n.results)))
    xx <- rbind(sites.results, n.results)

    #Relativization of the MultSE to the maximum for the minimum sampling effort

    max.mse.samples<- max(xx[xx$sv=="samples", 2])
    max.mse.sites<- max(xx[xx$sv=="sites", 2])
    xx$rel<-c((xx[xx$sv=="sites", 2]/max.mse.sites)*100, (xx[xx$sv=="samples", 2]/max.mse.samples)*100)
    xx$der<-c(rep(NA, nrow(xx)))
    mse.sites<-xx[xx$sv=="sites",c(1,6)]
    mse.residual<-xx[xx$sv=="samples",c(1,6)]

    for (i in 1:nrow(mse.sites)){
      mse.sites$der[i]<-(mse.sites$rel[i+1]-mse.sites$rel[i])/(mse.sites$samples[i]+1-mse.sites$samples[i])
    }

    for (i in 1:nrow(mse.residual)){
      mse.residual$der[i]<-(mse.residual$rel[i+1]-mse.residual$rel[i])/(mse.residual$samples[i]+1-mse.residual$samples[i])
    }

    xx$der<-abs(c(mse.sites$der, mse.residual$der))
    xx$der<-round(xx$der, 3)
    return(xx)
  }
  if (multi.site == FALSE) {
    #results$p.n <- factor(results$p.n)
    n.mse <- aggregate(mSE ~ dat.sim * p.n, data = results, mean)
    n.mean <- aggregate(mSE ~ p.n, data = n.mse, mean)
    colnames(n.mean) <- c("p.n", "mean")
    n.lower <- aggregate(mSE ~ p.n, data = n.mse, lower)
    colnames(n.lower) <- c("p.n", "lower")
    n.upper <- aggregate(mSE ~ p.n, data = n.mse, upper)
    colnames(n.upper) <- c("p.n", "upper")
    xx <- cbind(n.mean, n.upper[, 2], n.lower[, 2])
    colnames(xx) <- c("samples", "mean", "upper", "lower")

    #Relativization of the MultSE to the maximum for the minimum sampling effort
    xx$rel<-(xx$mean/xx$mean[1])*100
    xx$der<-c(rep(NA, nrow(xx)))

    for (i in 1:(nrow(xx)-1)){
     xx$der[i]<-(xx$rel[i+1]-xx$rel[i])/(xx$samples[i]+1-xx$samples[i])
    }
    xx$der<-abs(round(xx$der, 3))
    return(xx)
  }
}


####



