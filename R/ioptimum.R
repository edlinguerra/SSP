## ioptimum: Identification of the optimal sampling effort.
#'@export
ioptimum<-function(xx, multi.site = TRUE, c1 = 10, c2 =5, c3= 2.5){

  ##MultSE for multiple sites
  if (multi.site == TRUE) {
    #subsetting data in sites and samples
    mse.sites<- xx[xx$sv=="sites", c(1,7)]
    mse.samples<- xx[xx$sv=="samples", c(1,7)]

    #Find the minimun sampling effort for improvement of c1, c2 and c3 for each sampling unit.
    #matrix to store values
    sample.cut<-matrix(nrow=2, ncol = 3)
    rownames(sample.cut)<-c("sites", "Samples")
    colnames(sample.cut)<-c("c1", "c2", "c3")

    #cut points for sites
    s.c1<-mse.sites[mse.sites$der<=c1, 1]
    s.c2<-mse.sites[mse.sites$der<=c2, 1]
    s.c3<-mse.sites[mse.sites$der<=c3, 1]

    if (is.na(s.c1[2])){
      warning("Sites: the required optimizations were not achieved, instead, the maximum number of sites sampled with 'sampsd' is presented. Increase the 'm' value in 'sampsd' to achieve the specified optimization or relax the optimization values in arguments c1, c2 and c3", call. = FALSE)
      s.c3 <- s.c2 <- s.c1 <- max(mse.sites[,1])
    }

    if (is.na(s.c2[2])){
      warning("Sites: the required optimizations (c2, c3) were not achieved, instead, the maximum number of sites sampled with 'sampsd' is presented. Increase the 'm' value in 'sampsd' to achieve the specified optimization values or relax the optimization in arguments c2 and c3", call. = FALSE)
      s.c3 <- s.c2 <- max(mse.sites[,1])
    }

    if (is.na(s.c3[2])){
      warning("Sites: the required optimizations (c3) was not achieved, instead, the maximum number of sites sampled with 'sampsd' is presented. Increase the 'm' value in 'sampsd' to achieve the specified optimization or relax the optimization in argument c3", call. = FALSE)
      s.c3<-max(mse.sites[,1])
    }


    #cut points for samples
    n.c1<-mse.samples[mse.samples$der<=c1, 1]
    n.c2<-mse.samples[mse.samples$der<=c2, 1]
    n.c3<-mse.samples[mse.samples$der<=c3, 1]

    if (is.na(n.c1[2])){
      warning("Samples: the required optimizations were not achieved, instead, the maximum effort generated with 'sampsd' is presented. Increase the 'n' value in 'sampsd' to achieve the specified optimization or relax the optimization values in arguments c1, c2 and c3", call. = FALSE)
      n.c3 <- n.c2 <- n.c1 <- max(mse.samples[,1])
      }

    if (is.na(n.c2[2])){
      warning("Samples: the required optimizations (c2, c3) were not achieved, instead, the maximum effort generated with 'sampsd' is presented. Increase the 'n' value in 'sampsd' to achieve the specified optimization values or relax the optimization in arguments c2 and c3", call. = FALSE)
      n.c3 <- n.c2 <- max(mse.samples[,1])
    }

    if (is.na(n.c3[2])){
      warning("Samples: the required optimization was not achieved, instead, the maximum effort generated with 'sampsd' is presented. Increase the 'n' value in 'sampsd' to achieve the specified optimization or relax the optimization in argument c3", call. = FALSE)
      n.c3<-max(mse.samples[,1])
      }

    #add values to the matrix
    sample.cut[1,1]<- s.c1[2]
    sample.cut[1,2]<- s.c2[2]
    sample.cut[1,3]<- s.c3[2]
    sample.cut[2,1]<- n.c1[2]
    sample.cut[2,2]<- n.c2[2]
    sample.cut[2,3]<- n.c3[2]
    return(sample.cut)

  }
  ##MultSE for a single site
  if (multi.site == FALSE){
    n.c1<-xx[xx$der<=c1,1]
    n.c2<-xx[xx$der<=c2,1]
    n.c3<-xx[xx$der<=c3,1]

    if (is.na(n.c1[2])){
      warning("The required optimizations were not achieved, instead, the maximum effort generated with 'sampsd' is presented. Increase the 'n' value in 'sampsd' to achieve the specified optimization or relax the optimization values in arguments c1, c2 and c3", call. = FALSE)
      n.c3 <- n.c2 <- n.c1 <- max(xx[,1])
      sample.cut<-c(n.c1[1], n.c2[1], n.c3[1])
      names(sample.cut)<-c("c1", "c2", "c3")
      return(sample.cut)

    }

    if (is.na(n.c2[2])){
      warning("The required optimizations (c2, c3) were not achieved, instead, the maximum effort generated with 'sampsd' is presented. Increase the 'n' value in 'sampsd' to achieve the specified optimization values or relax the optimization in arguments c2 and c3", call. = FALSE)
      n.c3 <- n.c2 <- max(xx[,1])
      sample.cut<-c(n.c1[2], n.c2, n.c3)
      names(sample.cut)<-c("c1", "c2", "c3")
      return(sample.cut)
      }

    if (is.na(n.c3[2])){
      warning("The required optimization (c3) was not achieved, instead, the maximum effort generated with sampsd is presented. Increase the 'n' value in 'sampsd' to achieve the specified optimization or relax the optimization in argument c3", call. = FALSE)
      n.c3<-max(xx[,1])
      sample.cut<-c(n.c1[2], n.c2[2], n.c3)
      names(sample.cut)<-c("c1", "c2", "c3")
      return(sample.cut)

    } else {

    sample.cut<-c(n.c1[2], n.c2[2], n.c3[2])
    names(sample.cut)<-c("c1", "c2", "c3")
    return(sample.cut)
    }
  }
}
