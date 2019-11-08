## ioptimum: Identification of the optimal sampling effort.
#'@export
ioptimum<-function(xx, multi.site = TRUE, c1 = 10, c2 =5, c3= 1){
  xx<-xx

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

    s.c1<-mse.sites[mse.sites$der<=c1, 1]
    s.c2<-mse.sites[mse.sites$der<=c2, 1]
    s.c3<-mse.sites[mse.sites$der<=c3, 1]

    n.c1<-mse.samples[mse.samples$der<=c1, 1]
    n.c2<-mse.samples[mse.samples$der<=c2, 1]
    n.c3<-mse.samples[mse.samples$der<=c3, 1]

    #add values to the matrix
    sample.cut[1,1]<- s.c1[1]
    sample.cut[1,2]<- s.c2[1]
    sample.cut[1,3]<- s.c3[1]
    sample.cut[2,1]<- n.c1[1]
    sample.cut[2,2]<- n.c2[1]
    sample.cut[2,3]<- n.c3[1]
    return(sample.cut)

  }
  ##MultSE for a single site
  if (multi.site == FALSE){
    n.c1<-xx[xx$der<=c1,1]
    n.c2<-xx[xx$der<=c2,1]
    n.c3<-xx[xx$der<=c3,1]
    sample.cut<-c(n.c1[1], n.c2[1], n.c3[1])
    names(sample.cut)<-c("c1", "c2", "c3")
    return(sample.cut)

  }
}
