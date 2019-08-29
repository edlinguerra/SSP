## sampsd: Script to sample simulated data

# Multivariate standard error for a single site
mSE <- function (D) {
        n = dim(as.matrix(D))
        ss = sum(D^2)/n
        v = ss/(n-1)
        x = sqrt(v/n)
        return(x[1])
}

# Multivariate standard errors for a hierarchical model (multiple sites and sample replicates)
model.MSE <- function(D, y) {
                D = as.matrix(D)
                N=dim(D)[1]
                g=length(levels(factor(y$sites)))
                X=model.matrix(~factor(y$sites))   #model matrix
                H=X %*% solve(t(X)%*%X)%*%t(X)   #Hat matrix
                I=diag(N)                        #Identity matrix
                A=-0.5*D^2
                G=A-apply(A,1,mean)%o%rep(1,N)-rep(1,N)%o%apply(A,2,mean)+mean(A)
                MS1=sum(G * t(H))/(g-1) #Mean square of sites
                MS2=sum(G * t(I - H))/(N-g) #Mean square of residuals
                if (MS1>MS2){
                       CV1=(MS1-MS2)/(N/g)
                } else {CV1=0}
                MSE1=sqrt(CV1/g)
                MSE2=sqrt(MS2/(N/g))
                return(c(MSE1,MSE2))
                }

# Function to sample and estimate MSE

sampsd<-function(simdata, Sest, transformation, method, multi.site = TRUE, n, p.n, sites, p.s, k){
        names(simdata)<-sprintf("%i", 1:length(simdata))
        if (multi.site == TRUE){
                p.n<-seq(3,p.n, by=1)
                p.s<-seq(2,p.s, by=1)
                mse.results<-as.data.frame(matrix(nrow= length(p.n)*length(p.s)*length(simdata)*k, ncol=6))
                id.sd<-rep(1:length(simdata), each = length(p.n)*length(p.s)*k)
                id.k<-rep(1:k, times = length(simdata)*length(p.n)*length(p.s))
                id.p.s<-rep(p.s, times = length(simdata)*length(p.n), each = k)
                id.p.n<-rep(p.n, times = length(simdata), each = length(p.s)*k)
                mse.results[,1]<-id.sd
                mse.results[,2]<-id.k
                mse.results[,3]<-id.p.s
                mse.results[,4]<-id.p.n
                colnames(mse.results)<-c("simdata", "k", "p.s", "p.n", "MSE.sites", "MSE.n")
                Y <- cbind(1:(n*sites))
                YPU<-as.numeric(as.vector(gl(sites, n)))
                        for (i in 1:nrow(mse.results)){
                               for (j in 1:length(simdata)){
                                        if (mse.results[i,1] == names(simdata[j])){
                                        selection<-balancedtwostage(Y, selection = 1, m = mse.results[i,3], n = mse.results[i,3]*mse.results[i,4], PU = YPU, FALSE)
                                        for (l in 1:length(Y)){
                                                if (selection[l,1] <= -1){
                                                       selection[l,1]<- 0
                                                }
                                                if (selection[l,1] >= 2){
                                                        selection[l,1]<- 1
                                                }
                                        }
                                        y<-getdata(simdata[[j]], selection[,1])
                                        dat<-y[,2:(Sest+1)]
                                        if (transformation=="square root"){
                                                dat.t<-dat^0.5
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method)
                                        }
                                        if (transformation=="fourth root"){
                                                dat.t<-dat^0.25
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method)
                                        }
                                        if (transformation=="Log (X+1)"){
                                                dat.t<-log(dat+1)
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method)
                                        }
                                        if (transformation=="P/A"){
                                                dat.t<- 1 * (dat>0)
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method, binary = TRUE)
                                        }
                                        if (transformation=="none"){
                                                dat.t<- dat
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method)
                                        }
                                        #mse<-try(MSE(D=D, y=y), silent = T)#por ac? voy 09/04
                                        mse<-tryCatch(model.MSE(D=D, y=y), error = function(e) {return(c(NA,NA))})
                                        mse.results[i,5]<-mse[1]
                                        mse.results[i,6]<-mse[2]
                                }
                        }
                }
        return(mse.results)
        }
        if (multi.site == FALSE){
                p.n<-seq(3, p.n, by=1)
                mse.results<-as.data.frame(matrix(nrow= length(p.n)*length(simdata)*k, ncol=4))
                id.sd<-rep(1:length(simdata), each = length(p.n)*k)
                id.k<-rep(1:k, times = length(simdata)*length(p.n))
                id.p.n<-rep(p.n, times = length(simdata), each = k)
                mse.results[,1]<-id.sd
                mse.results[,2]<-id.k
                mse.results[,3]<-id.p.n
                colnames(mse.results)<-c("simdata", "k", "p.n", "mSE")
                for (i in 1:nrow(mse.results)){
                        for (j in 1:length(simdata)){
                                if (mse.results[i,1] == names(simdata[j])){
                                        y<-simdata[[j]][sample(nrow(simdata[[j]]),mse.results[i,3]), ]
                                        dat<-y[,1:Sest]
                                        if (transformation=="square root"){
                                                dat.t<-dat^0.5
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method)
                                        }
                                        if (transformation=="fourth root"){
                                                dat.t<-dat^0.25
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method)
                                        }
                                        if (transformation=="Log (X+1)"){
                                                dat.t<-log(dat+1)
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method)
                                        }
                                        if (transformation=="P/A"){
                                                dat.t<- 1 * (dat>0)
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method, binary = TRUE)
                                        }
                                        if (transformation=="none"){
                                                dat.t<- dat
                                                rm(dat)
                                                D<-vegdist(dat.t, method = method)
                                        }
                                        #mse<-try(MSE(D=D, y=y), silent = T)#por ac? voy 09/04
                                        mse<-tryCatch(mSE(D=D), error = function(e) {return(NA)})
                                        mse.results[i,4]<-mse
                                }
                        }
                }
                return(mse.results)
        }
}

##Description
# Function that iterates sampling (from the library 'sampling') with different effort for
# each of the simulated data sets. MultSE is estimated using pseudo-variance for single site evaluation
# or Mean Squares estimates in a distance-based model (MSresidual and MSsite) for multisite evaluations.

##Usage
#sampsd(simdata, Sest, transformation, method, multi.site = TRUE, n, p.n, sites, p.s, k)

## Arguments

#simdata            List of simulated data generated by SimDat.

#Sest               Number of species estimated with Chao2. See list generated by AssemPar.

#transformation     Mathematical function to reduce the weight of very dominant species: "square root", "fourth root", "Log (X+1)", "P/A", "none".

#method             The name of any method used in vegdist{vegan} to calculate pairwise distances

#multi.site         Logical argument indicating if several sites were simulated

#n                  Total number of samples simulated for each site.

#p.n                Maxinum number of samples to take at each site. Can be equal or less than n.

#sites              Total number of sites simulated for each data set.

#p.s                Maxinum number of sites to take at each data set. Can be equal or less than sites.

#k                  Number of repetitions of each combination between n and sites

