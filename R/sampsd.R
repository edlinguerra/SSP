## sampsd: Script to sample simulated data
#'@importFrom stats model.matrix
#'@importFrom vegan vegdist
#'@importFrom sampling balancedtwostage
#'@importFrom sampling getdata
#'@export

sampsd <- function(dat.sim, Par, transformation, method, multi.site = TRUE,n, p.n, sites, p.s, k) {
    names(dat.sim) <- sprintf("%i", 1:length(dat.sim))
    if (multi.site == TRUE) {
        # Function to estimate the multivariate standard errors (multiple sites and sample replicates)
        model.MSE <- function(D, y) {
            D = as.matrix(D)
            N = dim(D)[1]
            g = length(levels(factor(y$sites)))
            X = model.matrix(~factor(y$sites))  #model matrix
            H = X %*% solve(t(X) %*% X) %*% t(X)  #Hat matrix
            I = diag(N)  #Identity matrix
            A = -0.5 * D^2
            G = A - apply(A, 1, mean) %o% rep(1, N) - rep(1, N) %o% apply(A, 2, mean) + mean(A)
            MS1 = sum(G * t(H))/(g - 1)  #Mean square of sites
            MS2 = sum(G * t(I - H))/(N - g)  #Mean square of residuals

            #Components of variation and MultSE for each source of variation
            if (MS1 > MS2) {
                CV1 = (MS1 - MS2)/(N/g)
            } else {
                CV1 = 0
            }
            MSE1 = sqrt(CV1/g)
            MSE2 = sqrt(MS2/(N/g))
            return(c(MSE1, MSE2))
        }

        #Function to sample and estimate MultSE
        samp<-function (X = NULL){

            #Matrix to store the results
            p.n <- seq(2, p.n, by = 1)
            p.s <- seq(2, p.s, by = 1)
            mse.results <- matrix(nrow = length(p.n) * length(p.s) * k, ncol = 6)
            mse.results[, 2] <- base::rep(1:k, times = length(p.n) * length(p.s))
            mse.results[, 3] <- base::rep(p.s, times = length(p.n), each = k)
            mse.results[, 4] <- base::rep(p.n, times = 1, each = length(p.s) * k)
            colnames(mse.results) <- c("dat.sim", "k", "p.s", "p.n", "MSE.sites", "MSE.n")

            #Objects needed for loop in samp
            Y <- cbind(1:(n * sites))
            YPU <- as.numeric(as.vector(gl(sites, n)))
            NN<-nrow(mse.results)
            mm<-mse.results[,3]
            nn <- rep(NA, NN)
            for (i in seq_len(NN)){
                nn[i] <- mse.results[i, 3] * mse.results[i, 4]
            }
            Par<-Par
            Sest<-Par[[3]]

            #Sampling and estimation of MultSE for each sampling size
            for (i in seq_len(NN)){
                sel <- balancedtwostage(Y, selection = 1, m = mm[i], n = nn[i], PU = YPU, FALSE)

                #replace incorrect probabilities produced by balancetwostage
                sel[sel[,1]<= -1, 1] <- 0
                sel[sel[,1]>= 2, 1] <- 1

                #getting data
                dat.df<-X
                rownames(sel) <- c(1:nrow(dat.df))
                m<-sel[, 1]
                y <- dat.df[m==1,]
                dat <- as.matrix(y[ , 1:Sest])

                #Transformation of and estimatation of a dissimilarity matrix "D"

                if (transformation == "square root") {
                    dat.t <- sqrt(dat)
                    rm(dat)
                    D <- vegdist(dat.t, method = method)
                }
                if (transformation == "fourth root") {
                    dat.t <- sqrt(sqrt(dat))
                    rm(dat)
                    D <- vegdist(dat.t, method = method)
                }
                if (transformation == "Log (X+1)") {
                    dat.t <- log(dat + 1)
                    rm(dat)
                    D <- vegdist(dat.t, method = method)
                }
                if (transformation == "P/A") {
                    dat.t <- 1 * (dat > 0)
                    rm(dat)
                    D <- vegdist(dat.t, method = method, binary = TRUE)
                }
                if (transformation == "none") {
                    dat.t <- dat
                    rm(dat)
                    D <- vegdist(dat.t, method = method)
                }

                #estimation of MultSE for each source of variation
                mse <- tryCatch(model.MSE(D = D, y = y), error = function(e) {return(c(NA, NA))})

                #Store results
                mse.results[i,5] <- mse[1]
                mse.results[i,6] <- mse[2]
            }
            return(mse.results)
        }#end of samp

        #Apply samp to simulated data
        results<-lapply(dat.sim, samp)
        results<-do.call(rbind, results)

        #Identity of each simulated data
        p.n <- seq(2, p.n, by = 1)
        p.s <- seq(2, p.s, by = 1)
        results[, 1] <- base::rep(1:length(dat.sim), each = length(p.n) * length(p.s) * k)

        return(results)

    } # end of loop for multisites

    if (multi.site == FALSE) {

        # Function to estimate the multivariate standard error for a single site
        mSE <- function(D) {
            n = dim(as.matrix(D))
            ss = sum(D^2)/n
            v = ss/(n - 1)
            x = sqrt(v/n)
            return(x[1])
        }

        #Function to sample and estimate MultSE
        samp<-function (X = NULL){

            #Matrix to store the results
            p.n <- seq(2, p.n, by = 1)
            mse.results <- matrix(nrow = length(p.n) * k, ncol = 4)
            mse.results[, 2] <- rep(1:k, times = length(p.n))
            mse.results[, 3] <- rep(p.n, times = 1, each = k)
            colnames(mse.results) <- c("dat.sim", "k", "p.n", "mSE")

            #Objects needed for loop in samp
            NN<-nrow(mse.results)
            Par<-Par
            Sest<-Par[[3]]

            #estimation of MultSE for each sampling size
            for (i in seq_len(NN)) {
                y <- sample(nrow(X), mse.results[i,3])
                dat<-X[y,1:Sest]
                if (transformation == "square root") {
                    dat.t <- dat^0.5
                    rm(dat)
                    D <- vegdist(dat.t, method = method)
                }
                if (transformation == "fourth root") {
                    dat.t <- dat^0.25
                    rm(dat)
                    D <- vegdist(dat.t, method = method)
                }
                if (transformation == "Log (X+1)") {
                    dat.t <- log(dat + 1)
                    rm(dat)
                    D <- vegdist(dat.t, method = method)
                }
                if (transformation == "P/A") {
                    dat.t <- 1 * (dat > 0)
                    rm(dat)
                    D <- vegdist(dat.t, method = method, binary = TRUE)
                }
                if (transformation == "none") {
                    dat.t <- dat
                    rm(dat)
                    D <- vegdist(dat.t, method = method)
                }

                #estimation of MultSE
                mse <- tryCatch(mSE(D = D), error = function(e) {
                    return(NA)
                })

                #Store results
                mse.results[i, 4] <- mse
            }
            return(mse.results)
        }

        #Apply samp to simulated data
        results<-lapply(dat.sim, samp)
        results<-do.call(rbind, results)

        #Identity of each simulated data
        p.n <- seq(2, p.n, by = 1)
        results[, 1] <- base::rep(1:length(dat.sim), each = length(p.n) * k)

        return(results)
    } # end of loop for a single site

}



