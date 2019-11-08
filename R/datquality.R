## datquality: Script to compare the ecological parameters of
## simulated and original data
#'@importFrom vegan vegdist
#'@importFrom vegan diversity
#'@importFrom vegan specnumber
#'@export

datquality <- function(data, dat.sim, Par, transformation, method) {
    #Definition of preliminary functions

    MVD = function(Y, transformation, method) {
        if (transformation == "square root") {
            Y.t <- Y^0.5
            D <- vegdist(Y.t, method = method)
        }
        if (transformation == "fourth root") {
            Y.t <- Y^0.25
            D <- vegdist(Y.t, method = method)
        }
        if (transformation == "Log (X+1)") {
            Y.t <- log(Y + 1)
            D <- vegdist(Y.t, method = method)
        }
        if (transformation == "P/A") {
            Y.t <- 1 * (Y > 0)
            D <- vegdist(Y.t, method = method)
        }
        if (transformation == "none") {
            Y.t <- Y
            D <- vegdist(Y.t, method = method)
        }

        n = dim(as.matrix(D))
        ss = sum(D^2)/n
        v = ss/(n - 1)
        return(v[1])
    }

    alfa = function(Y) {
        mean = mean(specnumber(Y))
        sd = sd(specnumber(Y))
        x <- c(mean, sd)
        return(x)
    }

    simpson = function(Y) {
        mean = mean(diversity(Y, index = "simpson"))
        sd = sd(diversity(Y, index = "simpson"))
        x <- c(mean, sd)
        return(x)
    }

    #data frame to host results
    parameters <- as.data.frame(matrix(nrow = 2, ncol = 6))
    rownames(parameters) <- c("Pilot data", "Simulated data")
    colnames(parameters) <- c("S.mean", "S.sd", "1-D.mean", "1-D.sd", "MVDmin", "MVDmax")

    #Diversity metrics of pilot data
    data<-data[,2:length(data)]
    parameters[1, 1:2] <- alfa(data)
    parameters[1, 3:4] <- simpson(data)
    parameters[1, 5] <- MVD(data, transformation = transformation, method = method)
    parameters[1, 6] <- parameters[1, 5]

    #Diversity metrics of simulated data
    x <- matrix(nrow = length(dat.sim), ncol = 5)

    for (i in 1:length(dat.sim)) {
        x[i, 1:2] <- alfa(dat.sim[[i]][, 1:Par$Sest])
        x[i, 3:4] <- simpson(dat.sim[[i]][, 1:Par$Sest])
        x[i, 5] <- MVD(dat.sim[[i]][, 1:Par$Sest], transformation = transformation,
                       method = method)
    }
    parameters[2, 1:4] <- apply(x[, 1:4], 2, mean)
    parameters[2, 5:6] <- range(x[, 5])

    #End of estimations
    return(parameters)
}
