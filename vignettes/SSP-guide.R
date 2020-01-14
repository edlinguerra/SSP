## ----setup, include = FALSE----------------------------------------------
library(SSP)
library(ggplot2)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.retina=2,
  fig.align='center',
  fig.width = 7, 
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----eval=FALSE----------------------------------------------------------
#  ## install SSP package from CRAN
#  install.packages("SSP")
#  
#  ## install the latest version from github
#  install.packages('devtools')
#  library(devtools)
#  install_github('edlinguerra/SSP')
#  
#  ## import packages
#  library(SSP)
#  library(ggplot2)

## ----eval=FALSE----------------------------------------------------------
#  data(micromollusk)
#  
#  #Estimation of parameters
#  par.mic<-assempar(data = micromollusk, type = "P/A")
#  
#  #Simulation of data
#  sim.mic<-simdata(Par = par.mic, cases = 10, n = 50, site = 1)
#  
#  #Sampling and estimation of MultSE
#  samp.mic<-sampsd(sim.mic, par.mic,
#                          transformation = "P/A",
#                          method = "jaccard",
#                          multi.site = FALSE,
#                          n=50,
#                          p.n = 50,
#                          sites = 1,
#                          m = 1,
#                          k=10)
#  
#  #Summarizing results
#  sum.mic<-summary_ssp(results = samp.mic, multi.site = FALSE)
#  
#  #Identification of optimal effort
#  
#  opt.mic<-ioptimum(xx = sum.mic, multi.site = FALSE, c1=5, c2=3, c3=1)
#  

## ----eval=FALSE----------------------------------------------------------
#  fig.mic<-ggplot(sum.mic, aes(x=samples, y=mean))+
#    geom_point(size = 0.5)+
#    geom_errorbar(aes(ymin=lower, ymax=upper), size=0.1, width=.2)+
#    theme_bw(base_size=16) +
#    ylab ("Multivariate pseudo SE")+
#    xlab("Sampling effort (n)")+
#    scale_y_continuous(breaks=seq(0.0, 0.4, 0.025))+
#    scale_x_continuous(breaks=seq(2, 50, 2))+
#    theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
#          axis.text.y = element_text(colour="black", size=rel(0.7)),
#          axis.title.x = element_text(colour="black", size=rel(0.9)),
#          axis.title.y = element_text(colour="black", size=rel(0.9)),
#          panel.grid.major = element_blank(),
#          panel.grid.minor = element_blank(),
#          panel.border = element_rect(size=0.4),
#          axis.ticks= element_line(size=0.2))+
#    annotate("rect", xmin=opt.mic[2], xmax=opt.mic[3], ymin=min(sum.mic$lower),
#             ymax=max(sum.mic$upper), alpha=.1, fill="blue")+
#    annotate("text", x=12,  y =max(sum.mic$mean), label = "Optimal effort",
#             fontface = "bold", size = 3 )
#  fig.mic

## ---- echo = FALSE, out.width='85%', fig.align='center', fig.cap='Fig. 1. MultSE and smapling effort relationship using micromollusk simulated data'----
knitr::include_graphics('fig1.png')

## ----eval=FALSE----------------------------------------------------------
#  data(sponges)
#  
#  #Estimation of parameters
#  par.spo<-assempar(data = sponges, type = "counts")
#  
#  #Simulation of data
#  sim.spo<-simdata(Par = par.spo, cases = 10, N = 20, sites = 20)
#  
#  #Sampling and estimation of MultSE
#  samp.spo<-sampsd(sim.spo, par.spo,
#                          transformation = "square root",
#                          method = "bray",
#                          multi.site = TRUE,
#                          N = 20,
#                          n = 20,
#                          sites = 20,
#                          m = 20,
#                          k=10)
#  
#  #Summarizing results
#  sum.spo<-summary_ssp(results = samp.spo, multi.site = TRUE)
#  
#  #Identification of optimal effort
#  
#  opt.spo<-ioptimum(xx = sum.spo, multi.site = TRUE, c1=5, c2=3, c3=1)
#  
#  
#  #Data frame to generate shade rectangles for sub optimal and optimal efforts for each source of variation
#  
#  shade.opt<-data.frame("xmin" = c(opt.spo[1,2], opt.spo[2,2]),
#                        "xmax" = c(opt.spo[1,3], opt.spo[2,3]),
#                        "ymin" = c(min(sum.spo$lower), min(sum.spo$lower)),
#                        "ymax" = c(max(sum.spo$lower), max(sum.spo$lower)),
#                        "sv" = c("sites", "samples"),
#                        "labels" = c("Optimal improvement", "Optimal improvement"),
#                        "x" = c(13, 12))
#  
#  shade.sub<-data.frame("xmin" = c(opt.spo[1,1], opt.spo[2,1]),
#                        "xmax" = c(opt.spo[1,2], opt.spo[2,2]),
#                        "ymin" = c(min(sum.spo$lower), min(sum.spo$lower)),
#                        "ymax" = c(max(sum.spo$lower), max(sum.spo$lower)),
#                        "sv" = c("sites", "samples"),
#                        "labels" = c("Sub", "Sub"),
#                        "x" = c(7, 7))
#  
#  
#  fig.spo<-ggplot(sum.spo, aes(x=samples, y=mean))+
#          geom_point()+
#          geom_errorbar(aes(ymin=lower, ymax=upper), width=.1)+
#          facet_grid(.~sv)+
#          theme_bw(base_size=16) +
#          ylab ("Multivariate pseudo SE")+
#          xlab("Sampling effort")+
#          theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
#                axis.text.y = element_text(colour="black", size=rel(0.7)),
#                axis.title.x = element_text(colour="black", size=rel(0.9)),
#                axis.title.y = element_text(colour="black", size=rel(0.9)),
#                panel.grid.major = element_blank(),
#                panel.grid.minor = element_blank(),
#                panel.border = element_rect(size=0.4),
#                axis.ticks= element_line(size=0.2))+
#          geom_rect(data = shade.opt, aes(x = NULL,y = NULL,xmin=xmin, xmax=xmax,
#                                    ymin=ymin, ymax=ymax),
#                                    alpha=0.1, fill="blue")+
#          geom_text(data = shade.opt, y=0.25, aes(x = x, label=labels))+
#          geom_rect(data = shade.sub, aes(x = NULL,y = NULL,xmin=xmin, xmax=xmax,
#                                    ymin=ymin, ymax=ymax),
#                                    alpha=0.1, fill="green")+
#          geom_text(data = shade.sub, y=0.25, aes(x = x, label=labels))
#  
#  
#  fig.spo
#  
#  

## ---- echo = FALSE, out.width='85%', fig.align='center', fig.cap='Fig. 2. MultSE and smapling effort relationship using sponge simulated data'----
knitr::include_graphics('fig2.png')

## ------------------------------------------------------------------------
dat<-sponges[,2:length(sponges)]

#Square root transformation of abundances
dat.t<-sqrt(dat)

#Bray-Curtys
library(vegan)
bc<-vegdist(dat.t, method = "bray")

#function to estimate components of variation in PERMANOVA 
cv.permanova <- function(D, y) {
  D = as.matrix(D)
  N = dim(D)[1]
  g = length(levels(y[,1]))
  X = model.matrix(~y[,1])  #model matrix
  H = X %*% solve(t(X) %*% X) %*% t(X)  #Hat matrix
  I = diag(N)  #Identity matrix
  A = -0.5 * D^2
  G = A - apply(A, 1, mean) %o% rep(1, N) - rep(1, N) %o% apply(A, 2, mean) + mean(A)
  MS1 = sum(G * t(H))/(g - 1)  #Mean square of sites
  MS2 = sum(G * t(I - H))/(N - g)  #Mean square of residuals
  CV1 = (MS1 - MS2)/(N/g)# Components of variation of sites
  CV2 = MS2 # Components of variation of samples
  CV = c(CV1, CV2)
  sqrtCV = sqrt(CV)
  return(sqrtCV) #square root of components of variation
}

cv<-cv.permanova(D = bc, y = sponges)
cv


