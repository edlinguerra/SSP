### Pilot data Epibiont on mangrove roots from La Restinga National Park, Venezuela ####
### Data from Guerra-Castro et al., 2011


#Required libraries
library(vegan)

#Load pilot data
#estoy usando los datos pilotos de la tesis. Se convirtieron a un diseño balanceado por limitacions
#de adonis para la permutacion, se balancearon sitios sector externo
data(pilot)

factors<-pilot[,1:2]
dat<-pilot[,3:length(pilot)]

#PERMANOVA over fourth root transformation of abundances and Bray-Curtis similaritites using a
#nested model for Sector effect and site variability

dat.t<-dat^0.25
factors$Site<-as.factor(gl(6,10, 6*10*3))

perm <- how(within = Within("none"), plots = Plots(type = "free", strata = factors$Site), nperm = 999)
adonis(dat.t~Sector+Sector%in%Site, data = factors, permutations = perm)

PERMANOVA table of results
Unique
Source	 df	       SS	    MS	Pseudo-F	P(perm)	 perms
#Se   	  2	    83245	 41623	   7.641	  0.001	   999
#Si(Se)	 15	    83076	5538.4	  3.8119	  0.001	   997
#Res	  153	2.223E+05	1452.9
#Total	170	3.902E+05

####



D<-vegdist(dat.t)
y<-factors

perm=99

permanova<-function(data, y, method, perm){
  dat.t<-data^0.25
  D<-vegdist(dat.t, method= method)
  F.values<-rep(NA, perm+1)
  dbmanova <- function(D, y) {
  D = as.matrix(D)
  N = dim(D)[1]
  a = length(levels(y[,1]))
  b = length(levels(y[,2]))
  n = N/(a*b)
  Xfull = model.matrix(~Sector + Sector/Site + 0, data = y) #model matrix
  Hfull = Xfull %*% solve(t(Xfull) %*% Xfull) %*% t(Xfull)  #Hat matrix full model
  Xf1 = model.matrix(~Sector + 0, data = y) #model matrix factor 1
  Hf1 = Xf1 %*% solve(t(Xf1) %*% Xf1) %*% t(Xf1)  #Hat matrix factor 1
  I = diag(N)  #Identity matrix
  A = -0.5 * D^2 #centered
  G = A - apply(A, 1, mean) %o% rep(1, N) - rep(1, N) %o% apply(A, 2, mean) + mean(A) #Gower
  SSt = sum(G * t(I))
  SSres = sum(G * t(I - Hfull))
  SSf1 = sum(G * t(Hf1))
  SSf2 = SSt-(SSf1+SSres)
  glres = a*b*(n - 1)
  glf1 = a - 1
  glf2 = a*(b-1)
  glt = (a*b*n)-1
  MSres = SSres/ glres #Mean square of residuals
  MSf1 = SSf1/ glf1 #Mean square of sites
  MSf2 = SSf2/glf2
  f1 = MSf1/MSf2
  f2 = MSf2/MSres
  tab<-data.frame("Soure" = c("Sector", "Sites(Sector)", "residuals", "total"),
                  "df" = c(glf1, glf2, glres, glt),
                  "SS" = c(SSf1, SSf2, SSres, SSt),
                  "MS" = c(MSf1, MSf2, MSres, NA),
                  "pseudo-F" = c(f1,f2, NA, NA))

  return(tab) #square root of components of variation
}
  Fobs<-dbmanova(D, y)
  F.values[perm+1]<-Fobs$pseudo.F[1]

  #Permutation test for Sectors
  y.perm<-y
  a = length(levels(y[,1]))
  b = length(levels(y[,2]))
  n = N/(a*b)
  for (i in 1:perm){
    CTRL <- how(plots = Plots(gl(18,10), type = "free"), within = Within(type = "none")) #verifica Bloks
    pp<-shuffle(180, CTRL)
    dat.perm<-dat.t[pp,]
    D.perm<-vegdist(dat.perm)
    Fperm<-dbmanova(D.perm,y)
    F.values[i]<-Fperm$pseudo.F[1]
  }

  ## ahora a crear las permutaciones dentro de los sitios de cada sector

  F.ordered<-sort(F.values, decreasing = TRUE)
  t<-as.numeric(table(F.ordered >= Fobs$pseudo.F[1]))
  p.f1<-t[2]/(perm+1)
  Fobs$p.perm<-c(p.f1, NA, NA, NA)
  return(Fobs)
}

CTRL <- how(plots = gl(6*10, 3), within = Within(type = "free"))
shuffle(180, CTRL)






