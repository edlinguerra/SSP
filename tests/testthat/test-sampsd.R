context("Sampling and estimation of MultSE")

library(testthat)
library(SSP)

test_that("sampsd generates a matrix including all estimated MultSE for each simulated data, combination of sample replicates and sites for each k repetition", {
  #single site
  data("epibionts")
  x<-epibionts[1:8, 2:length(epibionts)]
  par.epi<-assempar(x, "P/A")
  cases<-3
  N<- n <-30
  sites <- m <-1
  k<-10
  sim.epi<-simdata(Par = par.epi, cases = cases, N = N, sites = sites)
  results<-sampsd(sim.epi, par.epi, transformation = "P/A", method = "jaccard", n = N, m = sites, k= 10)
  expect_is(results, "matrix")
  expect_equal(nrow(results), cases * (n-1) * k)

  #several sites
  data("epibionts")
  x<-epibionts[1:32, 2:length(epibionts)]
  par.epi<-assempar(x, "P/A")
  cases<-3
  N<- n <-30
  sites <- m <-3
  k<-10
  sim.epi<-simdata(Par = par.epi, cases = cases, N = N, sites = sites)
  results<-sampsd(sim.epi, par.epi, transformation = "P/A", method = "jaccard", n = N, m = sites, k= 10)
  expect_is(results, "matrix")
  expect_equal(nrow(results), cases * (n-1) * (sites-1) * k)
})
