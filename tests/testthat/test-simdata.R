context("Simulated data")

library(testthat)
library(SSP)

test_that("simdata generates a list of data frames with one or several sites", {
  #single site
  data("epibionts")
  x<-epibionts[1:8, 2:length(epibionts)]
  par.epi<-assempar(x, "P/A")
  cases<-3
  n<-30
  sites<-1
  sim.epi<-simdata(Par = par.epi, cases = cases, n = n, sites = sites)
  expect_is(sim.epi, "list")
  c1<-sim.epi[[1]]
  expect_equal(levels(c1$sites), "1")

  #several sites
  data("epibionts")
  x<-epibionts[1:32, 2:length(epibionts)]
  par.epi<-assempar(x, "P/A")
  cases<-3
  n<-30
  sites<-3
  sim.epi<-simdata(Par = par.epi, cases = cases, n = n, sites = sites)
  expect_is(sim.epi, "list")
  expect_equal(length(sim.epi), cases)
  c1<-sim.epi[[1]]
  expect_equal(levels(c1$sites), c("1","2","3"))
})

test_that("simdata generates data frames with rows n * sites", {
  #single site
  data("epibionts")
  x<-epibionts[1:8, 2:length(epibionts)]
  par.epi<-assempar(x, "P/A")
  cases<-3
  n<-30
  sites<-1
  sim.epi<-simdata(Par = par.epi, cases = cases, n = n, sites = sites)
  expect_is(sim.epi, "list")
  c1<-sim.epi[[1]]
  c2<-sim.epi[[2]]
  c3<-sim.epi[[3]]
  expect_is(c1, "data.frame")
  expect_is(c2, "data.frame")
  expect_is(c3, "data.frame")
  expect_equal(nrow(c1), n)
  expect_equal(nrow(c2), n)
  expect_equal(nrow(c3), n)

  #several sites
  data("epibionts")
  x<-epibionts[1:32, 2:length(epibionts)]
  par.epi<-assempar(x, "P/A")
  cases<-3
  n<-30
  sites<-3
  sim.epi<-simdata(Par = par.epi, cases = cases, n = n, sites = sites)
  expect_is(sim.epi, "list")
  c1<-sim.epi[[1]]
  c2<-sim.epi[[2]]
  c3<-sim.epi[[3]]
  expect_is(c1, "data.frame")
  expect_is(c2, "data.frame")
  expect_is(c3, "data.frame")
  expect_equal(nrow(c1), sites * n)
  expect_equal(nrow(c2), sites * n)
  expect_equal(nrow(c3), sites * n)
})

test_that("simdata do not generates negative abundances", {
  #data of type count
  data("sponges")
  par.spo<-assempar(sponges, "counts")
  cases<-3
  n<-40
  sites<-1
  sim.spo<-simdata(Par = par.spo, cases = cases, n = n, sites = sites)
  nz<-function(x){
    x1<-x[,1:(length(x)-2)]
    y<-apply(x1, 2, is.nan)
    ys<-sort(y, decreasing = F)
    ys[1]>=0
  }

  c1<-sim.spo[[1]]
  c2<-sim.spo[[2]]
  c3<-sim.spo[[3]]
  expect_true(nz(c1))
  expect_true(nz(c2))
  expect_true(nz(c3))

  #data of type cover
  data("epibionts")
  x<-epibionts[, 2:length(epibionts)]
  par.epi<-assempar(x, "cover")
  cases<-3
  n<-30
  sites<-1
  sim.epi<-simdata(Par = par.epi, cases = cases, n = n, sites = sites)
  nz<-function(x){
    x1<-x[,1:(length(x)-2)]
    y<-apply(x1, 2, is.nan)
    ys<-sort(y, decreasing = F)
    ys[1]>=0
  }

  c1<-sim.epi[[1]]
  c2<-sim.epi[[2]]
  c3<-sim.epi[[3]]
  expect_true(nz(c1))
  expect_true(nz(c2))
  expect_true(nz(c3))

})
