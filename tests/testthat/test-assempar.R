context("Assemblage paramaters")

library(testthat)
library(SSP)

test_that("assempar generates a data frame of lenght according to data type", {
  #P/A
  data("epibionts")
  x<-epibionts[1:8, 2:length(epibionts)]
  par.epi<-assempar(x, "P/A")
  expect_is(par.epi$par, "data.frame")
  expect_equal(length(par.epi$par), 3)

  #Counts
  data("sponges")
  x<-sponges[, 2:length(sponges)]
  par.spo<-assempar(x, "counts")
  expect_is(par.spo$par, "data.frame")
  expect_equal(length(par.spo$par), 8)

  #Cover
  data("epibionts")
  x<-epibionts[1:32, 2:length(epibionts)]
  par.epi<-assempar(x, "cover")
  expect_is(par.epi$par, "data.frame")
  expect_equal(length(par.epi$par), 5)

})

