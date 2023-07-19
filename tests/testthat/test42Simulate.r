#devtools::test("asremlPlus")
context("model_selection")
Sys.setenv("R_TESTS" = "")

cat("#### Test for parallel processing with asreml42\n")
test_that("parallel_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(parallel)
  library(foreach)
  library(doParallel)
  ### load a data set 
  data(Wheat.dat)
  
  #'## Analyze the wheat data n times using foreach
  n <- 100
  
  ## Testing of parallel processing
  (ncores <- parallel::detectCores())
  ncores <- 10
  (cl <- makeCluster(ncores))
  registerDoParallel(cl)
  fits <- foreach (i = 1:n, .packages = c("asreml","asremlPlus"))  %dopar%
    { 
      current.asr <- asreml(fixed    = yield ~ Rep + Variety, 
                            random   = ~ Row + units,
                            residual = ~ ar1(Row):ar1(Column), 
                            data = Wheat.dat)
    }
  stopCluster(cl)
  testthat::expect_equal(length(fits), 100)
})

cat("#### Test for simulate.asreml with asreml42\n")
test_that("simulate_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(Wheat.dat)
  
  # Fit initial model
  current.asr <- asreml(yield ~ Rep + Variety, 
                        random = ~ Row + units,
                        residual = ~ ar1(Row):ar1(Column), 
                        data=Wheat.dat, aom=TRUE)
  summary(current.asr)
  
  s2 <- current.asr$sigma2
  gamma.Row <- current.asr$vparameters[1]
  gamma.unit <- current.asr$vparameters[2]
  rho.r <- current.asr$vparameters[4]
  rho.c <- current.asr$vparameters[5]
  row.ar1 <- mat.ar1(order=10, rho=rho.r)
  col.ar1 <- mat.ar1(order=15, rho=rho.c)
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) + 
    gamma.unit * diag(1, nrow=150, ncol=150) + 
    mat.dirprod(row.ar1, col.ar1)
  V <- s2*V
  testthat::expect_equal(nrow(V), 150)
  testthat::expect_equal(ncol(V), 150)
  
  #Produce simulated and observed residuals
  resid <- simulate(current.asr, V=V, which="resid", nsim = 50)
  testthat::expect_equal(length(resid), 2)
  testthat::expect_equal(nrow(resid$observed), 150)
  testthat::expect_equal(ncol(resid$observed), 7)
  testthat::expect_equal(nrow(resid$residuals), 150)
  testthat::expect_equal(ncol(resid$residuals), 50)
  resid$residuals <- cbind(resid$observed[c("Row","Column")],
                           resid$residuals)
  #Plot variogram
  vario <- asreml::asr_varioGram(resid$observed[c("Row","Column","residuals")])
  testthat::expect_equal(nrow(vario), 150)
  testthat::expect_equal(ncol(vario), 4)
  lattice::wireframe(gamma ~ x*y, data=vario[1:3], 
                     xlab="Row differences", ylab="Column differences", zlab="")
  #Plot variofaces
  resid$residuals <- cbind(resid$observed[c("Row","Column")],
                           resid$residuals)
  plotVariofaces(data=resid$observed[c("Row","Column","residuals")],
                 residuals=resid$residuals, 
                 restype="Standardized conditional residuals")
  
  faces <- variofaces(current.asr, V=V, aom=TRUE, units="ignore", nsim = 50)
  testthat::expect_equal(nrow(faces$face1), 10)
  testthat::expect_equal(ncol(faces$face1), 5)
  testthat::expect_equal(nrow(faces$face2), 15)
  testthat::expect_equal(ncol(faces$face2), 5)
  plot.varioGram(varioGram.asreml(current.asr))

})
