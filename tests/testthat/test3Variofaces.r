#devtools::test("asremlPlus")
context("model_selection3")
asr3.lib <- "D:\\Analyses\\R oldpkg" 

cat("#### Test variofaces using Atieno with asreml3\n")
test_that("Variofaces_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)

  data(chkpeadat)
  testthat::expect_equal(nrow(chkpeadat), 1056)
  testthat::expect_equal(ncol(chkpeadat), 18)
  
  # Fit the final fitted asreml object
  current.asr <- asreml(fixed = Biomass.plant ~ Lines * TRT + Smarthouse/(vLanes + vPos), 
                        random = ~Smarthouse:Zone + Smarthouse:spl(vLanes), 
                        rcov = ~at(Smarthouse):ar1(Lane):ar1(Position), 
                        data = chkpeadat, trace = FALSE)
  summary(current.asr)$varcomp
  current.asrt <- as.asrtests(current.asr, denDF = "numeric")
  current.asrt$wald.tab
  recalcWaldTab(current.asrt, denDF="numeric", dDF.na = "maximum")
  testthat::expect_equal(length(current.asrt), 3)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 7)
  testthat::expect_equal(nrow(current.asrt$test.summary), 0)
  
  # Fit initial model
  current.asr <- asreml(Biomass.plant ~ Smarthouse + Lines*TRT , 
                        random = ~ Smarthouse:(Lane + Position),
                        rcov = ~ at(Smarthouse):ar1(Lane):ar1(Position), 
                        data=chkpeadat)
  summary(current.asr)$varcomp
  
  # Load current fit into an asrtests object
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  
  # Check for and remove any boundary terms
  current.asrt <- rmboundary(current.asrt)
  print(current.asrt)
  
  # Test Lanes autocorrelation
  current.asrt <- testresidual(current.asrt, "~ at(Smarthouse):Lane:ar1(Position)", 
                               label="Lane autocorrelation", simpler=TRUE)
  
  # Test Pos autocorrelation (depends on whether Lane autocorrelation retained)
  k <- match("Lane autocorrelation", current.asrt$test.summary$terms)
  p <- current.asrt$test.summary$p
  {if (p[k] <= 0.05)
    current.asrt <- testresidual(current.asrt, "~ at(Smarthouse):ar1(Lane):Position",  
                                 label="Pos autocorrelation", simpler=TRUE,
                                 update=FALSE)
    else
      current.asrt <- testresidual(current.asrt, "~ at(Smarthouse):Lane:Position", 
                                   label="Pos autocorrelation", simpler=TRUE,
                                   update=FALSE)
  }
  testthat::expect_equal(length(current.asrt), 3)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 5)
  testthat::expect_equal(nrow(current.asrt$test.summary), 3)
  print(current.asrt)
  
  # Get current fitted asreml object
  current.asr <- current.asrt$asreml.obj
  current.asr <- update(current.asr, aom=TRUE)
  
  
  # Form variance matrix based on estimated variance parameters
  gamma.Lane <- current.asr$gammas[1]
  s2.SW <- current.asr$gammas[2]
  rho.lSW <- current.asr$gammas[3]
  rho.pSW <- current.asr$gammas[4]
  lane.ar1SW <- mat.ar1(order=24, rho=rho.lSW) 
  pos.ar1SW <- mat.ar1(order=22, rho=rho.pSW)
  RSW <- s2.SW * mat.dirprod(lane.ar1SW, pos.ar1SW)
  s2.SE <- current.asr$gammas[5]
  rho.pSE <- current.asr$gammas[6]
  rho.lSE <- current.asr$gammas[7]
  lane.ar1SE <- mat.ar1(order=24, rho=rho.lSE) 
  pos.ar1SE <- mat.ar1(order=22, rho=rho.pSE)
  RSE <- s2.SE * mat.dirprod(lane.ar1SE, pos.ar1SE)
  V <- mat.dirsum(list(RSW, RSE))
  #2 Smarthouses = 1056
  V <- gamma.Lane * with(chkpeadat, fac.sumop(fac.combine(list(Smarthouse,Lane)))) + V 
  testthat::expect_equal(nrow(V), 1056)
  testthat::expect_equal(ncol(V), 1056)
  
  #Produce variogram faces plot (Stefanaova et al, 2009)
  faces <- variofaces(current.asr, V=V, nsim = 50, maxiter = 20, seed = 14522)
  testthat::expect_equal(nrow(faces$face1), 48)
  testthat::expect_equal(ncol(faces$face1), 6)
  testthat::expect_equal(nrow(faces$face2), 44)
  testthat::expect_equal(ncol(faces$face2), 6)
  
  #Get Variety predictions
  Var.pv <- predict(current.asr, classify="Lines")$predictions$pvals
  testthat::expect_equal(nrow(Var.pv), 247)
  testthat::expect_equal(ncol(Var.pv), 4)
})

