#devtools::test("asremlPlus")
context("model_selection")
asr41.lib <- "D:\\Analyses\\R ASReml4.1" 

cat("#### Test variofaces using Atieno with asreml41\n")
test_that("Variofaces_asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)

  data(chkpeadat)
  testthat::expect_equal(nrow(chkpeadat), 1056)
  testthat::expect_equal(ncol(chkpeadat), 18)

  # Fit a model and caluclate the Wald table
  current.asr <- asreml(fixed = Biomass.plant ~ Lines * TRT + Smarthouse/(vLanes + vPos), 
                        random = ~Smarthouse:Zone + Smarthouse:spl(vLanes), 
                        residual = ~dsum(~ar1(Lane):ar1(Position) | Smarthouse), 
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
                        residual = ~ dsum(~ar1(Lane):ar1(Position) | Smarthouse), 
                        data=chkpeadat)
  summary(current.asr)$varcomp
  
  # Load current fit into an asrtests object
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  
  # Check for and remove any boundary terms
  current.asrt <- rmboundary(current.asrt)
  print(current.asrt)
  
  # Test Lanes autocorrelation
  current.asrt <- testresidual(current.asrt, "~ dsum(~Lane:ar1(Position) | Smarthouse)", 
                               label="Lane autocorrelation", simpler=TRUE)
  
  # Test Pos autocorrelation (depends on whether Lane autocorrelation retained)
  k <- match("Lane autocorrelation", current.asrt$test.summary$terms)
  p <- current.asrt$test.summary$p
  {if (p[k] <= 0.05)
    current.asrt <- testresidual(current.asrt, "~ dsum(~ar1(Lane):Position | Smarthouse)", 
                                 label="Pos autocorrelation", simpler=TRUE,
                                 update=FALSE)
    else
      current.asrt <- testresidual(current.asrt, "~ dsum(~Lane:Position | Smarthouse)", 
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

  #Produce variogram faces plot (Stefanaova et al, 2009)
  faces <- variofaces(current.asr, nsim = 50, maxiter = 20, seed = 14522)
  testthat::expect_equal(nrow(faces$face1), 48)
  testthat::expect_equal(ncol(faces$face1), 6)
  testthat::expect_equal(nrow(faces$face2), 44)
  testthat::expect_equal(ncol(faces$face2), 6)
  
  #Get Variety predictions
  Var.pv <- predict(current.asr, classify="Lines")$pvals
  testthat::expect_equal(nrow(Var.pv), 247)
  testthat::expect_equal(ncol(Var.pv), 4)
})
