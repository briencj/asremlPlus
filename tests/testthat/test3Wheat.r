#devtools::test("asremlPlus")
context("model_selection3")
asr3.lib <- "D:\\Analyses\\R oldpkg" 
  
cat("#### Test using wheat example with asreml3\n")
test_that("Wheat_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asremlPlus)
  loadASRemlVersion(3, lib.loc = asr3.lib)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(Wheat.dat)
  
  # Fit initial model
  current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                        random = ~ Row + Column + units,
                        rcov = ~ ar1(Row):ar1(Column), 
                        data=Wheat.dat)
  summary(current.asr)
  info <- infoCriteria(current.asr)
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1346.76), 1e-02)
  
  # Load current fit into an asrtests object 
  # (Have to use REML because full likelihood not implemented for ASReml-R v3)
  current.asrt <- as.asrtests(current.asr, NULL, NULL, 
                              label = "Maximal model", IClikelihood = "REML")
  testthat::expect_lt(abs(current.asrt$test.summary$AIC - 1346.766), 0.10)
  
  # Check for and remove any boundary terms
  current.asrt <- rmboundary(current.asrt, IClikelihood = "REML")
  
  #Check term for within Column pairs
  current.asrt <- testranfix(current.asrt, term = "WithinColPairs", 
                             drop.fix.ns=TRUE, IClikelihood = "REML")
  
  # Test nugget term
  current.asrt <- testranfix(current.asrt, "units", positive=TRUE, IClikelihood = "REML")
  
  # Test Row autocorrelation
  current.asrt <- testresidual(current.asrt, "~ Row:ar1(Column)", 
                               label="Row autocorrelation", simpler=TRUE, 
                               IClikelihood = "REML")
  
  # Test Col autocorrelation (depends on whether Row autocorrelation retained)
  p <- getTestPvalue(current.asrt, label = "Row autocorrelation")
  testthat::expect_true((abs(p - 2.314881e-06) < 1e-05))
  { if (p <= 0.05)
    current.asrt <- testresidual(current.asrt, "~ ar1(Row):Column", 
                                 label="Col autocorrelation", simpler=TRUE,
                                 IClikelihood = "REML", update=FALSE)
    else
      current.asrt <- testresidual(current.asrt, "~ Row:Column", 
                                   label="Col autocorrelation", simpler=TRUE,
                                   IClikelihood = "REML", update=FALSE)
  }
  print(current.asrt)
  testthat::expect_equal(length(current.asrt), 3)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 3)
  testthat::expect_equal(nrow(current.asrt$test.summary), 6)
  testthat::expect_lt(abs(current.asrt$test.summary$AIC[6] - 1353.762), 0.10)
  testthat::expect_lt(abs(current.asrt$test.summary$BIC[6] - 1367.700), 0.10)
  info <- infoCriteria(current.asrt$asreml.obj)
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1353.762), 1e-03)

  # Get current fitted asreml object
  current.asr <- current.asrt$asreml.obj
  current.asr <- update(current.asr, aom=TRUE)
  
  
  # Do residuals-versus-fitted values plot
  plot(fitted(current.asr),residuals(current.asr))

  # Form variance matrix based on estimated variance parameters
  s2 <- current.asr$sigma2
  gamma.Row <- current.asr$gammas[1]
  gamma.unit <- current.asr$gammas[2]
  rho.r <- current.asr$gammas[4]
  rho.c <- current.asr$gammas[5]
  row.ar1 <- mat.ar1(order=10, rho=rho.r)
  col.ar1 <- mat.ar1(order=15, rho=rho.c)
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) + 
    gamma.unit * diag(1, nrow=150, ncol=150) + 
    mat.dirprod(row.ar1, col.ar1)
  testthat::expect_equal(nrow(V), 150)
  testthat::expect_equal(ncol(V), 150)
  V <- s2*V
  
  #Produce variogram and variogram faces plot (Stefanaova et al, 2009)
  if (requireNamespace("lattice"))
  {
    plot.asrVariogram(variogram.asreml(current.asr))
  }
  faces <- variofaces(current.asr, V=V, maxiter=50, units="addtores")
  testthat::expect_equal(nrow(faces$face1), 10)
  testthat::expect_equal(nrow(faces$face2), 15)
  
  #Get Variety predictions and all pairwise prediction differences and p-values
  Var.diffs <- predictPlus(classify = "Variety", 
                           asreml.obj=current.asr, 
                           error.intervals="halfLeast",
                           wald.tab=current.asrt$wald.tab,
                           tables = "predictions", 
                           sortFactor = "Variety")
  testthat::expect_equal(nrow(Var.diffs$predictions), 25)
  testthat::expect_equal(ncol(Var.diffs$predictions), 6)
  testthat::expect_equal(nrow(Var.diffs$differences), 25)
  testthat::expect_equal(ncol(Var.diffs$differences), 25)
  testthat::expect_equal(nrow(Var.diffs$p.differences), 25)
  testthat::expect_equal(ncol(Var.diffs$p.differences), 25)
  testthat::expect_equal(length(Var.diffs$LSD), 3)
  testthat::expect_true("lower.halfLeastSignificant.limit" %in% names(Var.diffs$predictions))
  testthat::expect_equal(Var.diffs$backtransforms, NULL)
  testthat::expect_equal(as.character(Var.diffs$predictions$Variety[[1]]),"10")
  testthat::expect_silent(plotPvalues(Var.diffs))
  testthat::expect_silent(plotPvalues(Var.diffs, show.sig = TRUE, alpha = 0.05))
  
  #Test for single-value LSDs
  diffs <- predictPlus(classify = "Variety", 
                           asreml.obj=current.asr, 
                           error.intervals="halfLeast",
                           meanLSD.type = "fact", LSDby = "Variety",
                           wald.tab=current.asrt$wald.tab,
                           tables = "predictions", 
                           sortFactor = "Variety")
  testthat::expect_equal(nrow(diffs$LSD), 25)
  testthat::expect_true("lower.halfLeastSignificant.limit" %in% names(diffs$predictions))
  testthat::expect_true(all((diffs$predictions$upper.halfLeastSignificant.limit - 
                               diffs$predictions$lower.halfLeastSignificant.limit - 
                               diffs$LSD$meanLSD[as.numfac(diffs$predictions$Variety)]) < 1e-05))
  diffs$predictions$upper.halfLeastSignificant.limit - diffs$predictions$lower.halfLeastSignificant.limit
})

