#devtools::test("asremlPlus")
context("model_selection")

cat("#### Test for wheat76 example with asreml4\n")
test_that("Wheat_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(Wheat.dat)
  
  # Fit initial model
  current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                        random = ~ Row + Column + units,
                        residual = ~ ar1(Row):ar1(Column), 
                        data=Wheat.dat)
  summary(current.asr)
  info <- infoCriteria(current.asr)
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1346.768), 1e-03)
  
  # Load current fit into an asrtests object
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  
  # Check for and remove any boundary terms
  current.asrt <- rmboundary(current.asrt)
  
  #Check term for within Column pairs
  current.asrt <- testranfix(current.asrt, term = "WithinColPairs", 
                             drop.fix.ns=TRUE)
  
  # Test nugget term
  current.asrt <- testranfix(current.asrt, "units", positive=TRUE)
  
  # Test Row autocorrelation
  current.asrt <- testresidual(current.asrt, "~ Row:ar1(Column)", 
                               label="Row autocorrelation", simpler=TRUE)
  
  # Test Col autocorrelation (depends on whether Row autocorrelation retained)
  p <- getTestPvalue(current.asrt, label = "Row autocorrelation")
  testthat::expect_true((abs(p - 2.314881e-06) < 1e-05))
  { if (p <= 0.05)
    current.asrt <- testresidual(current.asrt, "~ ar1(Row):Column", 
                                 label="Col autocorrelation", simpler=TRUE,
                                 update=FALSE)
    else
      current.asrt <- testresidual(current.asrt, "~ Row:Column", 
                                   label="Col autocorrelation", simpler=TRUE,
                                   update=FALSE)
  }
  print(current.asrt)
  testthat::expect_equal(length(current.asrt), 3)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 3)
  testthat::expect_equal(nrow(current.asrt$test.summary), 5)
  info <- infoCriteria(current.asrt$asreml.obj)
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1353.762), 1e-03)
  
  # Get current fitted asreml object
  current.asr <- current.asrt$asreml.obj
  current.asr <- update(current.asr, aom=TRUE)
  
  
  # Do residuals-versus-fitted values plot
  plot(fitted(current.asr),residuals(current.asr))
  
  #Produce variogram and variogram faces plot (Stefanaova et al, 2009)
  if (requireNamespace("lattice"))
  {
    plot.varioGram(varioGram(current.asr))
  }
  V <- estimateV(current.asr)
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
