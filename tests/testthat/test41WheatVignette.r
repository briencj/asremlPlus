#devtools::test("asremlPlus")
context("model_selection")
asr41.lib <- "D:\\Analyses\\R ASReml4.1" 

cat("#### Test for wheat76 example with asreml41\n")
test_that("Wheat_asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml, lib.loc = asr41.lib)
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
  testthat::expect_lt(abs(info$AIC - 1346.768), 0.10)
  
  # Load current fit into an asrtests object
  current.asrt <- as.asrtests(current.asr, NULL, NULL, 
                              label = "Maximal model", IClikelihood = "full")
  testthat::expect_lt(abs(current.asrt$test.summary$AIC - 1653.1), 0.50)
  
  
  # Check for and remove any boundary terms
  current.asrt <- rmboundary(current.asrt, IClikelihood = "full")
  
  #Check term for within Column pairs
  current.asrt <- testranfix(current.asrt, term = "WithinColPairs", 
                             drop.fix.ns=TRUE, IClikelihood = "full")
  
  # Test nugget term
  current.asrt <- testranfix(current.asrt, "units", positive=TRUE, IClikelihood = "full")
  
  # Test Row autocorrelation
  current.asrt <- testresidual(current.asrt, "~ Row:ar1(Column)", 
                               label="Row autocorrelation", simpler=TRUE, 
                               IClikelihood = "full")
  
  # Test Col autocorrelation (depends on whether Row autocorrelation retained)
  p <- getTestPvalue(current.asrt, label = "Row autocorrelation")
  testthat::expect_true((abs(p - 2.314881e-06) < 1e-05))
  { if (p <= 0.05)
    current.asrt <- testresidual(current.asrt, "~ ar1(Row):Column", 
                                 label="Col autocorrelation", simpler=TRUE,
                                 IClikelihood = "full", update=FALSE)
    else
      current.asrt <- testresidual(current.asrt, "~ Row:Column", 
                                   label="Col autocorrelation", simpler=TRUE,
                                   IClikelihood = "full", update=FALSE)
  }
  print(current.asrt)
  testthat::expect_equal(length(current.asrt), 3)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 3)
  testthat::expect_equal(nrow(current.asrt$test.summary), 6)
  testthat::expect_lt(abs(current.asrt$test.summary$AIC[6] - 1651.329), 0.10)
  testthat::expect_lt(abs(current.asrt$test.summary$BIC[6] - 1756.701), 0.10)
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
    plot.varioGram(varioGram.asreml(current.asr))
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
  testthat::expect_equal(length(Var.diffs$LSD), 8)
  testthat::expect_true("lower.halfLeastSignificant.limit" %in% names(Var.diffs$predictions))
  testthat::expect_equal(Var.diffs$backtransforms, NULL)
  testthat::expect_equal(as.character(Var.diffs$predictions$Variety[[1]]),"10")
  testthat::expect_silent(plotPvalues(Var.diffs))
  testthat::expect_silent(plotPvalues(Var.diffs, show.sig = TRUE, alpha = 0.05))
  
  #Test for single-value LSDs
  diffs <- predictPlus(classify = "Variety", 
                           asreml.obj=current.asr, 
                           error.intervals="halfLeast",
                           LSDtype = "fact", LSDby = "Variety",
                           wald.tab=current.asrt$wald.tab,
                           tables = "predictions", 
                           sortFactor = "Variety")
  testthat::expect_equal(nrow(diffs$LSD), 25)
  testthat::expect_true("lower.halfLeastSignificant.limit" %in% names(diffs$predictions))
  #Are the LSDs sorted?
  testthat::expect_true(all((diffs$predictions$upper.halfLeastSignificant.limit - 
                               diffs$predictions$lower.halfLeastSignificant.limit - 
                               diffs$LSD$meanLSD) < 1e-05))
  diffs$predictions$upper.halfLeastSignificant.limit - diffs$predictions$lower.halfLeastSignificant.limit
})



cat("#### Test for wheat76 addtoSummary with asreml41\n")
test_that("Wheat_addtoSummary_asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  data(Wheat.dat)
  
  ## Fit an autocorrelation model
  ar1.asr <- do.call(asreml, 
                     list(yield ~ Rep + WithinColPairs + Variety, 
                          random = ~ Row + Column + units,
                          residual = ~ ar1(Row):ar1(Column), 
                          data=Wheat.dat))
  ar1.asrt <- as.asrtests(ar1.asr, NULL, NULL, 
                          label = "Autocorrelation model")
  ar1.asrt <- rmboundary.asrtests(ar1.asrt)
  testthat::expect_equal(nrow(ar1.asrt$test.summary),2)
  
  ## Fit a tensor spline
  Wheat.dat <- within(Wheat.dat, 
                      {
                        cRow <- dae::as.numfac(Row)
                        cRow <- cRow - mean(unique(cRow))
                        cColumn <- dae::as.numfac(Column)
                        cColumn <- cColumn - mean(unique(cColumn))
                      })
  ts.asr <- do.call(asreml,
                    list(yield ~ Rep + cRow + cColumn + WithinColPairs + 
                           Variety, 
                         random = ~ spl(cRow) + spl(cColumn) + 
                           dev(cRow) + dev(cColumn) + 
                           spl(cRow):cColumn + cRow:spl(cColumn) + 
                           spl(cRow):spl(cColumn),
                         residual = ~ Row:Column, 
                         data=Wheat.dat))
  ts.asrt <- as.asrtests(ts.asr, NULL, NULL, 
                         label = "Tensor spline model")
  ts.asrt <- rmboundary.asrtests(ts.asrt)
  testthat::expect_equal(nrow(ts.asrt$test.summary),3)

  ar1.ic <- infoCriteria(ar1.asrt$asreml.obj)
  ts.ic <- infoCriteria(ts.asrt$asreml.obj)
  if (ar1.ic$AIC < ts.ic$AIC)
  {
    ic.diff <- ar1.ic - ts.ic
    new.asrt <- ar1.asrt 
    new.asrt$test.summary <- addto.test.summary(ar1.asrt$test.summary, 
                                                terms = "Compare ar1 to ts", 
                                                DF = ic.diff$varDF, 
                                                AIC = ic.diff$AIC, BIC = ic.diff$BIC, 
                                                action = "Chose ar1")
  } else
  {
    ic.diff <- ts.ic - ar1.ic
    new.asrt <- ts.asrt
    new.asrt$test.summary <- addto.test.summary(ts.asrt$test.summary, 
                                                terms = "Compare ar1 to ts", 
                                                DF = ic.diff$varDF, 
                                                AIC = ic.diff$AIC, BIC = ic.diff$BIC, 
                                                action = "Chose ts")
  }
  testthat::expect_equal(nrow(new.asrt$test.summary),3)
})
