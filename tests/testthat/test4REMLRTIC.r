#devtools::test("asremlPlus")
context("model_selection")

cat("#### Test for REMLRT with asreml4\n")
test_that("REMLRT_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(Wheat.dat)
  
  # Fit initial model
  m1.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column + units,
                   residual = ~ ar1(Row):ar1(Column), 
                   data=Wheat.dat)
  summary(m1.asr)$varcomp
  info <- infoCriteria(m1.asr)
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1346.76764), 1e-02)
  
  #Fit model without the units term
  m2.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column,
                   residual = ~ ar1(Row):ar1(Column), 
                   data=Wheat.dat)
  summary(m2.asr)$varcomp
  info <- infoCriteria(m2.asr)
  testthat::expect_equal(info$varDF, 4)
  testthat::expect_lt(abs(info$AIC - 1352.941), 1e-03)
  test <- REMLRT(m2.asr, m1.asr)
  testthat::expect_lt(abs(test$p - 0.004232946), 1e-03)
  testthat::expect_equal(test$DF, 1)
  test <- REMLRT(m2.asr, m1.asr, DF = 1)
  testthat::expect_lt(abs(test$p - 0.004232946), 1e-03)
  testthat::expect_equal(test$DF, 1)

  m3.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column,
                   residual = ~ Row:Column, 
                   data=Wheat.dat)
  summary(m3.asr)$varcomp
  test3 <- REMLRT(m3.asr, m1.asr)
  testthat::expect_lt(abs(test3$p - 2.596812e-13), 1e-03)
  testthat::expect_equal(test3$DF, 2)
  test3 <- REMLRT(m3.asr, m1.asr, DF = 3)
  testthat::expect_lt(abs(test3$p - 1.603828e-12), 1e-03)
  testthat::expect_equal(test3$DF, 3)
  
  info <- infoCriteria(m3.asr, IClikelihood = "full")
  testthat::expect_equal(info$fixedDF, 31)
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1720.888), 5e-03)
  testthat::expect_lt(abs(info$BIC - 1823.25), 5e-03)
  testthat::expect_lt(abs(info$loglik - m3.asr$loglik), 130)
})

cat("#### Test for wheat76 example with asreml4\n")
test_that("Wheat_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  ## Fit several models to the wheat data and caclulate their ICs
  data(Wheat.dat)
  
  # Fit initial model
  m.max <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                  random = ~ Row + Column + units,
                  residual = ~ ar1(Row):ar1(Column), 
                  data=Wheat.dat)

  #Drop term for within Column pairs
  m1 <- asreml(yield ~ Rep + Variety, 
               random = ~ Row + Column + units,
               residual = ~ ar1(Row):ar1(Column), 
               data=Wheat.dat)
  
  #Drop nugget term
  m2 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
               random = ~ Row + Column,
               residual = ~ ar1(Row):ar1(Column), 
               data=Wheat.dat)

  #Drop Row autocorrelation
  m3 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                  random = ~ Row + Column + units,
                  residual = ~ Row:ar1(Column), 
                  data=Wheat.dat)

  #Drop Col autocorrelation
  m4 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
               random = ~ Row + Column + units,
               residual = ~ ar1(Row):Column, 
               data=Wheat.dat)

  mods.asr <- list(m.max, m1, m2, m3, m4)
  ic <- infoCriteria(mods.asr, IClikelihood = "full")
  testthat::expect_equal(nrow(ic), 5)
  testthat::expect_true(all(ic$fixedDF == c(31, 30, 31, 31, 31)))
  testthat::expect_true(all(ic$varDF == c(5, 5, 4, 4, 5)))
  testthat::expect_true(all(abs(ic$AIC - c(1653.100,1651.294,1654.613,1669.928,1708.997)) < 1e-01))
  testthat::expect_true(abs(ic$BIC[1] - 1761.483) < 1)
  testthat::expect_true(abs(ic$loglik[1] - (-790.5502)) < 1e-01)

})

cat("#### Test for IC with wheat94 using asreml4\n")
test_that("IC_wheat94_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(wheat94.dat)
  
  ### Start with a simple model
  fm0 <- asreml(yield ~ 1,
                random = ~ Variety + Block,
                data = wheat94.dat)
  
  fm0 <- update(fm0)
  
  current.asrt <- as.asrtests(fm0, NULL, NULL, 
                              label = "Simple model", IClikelihood = "full")
  testthat::expect_equal(nrow(current.asrt$wald.tab), 1)
  
  #Add autocorrelation
  current.asrt <- changeTerms(current.asrt, newResidual = "ar1(Col):ar1(Row)", 
                              label = "Add autocorrelation", IClikelihood = "full")
  testthat::expect_true(abs(diff(current.asrt$test.summary$AIC)) - 312.1018 < 1e-03)
  
  #Add units term
  current.asrt <- changeTerms(current.asrt, addRandom = "units", 
                              label = "Add units", IClikelihood = "full")
  
  vpar3 <- current.asrt$asreml.obj$vparameters[1:3]
  current.asrt <- iterate(current.asrt)
  testthat::expect_true(all(abs(current.asrt$asreml.obj$vparameters[1:3] - vpar3) < 1e-03))
  testthat::expect_true(abs(current.asrt$asreml.obj$loglik - -1563.459) < 1)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 1)
  
  #Add random Row and Col terms
  current.asrt <- changeTerms(current.asrt, addRandom = "Row + Col", 
                              label = "Add Row + Col", IClikelihood = "full")
  
  current.asrt <- iterate(current.asrt)
  #check that denDf for current model is the same as the number of variance parameters
  testthat::expect_true(
    length(vpc.char(current.asrt$asreml.obj)[!(vpc.char(current.asrt$asreml.obj) 
                                               %in% c("F","S","B"))])  == 
      current.asrt$test.summary$denDF[current.asrt$test.summary$terms == "Add Row + Col"])
  
  #Add fixed lin(Row) and lin(Col) terms
  current.asrt <- changeTerms(current.asrt, addFixed = "lin(Row) + lin(Col)", 
                              label = "Add lin(Row) + lin(Col)", IClikelihood = "full")
  #three fixed parameters?
  testthat::expect_equal(nrow(current.asrt$wald.tab), 
                         current.asrt$test.summary$DF[current.asrt$test.summary$terms == 
                                                        "Add lin(Row) + lin(Col)"])
  
  #Add random spl(Col) term
  current.asrt <- changeTerms(current.asrt, addRandom = "spl(Col)", 
                              label = "Add spl(Col)", 
                              IClikelihood = "full")
  testthat::expect_equal(
    length(vpc.char(current.asrt$asreml.obj)[!(vpc.char(current.asrt$asreml.obj) 
                                               %in% c("F","S","B"))]),  
    current.asrt$test.summary$denDF[current.asrt$test.summary$terms == "Add spl(Col)"][1])
  
  #Restart with fixed Rowcode and Colcode covariates, units and autocorrelation
  fm6 <- asreml(yield ~ Rowcode + Colcode,
                random = ~ Variety + Block + units,
                residual = ~ ar1(Col):ar1(Row),
                data = wheat94.dat)
  
  fm6 <- update(fm6)
  
  current.asrt <- as.asrtests(fm6, wald.tab = NULL, 
                              test.summary = current.asrt$test.summary,
                              label = "Basic + Row/Col covariates", IClikelihood = "full")
  testthat::expect_true(tail(current.asrt$test.summary$action,1) == "Starting model")
  
  #Add random Row and Col terms
  current.asrt <- changeTerms(current.asrt, addRandom = "Row + Col", 
                              label = "Add Row + Col", IClikelihood = "full")
  
  current.asrt <- iterate(current.asrt)
  
  #Add fixed lin(Row) and lin(Col) terms
  current.asrt <- changeTerms(current.asrt, addFixed = "lin(Row) + lin(Col)", 
                              label = "Add lin(Row) + lin(Col)", IClikelihood = "full")
  
  #Add random spl(Col) term
  current.asrt <- changeTerms(current.asrt, addRandom = "spl(Col)", 
                              label = "Add spl(Col)", IClikelihood = "full")
  
  current.asrt <- iterate(current.asrt)
  testthat::expect_true(nrow(current.asrt$test.summary) %in% c(13,14))
  print(current.asrt$test.summary, omit.columns = "p")
  
  
  #Start with Maximal model
  fm.max <- asreml(yield ~ lin(Row) + lin(Col) + Rowcode + Colcode,
                   random = ~ Variety + Block + Row + spl(Col) + Col + units,
                   residual = ~ ar1(Col):ar1(Row),
                   data = wheat94.dat)
  
  current.asrt <- as.asrtests(fm.max, NULL, NULL, 
                              label = "Maximal model", IClikelihood = "full")
  current.asrt <- iterate(current.asrt)
  testthat::expect_true(tail(current.asrt$test.summary$action,1) == "Starting model")
  testthat::expect_equal(current.asrt$test.summary$DF, 7)
  testthat::expect_equal(current.asrt$test.summary$denDF, 8)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 9) #includes bound Block
  
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 
                         current.asrt$test.summary$denDF[1])
  
  #Drop random Row and Col terms
  current.asrt <- changeTerms(current.asrt, dropRandom = "Row + Col", 
                              label = "Drop Row + Col", IClikelihood = "full")
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 
                         current.asrt$test.summary$denDF[3])
  
  #Drop random spl(Col) term
  current.asrt <- changeTerms(current.asrt, dropRandom = "spl(Col)", 
                              label = "Drop spl(Col)", IClikelihood = "full")
  testthat::expect_equal(
    length(vpc.char(current.asrt$asreml.obj)[!(vpc.char(current.asrt$asreml.obj) 
                                               %in% c("F","S","B"))]), 
    current.asrt$test.summary$denDF[4])
  testthat::expect_true((abs(diff(current.asrt$test.summary$BIC[3:4])) - 4.062308) < 1e-05)
  
  #Use hypothesis testing with the maximal model
  current.asrt <- as.asrtests(fm.max, NULL, test.summary = current.asrt$test.summary, 
                              label = "Maximal model", IClikelihood = "full")
  current.asrt <- iterate(current.asrt)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(current.asrt$test.summary), 6)
  
  #Test random Row term
  current.asrt <- testranfix(current.asrt, term = "Row", alpha = 0.20)
  
  #Test random Col term
  current.asrt <- testranfix(current.asrt, term = "Col", alpha = 0.20)
  current.asrt <- iterate(current.asrt)
  
  #test random spl(Col) term
  if (getTestPvalue(current.asrt, label = "Col") > 0.05)
    current.asrt <- testranfix(current.asrt, term = "spl(Col)", alpha = 0.20)
  
  #tests for getTestPvalue
  testthat::expect_true(abs(getTestPvalue(current.asrt, label = "Col") - .5944761) < 1e-05)
  testthat::expect_error(getTestPvalue(current.asrt, label = "Co"))
  
  #Test units term
  current.asrt <- testranfix(current.asrt, term = "units", alpha = 0.20)
  testthat::expect_equal(nrow(current.asrt$test.summary), 10)
  
})

cat("#### Test for getFormulae with wheat94 using asreml4\n")
test_that("Formulae_wheat94_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(wheat94.dat)

  fm.max <- asreml(yield ~ lin(Row) + lin(Col) + Rowcode + Colcode,
                   random = ~ Variety + Block + Row + spl(Col) + Col + units,
                   residual = ~ ar1(Col):ar1(Row),
                   data = wheat94.dat)
  
  mod <- getFormulae(fm.max, which = "all", envir = sys.frame(sys.nframe()))
  testthat::expect_true(all(unlist(lapply(mod, function(form) is.null(form) | 
                                            inherits(form, what = "formula")))))
  testthat::expect_true(all(names(mod) == c("fixed", "random", "residual", "sparse")))
  
  #Print fitted model
  testthat::expect_equal(length(mod), 4)
  testthat::expect_true(is.null(mod$sparse))
  
  p <- printFormulae(fm.max, which = "all")
  testthat::expect_true(all(nchar(p) > 11))
  testthat::expect_equal(length(p), 4)
  p <- printFormulae(fm.max, expanded = TRUE)
  testthat::expect_equal(length(p), 3)
  p <- printFormulae(fm.max, which = c("fixed", "random"))
  testthat::expect_equal(length(p), 2)
  p <- printFormulae(fm.max, which = "fixed")
  testthat::expect_equal(length(p), 1)
  
  #Test when have formulae in a character or list
  fix.mod <- mod$fixed
  fm.max <- asreml(fixed = fix.mod,
                   random = mod$random,
                   residual = mod$residual,
                   data = wheat94.dat)
  mod.fm <- getFormulae(fm.max, which = "all", envir = sys.frame(sys.nframe()))
  testthat::expect_true(all(unlist(lapply(mod.fm, function(form) is.null(form) | 
                                            inherits(form, what = "formula")))))
  testthat::expect_equal(length(mod), 4)
  testthat::expect_true(is.null(mod$sparse))
  p <- printFormulae(fm.max, which = "all")
  testthat::expect_true(all(nchar(p) > 11))
  
})

  