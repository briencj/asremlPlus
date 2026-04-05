# Extracted from test42WheatVignette.r:259

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("model_selection")
cat("#### Test for wheat76 example with asreml42\n")
cat("#### Test for wheat76 example using AIC with asreml42\n")
cat("#### Test for wheat76 addtoSummary with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
data(Wheat.dat)
ar1.asr <- do.call(asreml, 
                     list(yield ~ Rep + WithinColPairs + Variety, 
                          random = ~ Row + Column + units,
                          residual = ~ ar1(Row):ar1(Column), 
                          data=Wheat.dat))
ar1.asrt <- as.asrtests(ar1.asr, NULL, NULL, 
                          label = "Autocorrelation model")
ar1.asrt <- rmboundary.asrtests(ar1.asrt)
testthat::expect_equal(nrow(ar1.asrt$test.summary),2)
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
