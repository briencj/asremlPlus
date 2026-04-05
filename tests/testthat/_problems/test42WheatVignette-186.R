# Extracted from test42WheatVignette.r:186

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("model_selection")
cat("#### Test for wheat76 example with asreml42\n")
cat("#### Test for wheat76 example using AIC with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                        random = ~ Row + Column + units,
                        residual = ~ ar1(Row):ar1(Column), 
                        maxit = 30, data=Wheat.dat)
summary(current.asr)
info <- infoCriteria(current.asr)
testthat::expect_equal(info$varDF, 5)
testthat::expect_lt(abs(info$AIC - 1346.76), 0.10)
current.asrt <- as.asrtests(current.asr, NULL, NULL, 
                              label = "Maximal model", IClikelihood = "full")
testthat::expect_lt(abs(current.asrt$test.summary$AIC - 1653.098), 0.10)
current.asrt <- rmboundary(current.asrt, IClikelihood = "full")
current.asrt <- testranfix(current.asrt, term = "WithinColPairs", 
                             drop.fix.ns=TRUE, IClikelihood = "full")
current.asrt <- changeModelOnIC(current.asrt, addRandom = "units", 
                                  label = "units", IClikelihood = "full")
current.asrt <- changeModelOnIC(current.asrt, newResidual = "Row:ar1(Column)", 
                                  label="Row autocorrelation", 
                                  IClikelihood = "full")
result <- getTestEntry(current.asrt, label = "Row autocorrelation")$action
testthat::expect_true(grepl("Unswapped", result))
{ 
    if (grepl("Unswapped", result) || grepl("Unchanged", result))
      current.asrt <- changeModelOnIC(current.asrt, newResidual = "ar1(Row):Column", 
                                      label="Col autocorrelation", 
                                      IClikelihood = "full", update=FALSE)
    else
      current.asrt <- changeModelOnIC(current.asrt, newResidual = "Row:Column", 
                                      label="Col autocorrelation", 
                                      IClikelihood = "full", update=FALSE)
  }
print(current.asrt)
testthat::expect_equal(length(current.asrt), 3)
testthat::expect_equal(nrow(current.asrt$wald.tab), 3)
testthat::expect_equal(nrow(current.asrt$test.summary), 8)
testthat::expect_lt(abs(current.asrt$test.summary$AIC[8] - 54.19476), 0.10)
testthat::expect_lt(abs(current.asrt$test.summary$BIC[8] - 51.1841258), 0.10)
