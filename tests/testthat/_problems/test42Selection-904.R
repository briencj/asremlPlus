# Extracted from test42Selection.r:904

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("model_selection")
cat("#### Test for chooseModel.data.frame with asreml42\n")
cat("#### Test for chooseModel.asrtests with asreml42\n")
cat("#### Test for testing fixed at terms with asreml42\n")
cat("#### Test for changeTerms with random at terms with asreml42\n")
cat("#### Test for testing MET at terms with asreml42\n")
cat("#### Test for at terms in testswapran with asreml42\n")
cat("#### Test for spline testing with asreml42\n")
cat("#### Test for reparamSigDevn.asrtests with asreml42\n")
cat("#### Test for changeModelOnIC with wheat94 using asreml42\n")
cat("#### Test for changeModelOnIC example using asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
data(Wheat.dat)
current.asr <- do.call(asreml,
                         list(yield ~ Rep + WithinColPairs + Variety, 
                              random = ~ Row + Column + units,
                              residual = ~ ar1(Row):ar1(Column), 
                              data=Wheat.dat))
current.asr <- update(current.asr)
current.asrt <- as.asrtests(current.asr, NULL, NULL, 
                              label = "Maximal model", IClikelihood = "full")
testthat::expect_true(current.asrt$asreml.obj$converge)
testthat::expect_true(current.asrt$test.summary$action[1] == "Starting model")
testthat::expect_equal(current.asrt$test.summary$DF[1], 31)
testthat::expect_equal(current.asrt$test.summary$denDF[1], 5)
testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 6)
current.asrt <- changeModelOnIC(current.asrt, 
                                  dropRandom = "Row + Column", label = "Drop Row + Column",
                                  checkboundaryonly = TRUE,
                                  which.IC = "AIC", IClikelihood = "full")
testthat::expect_true(current.asrt$asreml.obj$converge)
testthat::expect_equal(current.asrt$test.summary$denDF[2], -1)
current.asrt <- changeModelOnIC(current.asrt, 
                                  newResidual = "Row:ar1(Column)", 
                                  label="Row autocorrelation",
                                  IClikelihood = "full")
testthat::expect_true(current.asrt$asreml.obj$converge)
testthat::expect_equal(current.asrt$test.summary$denDF[3], -2)
testthat::expect_true((abs(current.asrt$test.summary$AIC[3]) - 21.72047) < 1e-02)
mod <- printFormulae(current.asrt$asreml.obj)
testthat::expect_equal(length(mod), 3)
testthat::expect_true(grepl("units", mod[2], fixed = TRUE))
