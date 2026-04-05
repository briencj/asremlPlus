# Extracted from test42WheatVignette.r:21

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("model_selection")
cat("#### Test for wheat76 example with asreml42\n")

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
