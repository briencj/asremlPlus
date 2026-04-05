# Extracted from test42REMLRTIC.r:36

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("model_selection")
cat("#### Test for REMLRT with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
data(Wheat.dat)
asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
m1.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column + units,
                   residual = ~ ar1(Row):ar1(Column), 
                   maxit = 30, data=Wheat.dat)
summary(m1.asr)$varcomp
info <- infoCriteria(m1.asr)
testthat::expect_equal(info$varDF, 5)
testthat::expect_lt(abs(info$AIC - 1346.76764), 1e-02)
m2.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column,
                   residual = ~ ar1(Row):ar1(Column), 
                   maxit = 30, data=Wheat.dat)
summary(m2.asr)$varcomp
info <- infoCriteria(m2.asr)
testthat::expect_equal(info$varDF, 4)
testthat::expect_lt(abs(info$AIC - 1352.941), 1e-03)
testthat::expect_warning(
    test <- REMLRT(m2.asr, m1.asr), 
    regexp = "There were a total of 1 bound terms. These bound terms occur in both models")
