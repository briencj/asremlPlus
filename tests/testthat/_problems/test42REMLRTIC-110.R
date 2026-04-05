# Extracted from test42REMLRTIC.r:110

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("model_selection")
cat("#### Test for REMLRT with asreml42\n")
cat("#### Test for wheat76 example with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
data(Wheat.dat)
asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
m.max <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                  random = ~ Row + Column + units,
                  residual = ~ ar1(Row):ar1(Column), 
                  maxit = 30, data=Wheat.dat)
m1 <- asreml(yield ~ Rep + Variety, 
               random = ~ Row + Column + units,
               residual = ~ ar1(Row):ar1(Column), 
               maxit = 30, data=Wheat.dat)
m2 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
               random = ~ Row + Column,
               residual = ~ ar1(Row):ar1(Column), 
               maxit = 30, data=Wheat.dat)
m3 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                  random = ~ Row + Column + units,
                  residual = ~ Row:ar1(Column), 
                  data=Wheat.dat)
m4 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
               random = ~ Row + Column + units,
               residual = ~ ar1(Row):Column, 
               maxit = 30, data=Wheat.dat)
mods.asr <- list(m.max, m1, m2, m3, m4)
ic <- infoCriteria(mods.asr, IClikelihood = "full")
testthat::expect_equal(nrow(ic), 5)
testthat::expect_true(all(ic$fixedDF == c(31, 30, 31, 31, 31)))
testthat::expect_true(all(ic$varDF == c(5, 5, 4, 4, 5)))
testthat::expect_true(all(abs(ic$AIC - c(1653.100,1651.294,1654.613,1669.928,1708.997)) < 1e-01))
