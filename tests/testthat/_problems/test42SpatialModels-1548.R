# Extracted from test42SpatialModels.r:1548

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("spatial_modelling")
cat("#### Test for makeTPPSplineMats with both wheat datasets with asreml42\n")
cat("#### Test makeTPPSplineMats with chick pea example with asreml42\n")
cat("#### Test for wheat76 spatial models with asreml42\n")
cat("#### Test for wheat76 spatial models using mbf with asreml42\n")
cat("#### Test for wheat703 corr spatial models with asreml42\n")
cat("#### Test for wheat76 corr spatial models with asreml42\n")
cat("#### Test for PSA_NW corb spatial models with asreml42\n")
cat("#### Test for PSA_NW with fixed correlations with asreml42\n")
cat("#### Test for spatial models with asreml42\n")
cat("#### Test for nonfitting spatial models with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
data("gw.dat")
asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
gw.dat <- within(gw.dat, 
                   {
                     cRow <- as.numfac(Row, center = TRUE)
                     cCol <- as.numfac(Column, center = TRUE)
                   })
current.asr <- do.call(asreml, 
                         args = list(y ~ Species:Substrate:Irrigation + cRow +cCol, 
                                     data = gw.dat, maxit = 50))
init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Row and Column trends")
init.asrt <- rmboundary(init.asrt)
spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = "none")
testthat::expect_true(all(names(spatial.asrts) == 
                              c("asrts","spatial.IC","best.spatial.mod","best.spatial.IC")))
testthat::expect_equal(names(spatial.asrts$asrts), "nonspatial")
testthat::expect_equal(spatial.asrts$best.spatial.mod, "nonspatial")
testthat::expect_true(abs(spatial.asrts$best.spatial.IC - 892.861) < 1e-04)
testthat::expect_true(abs(spatial.asrts$spatial.IC$AIC - 892.861) < 1e-04)
spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("TPN", "TPPSC"), 
                                          row.covar = "cRow", col.covar = "cCol",
                                          row.factor = "Row", col.factor = "Column", 
                                          dropRandom = c("Row + Column"),
                                          asreml.option = "grp", return.asrts = "all")
testthat::expect_equal(length(spatial.asrts$asrts), 2)
testthat::expect_equal(names(spatial.asrts$asrts), c("TPNCSS", "TPPSC2"))
testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "TPNCSS", "TPPSC2")))
testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(892.861, 892.861, 892.861)) < 0.10))
spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr", "TPN", "TPPSC"), 
                                          row.covar = "cRow", col.covar = "cCol",
                                          row.factor = "Row", col.factor = "Column", 
                                          dropRandom = c("Row + Column"),
                                          asreml.option = "grp", return.asrts = "all")
testthat::expect_equal(length(spatial.asrts$asrts), 3)
testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPSC2"))
testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "corr", "TPNCSS", "TPPSC2")))
testthat::expect_true(all(abs(na.omit(spatial.asrts$spatial.IC$AIC) - 
                                  c(892.861, 887.718, 892.861, 892.861)) < 0.10))
testthat::expect_equal(spatial.asrts$best.spatial.mod, "corr")
testthat::expect_true(all(spatial.asrts$asrts$corr$asreml.obj$vparameters.con == c("F","U","P")))
