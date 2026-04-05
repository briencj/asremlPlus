# Extracted from test42WheatSpatialVignette.r:119

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("model_selection")
if (Sys.getenv("NOT_CRAN") == "true") require(asreml)
library(asremlPlus)
cat("#### Test for wheat76 spatial example with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(qqplotr)
data(Wheat.dat)
asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
tmp.dat <- within(Wheat.dat, 
                    {
                      cColumn <- dae::as.numfac(Column)
                      cColumn <- cColumn  - mean(unique(cColumn))
                      cRow <- dae::as.numfac(Row)
                      cRow <- cRow - mean(unique(cRow))
                    })
current.asr <- do.call(asreml, 
                         list(yield ~ Rep + WithinColPairs + Variety, 
                              random = ~ Row + Column,
                              residual = ~ Row:Column,
                              data=tmp.dat, maxit = 10))
summary(current.asr)$varcomp
info <- infoCriteria(current.asr, IClikelihood = "full")
testthat::expect_equal(info$varDF, 3)
testthat::expect_lt(abs(info$AIC - 1720.891), 0.10)
current.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                              label = "Initial model")
testthat::expect_lt(abs(current.asrt$test.summary$AIC - 1720.891), 0.50)
current.asrt <- rmboundary(current.asrt, IClikelihood = "full")
current.asrt <- changeModelOnIC(current.asrt, dropFixed = "WithinColPairs", 
                                  label = "Try dropping withinColPairs", IClikelihood = "full")
print(current.asrt)
corb.asrt <- addSpatialModelOnIC(current.asrt, spatial.model = "corr", 
                                   row.covar = "cRow", col.covar = "cColumn", 
                                   row.factor = "Row", col.factor = "Column", 
                                   corr.funcs = c("corb", "corb"), corr.orders = c(0,0),
                                   nugget.variance = TRUE, allow.corrsJointFit = TRUE, 
                                   IClikelihood = "full")
corb.asrt <- rmboundary(corb.asrt, IClikelihood = "full")
inf <- infoCriteria(corb.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(inf$varDF, 6)
testthat::expect_true(abs(inf$AIC - 1666.329) < 0.1)
spatialEach.asrts <- list()
spatialEach.asrts[["corr"]] <- addSpatialModelOnIC(current.asrt, spatial.model = "corr", 
                                                     row.covar = "cRow", col.covar = "cColumn", 
                                                     row.factor = "Row", col.factor = "Column", 
                                                     IClikelihood = "full")
spatialEach.asrts[["corr"]] <- rmboundary(spatialEach.asrts[["corr"]], IClikelihood = "full")
spatialEach.asrts[["TPNCSS"]] <- addSpatialModelOnIC(current.asrt, spatial.model = "TPNCSS", 
                                                       row.covar = "cRow", col.covar = "cColumn", 
                                                       row.factor = "Row", col.factor = "Column", 
                                                       dropRandom = "Row + Column",
                                                       IClikelihood = "full")
spatialEach.asrts[["TPNCSS"]] <- rmboundary(spatialEach.asrts[["TPNCSS"]], IClikelihood = "full")
spatialEach.asrts[["TPPSC2"]] <- addSpatialModelOnIC(current.asrt, spatial.model = "TPPS", 
                                                      row.covar = "cRow", col.covar = "cColumn", 
                                                      row.factor = "Row", col.factor = "Column", 
                                                      dropRandom = "Row + Column",
                                                      degree = c(3,3), difforder = c(2,2), 
                                                      rotateX = TRUE, ngridangles = NULL, 
                                                      asreml.option = "grp", 
                                                      IClikelihood = "full")
spatialEach.asrts[["TPPSC2"]] <- rmboundary(spatialEach.asrts[["TPPSC2"]], IClikelihood = "full")
spatialEach.asrts[["TPPSL1"]] <- addSpatialModelOnIC(current.asrt, spatial.model = "TPPS", 
                                                      row.covar = "cRow", col.covar = "cColumn", 
                                                      row.factor = "Row", col.factor = "Column", 
                                                      dropRandom = "Row + Column",
                                                      degree = c(1,1), difforder = c(1,1),
                                                      asreml.option = "grp", 
                                                      IClikelihood = "full")
spatialEach.asrts[["TPPSL1"]] <- rmboundary(spatialEach.asrts[["TPPSL1"]], IClikelihood = "full")
infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
(infoEach)
spatial.asrts <- chooseSpatialModelOnIC(current.asrt, 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          dropRandom = "Row + Column",
                                          rotateX = TRUE, ngridangles = NULL, 
                                          asreml.option = "mbf", return.asrts = "all")
print(spatial.asrts$spatial.IC)
print(spatial.asrts$asrts$TPNCSS)
testthat::expect_equal(length(spatial.asrts$asrts), 4)
testthat::expect_equal(spatial.asrts$spatial.IC$varDF, c(3,5,6,7,3))
testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1718.609, 1651.314, 1639.489, 1642.838, 1710.225) ) < 1e-02))
