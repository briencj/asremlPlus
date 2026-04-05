# Extracted from test42SpatialModels.r:603

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("spatial_modelling")
cat("#### Test for makeTPPSplineMats with both wheat datasets with asreml42\n")
cat("#### Test makeTPPSplineMats with chick pea example with asreml42\n")
cat("#### Test for wheat76 spatial models with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
data(Wheat.dat)
asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
tmp.dat <- within(Wheat.dat, 
                    {
                      cColumn <- dae::as.numfac(Column, center = TRUE)
                      cRow <- dae::as.numfac(Row, center = TRUE)
                    })
asreml.options(design = TRUE)
current.asr <- do.call(asreml, 
                         list(yield ~ Rep + WithinColPairs + Variety, 
                              random = ~ Row + Column,
                              data=tmp.dat, maxit = 50))
info <- infoCriteria(current.asr, IClikelihood = "full")
testthat::expect_equal(info$varDF, 3)
testthat::expect_lt(abs(info$AIC - 1720.891), 0.10)
init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Row and Column effects")
init.asrt <- rmboundary(init.asrt, IClikelihood = "full")
testthat::expect_lt(abs(init.asrt$test.summary$AIC - 1720.891), 0.50)
testthat::expect_error(
    current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                        row.covar = "cRow", col.covar = "cColumn",
                                        dropRandom = c("Row + Column"), 
                                        nsect = 2,
                                        asreml.option = "grp"), 
    regexp = "the argument\\(s\\) nsect are not legal arguments for 'changeModelOnIC.asrtests', 'asreml'")
grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  dropRandom = c("Row + Column"), 
                                  asreml.option = "grp")
info <- infoCriteria(grp.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 7)
testthat::expect_lt(abs(info$AIC - 1643.467), 0.10)
testthat::expect_equal(rownames(summary(grp.asrt$asreml.obj)$varcomp), 
                         c("grp(TP.C.2_frow)", "dev(cRow)", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)!cor", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)_1", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)_2", 
                           "grp(TP_fcol_frow)", "units!R"))
testthat::expect_equal(rownames(grp.asrt$wald.tab), c("(Intercept)", "Rep", "WithinColPairs", 
                                                        "Variety", 
                                                        "TP.CR.2", "TP.CR.3", "TP.CR.4"))
current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRandom = c("Row + Column"), 
                                      asreml.option = "grp")
info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 7)
testthat::expect_lt(abs(info$AIC - 1643.467), 0.10)
mbf.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  dropRandom = c("Row + Column"), 
                                  asreml.option = "mbf")
info <- infoCriteria(list(grp.asrt$asreml.obj, mbf.asrt$asreml.obj), IClikelihood = "full")
testthat::expect_true(all.equal(info[1,], info[2,], tolerance = 1e-06, check.attributes = FALSE ))
mbf.logl.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                              row.covar = "cRow", col.covar = "cColumn",
                              dropRandom = c("Row + Column"), 
                              rotateX = TRUE, ngridangles = NULL, 
                              which.rotacriterion = "likelihood",
                              asreml.option = "mbf")
info <- infoCriteria(mbf.logl.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 7)
testthat::expect_lt(abs(info$AIC - 1650.335), 0.10)
testthat::expect_lt(abs(info$loglik - -784.1677), 0.10)
testthat::expect_true(all(abs(attr(mbf.logl.asrt$asreml.obj, which = "theta.opt")[[1]] - c(20.19972, 64.98770)) < 0.001))
grp.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                              row.covar = "cRow", col.covar = "cColumn",
                              dropRandom = c("Row + Column"), 
                              rotateX = TRUE, ngridangles = c(9,9), 
                              asreml.option = "grp")
info <- infoCriteria(grp.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 7)
testthat::expect_false(info$AIC < 1643.467)
mbf.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                              row.covar = "cRow", col.covar = "cColumn",
                              dropRandom = c("Row + Column"), 
                              rotateX = TRUE, ngridangles = c(9,9), 
                              asreml.option = "mbf")
info <- infoCriteria(list(grp.asrt$asreml.obj, mbf.asrt$asreml.obj), IClikelihood = "full")
testthat::expect_true(all.equal(info[1,], info[2,], check.attributes = FALSE))
info <- infoCriteria(mbf.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 7)
testthat::expect_false(info$AIC < 1643.467)
testthat::expect_equal(rownames(summary(mbf.asrt$asreml.obj)$varcomp), 
                         c("dev(cRow)", "mbf(TP.row):TP.C.2", 
                           "TP.R.1:mbf(TP.col)+mbf(TP.col):TP.R.2!corh(2)!cor", 
                           "TP.R.1:mbf(TP.col)+mbf(TP.col):TP.R.2!corh(2)_1", 
                           "TP.R.1:mbf(TP.col)+mbf(TP.col):TP.R.2!corh(2)_2", 
                           "mbf(TP.CxR)", "units!R"))
testthat::expect_equal(rownames(mbf.asrt$wald.tab), c("(Intercept)", "Rep", "WithinColPairs", 
                                                        "Variety", 
                                                        "TP.CR.2", "TP.CR.3", "TP.CR.4"))
testthat::expect_true(all(attr(mbf.asrt, which = "theta.opt")[[1]] == c(20,60)))
grp.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                              row.covar = "cRow", col.covar = "cColumn",
                              dropRandom = c("Row + Column"), 
                              rotateX = TRUE, ngridangles = NULL, 
                              asreml.option = "grp")
info <- infoCriteria(grp.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 7)
testthat::expect_lt(abs(info$AIC - 1650.335), 0.10)
testthat::expect_equal(rownames(summary(grp.asrt$asreml.obj)$varcomp), 
                         c("grp(TP.C.2_frow)", "dev(cRow)", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)!cor", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)_1", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)_2", 
                           "grp(TP_fcol_frow)", "units!R"))
testthat::expect_equal(rownames(grp.asrt$wald.tab), c("(Intercept)", "Rep", "WithinColPairs", 
                                                        "Variety", 
                                                        "TP.CR.2", "TP.CR.3", "TP.CR.4"))
testthat::expect_true(all(abs(attr(grp.asrt$asreml.obj, which = "theta.opt")[[1]] - 
                                  c(20.19973, 64.98769)) < 0.01))
mbf.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                              row.covar = "cRow", col.covar = "cColumn",
                              dropRandom = c("Row + Column"), 
                              rotateX = TRUE, ngridangles = NULL, 
                              asreml.option = "mbf")
info <- infoCriteria(list(grp.asrt$asreml.obj, mbf.asrt$asreml.obj), IClikelihood = "full")
info <- infoCriteria(mbf.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 7)
testthat::expect_lt(abs(info$AIC - 1650.335), 0.10)
testthat::expect_equal(rownames(summary(mbf.asrt$asreml.obj)$varcomp), 
                         c("dev(cRow)", "mbf(TP.row):TP.C.2", 
                           "TP.R.1:mbf(TP.col)+mbf(TP.col):TP.R.2!corh(2)!cor", 
                           "TP.R.1:mbf(TP.col)+mbf(TP.col):TP.R.2!corh(2)_1", 
                           "TP.R.1:mbf(TP.col)+mbf(TP.col):TP.R.2!corh(2)_2", 
                           "mbf(TP.CxR)", "units!R"))
testthat::expect_equal(rownames(mbf.asrt$wald.tab), c("(Intercept)", "Rep", "WithinColPairs", 
                                                        "Variety", 
                                                        "TP.CR.2", "TP.CR.3", "TP.CR.4"))
testthat::expect_true(all(abs(attr(mbf.asrt, which = "theta.opt")[[1]] - c(20.1997, 64.9876)) < 0.001))
testthat::expect_true(all(abs(attr(grp.asrt, which = "theta.opt")[[1]] - 
                                  attr(mbf.asrt, which = "theta.opt")[[1]]) < 0.001))
grp.dat <- grp.asrt$asreml.obj$call$data
mbf.dat <- mbf.asrt$asreml.obj$call$data
testthat::expect_true(all.equal(grp.dat[c("TP.col","TP.row","TP.CxR",
                                            "TP.C.1","TP.C.2","TP.R.1","TP.R.2", 
                                            "TP.CR.1","TP.CR.2","TP.CR.3","TP.CR.4")], 
                                  mbf.dat[c("TP.col","TP.row","TP.CxR",
                                            "TP.C.1","TP.C.2","TP.R.1","TP.R.2", 
                                            "TP.CR.1","TP.CR.2","TP.CR.3","TP.CR.4")], 
                                  tolerance = 1e-05))
current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPNCSS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRandom = c("Row + Column"), 
                                      asreml.option = "grp")
info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 6)
testthat::expect_lt(abs(info$AIC - 1639.792), 0.10)
current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "mbf", IClikelihood = "full")
info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 5)
testthat::expect_lt(abs(info$AIC - 1653.096), 0.10)
current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn", 
                                      row.factor = "Row", col.factor = "Column", 
                                      row.corrFitfirst = FALSE,
                                      asreml.option = "mbf", IClikelihood = "full")
info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 5)
testthat::expect_lt(abs(info$AIC - 1653.096), 0.10)
current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn", 
                                      row.factor = "Row", col.factor = "Column", 
                                      dropRandom = "Row + Column", 
                                      difforder = c(1,1), degree = c(1,1),
                                      asreml.option = "mbf", IClikelihood = "full")
info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
testthat::expect_equal(info$varDF, 3)
