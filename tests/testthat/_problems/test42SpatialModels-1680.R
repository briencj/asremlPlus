# Extracted from test42SpatialModels.r:1680

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
cat("#### Test spatial modelling for chick pea example with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
data(chkpeadat)
tmp.dat <- within(chkpeadat, 
                    {
                      vMPosn <- as.numfac(fac.recast(Mainplot, newlevels = rep(1:11, times = 4)))
                      vMPosn <- vMPosn - mean(unique(vMPosn))
                    })
asreml.options(design = TRUE)
current.asr <- do.call(asreml, 
                         list(fixed  = Biomass.plant ~ Smarthouse + Lines * TRT, 
                              random = ~ (at(Smarthouse, "SW") + at(Smarthouse, "SE")):Zone + 
                                (at(Smarthouse, "SW") + at(Smarthouse, "SE")):Zone:Mainplot, #nugget terms
                              data = tmp.dat, maxit = 50))
init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Zone effects")
init.asrt <- rmboundary(init.asrt)
TPPS.Main.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                            sections = "Smarthouse", 
                                            row.covar = "vLanes", col.covar = "vMPosn",
                                            dropRandom = c('at(Smarthouse, "SW"):Zone', 
                                                           'at(Smarthouse, "SE"):Zone'),
                                            asreml.option = "grp")
info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.Main.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
testthat::expect_true(all(info$varDF == c(5,12)))
testthat::expect_true(all(abs(info$AIC - c(4263.948, 4004.857)) < 0.10))
testthat::expect_false(all(c("at(Smarthouse, 'SW'):Zone", "at(Smarthouse, 'SE'):Zone") %in% 
                               names(TPPS.Main.grp.asrt$asreml.obj$vparameters)))
TPPS.LP.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                          sections = "Smarthouse", 
                                          row.covar = "vLanes", col.covar = "vPos",
                                          dropFixed = NULL, dropRandom = NULL, 
                                          asreml.option = "grp")
info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.LP.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
testthat::expect_true(all(info$varDF == c(5,11)))
testthat::expect_true(all(abs(info$AIC - c(4263.948, 3995.805)) < 0.10))
library(tictoc)
tic.clearlog()
tic("grid search")
TPPSRot.Main.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                               sections = "Smarthouse", 
                                               row.covar = "vLanes", col.covar = "vMPosn",
                                               dropFixed = NULL, dropRandom = NULL, 
                                               rotateX = TRUE, ngridangles = c(3,3),
                                               asreml.option = "grp")
toc(log = TRUE)
info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPSRot.Main.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
testthat::expect_true(all(info$varDF == c(5,9)))
testthat::expect_true(all(abs(info$AIC - c(4263.948, 4000.055)) < 0.10))
theta.opt <- attr(TPPSRot.Main.grp.asrt$asreml.obj, which = "theta.opt")
testthat::expect_true(all(theta.opt$SW == c(60,90)))
testthat::expect_true(all(theta.opt$SE == c(30,30)))
tic("optimize")
TPPSRot.Main.grp.opt.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                                   sections = "Smarthouse", 
                                                   row.covar = "vLanes", col.covar = "vMPosn",
                                                   dropFixed = NULL, dropRandom = NULL, 
                                                   rotateX = TRUE, ngridangles = NULL,
                                                   asreml.option = "grp")
toc(log = TRUE)
info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPSRot.Main.grp.opt.asrt$asreml.obj), 
                       IClikelihood = "full")
testthat::expect_true(all(info$varDF == c(5,9)))
testthat::expect_true(all(abs(info$AIC - c(4263.948, 3984.119)) < 0.10))
theta.opt <- attr(TPPSRot.Main.grp.opt.asrt$asreml.obj, which = "theta.opt")
testthat::expect_true(all(abs(theta.opt$SW - c( 64.022114, 49.92560)) < 0.001))
testthat::expect_true(all(abs(theta.opt$SE - c(34.313518, 15.6094129)) < 0.001))
library(tictoc)
tic.clearlog()
tic("optimize LP")
TPPS.LP.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                          sections = "Smarthouse", 
                                          row.covar = "vLanes", col.covar = "vPos",
                                          dropFixed = NULL, dropRandom = NULL, 
                                          rotateX = TRUE, ngridangles = NULL,
                                          asreml.option = "grp")
toc(log = TRUE)
info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.LP.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
testthat::expect_true(all(info$varDF == c(5,9)))
testthat::expect_true(all(abs(info$AIC - c(4263.948, 3978.829)) < 0.10))
