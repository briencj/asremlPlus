# Extracted from test42SpatialModels.r:865

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

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
data(indv703.dat)
data(summ703)
for (kresp in responses.test)
  { 
    mod.ch <- paste(kresp, "~ Block + Line")
    mod <- as.formula(mod.ch)
    cat("\n\n#### ",mod.ch,"\n\n")
    asreml.options(keep.order = TRUE)
    current.asr <- do.call(asreml, 
                           args=list(fixed = mod,
                                     random = ~ SubBlock/Block, 
                                     data = indv703.dat, maxiter=50))
    current.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                                label = "Initial model")
    current.asrt <- rmboundary(current.asrt)
    corr.asrt <- chooseSpatialModelOnIC(current.asrt, 
                                        row.covar = "cLane", col.covar = "cPosn", 
                                        row.factor = "Lane", col.factor = "Position", 
                                        trySpatial = "corr")
    ksumm <- summary(corr.asrt$asrts[[1]]$asreml.obj)$varcomp
#    print(all.equal(ksumm, summ[[kresp]], tolerance = 1e-05))
    testthat::expect_true(all.equal(ksumm, summ[[kresp]], tolerance = 1e-05))
#    summ <- c(summ, list(summary(corr.asrt$asrts[[1]]$asreml.obj)$varcomp))
  }
