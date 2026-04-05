# Extracted from test42SpatialModels.r:1477

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

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
data("barley.dat")
current.asr <- do.call(asreml, 
                         list(yield ~ rep + gen, 
                              random = ~ row + col,
                              data=barley.dat, maxit = 50))
info <- infoCriteria(current.asr, IClikelihood = "full")
testthat::expect_equal(info$varDF, 3)
testthat::expect_lt(abs(info$AIC - -484.1135), 0.10)
init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random row and col effects")
init.asrt <- rmboundary(init.asrt)
spatialEach.asrts <- list()
spatialEach.asrts[["corr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                                 row.covar = "crow", col.covar = "ccol",
                                                 row.factor = "row", col.factor = "col")
spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "crow", col.covar = "ccol",
                                                   dropRandom = c("row + col"))
spatialEach.asrts[["TPPSC2"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "crow", col.covar = "ccol",
                                                  dropRandom = c("row + col"),
                                                  asreml.option = "grp")
spatialEach.asrts[["TPPSL1"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                   row.covar = "crow", col.covar = "ccol",
                                                   dropRandom = c("row + col"),
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
infoEach <- lapply(spatialEach.asrts, function(asrt) infoCriteria(asrt$asreml.obj, 
                                                                    IClikelihood = "full"))
(infoEach <- do.call(rbind, infoEach))
testthat::expect_true(all.equal(infoEach$AIC, c(-641.2598, -611.8811, -616.8260, -646.7571), 
                                  tolerance = 1e-02))
infoEach <- lapply(spatialEach.asrts, function(asrt) infoCriteria(asrt$asreml.obj, 
                                                                    IClikelihood = "REML"))
(infoEach <- do.call(rbind, infoEach))
testthat::expect_true(all.equal(infoEach$AIC, c(-230.4462, -191.8063, -226.5424, -230.1942), 
                                  tolerance = 1e-02))
R2adj <- sapply(spatialEach.asrts, function(asrt) R2adj(asrt$asreml.obj, 
                                                          include.which.fixed = ~ gen,
                                                          orthogonalize = "eigen"))
testthat::expect_true(all(abs(R2adj - c(33.42414, 18.13152, 16.2514, 29.30381)) < 0.001))
R2adj.corr <- R2adj(spatialEach.asrts[["corr"]]$asreml.obj, 
                      include.which.fixed = ~ ., include.which.random = ~ .)
testthat::expect_true(abs(R2adj.corr - 79.27872) < 1e-03)
R2adj.TPNCSS <- R2adj(spatialEach.asrts[["TPNCSS"]]$asreml.obj, 
                        include.which.fixed = ~ ., include.which.random = ~ .)
testthat::expect_true(abs(R2adj.TPNCSS - 85.47652) < 1e-03)
