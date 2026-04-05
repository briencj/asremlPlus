# Extracted from test42SpatialModels.r:1874

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
cat("#### Test hetero variances for HEB25 with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
data(cart.dat)
tmp.dat <- within(cart.dat, 
                    { 
                      Smarthouse.Treat <- fac.combine(list(Smarthouse, Treatment.1))
                      Lanes <- factor(Lanes)
                      xPosition <- dae::as.numfac(Positions, center = TRUE)
                    })
tmp.dat <- tmp.dat[c("Snapshot.ID.Tag", "Smarthouses", "Lanes", "Positions", 
                       "Genotype.ID", "Lines.nos", "Check", "Treatment.1", "Conditions", 
                       "Smarthouse", "Treat.Smarthouse", "Smarthouse.Treat", 
                       "Zones", "Rows", "Mainplots", "Subplots", 
                       "xLane", "xPosition", "xMainPosn", "MainCol", 
                       "Fresh.Weight", "Dry.Weight", "Number.Tillers.correct" , 
                       "Plant.Length", "ratio", "Water_Amount", "SSA", "ASA", 
                       "Caliper.Length", "Convex.Hull.Area", "Height", "SCR", 
                       "WUE.2", "agrOST", "rgrOST", "linOST", "logOST", "linm4OST", 
                       "logm4OST", "linm5OST", "logm5OST", 
                       "agrm4OST", "rgrm4OST", "agrm5OST", "rgrm5OST", 
                       "WUE100", "CHA10000", "ASA10000", "SSA10000", 
                       "ShootArea_sm", "AGR_sm_32_42", "RGR_sm_32_42", 
                       "AGR_sm_42_50", "RGR_sm_42_50", "AGR_sm_50_59", "RGR_sm_50_59")]
names(tmp.dat)[match(c("xLane", "xPosition", "xMainPosn", "MainCol"), names(tmp.dat))] <- 
    c("cLane", "cPosition", "cMainPosn", "MainPosn")
tmp.dat <- with(tmp.dat, tmp.dat[order(Treat.Smarthouse, Zones, Mainplots), ])
asreml.options(keep.order = TRUE)
HEB25.asr <- do.call(asreml, 
                       list(fixed = Dry.Weight ~ Smarthouse + Check + Treatment.1 + 
                              Check:Treatment.1, 
                            random = ~ us(Treatment.1):Genotype.ID + 
                              (at(Smarthouse, 'NW') + at(Smarthouse, 'NE')):Zones:Mainplots,  #nugget terms
                            residual = ~idh(Treat.Smarthouse):Zones:Mainplots, 
                            data = tmp.dat, na.action=na.method(y="include", x="include"), 
                            maxit = 100, trace = FALSE))
summ <- summary(HEB25.asr)$varcomp
testthat::expect_equal(nrow(summ), 10)
testthat::expect_equal(summ$bound, c("P","P","P","P","P","F","P","P","P","P"))
HEB25.idh.asrt <- as.asrtests(HEB25.asr, NULL, NULL, label = "Nonspatial model", 
                                IClikelihood = "full")
suppressWarnings(
    testthat::expect_true(all(abs(infoCriteria(HEB25.idh.asrt$asreml.obj)[c("AIC","BIC")] - 
                                    c(539.034, 583.5741)) < 1e-03)))
tpsLM.mat <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                                 row.covar = "cLane", col.covar = "cMainPosn",
                                 asreml.option = "grp")
testthat::expect_equal(names(tpsLM.mat), c("NW", "NE"))
testthat::expect_equal(c(nrow(tpsLM.mat$NW$data), nrow(tpsLM.mat$NE$data)), c(264, 264))
testthat::expect_equal(c(ncol(tpsLM.mat$NW$data), ncol(tpsLM.mat$NE$data)), c(348, 348))
testthat::expect_equal(c(nrow(tpsLM.mat$NW$data.plus), nrow(tpsLM.mat$NE$data.plus)), 
                         c(1056, 1056))
testthat::expect_equal(c(ncol(tpsLM.mat$NW$data.plus), ncol(tpsLM.mat$NE$data.plus)), 
                         c(401, 401))
testthat::expect_equal(tpsLM.mat$NW$data.plus$Snapshot.ID.Tag, tmp.dat$Snapshot.ID.Tag)
HEB25.spatialLM.asrts <- 
    chooseSpatialModelOnIC(HEB25.idh.asrt, 
                           sections = "Smarthouse", 
                           row.covar = "cLane", col.covar = "cMainPosn",
                           row.factor = "Lanes", col.factor = "MainPosn",
                           allow.fixedcorrelation = FALSE, 
                           asreml.option = "grp", return.asrts = "all")
testthat::expect_true(all(abs(HEB25.spatialLM.asrts$spatial.IC$AIC - 
                                  c(525.5955, 489.5259, 470.8116, 473.2412, 479.7448) < 1e-03)))
testthat::expect_equal(names(HEB25.spatialLM.asrts$asrts), 
                         c("corr",  "TPNCSS", "TPPSC2",  "TPPSL1"))
summ <- summary(HEB25.spatialLM.asrts$asrts$TPPSC2$asreml.obj)$varcomp
summ$bound[summ$bound == " "] <- "P"
testthat::expect_equal(nrow(summ), 18)
testthat::expect_true(all((summ$bound[-14] == "P")))
testthat::expect_true(all((summ$bound[14] == "F")))
summ <- summary(HEB25.spatialLM.asrts$asrts$corr$asreml.obj)$varcomp
testthat::expect_equal(nrow(summ), 13)
