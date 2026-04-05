# Extracted from test42SpatialModels.r:1987

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
testthat::expect_equal(summ$bound, c("P","U","U","P","U","P","P","P",
                                       "F","P","P","P","P"))
spatialEach.asrts <- list()
spatialEach.asrts[["corr"]] <- 
    addSpatialModelOnIC(HEB25.idh.asrt, spatial.model = "corr", 
                        sections = "Smarthouse", 
                        row.covar = "cLane", col.covar = "cMainPosn",
                        row.factor = "Lanes", col.factor = "MainPosn", 
                        allow.fixedcorrelation = FALSE)
testthat::expect_true(any(grepl("NW", spatialEach.asrts[["corr"]]$test.summary$terms)) & 
                          any(grepl("NE", spatialEach.asrts[["corr"]]$test.summary$terms)))
spatialEach.asrts[["TPNCSS"]] <- 
    addSpatialModel(HEB25.idh.asrt, spatial.model = "TPN", 
                    sections = "Smarthouse", 
                    row.covar = "cLane", col.covar = "cMainPosn")
testthat::expect_true(any(grepl("NW", spatialEach.asrts[["TPNCSS"]]$test.summary$terms)) & 
                          any(grepl("NE", spatialEach.asrts[["TPNCSS"]]$test.summary$terms)))
spatialEach.asrts[["TPPSC2"]] <- 
    addSpatialModel(HEB25.idh.asrt, spatial.model = "TPPS", 
                    sections = "Smarthouse", 
                    row.covar = "cLane", col.covar = "cMainPosn",
                    asreml.option = "grp")
testthat::expect_true(any(grepl("NW", spatialEach.asrts[["TPPSC2"]]$test.summary$terms)) & 
                          any(grepl("NE", spatialEach.asrts[["TPPSC2"]]$test.summary$terms)))
spatialEach.asrts[["TPPSL1"]] <- 
    addSpatialModel(HEB25.idh.asrt, spatial.model = "TPPS",
                    sections = "Smarthouse", 
                    row.covar = "cLane", col.covar = "cMainPosn",
                    degree = c(1,1), difforder = c(1,1),
                    asreml.option = "grp")
testthat::expect_true(any(grepl("NW", spatialEach.asrts[["TPPSL1"]]$test.summary$terms)) & 
                          any(grepl("NE", spatialEach.asrts[["TPPSL1"]]$test.summary$terms)))
infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, 
                                                         IClikelihood = "full")))
testthat::expect_true(all.equal(HEB25.spatialLM.asrts$spatial.IC[c(2,4:5),], 
                                  infoEach[c(1,3:4), -3], 
                                  tolerance = 1e-05))
asreml.obj <- spatialEach.asrts[["corr"]]$asreml.obj
asreml.obj <- update(asreml.obj, aom = TRUE)
G.dat <- lapply(c("at(Smarthouse, 'NW'):Lanes:MainPosn", 
                    "at(Smarthouse, 'NE'):Lanes:MainPosn"), 
                  function(term, asreml.obj, use)
                    t <- convEffectNames2DataFrame.asreml(asreml.obj, term = term, 
                                                          use = use),
                  asreml.obj = asreml.obj, use = "G.aom")
testthat::expect_true(all(lapply(G.dat, nrow) == 264))
G.dat <- do.call(rbind, G.dat)
testthat::expect_true(nrow(G.dat) == 528)
testthat::expect_true(all(sapply(G.dat, is.factor)))
tpsLP.mat <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                                 row.covar = "cLane", col.covar = "cPosition",
                                 asreml.option = "grp")
testthat::expect_equal(names(tpsLP.mat), c("NW", "NE"))
testthat::expect_equal(c(nrow(tpsLP.mat$NW$data), nrow(tpsLP.mat$NE$data)), c(528, 528))
testthat::expect_equal(c(ncol(tpsLP.mat$NW$data), ncol(tpsLP.mat$NE$data)), c(634, 634))
testthat::expect_equal(c(nrow(tpsLP.mat$NW$data.plus), nrow(tpsLP.mat$NE$data.plus)), c(1056, 1056))
testthat::expect_equal(c(ncol(tpsLP.mat$NW$data.plus), ncol(tpsLP.mat$NE$data.plus)), 
                         c(687, 687))
HEB25.spatialLP.asrts <- 
    chooseSpatialModelOnIC(HEB25.idh.asrt,  
                           sections = "Smarthouse", 
                           row.covar = "cLane", col.covar = "cPosition",
                           row.factor = "Lanes", col.factor = "Positions",
                           asreml.option = "grp", return.asrts = "all")
testthat::expect_true(all(abs(HEB25.spatialLP.asrts$spatial.IC$AIC - 
                                  c(525.5955, 513.2121, 471.5088, 472.8215, 476.6325) < 0.1)))
testthat::expect_equal(names(HEB25.spatialLP.asrts$asrts), 
                         c("corr",  "TPNCSS", "TPPSC2",  "TPPSL1"))
summ <- summary(HEB25.spatialLP.asrts$asrts$TPPSC2$asreml.obj)$varcomp
summ$bound[summ$bound == " "] <- "P"
testthat::expect_equal(nrow(summ), 19)
testthat::expect_true(all((summ$bound[-15] == "P")))
testthat::expect_true(all((summ$bound[15] == "F")))
summ <- summary(HEB25.spatialLP.asrts$asrts$corr$asreml.obj)$varcomp
testthat::expect_equal(nrow(summ), 15)
testthat::expect_equal(summ$bound, c("P","P","U","U","P","U","U","P","P","P",
                                       "F","P","B","P","P"))
HEB25Rot.spatialLP.asrts <- 
    chooseSpatialModelOnIC(HEB25.idh.asrt,  trySpatial = c("TPPSC2", "TPPSL1"),
                           sections = "Smarthouse", 
                           row.covar = "cLane", col.covar = "cPosition",
                           row.factor = "Lanes", col.factor = "Positions",
                           allow.fixedcorrelation = FALSE,
                           rotateX = TRUE, ngridangles = NULL,
                           asreml.option = "grp", return.asrts = "all")
testthat::expect_true(all(abs(HEB25Rot.spatialLP.asrts$spatial.IC$AIC - 
                                  c(525.5955, 469.2255, 476.6325) < 0.1)))
testthat::expect_equal(names(HEB25Rot.spatialLP.asrts$asrts), c("TPPSC2",  "TPPSL1"))
summ <- summary(HEB25Rot.spatialLP.asrts$asrts$TPPSC2$asreml.obj)$varcomp
summ$bound[summ$bound == " "] <- "P"
testthat::expect_equal(nrow(summ), 16)
testthat::expect_true(all((summ$bound[-c(12)] == "P")))
testthat::expect_true(all((summ$bound[12] == "F")))
summ <- summary(HEB25Rot.spatialLP.asrts$asrts$TPPSL1$asreml.obj)$varcomp
summ$bound[summ$bound == " "] <- "P"
testthat::expect_equal(nrow(summ), 17)
testthat::expect_true(all((summ$bound[-13] == "P")))
testthat::expect_true(all((summ$bound[13] == "F")))
theta.opt <- attr(HEB25Rot.spatialLP.asrts$asrts$TPPSC2$asreml.obj, which = "theta.opt")
testthat::expect_true(all(abs(theta.opt$NW - c(47.54955, 89.92521)) < 0.001))
