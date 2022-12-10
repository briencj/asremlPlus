#devtools::test("asremlPlus")
context("spatial_modelling")



cat("#### Test for wheat76 spatial models with asreml4\n")
test_that("Wheat_spatial_models_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(Wheat.dat)
  
  #Add row and column covariates
  Wheat.dat <- within(Wheat.dat, 
                      {
                        cColumn <- dae::as.numfac(Column)
                        cColumn <- cColumn  - mean(unique(cColumn))
                        cRow <- dae::as.numfac(Row)
                        cRow <- cRow - mean(unique(cRow))
                      })
  
  #Fit initial model - Row and column random
  current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                        random = ~ Row + Column,
                        data=Wheat.dat)
  info <- infoCriteria(current.asr)
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1400.719), 0.10)
  
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, 
                           label = "Random Row and Column effects")
  init.asrt <- rmboundary(init.asrt)
  
  # Try TPPS model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj)
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1302.258), 0.10)
  
  #Repeat to make sure no carry-over effects non-NULL for factors
  current.asrt <- addSpatialModelOnIC(current.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj)
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1302.258), 0.10)
  
  # Try TPPS model using mbf
  tps <- makeTPSPlineXZMats(Wheat.dat, row.covar = "cRow", col.covar = "cColumn")
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "mbf", tpps4mbf.obj = tps, 
                                      update = FALSE)
  info <- infoCriteria(current.asrt$asreml.obj)
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1302.258), 0.10)
  
  # Try TPNCSS model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPNCSS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj)
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1329.024), 0.10)
  
  # Try corr model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn")
  info <- infoCriteria(current.asrt$asreml.obj)
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1399.628), 0.10)
  
  #Choose the best model
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          asreml.option = "grp")
  testthat::expect_equal(length(spatial.asrts$asrts), 1)
  testthat::expect_equal(names(spatial.asrts$asrts), "TPPS")
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "corr", "TPNCSS", "TPPS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - c(1400.719, 1399.628, 1329.024, 1302.258)) < 0.10))
  
  #Fit two models and return both
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("TPN", "TPP"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 2)
  testthat::expect_equal(names(spatial.asrts$asrts), c("TPNCSS", "TPPS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "TPNCSS", "TPPS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - c(1400.719, 1329.024, 1302.258)) < 0.10))
  
  #Fit initial model - Row and column fixed
  current.asr <- asreml(yield ~ Rep + WithinColPairs + Row + Column + Variety, 
                        data=Wheat.dat)
  info <- infoCriteria(current.asr)
  testthat::expect_equal(info$varDF, 1)
  testthat::expect_lt(abs(info$AIC - 1191.179), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Row and Column effects")
  init.asrt <- rmboundary(init.asrt)
  
  # Try at TPNCSS model with fixed Row and Column
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPNCSS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj)
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1329.024), 0.10)
  facs <- c("Row", "Column")
  #Check Row and COlumn terms not in model
  testthat::expect_false(any(facs %in% rownames(current.asrt$wald.tab)) &&
                           any(facs %in% names(current.asrt$asreml.obj$vparameters)))
  
  # Try TPPS model with fixed Row and Column
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj)
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1302.258), 0.10)
  #Check Row and COlumn terms not in model
  facs <- c("Row", "Column")
  testthat::expect_false(any(facs %in% rownames(current.asrt$wald.tab)) &&
                           any(facs %in% names(current.asrt$asreml.obj$vparameters)))
  
})

cat("#### Test soatial modelling for chick pea example with asreml4\n")
test_that("chickpea_spatial_mod_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(chkpeadat)
  chkpeadat$vMPosn <- as.numfac(fac.recast(chkpeadat$Mainplot, newlevels = rep(1:11, times = 4)))
  chkpeadat$vMPosn <- with(chkpeadat, vMPosn - mean(unique(vMPosn)))
  asreml.options(design = TRUE)
  current.asr <- asreml(fixed = Biomass.plant ~ Smarthouse + Lines * TRT, 
                        random = ~Smarthouse:Zone/Mainplot, 
                        data = chkpeadat)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Lane and Position effects")
  init.asrt <- rmboundary(init.asrt)
  
  # Try TPPS model with Mainplots and two Smarthouses
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      sections = "Smarthouse", 
                                      row.covar = "vLanes", col.covar = "vMPosn",
                                      row.factor = "Lane", col.factor = NULL,
                                      asreml.option = "grp")
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = current.asrt$asreml.obj))
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(2375.035, 2217.902)) < 0.10))

  # Try TPPS model with Lanes x Positions and two Smarthouses
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      sections = "Smarthouse", 
                                      row.covar = "vLanes", col.covar = "vPos",
                                      row.factor = NULL, col.factor = NULL,
                                      asreml.option = "grp")
  
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = current.asrt$asreml.obj))
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(2375.035, 2213.375)) < 0.10))
  
})