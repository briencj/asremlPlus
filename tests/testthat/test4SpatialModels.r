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
  tmp.dat <- within(Wheat.dat, 
                    {
                      cColumn <- dae::as.numfac(Column)
                      cColumn <- cColumn  - mean(unique(cColumn))
                      cRow <- dae::as.numfac(Row)
                      cRow <- cRow - mean(unique(cRow))
                    })
  
  #Fit initial model - Row and column random
  current.asr <- do.call(asreml, 
                         list(yield ~ Rep + WithinColPairs + Variety, 
                              random = ~ Row + Column,
                              data=tmp.dat))
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1720.891), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Row and Column effects")
  init.asrt <- rmboundary(init.asrt)

  # Try TPPS model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1644.007), 0.10)
  
  #Repeat to make sure no carry-over effects non-NULL for factors
  current.asrt <- addSpatialModelOnIC(current.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1644.007), 0.10)
  
  # Try TPPS model using mbf
  tps <- makeTPSPlineXZMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn")
  testthat::expect_error(
    current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                        row.covar = "cRow", col.covar = "cColumn",
                                        row.factor = "Row", col.factor = "Column",
                                        asreml.option = "mbf", tpps4mbf.obj = tps, 
                                        update = FALSE), 
    regexp = 'Sorry, but the mbf setting of asreml.opt is not functioning yet')
#  info <- infoCriteria(current.asrt$asreml.obj)
#  testthat::expect_equal(info$varDF, 6)
#  testthat::expect_lt(abs(info$AIC - 1302.258), 0.10)
  
  # Try TPNCSS model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPNCSS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1639.792), 0.10)
  
  # Try corr model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1719.136), 0.10)
  
  #Choose the best model
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          asreml.option = "grp")
  testthat::expect_equal(length(spatial.asrts$asrts), 1)
  testthat::expect_equal(names(spatial.asrts$asrts), "TPNCSS")
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "corr", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - c(1720.891, 1719.136, 1639.792, 1644.007)) < 0.10))
  
  #Fit two models and return both
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 2)
  testthat::expect_equal(names(spatial.asrts$asrts), c("TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - c(1720.891, 1639.792, 1644.007)) < 0.10))
  
  #Fit all models with Row and Column random and return all
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr", "TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 3)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "corr", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - c(1720.891, 1719.136, 1639.792, 1644.007)) < 0.10))
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                                 row.covar = "cRow", col.covar = "cColumn")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   row.factor = "Row", col.factor = "Column")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "cRow", col.covar = "cColumn",
                                                  row.factor = "Row", col.factor = "Column",
                                                  asreml.option = "grp")
  # spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
  #                                                  row.covar = "cRow", col.covar = "cColumn",
  #                                                  row.factor = "Row", col.factor = "Column",
  #                                                  degree = 1, difforder = 1, 
  #                                                  asreml.option = "grp")
  infoEach <- lapply(spatialEach.asrts, function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full"))
  testthat::expect_true(all.equal(spatial.asrts$spatial.IC[-1,], do.call(rbind, infoEach)[,-3], 
                                  tolerance = 1e-05))
  
  #Fit initial model - Row and column fixed
  current.asr <- do.call(asreml, 
                         list(yield ~ Rep + WithinColPairs + Row + Column + Variety, 
                              data=tmp.dat))
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 1)
  testthat::expect_lt(abs(info$AIC - 1690.964), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Row and Column effects")
  init.asrt <- rmboundary(init.asrt)
  
  # Try a TPNCSS model with fixed Row and Column
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPNCSS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1639.792), 0.10)
  facs <- c("Row", "Column")
  #Check Row and COlumn terms not in model
  testthat::expect_false(any(facs %in% rownames(current.asrt$wald.tab)) &&
                           any(facs %in% names(current.asrt$asreml.obj$vparameters)))
  
  # Try TPPS model with fixed Row and Column
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1644.007), 0.10)
  #Check Row and COlumn terms not in model
  facs <- c("Row", "Column")
  testthat::expect_false(any(facs %in% rownames(current.asrt$wald.tab)) &&
                           any(facs %in% names(current.asrt$asreml.obj$vparameters)))

  #Fit all models with Row and Column fixed and return all
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr", "TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 3)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "corr", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - c(1690.964, 1687.705, 1639.792, 1644.007)) < 0.10))
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                                 row.covar = "cRow", col.covar = "cColumn")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   row.factor = "Row", col.factor = "Column")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "cRow", col.covar = "cColumn",
                                                  row.factor = "Row", col.factor = "Column",
                                                  asreml.option = "grp")
  # spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
  #                                                  row.covar = "cRow", col.covar = "cColumn",
  #                                                  row.factor = "Row", col.factor = "Column",
  #                                                  degree = 1, difforder = 1, 
  #                                                  asreml.option = "grp")
  infoEach <- lapply(spatialEach.asrts, function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full"))
  testthat::expect_true(all.equal(spatial.asrts$spatial.IC[-1,], do.call(rbind, infoEach)[,-3]))
})

cat("#### Test for nonfitting spatial models with asreml4\n")
test_that("nonfit_spatial_models_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data("gw.dat")  
  
  gw.dat <- within(gw.dat, 
                   {
                     cRow <- as.numfac(Row)
                     cRow <- cRow - mean(unique(cRow))
                     cCol <- as.numfac(Column)
                     cCol <- cCol - mean(unique(cCol))
                   })
  
  #Fit initial model
  current.asr <- do.call(asreml, 
                         args = list(y ~ Species:Substrate:Irrigation + cRow +cCol, 
                                     data = gw.dat))
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Row and Column trends")
  init.asrt <- rmboundary(init.asrt)
  
  #Fit two models and return both - neither fits
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cCol",
                                          row.factor = "Row", col.factor = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 2)
  testthat::expect_equal(names(spatial.asrts$asrts), c("TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(892.861, 897.436, 899.239)) < 0.10))
  
  #Fit all models and return all - none fits
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr", "TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cCol",
                                          row.factor = "Row", col.factor = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 3)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "corr", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(na.omit(spatial.asrts$spatial.IC$AIC) - 
                                  c(892.861, 897.436, 899.239)) < 0.10))
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                                 row.covar = "cRow", col.covar = "cCol")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "cRow", col.covar = "cCol")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "cRow", col.covar = "cCol",
                                                  row.factor = "Row", col.factor = "Column",
                                                  asreml.option = "grp")
  # spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
  #                                                  row.covar = "cRow", col.covar = "cCol",
  #                                                  row.factor = "Row", col.factor = "Col",
  #                                                  degree = 1, difforder = 1, 
  #                                                  asreml.option = "grp")
  infoEach <- lapply(spatialEach.asrts, function(asrt) infoCriteria(asrt$asreml.obj, , IClikelihood = "full"))
  testthat::expect_true(all.equal(spatial.asrts$spatial.IC[3:4,], do.call(rbind, infoEach)[-1,-3]))
})

cat("#### Test spatial modelling for chick pea example with asreml4\n")
test_that("chickpea_spatial_mod_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(chkpeadat)
  tmp.dat <- within(chkpeadat, 
                    {
                      vMPosn <- as.numfac(fac.recast(Mainplot, newlevels = rep(1:11, times = 4)))
                      vMPosn <- vMPosn - mean(unique(vMPosn))
                    })
  asreml.options(design = TRUE)
  current.asr <- do.call(asreml, 
                         list(fixed = Biomass.plant ~ Smarthouse + Lines * TRT, 
                              random = ~Smarthouse:Zone/Mainplot, 
                              data = tmp.dat))
  
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
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = current.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 4001.819)) < 0.10))

  # Try TPPS model with Lanes x Positions and two Smarthouses
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      sections = "Smarthouse", 
                                      row.covar = "vLanes", col.covar = "vPos",
                                      row.factor = NULL, col.factor = NULL,
                                      asreml.option = "grp")
  
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = current.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 3999.176)) < 0.10))
  
})
