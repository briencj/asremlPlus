#devtools::test("asremlPlus")
context("spatial_modelling")


cat("#### Test for wheat76 spatial models with asreml42\n")
test_that("Wheat_spatial_models_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(Wheat.dat)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
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
  
  
  # Check for and remove any boundary terms
  init.asrt <- rmboundary(init.asrt, IClikelihood = "full")
  testthat::expect_lt(abs(init.asrt$test.summary$AIC - 1720.891), 0.50)

  # Try call with illegal argument
  testthat::expect_error(
    current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                        row.covar = "cRow", col.covar = "cColumn",
                                        dropRowterm = "Row", dropColterm = "Column",
                                        nsect = 2,
                                        asreml.option = "grp"), 
    regexp = "the argument\\(s\\) nsect are not legal arguments for 'changeModelOnIC', 'asreml'")
  
  # Try TPPS model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1643.467), 0.10)
  
  #Repeat to make sure no carry-over effects 
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1643.467), 0.10)
  
  #Rotate the penalty matrix
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column", 
                                      rotateX = TRUE, ngridangles = 9, 
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1641.697), 0.10)
  testthat::expect_equal(rownames(summary(current.asrt$asreml.obj)$varcomp), 
                         c("grp(TP.C.2_frow)", "dev(cRow)", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!us(2)_1:1", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!us(2)_2:1", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!us(2)_2:2", 
                           "grp(TP_fcol_frow)", "units!R"))
  testthat::expect_equal(attr(current.asrt, which = "theta.opt"), c(10,90))
  

  
  #Test makeTPPSplineMats with grp
  tps <- makeTPPSplineMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn", 
                             asreml.option = "grp")
  testthat::expect_true(all(names(tps[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                 "BcrZ.df","dim","trace","grp","data.plus")))
  testthat::expect_true(all(names(tps[[1]]$data.plus[,1:19]) == 
                              c("cRow","cColumn","Rep","Row","Column", 
                                "WithinColPairs","Variety","yield","TP.col","TP.row",
                                "TP.CxR","TP.C.1","TP.C.2","TP.R.1","TP.R.2", 
                                "TP.CR.1","TP.CR.2","TP.CR.3","TP.CR.4")))
  testthat::expect_true(all(grepl("TP\\.",names(tps[[1]]$data.plus[,20:50]))))
  testthat::expect_true(all(grepl("TP\\_",names(tps[[1]]$data.plus)[81:ncol(tps[[1]]$data.plus)])))
  
  #Test trapping of illegal nsect argument
  testthat::expect_error(
    tps <- makeTPPSplineMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn", nsect = 2, 
                               asreml.option = "grp"),
    regexp = "the argument\\(s\\) nsect are not legal arguments for 'tpsmmb'")
  
  
  # Try TPPS model using mbf - does not fit the same model as grp because the model is unconverged
  testthat::expect_silent(
    tps <- makeTPPSplineMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn"))
  testthat::expect_true(all(names(tps[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                 "BcrZ.df","dim","trace","data.plus")))
  testthat::expect_true(all(names(tps[[1]]$data.plus) == 
                              c("cRow","cColumn","Rep","Row","Column", 
                                "WithinColPairs","Variety","yield","TP.col","TP.row",
                                "TP.CxR","TP.C.1","TP.C.2","TP.R.1","TP.R.2", 
                                "TP.CR.1","TP.CR.2","TP.CR.3","TP.CR.4")))
  
  testthat::expect_error(
    current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                        row.covar = "cRow", col.covar = "cColumn",
                                        dropRowterm = "Row", dropColterm = "Column",
                                        asreml.option = "mbf", tpps4mbf.obj = tps, 
                                        update = FALSE), 
    regexp = "Can't find BcZxx.df")
  # info <- infoCriteria(current.asrt$asreml.obj)
  #  testthat::expect_equal(info$varDF, 6)
  #  testthat::expect_lt(abs(info$AIC - 1331.759), 0.10)
  
  # Try TPNCSS model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPNCSS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1639.792), 0.10)
  
  # Try corr model - at the moment this fit fails because the addition of a unit terms 
  # clashes with the addition of a variance term for cRow:exp(cColumn)
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column", 
                                      IClikelihood = "full")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1653.096), 0.10)
  
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn", 
                                      row.factor = "Row", col.factor = "Column", 
                                      row.corrFitfirst = FALSE,
                                      IClikelihood = "full")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1653.096), 0.10)

  #Return all of the models
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 4)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS", "TPP1LS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == 
                              c("nonspatial", "corr", "TPNCSS", "TPPCS", "TPP1LS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1720.891, 1653.094, 1639.792, 1644.007, 1710.282)) < 0.10))
  testthat::expect_equal(spatial.asrts$best.spatial.mod, "TPNCSS")
  
  #Fit two models and return both
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 2)
  testthat::expect_equal(names(spatial.asrts$asrts), c("TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - c(1720.891, 1639.792, 1644.007)) < 0.10))
  
  #Fit all models with Row and Column random and return all
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr", "TPN", "TPPC", "TPP1"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 4)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS", "TPP1LS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == 
                              c("nonspatial", "corr", "TPNCSS", "TPPCS", "TPP1LS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1720.891, 1653.094, 1639.792, 1644.007, 1710.282)) < 0.10))
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                                 row.covar = "cRow", col.covar = "cColumn",
                                                 row.factor = "Row", col.factor = "Column")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   dropRowterm = "Row", dropColterm = "Column")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "cRow", col.covar = "cColumn",
                                                  dropRowterm = "Row", dropColterm = "Column",
                                                  asreml.option = "grp")
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   dropRowterm = "Row", dropColterm = "Column",
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
  infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
  testthat::expect_true(all.equal(spatial.asrts$spatial.IC[-1,], infoEach[1:4,-3], 
                                  tolerance = 0.5))
  testthat::expect_true(abs(infoEach$AIC[rownames(infoEach) == 'TPP1LS'] - 1710.282) < 1e-03)
  
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
                                      dropRowterm = "Row", dropColterm = "Column",
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
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1644.007), 0.10)
  #Check Row and COlumn terms not in model
  facs <- c("Row", "Column")
  testthat::expect_false(any(facs %in% rownames(current.asrt$wald.tab)) &&
                           any(facs %in% names(current.asrt$asreml.obj$vparameters)))
  
  #Fit all models with Row and Column fixed and return all
  #NB chooseSpatialModel returns that spatialICs for the best fitting correlation model
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = "all", 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 4)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS", "TPP1LS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == 
                              c("nonspatial", "corr", "TPNCSS", "TPPCS", "TPP1LS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1690.964, 1653.978, 1639.792, 1644.007, 1653.111)) < 0.10))
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                                     row.covar = "cRow", col.covar = "cColumn",
                                                     row.factor = "Row", col.factor = "Column")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   dropRowterm = "Row", dropColterm = "Column")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "cRow", col.covar = "cColumn",
                                                  dropRowterm = "Row", dropColterm = "Column",
                                                  asreml.option = "grp")
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS",
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   dropRowterm = "Row", dropColterm = "Column",
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
  infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
  testthat::expect_true(all.equal(spatial.asrts$spatial.IC[2:5,], infoEach[ ,-3], 
                                  tolerance = 0.5))
})

cat("#### Test for wheat76 corr spatial models with asreml42\n")
test_that("Wheat_corr_models_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(Wheat.dat)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
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
  
  # Try Row:ar1(Column) model
  current.asrt <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  row.factor = "Row", col.factor = "Column",
                                  corr.funcs = c("", "ar1"))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 4)
  testthat::expect_lt(abs(info$AIC - 1669.928), 0.10)
  
  # Try exp(cRow):Column model
  current.asrt <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  row.factor = "Row", col.factor = "Column",
                                  corr.funcs = c("exp", ""))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1714.379), 0.10)
  
  # Try exp(cRow):ar1(Column) model
  current.asrt <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  row.factor = "Row", col.factor = "Column",
                                  corr.funcs = c("exp", "ar1"))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1714.379), 0.10)
  testthat::expect_equal(names(current.asrt$asreml.obj$vparameters), 
                         c("Row", "Column", "Column:cRow", "Column:cRow!cRow!pow", "units!R"))
  
  #Compare lvr and TPP1LS models
  spatial.asrts <- list()
  spatial.asrts[["lvr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                            row.covar = "cRow", col.covar = "cColumn",
                                            row.factor = "Row", col.factor = "Column",
                                            corr.funcs = c("lvr", "lvr"))
  spatial.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS",
                                               row.covar = "cRow", col.covar = "cColumn",
                                               dropRowterm = "Row", dropColterm = "Column",
                                               degree = c(1,1), difforder = c(1,1),
                                               asreml.option = "grp")
  infoEach <- do.call(rbind, 
                      lapply(spatial.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
  
  testthat::expect_true(all.equal(infoEach$AIC, c(1714.861, 1653.111), tolerance = 1e-05))

  #Check trap for all id 
  testthat::expect_error(
    current.asrt <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                    row.covar = "cRow", col.covar = "cColumn",
                                    row.factor = "Row", col.factor = "Column",
                                    corr.funcs = c("idv", "")), 
    regexp = "Both correlation functions are id or equivalent")
  
  # Try id(Row):ar1(Column) model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      corr.funcs = c("id", "ar1"))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1667.54), 0.10)
  testthat::expect_equal(names(current.asrt$asreml.obj$vparameters), 
                         c("Column", "Row:Column", "Row:Column!Column!cor", "units!R"))
  tests <- current.asrt$test.summary
  testthat::expect_equal(nrow(tests), 4)
  testthat::expect_false(all(grepl("Try row", tests$terms)))
  
  # Try ar1(Row):id(Column) model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      corr.funcs = c("ar1", "id"))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1709.001), 0.10)
  testthat::expect_equal(names(current.asrt$asreml.obj$vparameters), 
                         c("Row", "Column", "Column:Row", "Column:Row!Row!cor", "units!R"))
  
  #Fit a correlation model with id to check spatial.IC
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr",  "TPPC"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column", 
                                          corr.funcs = c("id", "ar1"),
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 2)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPPCS"))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1720.891, 1667.540, 1644.007)) < 0.10))
})  

cat("#### Test for barley03 spatial models with asreml42\n")
test_that("barely_spatial_models_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  #This data is for the Durban 2003 barley data cited in Piepho, Boer and Williams (2022)
  data("barley.dat")
  
  #Fit initial model - Row and column random
  current.asr <- do.call(asreml, 
                         list(yield ~ rep + gen, 
                              random = ~ row + col,
                              data=barley.dat))
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - -484.1135), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random row and col effects")
  init.asrt <- rmboundary(init.asrt)
  
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                                 row.covar = "crow", col.covar = "ccol",
                                                 row.factor = "row", col.factor = "col")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "crow", col.covar = "ccol",
                                                   dropRowterm = "row", dropColterm = "col")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "crow", col.covar = "ccol",
                                                  dropRowterm = "row", dropColterm = "col",
                                                  asreml.option = "grp")
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                   row.covar = "crow", col.covar = "ccol",
                                                   dropRowterm = "row", dropColterm = "col",
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
})

cat("#### Test for nonfitting spatial models with asreml42\n")
test_that("nonfit_spatial_models_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data("gw.dat")  
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
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
  
  #Test for trySpatial = "none"
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = "none")
  testthat::expect_true(all(names(spatial.asrts) == 
                              c("asrts","spatial.IC","best.spatial.mod","best.spatial.IC")))
  testthat::expect_equal(names(spatial.asrts$asrts), "nonspatial")
  testthat::expect_equal(spatial.asrts$best.spatial.mod, "nonspatial")
  testthat::expect_true(abs(spatial.asrts$best.spatial.IC - 892.861) < 1e-04)
  testthat::expect_true(abs(spatial.asrts$spatial.IC$AIC - 892.861) < 1e-04)
  
  #Fit two models and return both - neither fits
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cCol",
                                          row.factor = "Row", col.factor = "Column", 
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 2)
  testthat::expect_equal(names(spatial.asrts$asrts), c("TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(892.861, 897.436, 899.239)) < 0.10))
  
  #Fit all models and return all - corr fits
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr", "TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cCol",
                                          row.factor = "Row", col.factor = "Column", 
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 3)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "corr", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(na.omit(spatial.asrts$spatial.IC$AIC) - 
                                  c(892.861, 887.718, 897.436, 899.239)) < 0.10))
  testthat::expect_equal(spatial.asrts$best.spatial.mod, "corr")
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                                     row.covar = "cRow", col.covar = "cCol",
                                                     row.factor = "Row", col.factor = "Column")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "cRow", col.covar = "cCol")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "cRow", col.covar = "cCol",
                                                  dropRowterm = "Row", dropColterm = "Column",
                                                  asreml.option = "grp")
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS",
                                                   row.covar = "cRow", col.covar = "cCol",
                                                   dropRowterm = "Row", dropColterm = "Column",
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
  infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, , IClikelihood = "full")))
  testthat::expect_true(all.equal(spatial.asrts$spatial.IC[2:4,], infoEach[-4,-3], tolerance = 1e-01))
})

cat("#### Test spatial modelling for chick pea example with asreml42\n")
test_that("chickpea_spatial_mod_asreml42", {
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
                         list(fixed = Biomass.plant ~ Smarthouse + Lines * TRT, 
                              random = ~Smarthouse:Zone/Mainplot, 
                              data = tmp.dat))
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Lane and Position effects")
  init.asrt <- rmboundary(init.asrt)
  
  #Test makeTPPSplineMats with sections and grp
  tps <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                             row.covar = "vLanes", col.covar = "vMPosn",
                             asreml.option = "grp")
  testthat::expect_true(all(names(tps) == c("SW","SE")))
  testthat::expect_true(all(names(tps[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                 "BcrZ.df","dim","trace","grp","data.plus")))
  testthat::expect_true(all(names(tps[[1]]$data.plus[,1:19]) == 
                              c("Smarthouse", "vLanes","vMPosn","Lane","Position","Zone","vPos",
                                "Mainplot","Subplot","Lines","TRT","Rep",
                                "X100.SW","Biomass.plant","Pods.plant","Filled.pods.plant", 
                                "Empty.pods.plant","Seed.No.plant","Seed.weight.plant")))
  testthat::expect_true(all(grepl("TP\\.",names(tps[[1]]$data.plus[,20:100]))))
  testthat::expect_true(all(grepl("TP\\_",names(tps[[1]]$data.plus)[101:ncol(tps[[1]]$data.plus)])))
  testthat::expect_equal(tps[[1]]$grp$TP.C.1_frow[1], tps[[1]]$grp$All[1])
  testthat::expect_equal(length(tps[[1]]$grp$All), 334)
  
  # Try TPPS model with Mainplots and two Smarthouses
  TPPS.Main.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                            sections = "Smarthouse", 
                                            row.covar = "vLanes", col.covar = "vMPosn",
                                            dropRowterm = "Lane", dropColterm = NULL,
                                            asreml.option = "grp")
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.Main.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 4013.592)) < 0.10))
  
  # Try TPPS model with Lanes x Positions and two Smarthouses
  TPPS.LP.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                          sections = "Smarthouse", 
                                          row.covar = "vLanes", col.covar = "vPos",
                                          dropRowterm = NULL, dropColterm = NULL,
                                          asreml.option = "grp")
  
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.LP.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 3999.176)) < 0.10))
  
  # Try TPPS model with rotation for Mainplots and two Smarthouses
  TPPSRot.Main.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                               sections = "Smarthouse", 
                                               row.covar = "vLanes", col.covar = "vMPosn",
                                               dropRowterm = "Lane", dropColterm = NULL,
                                               rotateX = TRUE, ngridangles = 3,
                                               asreml.option = "grp")
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPSRot.Main.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,7)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 3997.087)) < 0.10))
  theta.opt <- attr(TPPSRot.Main.grp.asrt$asreml.obj, which = "theta.opt")
  testthat::expect_true(all(theta.opt$SW == 0))
  testthat::expect_true(all(theta.opt$SE == c(30,60)))
  
  # Try TPPS model with rotation for Lanes x Positions and two Smarthouses
  TPPS.LP.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                          sections = "Smarthouse", 
                                          row.covar = "vLanes", col.covar = "vPos",
                                          dropRowterm = NULL, dropColterm = NULL,
                                          rotateX = TRUE, ngridangles = 3,
                                          asreml.option = "grp")
  
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.LP.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 3998.683)) < 0.10))
  
  #Test makeTPPSplineMats with sections with Lanes x Positions  and mbf
  tps.mbf <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                               row.covar = "vLanes", col.covar = "vPos",
                               asreml.option = "mbf")
  testthat::expect_true(all(names(tps.mbf) == c("SW","SE")))
  testthat::expect_true(all(names(tps.mbf[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                     "BcrZ.df","dim","trace","data.plus")))
  testthat::expect_true(all(names(tps.mbf[[1]]$data.plus[,1:19]) == 
                              c("Smarthouse", "vLanes","vPos","Lane","Position","Zone",
                                "Mainplot","Subplot","Lines","TRT","Rep",
                                "X100.SW","Biomass.plant","Pods.plant","Filled.pods.plant", 
                                "Empty.pods.plant","Seed.No.plant","Seed.weight.plant","vMPosn")))
  testthat::expect_true(all(grepl("TP\\.",names(tps.mbf[[1]]$data.plus[,20:30]))))
  testthat::expect_equal(nrow(tps.mbf[[1]]$data.plus), 1056)
  testthat::expect_equal(ncol(tps.mbf[[1]]$data.plus), 30)
  testthat::expect_equal(sapply(list(BcZSE.df,BrZSE.df,BcrZSE.df), nrow), c(22, 24, 528))
  testthat::expect_equal(sapply(list(BcZSE.df,BrZSE.df,BcrZSE.df), ncol), c(4, 4, 10))
  
  # Try TPPS model with Lanes x Positions and two Smarthouses and mbf
  TPPS.LP.mbf.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                          sections = "Smarthouse", 
                                          row.covar = "vLanes", col.covar = "vPos",
                                          dropRowterm = NULL, dropColterm = NULL,
                                          asreml.option = "mbf", tpps4mbf.obj = tps.mbf)
  
  print(info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.LP.mbf.asrt$asreml.obj), 
                       IClikelihood = "full"))
  testthat::expect_true(all(info$varDF == c(3,7)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 4088.381 )) < 0.10))
  
 
  #Test makeTPPSplineMats with sections with Lanes x MainPosns  and mbf
  tps.mbf <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                               row.covar = "vLanes", col.covar = "vMPosn",
                               asreml.option = "mbf")
  testthat::expect_true(all(names(tps.mbf) == c("SW","SE")))
  testthat::expect_true(all(names(tps.mbf[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                     "BcrZ.df","dim","trace","data.plus")))
  testthat::expect_true(all(names(tps.mbf[[1]]$data.plus[,1:19]) == 
                              c("Smarthouse", "vLanes","vMPosn","Lane","Position","Zone",
                                "vPos", "Mainplot","Subplot","Lines","TRT","Rep",
                                "X100.SW","Biomass.plant","Pods.plant","Filled.pods.plant", 
                                "Empty.pods.plant","Seed.No.plant","Seed.weight.plant")))
  testthat::expect_true(all(grepl("TP\\.",names(tps.mbf[[1]]$data.plus[,20:30]))))
  testthat::expect_equal(nrow(tps.mbf[[1]]$data.plus), 1056)
  testthat::expect_equal(ncol(tps.mbf[[1]]$data.plus), 30)
  testthat::expect_equal(sapply(list(BcZSE.df,BrZSE.df,BcrZSE.df), nrow), c(11, 24, 264))
  testthat::expect_equal(sapply(list(BcZSE.df,BrZSE.df,BcrZSE.df), ncol), c(4, 4, 10))
  
  # Try TPPS model with Lanes x Positions and two Smarthouses
  TPPS.Main.mbf.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                            sections = "Smarthouse", 
                                            row.covar = "vLanes", col.covar = "vMPosn",
                                            dropRowterm = NULL, dropColterm = NULL,
                                            asreml.option = "mbf", tpps4mbf.obj = tps.mbf)
  
  print(info <- infoCriteria(list(split = init.asrt$asreml.obj, 
                                  TPPS = TPPS.Main.mbf.asrt$asreml.obj), 
                             IClikelihood = "full"))
  testthat::expect_true(all(info$varDF == c(3,8)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 4090.482)) < 0.10))
})

cat("#### Test hetero variances for HEB25 with asreml42\n")
test_that("HEB25_heterovar_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  #Re-arrange and re-order the data.frame and add factors 
  data(cart.dat)
  tmp.dat <- within(cart.dat, 
                    { 
                      Smarthouse.Treat <- fac.combine(list(Smarthouse, Treatment.1))
                      Lanes <- factor(Lanes)
                      xPosition <- dae::as.numfac(Positions)
                      xPosition <- xPosition - mean(unique(xPosition))
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
  
  
  #Fit an initial model that includes the random term us(Treatment):Genotype
  asreml.options(keep.order = TRUE) #required for asreml4 only
  HEB25.asr <- do.call(asreml, 
                       list(fixed = Dry.Weight ~ Smarthouse + Check + Treatment.1 + 
                              Check:Treatment.1, 
                            random = ~ us(Treatment.1):Genotype.ID + Smarthouse:Zones:Mainplots, 
                            residual = ~idh(Treat.Smarthouse):Zones:Mainplots, 
                            data = tmp.dat, na.action=na.method(y="include", x="include"), 
                            maxit = 100, trace = FALSE))
  summ <- summary(HEB25.asr)$varcomp
  testthat::expect_equal(nrow(summ), 9)
  testthat::expect_equal(summ$bound, c("P","P","P","P","F","P","P","P","P"))
  
  HEB25.idh.asrt <- as.asrtests(HEB25.asr, NULL, NULL, label = "Nonspatial model", 
                                IClikelihood = "full")
  suppressWarnings(
    testthat::expect_true(all(abs(infoCriteria(HEB25.idh.asrt$asreml.obj)[c("AIC","BIC")] - 
                                    c(537.3956, 576.9867)) < 1e-03)))
  #print(HEB25.idh.asrt)
  
  #Test spatial models on Lanes x MainPosn
  #Check makeTPPSplineMats - must be ordered for Smarthouse then Treatment.1
  tpsLM.mat <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                                 row.covar = "cLane", col.covar = "cMainPosn",
                                 asreml.option = "grp")
  testthat::expect_equal(names(tpsLM.mat), c("NW", "NE"))
  testthat::expect_equal(c(nrow(tpsLM.mat$NW$data), nrow(tpsLM.mat$NE$data)), c(264, 264))
  testthat::expect_equal(c(ncol(tpsLM.mat$NW$data), ncol(tpsLM.mat$NE$data)), c(348, 348))
  testthat::expect_equal(c(nrow(tpsLM.mat$NW$data.plus), nrow(tpsLM.mat$NE$data.plus)), c(1056, 1056))
  testthat::expect_equal(c(ncol(tpsLM.mat$NW$data.plus), ncol(tpsLM.mat$NE$data.plus)), c(401, 401))
  #Check that order of dat.plus is the same as in the original tmp.dat
  testthat::expect_equal(tpsLM.mat$NW$data.plus$Snapshot.ID.Tag, tmp.dat$Snapshot.ID.Tag)
  
  #Choose best model for L x M spatial variation
  HEB25.spatialLM.asrts <- chooseSpatialModelOnIC(HEB25.idh.asrt, 
                                                  sections = "Smarthouse", 
                                                  row.covar = "cLane", col.covar = "cMainPosn",
                                                  row.factor = "Lanes", col.factor = "MainPosn",
                                                  asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25.spatialLM.asrts$spatial.IC$AIC - 
                                  c(524.1956, 483.8668, 482.3137, 475.9464, 479.9971) < 1e-03)))
  testthat::expect_equal(names(HEB25.spatialLM.asrts$asrts), c("corr",  "TPNCSS", "TPPCS",  "TPP1LS"))
  summ <- summary(HEB25.spatialLM.asrts$asrts$TPPCS$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ), 19)
  testthat::expect_true(all((summ$bound[-15] == "P")))
  testthat::expect_true(all((summ$bound[15] == "F")))
  summ <- summary(HEB25.spatialLM.asrts$asrts$corr$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ), 15)
  testthat::expect_equal(summ$bound, c("P","U","U","P","U","U","P",
                                       "P","P","P","F","P","P","P","P"))
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModelOnIC(HEB25.idh.asrt, spatial.model = "corr", 
                                                     sections = "Smarthouse", 
                                                     row.covar = "cLane", col.covar = "cMainPosn",
                                                     row.factor = "Lanes", col.factor = "MainPosn")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(HEB25.idh.asrt, spatial.model = "TPN", 
                                                   sections = "Smarthouse", 
                                                   row.covar = "cLane", col.covar = "cMainPosn")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(HEB25.idh.asrt, spatial.model = "TPPS", 
                                                  sections = "Smarthouse", 
                                                  row.covar = "cLane", col.covar = "cMainPosn",
                                                  asreml.option = "grp")
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(HEB25.idh.asrt, spatial.model = "TPPS",
                                                   sections = "Smarthouse", 
                                                   row.covar = "cLane", col.covar = "cMainPosn",
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
  infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
  testthat::expect_true(all.equal(HEB25.spatialLM.asrts$spatial.IC[2:5,], infoEach[ ,-3], 
                                  tolerance = 1e-05))
  
  #Test spatial models on Lanes x Positions
  #Check makeTPPSplineMats - must be ordered for Smarthouse then Treatment.1
  tpsLP.mat <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                                 row.covar = "cLane", col.covar = "cPosition",
                                 asreml.option = "grp")
  testthat::expect_equal(names(tpsLP.mat), c("NW", "NE"))
  testthat::expect_equal(c(nrow(tpsLP.mat$NW$data), nrow(tpsLP.mat$NE$data)), c(528, 528))
  testthat::expect_equal(c(ncol(tpsLP.mat$NW$data), ncol(tpsLP.mat$NE$data)), c(634, 634))
  testthat::expect_equal(c(nrow(tpsLP.mat$NW$data.plus), nrow(tpsLP.mat$NE$data.plus)), c(1056, 1056))
  testthat::expect_equal(c(ncol(tpsLP.mat$NW$data.plus), ncol(tpsLP.mat$NE$data.plus)), c(687, 687))
  
  #Choose best model for L x P spatial variation
  HEB25.spatialLP.asrts <- chooseSpatialModelOnIC(HEB25.idh.asrt,  
                                                  sections = "Smarthouse", 
                                                  row.covar = "cLane", col.covar = "cPosition",
                                                  row.factor = "Lanes", col.factor = "Positions",
                                                  asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25.spatialLP.asrts$spatial.IC$AIC - 
                                  c(524.1825, 480.9052, 477.0927, 471.3009, 476.9187) < 0.1)))
  testthat::expect_equal(names(HEB25.spatialLP.asrts$asrts), c("corr",  "TPNCSS", "TPPCS",  "TPP1LS"))
  summ <- summary(HEB25.spatialLP.asrts$asrts$TPPCS$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ), 18)
  testthat::expect_true(all((summ$bound[-14] == "P")))
  testthat::expect_true(all((summ$bound[14] == "F")))
  summ <- summary(HEB25.spatialLP.asrts$asrts$corr$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ), 15)
  testthat::expect_equal(summ$bound, c("P","P","U","U","P","U","U",
                                       "P","P","P","F","P","P","P","P"))
  
  #Return two P-spline models with rotation  for L x P spatial variation
  HEB25Rot.spatialLP.asrts <- chooseSpatialModelOnIC(HEB25.idh.asrt,  trySpatial = c("TPPCS", "TPP1LS"),
                                                     sections = "Smarthouse", 
                                                     row.covar = "cLane", col.covar = "cPosition",
                                                     row.factor = "Lanes", col.factor = "Positions",
                                                     rotateX = TRUE, ngridangles = 2,
                                                     asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25Rot.spatialLP.asrts$spatial.IC$AIC - 
                                  c(524.1957, 466.9671, 476.9318) < 0.1)))
  testthat::expect_equal(names(HEB25Rot.spatialLP.asrts$asrts), c("TPPCS",  "TPP1LS"))
  summ <- summary(HEB25Rot.spatialLP.asrts$asrts$TPPCS$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ), 15)
  testthat::expect_true(all((summ$bound[-c(7,11)] == "P")))
  testthat::expect_true(all((summ$bound[7] == " ")))
  testthat::expect_true(all((summ$bound[11] == "F")))
  summ <- summary(HEB25Rot.spatialLP.asrts$asrts$TPP1LS$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ), 17)
  testthat::expect_true(all((summ$bound[-13] == "P")))
  testthat::expect_true(all((summ$bound[13] == "F")))
  theta.opt <- attr(HEB25Rot.spatialLP.asrts$asrts$TPPCS$asreml.obj, which = "theta.opt")
  testthat::expect_true(all(theta.opt$NW == c(45,90)))
  testthat::expect_true(all(theta.opt$NE == c(45,45)))
  
  #Test dsum 
  HEB25.asr <- do.call(asreml,
                       list(fixed = Dry.Weight ~ Smarthouse + Check + Treatment.1 + 
                              Check:Treatment.1, 
                            random = ~ us(Treatment.1):Genotype.ID + Smarthouse:Zones:Mainplots, 
                            residual = ~ dsum(~ Zones:Mainplots | Treat.Smarthouse), 
                            data = tmp.dat, na.action=na.method(y="include", x="include"), 
                            maxit = 100, trace = FALSE))
  
  HEB25.ds.asrt <- as.asrtests(HEB25.asr, NULL, NULL, label = "Nonspatial model", 
                               IClikelihood = "full")
  suppressWarnings(
    testthat::expect_true(all(abs(infoCriteria(HEB25.ds.asrt$asreml.obj)[c("AIC","BIC")] - 
                                    c(537.3956, 576.9867)) < 1e-03)))
  summ.idh <- summary(HEB25.idh.asrt$asreml.obj)$varcomp
  summ.ds <- summary(HEB25.ds.asrt$asreml.obj)$varcomp
  #Check that varcomp is the same for idh and dsum
  testthat::expect_true(all.equal(summ.idh[-5, ], summ.ds, tolerance = 1e-03, 
                                  check.attributes = FALSE))
  #print(HEB25.ds.asrt)
  
  #Choose best model for L x M spatial variation when variance specified using dsum
  HEB25.spatialLM.ds.asrts <- chooseSpatialModelOnIC(HEB25.ds.asrt, 
                                                     sections = "Smarthouse", 
                                                     row.covar = "cLane", col.covar = "cMainPosn",
                                                     row.factor = "Lanes", col.factor = "MainPosn",
                                                     asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25.spatialLM.ds.asrts$spatial.IC$AIC - 
                                  c(524.1956, 488.4492, 482.3137, 475.9464, 479.9971) < 1e-03)))
  testthat::expect_equal(names(HEB25.spatialLM.ds.asrts$asrts), c("corr",  "TPNCSS", "TPPCS",  "TPP1LS"))
  #Check TPPCS
  summ.idh <- summary(HEB25.spatialLM.asrts$asrts$TPPCS$asreml.obj)$varcomp
  summ.ds <- summary(HEB25.spatialLM.ds.asrts$asrts$TPPCS$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ.ds), 18)
  testthat::expect_true(all((summ.idh$bound[-15] == "P")))
  testthat::expect_true(all((summ.idh$bound[15] == "F")))
  testthat::expect_equal(rownames(summ.idh)[1:14], rownames(summ.ds)[1:14])
  testthat::expect_true(all.equal(summ.idh[-15,-5], summ.ds[,-5], tolerance = 1e-02, 
                                  check.attributes = FALSE))
  #Check corr
  summ.idh <- summary(HEB25.spatialLM.asrts$asrts$corr$asreml.obj)$varcomp
  summ.ds <- summary(HEB25.spatialLM.ds.asrts$asrts$corr$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ.ds), 12)
  testthat::expect_equal(summ.ds$bound, c("P","U","U","P","U","P",
                                          "P","P","P","P","P","P"))
  #Two components are missing from dsum
  testthat::expect_equal(rownames(summ.idh)[c(1:5,8:10)], rownames(summ.ds)[1:8])
  # testthat::expect_true(all.equal(summ.idh[-c(6,7,11), 1:2], summ.ds[, 1:2], tolerance = 0.1, 
  #                                 check.attributes = FALSE))
  
  #Choose best model for L x M spatial variation when variance specified using dsum
  HEB25.spatialLP.ds.asrts <- chooseSpatialModelOnIC(HEB25.ds.asrt,  
                                                     sections = "Smarthouse", 
                                                     row.covar = "cLane", col.covar = "cPosition",
                                                     row.factor = "Lanes", col.factor = "Positions",
                                                     asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25.spatialLP.ds.asrts$spatial.IC$AIC - 
                                  c(524.1825, 480.9052, 477.0927, 471.3009, 476.9187) < 0.1)))
  testthat::expect_equal(names(HEB25.spatialLP.ds.asrts$asrts), c("corr",  "TPNCSS", "TPPCS",  "TPP1LS"))
  #Check TPPCS
  summ.idh <- summary(HEB25.spatialLP.asrts$asrts$TPPCS$asreml.obj)$varcomp
  summ.ds <- summary(HEB25.spatialLP.ds.asrts$asrts$TPPCS$asreml.obj)$varcomp
  testthat::expect_equal(rownames(summ.idh)[1:13], rownames(summ.ds)[1:13])
  testthat::expect_true(all.equal(summ.idh[-14,], summ.ds, tolerance = 1e-03, 
                                  check.attributes = FALSE))
  #Check corr
  summ.idh <- summary(HEB25.spatialLP.asrts$asrts$corr$asreml.obj)$varcomp
  summ.ds <- summary(HEB25.spatialLP.ds.asrts$asrts$corr$asreml.obj)$varcomp
  testthat::expect_equal(rownames(summ.idh)[1:10], rownames(summ.ds)[1:10])
  testthat::expect_true(all.equal(summ.idh[-11,], summ.ds, tolerance = 1e-03, 
                                  check.attributes = FALSE))
})
