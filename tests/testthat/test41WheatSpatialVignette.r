#devtools::test("asremlPlus")
context("model_selection")
asr41.lib <- "D:\\Analyses\\R ASReml4.1" 

cat("#### Test for wheat76 spatial example with asreml41\n")
test_that("Wheat_spatial_asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(qqplotr)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(Wheat.dat)
  
  #Add row and column covariates
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
                              data=tmp.dat))
  summary(current.asr)$varcomp
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1720.891), 0.10)

  # Load init fit into an asrtests object
  current.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                              label = "Initial model")
  testthat::expect_lt(abs(current.asrt$test.summary$AIC - 1720.891), 0.50)


  # Check for and remove any boundary terms
  current.asrt <- rmboundary(current.asrt, IClikelihood = "full")


  #Check term for within Column pairs
  current.asrt <- changeModelOnIC(current.asrt, dropFixed = "WithinColPairs", 
                                  label = "Try dropping withinColPairs", IClikelihood = "full")
  print(current.asrt)

  #Fit autocorrelation model
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModelOnIC(current.asrt, spatial.model = "corr", 
                                                     row.covar = "cRow", col.covar = "cColumn", 
                                                     row.factor = "Row", col.factor = "Column", 
                                                     allow.fixedcorrelation = TRUE,
                                                     checkboundaryonly = TRUE, IClikelihood = "full")
  spatialEach.asrts[["corr"]] <- rmboundary(spatialEach.asrts[["corr"]], IClikelihood = "full")

  spatialEach.asrts[["TPNCSS"]] <- addSpatialModelOnIC(current.asrt, spatial.model = "TPNCSS", 
                                                       row.covar = "cRow", col.covar = "cColumn", 
                                                       row.factor = "Row", col.factor = "Column", 
                                                       allow.fixedcorrelation = TRUE,
                                                       checkboundaryonly = TRUE, IClikelihood = "full")
  spatialEach.asrts[["TPNCSS"]] <- rmboundary(spatialEach.asrts[["TPNCSS"]], IClikelihood = "full")
  
  spatialEach.asrts[["TPPCS"]] <- addSpatialModelOnIC(current.asrt, spatial.model = "TPPS", 
                                                      row.covar = "cRow", col.covar = "cColumn", 
                                                      row.factor = "Row", col.factor = "Column", 
                                                      degree = c(3,3), difforder = c(2,2),
                                                      asreml.option = "grp", allow.fixedcorrelation = TRUE,
                                                      checkboundaryonly = TRUE, IClikelihood = "full")
  spatialEach.asrts[["TPPCS"]] <- rmboundary(spatialEach.asrts[["TPPCS"]], IClikelihood = "full")
  
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModelOnIC(current.asrt, spatial.model = "TPPS", 
                                                      row.covar = "cRow", col.covar = "cColumn", 
                                                      row.factor = "Row", col.factor = "Column", 
                                                      degree = c(1,1), difforder = c(1,1),
                                                      asreml.option = "grp", allow.fixedcorrelation = TRUE,
                                                      checkboundaryonly = TRUE, IClikelihood = "full")
  spatialEach.asrts[["TPP1LS"]] <- rmboundary(spatialEach.asrts[["TPP1LS"]], IClikelihood = "full")
  
  infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
  (infoEach)
  
  #Choose  spatial model
  spatial.asrts <- chooseSpatialModelOnIC(current.asrt, 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  
  print(spatial.asrts$spatial.IC)
  print(spatial.asrts$asrts$TPNCSS)
  testthat::expect_equal(length(spatial.asrts$asrts), 4)
  testthat::expect_equal(spatial.asrts$spatial.IC$varDF, c(3,5,6,6,3))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1718.609, 1651.314, 1639.489, 1644.190, 1708.443) ) < 1e-03))
  testthat::expect_true(all.equal(spatial.asrts$spatial.IC[2:4,], infoEach[1:3 ,-3], 
                                  tolerance = 1e-05))
  
  current.asr <- spatial.asrts$asrts$TPNCSS$asreml.obj
  printFormulae(current.asr)
  
  ## Get current fitted asreml object and update to include standardized residuals
  
  current.asr <- update(current.asr, aom=TRUE)
  Wheat.dat$res <- residuals(current.asr, type = "stdCond")
  Wheat.dat$fit <- fitted(current.asr)

  ## Do diagnostic checking
  
  ### Do residuals-versus-fitted values plot
  
  with(Wheat.dat, plot(fit, res))

  ### Plot variofaces
  
  variofaces(current.asr, V=NULL, units="addtores", 
             maxit=50, update = FALSE)

  ### Plot normal quantile plot
  
  ggplot(data = Wheat.dat, mapping = aes(sample = res)) +
    qqplotr::stat_qq_band(bandType = "ts") + 
    qqplotr::stat_qq_line() + 
    qqplotr::stat_qq_point() +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles",
         title = "Normal probability plot") +
    theme(plot.title = element_text(size = 12, face = "bold")) + theme_bw()

  ## Get Variety predictions and all pairwise prediction differences and p-values
  Var.diffs <- predictPlus(classify = "Variety", 
                           asreml.obj=current.asr, 
                           error.intervals="halfLeast",
                           wald.tab=current.asrt$wald.tab, 
                           sortFactor = "Variety",
                           tables = "predictions")
  
  ## Plot the Variety predictions, with halfLSD intervals, and the p-values
  
  plotPredictions(Var.diffs$predictions, 
                  classify = "Variety", y = "predicted.value", 
                  error.intervals = "half")
  plotPvalues(Var.diffs)
})
