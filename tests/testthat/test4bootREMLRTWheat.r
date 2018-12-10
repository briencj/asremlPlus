#devtools::test("asremlPlus")
context("model_selection")

cat("#### Test for boot using wheat example with asreml4\n")
test_that("Wheatboot_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(Wheat.dat)
  #'## Add cubic trend to Row so that spline is not bound
  Wheat.dat <- within(Wheat.dat, 
                      {
                        vRow <- as.numeric(Row)
                        vRow <- vRow - mean(unique(vRow))
                        yield <- yield + 10*vRow + 10 * (vRow^2) + 5 * (vRow^3)
                      })
  
  # Fit initial model
  asreml.options(design = TRUE)
  current.asr <- do.call("asreml", 
                         list(yield ~ Rep + WithinColPairs + Variety, 
                              random = ~ spl(vRow) + Row + Column + units,
                              residual = ~ ar1(Row):ar1(Column), 
                              maxit=100, data=Wheat.dat))
  summary(current.asr)$varcomp
  info <- infoCriteria(current.asr)
  testthat::expect_equal(info$DF, 5)
  testthat::expect_lt(abs(info$AIC - 1357.118), 1e-03)
  
  # Load current fit into an asrtests object
  current.asrt <- asrtests(current.asr, NULL, NULL)
  
  # Check for and remove any boundary terms
  current.asrt <- rmboundary(current.asrt)
  
  #Check term for within Column pairs
  current.asrt <- testranfix(current.asrt, "WithinColPairs", drop.fix.ns=TRUE)
  
  # Test nugget term
  current.asrt <- testranfix(current.asrt, "units", positive=TRUE)
  current.asrt$test.summary
  
  # Use bootstrap to test units
  full.asreml.obj <- current.asrt$asreml.obj
  reduced.asreml.obj <- update(full.asreml.obj, random = ~ . - units)
  units.lrt <- REMLRT(h0.asreml.obj = reduced.asreml.obj, 
                      h1.asreml.obj = full.asreml.obj)
  boot.units <- bootREMLRT(h0.asreml.obj = reduced.asreml.obj, 
                           h1.asreml.obj = full.asreml.obj, 
                           nboot = 100, seed = 6250,
                           fixed.spline.terms = "spl(vRow)")
  testthat::expect_equal(length(boot.units), 6)
  testthat::expect_equal(boot.units$DF, 1)
  testthat::expect_equal(length(boot.units$REMLRT.sim) + 
                           length(attr(boot.units$REMLRT.sim, which = "na.action")), 100)
  testthat::expect_equal(length(boot.units$nunconverged), 100)
  testthat::expect_lt(abs(boot.units$REMLRT - 10.9), 1e-01)
  testthat::expect_lt(abs(boot.units$REMLRT - units.lrt$REMLRT), 1e-01)
  
  # Use bootstrap to test Row autocorrelation
  full.asreml.obj <- current.asrt$asreml.obj
  reduced.asreml.obj <- update(full.asreml.obj, residual. = ~ Row:ar1(Column))
  arR.lrt <- REMLRT(h0.asreml.obj = reduced.asreml.obj, 
                    h1.asreml.obj = full.asreml.obj)
  REMLRT(h0.asreml.obj = reduced.asreml.obj, 
         h1.asreml.obj = full.asreml.obj, bound.exclusions = "B")
  boot.units <- bootREMLRT(h0.asreml.obj = reduced.asreml.obj, 
                           h1.asreml.obj = full.asreml.obj, 
                           nboot = 100, seed = 6250,
                           fixed.spline.terms = "spl(vRow)")
  testthat::expect_equal(length(boot.units), 6)
  testthat::expect_equal(boot.units$DF, 2)
  testthat::expect_equal(length(boot.units$REMLRT.sim) + 
                           length(attr(boot.units$REMLRT.sim, which = "na.action")), 100)
  testthat::expect_equal(length(boot.units$nunconverged), 100)
  testthat::expect_lt(abs(boot.units$REMLRT - 29.8), 1e-01)
  testthat::expect_lt(abs(boot.units$REMLRT - arR.lrt$REMLRT), 1e-01)

  # Use bootstrap to test Col autocorrelation
  reduced.asreml.obj <- update(full.asreml.obj, residual. = ~ ar1(Row):Column)
  arC.lrt <- REMLRT(h0.asreml.obj = reduced.asreml.obj, 
                    h1.asreml.obj = full.asreml.obj)
  boot.units <- bootREMLRT(h0.asreml.obj = reduced.asreml.obj, 
                           h1.asreml.obj = full.asreml.obj, 
                           nboot = 100, seed = 6250,
                           fixed.spline.terms = "spl(vRow)")
  testthat::expect_equal(length(boot.units), 6)
  testthat::expect_equal(boot.units$DF, 1)
  testthat::expect_equal(length(boot.units$REMLRT.sim) + 
                           length(attr(boot.units$REMLRT.sim, which = "na.action")), 100)
  testthat::expect_equal(length(boot.units$nunconverged), 100)
  testthat::expect_lt(abs(boot.units$REMLRT - 56.1), 1e-01)
  testthat::expect_lt(abs(boot.units$REMLRT - arC.lrt$REMLRT), 1e-01)
  
  asreml.options(design = FALSE) 
})
