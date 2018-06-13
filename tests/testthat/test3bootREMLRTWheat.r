#devtools::test("asremlPlus")
context("model_selection3")

#' cat("#### Test for boot using wheat example with asreml3\n")
#' test_that("Wheatboot_asreml3", {
#'   skip_if_not_installed("asreml")
#'   skip_on_cran()
#'   library(dae)
#'   library(asreml)
#'   library(asremlPlus)
#'   ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
#'   data(Wheat.dat)
#'   #'## Add cubic trend to Row so that spline is not bound
#'   Wheat.dat <- within(Wheat.dat, 
#'                       {
#'                         vRow <- as.numeric(Row)
#'                         vRow <- vRow - mean(unique(vRow))
#'                         yield <- yield + 10*vRow + 10 * (vRow^2) + 5 * (vRow^3)
#'                       })
#'   
#'   # Fit initial model
#'   current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
#'                         random = ~ spl(vRow) + Row + Column + units,
#'                         rcov = ~ ar1(Row):ar1(Column), 
#'                         data=Wheat.dat)
#'   summary(current.asr)
#'   info <- infoCriteria(current.asr)
#'   testthat::expect_equal(info$DF, 5)
#'   testthat::expect_lt(abs(info$AIC - 1357.65), 1e-02)
#'   
#'   # Load current fit into an asrtests object
#'   current.asrt <- asrtests(current.asr, NULL, NULL)
#'   
#'   # Check for and remove any boundary terms
#'   current.asrt <- rmboundary(current.asrt)
#'   
#'   #Check term for within Column pairs
#'   current.asrt <- testranfix(current.asrt, "WithinColPairs", drop.fix.ns=TRUE)
#'   
#'   # Test nugget term
#'   current.asrt <- testranfix(current.asrt, "units", positive=TRUE)
#'   
#'   # Test Row autocorrelation
#'   current.asrt <- testresidual(current.asrt, "~ Row:ar1(Column)", 
#'                                label="Row autocorrelation", simpler=TRUE)
#'   
#'   # Test Col autocorrelation (depends on whether Row autocorrelation retained)
#'   k <- match("Row autocorrelation", current.asrt$test.summary$terms)
#'   p <- current.asrt$test.summary$p
#'   { if (p[k] <= 0.05)
#'     current.asrt <- testresidual(current.asrt, "~ ar1(Row):Column", 
#'                                  label="Col autocorrelation", simpler=TRUE,
#'                                  update=FALSE)
#'     else
#'       current.asrt <- testresidual(current.asrt, "~ Row:Column", 
#'                                    label="Col autocorrelation", simpler=TRUE,
#'                                    update=FALSE)
#'   }
#'   print(current.asrt)
#'   testthat::expect_equal(length(current.asrt), 3)
#'   testthat::expect_equal(nrow(current.asrt$wald.tab), 3)
#'   testthat::expect_equal(nrow(current.asrt$test.summary), 6)
#'   info <- infoCriteria(current.asrt$asreml.obj)
#'   testthat::expect_equal(info$DF, 5)
#'   testthat::expect_lt(abs(info$AIC - 1364.143), 1e-03)
#'   
#'   # Get current fitted asreml object
#'   current.asr <- current.asrt$asreml.obj
#'   current.asr <- update(current.asr, aom=TRUE)
#'   
#'   
#'   # Form variance matrix based on estimated variance parameters
#'   s2 <- current.asr$sigma2
#'   gamma.Row <- current.asr$gammas[1]
#'   gamma.unit <- current.asr$gammas[2]
#'   rho.r <- current.asr$gammas[4]
#'   rho.c <- current.asr$gammas[5]
#'   row.ar1 <- mat.ar1(order=10, rho=rho.r)
#'   col.ar1 <- mat.ar1(order=15, rho=rho.c)
#'   V <- fac.vcmat(Wheat.dat$Row, gamma.Row) + 
#'     gamma.unit * diag(1, nrow=150, ncol=150) + 
#'     mat.dirprod(row.ar1, col.ar1)
#'   testthat::expect_equal(nrow(V), 150)
#'   testthat::expect_equal(ncol(V), 150)
#'   V <- s2*V
#'   
#'   # Use bootstrap to test units
#'   full.asreml.obj <- current.asrt$asreml.obj
#'   reduced.asreml.obj <- update(full.asreml.obj, random = ~ . - units)
#'   expect_error(boot.units <- bootREMLRT(h1 = full.asreml.obj, h0 = reduced.asreml.obj, 
#'                            nboot = 100, seed = 6250,
#'                            V = V, fixed.spline.terms = "spl(vRow)"))
#' 
#' })

