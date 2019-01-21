#devtools::test("asremlPlus")
context("model_selection3")
asr3.lib <- "D:\\Analyses\\R oldpkg" 

cat("#### Test for changeTerms using wheat example with asreml3\n")
test_that("Wheatchange_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(Wheat.dat)
  #'## Add cubic trend to Row so that spline is not bound
  Wheat.dat <- within(Wheat.dat, 
                      {
                        vRow <- as.numeric(Row)
                        vRow <- vRow - mean(unique(vRow))
                        yield <- yield + 10*vRow + 20 * (vRow^2) + 5 * (vRow^3)
                      })
  
  #Function to check with models have been changes
  chk.changes <- function(test.summary)
  {
    changes <- tail(test.summary[test.summary$action != "Boundary",], 1)
    changes <- unlist(lapply(c("fixed", "random", "residual"), 
                             function(x, changes){grepl(x, changes, fixed = TRUE)}, 
                             changes = changes$action))
    return(changes)
  }
  
  # Fit initial model
  
  current.asr <- do.call("asreml", 
                         args = list(yield ~ Rep + WithinColPairs + Variety, 
                                     random = ~ Row + Column,
                                     rcov = ~ ar1(Row):ar1(Column), 
                                     data=Wheat.dat))
  summary(current.asr)$varcomp
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 5)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  
  # Add and drop both fixed and random terms (No units term in model)
  current.asrt <- changeTerms(current.asrt, addFixed = "vRow", dropFixed = "WithinColPairs", 
                              addRandom = "spl(vRow)", dropRandom = "units+Row",
                              checkboundaryonly = TRUE, data = Wheat.dat)
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 2)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 5)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  
  # Drop all fixed terms
  current.asrt <- changeTerms(current.asrt, dropFixed = "Rep + Variety + vRow")
  print(summary(current.asrt$asreml.obj)$varcomp)
  print(current.asrt$wald.tab)
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 1)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 4)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 1)
  
  # Add back fixed terms
  current.asrt <- changeTerms(current.asrt, addFixed = "Rep + Variety + vRow")
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 1)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 4)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  
  # Restart with residual, remove it and then return it
  current.asr <- do.call("asreml", 
                         args = list(yield ~ Rep + WithinColPairs + Variety, 
                                     random = ~ Row + Column,
                                     rcov = ~ ar1(Row):ar1(Column), 
                                     data=Wheat.dat))
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- changeTerms(current.asrt, newResidual = "-(.)")
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 1)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 2)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  current.asrt <- changeTerms(current.asrt, newResidual = "ar1(Row):ar1(Column)")
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 1)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 4)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  
  # Restart with no random and add them back in
  current.asr <- do.call("asreml", 
                         args = list(yield ~ Rep + WithinColPairs + Variety, 
                                     rcov = ~ ar1(Row):ar1(Column), 
                                     data=Wheat.dat))
  summary(current.asr)$varcomp
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- changeTerms(current.asrt, addRandom = "Column")
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 1)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 3)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  
  #Change residual with no random terms  
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- changeTerms(current.asrt, 
                              newResidual = "Row:ar1(Column)", 
                              label="Row autocorrelation")
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 1)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 2)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
})

