#devtools::test("asremlPlus")
context("model_selection")
asr41.lib <- "D:\\Analyses\\R ASReml4.1" 

cat("#### Test for changeTerms using wheat example with asreml41\n")
test_that("Wheatchange_asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml, lib.loc = asr41.lib)
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
  
  #Function to check which models have been changed
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
                                     random = ~ Row + Column + units,
                                     residual = ~ ar1(Row):ar1(Column), 
                                     data=Wheat.dat))
  summary(current.asr)$varcomp
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 6)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  
  # Add and drop both fixed and random terms
  current.asrt <- changeTerms(current.asrt, addFixed = "vRow", dropFixed = "WithinColPairs", 
                              addRandom = "spl(vRow)", dropRandom = "units", 
                              checkboundaryonly = TRUE)
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 2)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 6)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  
  # Drop all fixed terms
  current.asrt <- changeTerms(current.asrt, dropFixed = "Rep + Variety + vRow", 
                              checkboundaryonly = TRUE)
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 1)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 6)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 1)

  # Add back fixed terms
  current.asrt <- changeTerms(current.asrt, addFixed = "Rep + Variety + vRow", 
                              checkboundaryonly = TRUE)
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 1)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 6)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  
  # Restart with residual, remove it and then return it
  current.asr <- do.call("asreml", 
                         args = list(yield ~ Rep + WithinColPairs + Variety, 
                                     random = ~ Row + Column + units,
                                     residual = ~ ar1(Row):ar1(Column), 
                                     data=Wheat.dat))
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- changeTerms(current.asrt, dropRandom = "units",newResidual = "-(1)")
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 2)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 1)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)
  current.asrt <- changeTerms(current.asrt, newResidual = "ar1(Row):ar1(Column)")
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 1)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 3)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 4)

  # Restart with no random and add them back in
  current.asr <- do.call("asreml", 
                         args = list(yield ~ Rep + WithinColPairs + Variety, 
                                     residual = ~ ar1(Row):ar1(Column), 
                                     data=Wheat.dat))
  summary(current.asr)$varcomp
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- changeTerms(current.asrt, addRandom = "Row + Column + units", 
                              newResidual = "Row:ar1(Column)")
  testthat::expect_equal(sum(chk.changes(current.asrt$test.summary)), 2)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 5)
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

cat("#### Test for changing the residual model with asreml41\n")
test_that("residual_changeTerms_asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  print(packageVersion("asreml"))
  #'## Load data
  data("Exp355.Control.dat")
  
  dat <- with(dat, dat[order(Genotype, Salt, NP_AMF, InTreat), ])
  
  #Fit model with quotes around AMF_plus 
  # NB must have quotes for character levels, 
  # but cannot have in testranfix term because not in wald.tab or varcomp rownames
  current.asr <- asreml(fixed = TSP ~ Lane + xPosn + AMF*Genotype*NP + 
                          at(AMF, "AMF_plus"):per.col + (Genotype*NP):at(AMF, "AMF_plus"):per.col,
                        random = ~ spl(xPosn) + Position ,
                        residual = ~ Genotype:idh(NP_AMF):InTreat,
                        keep.order=TRUE, data = dat, 
                        maxiter=50, na.action = na.method(x="include"))
  
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(current.asrt$wald.tab),14)
  testthat::expect_true(all(c("AMF:Genotype:NP", "at(AMF, AMF_plus):per.col", "Genotype:at(AMF, AMF_plus):per.col", 
                              "NP:at(AMF, AMF_plus):per.col", "Genotype:NP:at(AMF, AMF_plus):per.col") %in% 
                              rownames(current.asrt$wald.tab)))
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp),7)
  
  t.asrt <- changeTerms(current.asrt, newResidual = "Genotype:NP_AMF:InTreat", 
                        set.terms = "Genotype:NP_AMF:InTreat!R", 
                        initial.values = 1, bounds = "P", ignore.suffices = FALSE)
  testthat::expect_equal(nrow(summary(t.asrt$asreml.obj)$varcomp),1)
  testthat::expect_true(vpc.char(t.asrt$asreml.obj)["Genotype:NP_AMF:InTreat!R"] == "P")
})
