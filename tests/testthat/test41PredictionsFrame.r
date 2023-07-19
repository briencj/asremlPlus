#devtools::test("asremlPlus")
context("predictions_alldiffs")
asr41.lib <- "D:\\Analyses\\R ASReml4.1" 

cat("#### Test for predictions.frame on Oats with asreml41\n")
test_that("PredictionsFrame_asreml41", {
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  ## Use asreml to get predictions and associated statistics
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  current.asrt <- as.asrtests(m1.asr)
  Var.pred <- asreml::predict.asreml(m1.asr, classify="Nitrogen:Variety", 
                                     sed=TRUE)
  if (getASRemlVersionLoaded(nchar = 1) == "3")
    Var.pred <- Var.pred$predictions
  Var.preds <- as.predictions.frame(Var.pred$pvals, se = "std.error", 
                                    est.status = "status")
  
  ## Check the class and validity of the predictions.frame
  testthat::expect_true(is.predictions.frame(Var.preds))
  testthat::expect_true(validPredictionsFrame(Var.preds))
})
