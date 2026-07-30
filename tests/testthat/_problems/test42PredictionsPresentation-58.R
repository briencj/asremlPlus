# Extracted from test42PredictionsPresentation.r:58

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("prediction_presentation")
cat("#### Test for Intercept prediction on Oats with asreml42\n")
cat("#### Test for NA on Oats with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(dae)
data(Oats.dat)
Oats.dat$Nitrogen[8] <- NA
Oats.dat$xNitrogen <- as.numfac(Oats.dat$Nitrogen)
m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   na.action = na.method(x = "include"),
                   data=Oats.dat)
testthat::expect_equal(length(m1.asr$vparameters),3)
Trt.pred <- predictPlus(m1.asr, classify="(Variety:Nitrogen)", tables = "none")$predictions
testthat::expect_equal(nrow(Trt.pred), 12)
testthat::expect_true(abs( Trt.pred$predicted.value[3] - 109.07892) < 1e-04)
testthat::expect_true(abs( Trt.pred$standard.error[3] - 9.414492) < 1e-04)
