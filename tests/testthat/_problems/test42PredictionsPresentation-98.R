# Extracted from test42PredictionsPresentation.r:98

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
m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   na.action = na.method(x = "omit"),
                   data=Oats.dat)
testthat::expect_equal(length(m1.asr$vparameters),3)
Trt.pred <- predictPlus(m1.asr, classify="(Variety:Nitrogen)", tables = "none")$predictions
testthat::expect_equal(nrow(Trt.pred), 12)
testthat::expect_true(abs( Trt.pred$predicted.value[3] - 109.07892) < 1e-04)
testthat::expect_true(abs( Trt.pred$standard.error[3] - 9.414492) < 1e-04)
m1.asr <- asreml(Yield ~ xNitrogen*Variety, 
                   random=~Blocks/Wplots,
                   na.action = na.method(x = "include"),
                   data=Oats.dat)
testthat::expect_equal(length(m1.asr$vparameters),3)
Trt.pred <- predictPlus(m1.asr, classify="Variety:xNitrogen", 
                          levels = list(xNitrogen =  c(0, 0.2, 0.4, 0.6)),
                          tables = "none")$predictions
testthat::expect_equal(nrow(Trt.pred), 12)
testthat::expect_true(abs( Trt.pred$predicted.value[3] - 105.99910) < 1e-04)
testthat::expect_true(abs( Trt.pred$standard.error[3] - 8.403037) < 1e-04)
Oats.dat$Nitrogen <- factor(Oats.dat$Nitrogen, exclude = NULL)
m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   na.action = na.method(x = "include"),
                   data=Oats.dat)
testthat::expect_equal(length(m1.asr$vparameters),3)
Trt.pred <- predictPlus(m1.asr, classify="(Variety:Nitrogen)", tables = "none")$predictions
testthat::expect_equal(nrow(Trt.pred), 13)
testthat::expect_true(abs( Trt.pred$predicted.value[3] - 109.07892) < 1e-04)
