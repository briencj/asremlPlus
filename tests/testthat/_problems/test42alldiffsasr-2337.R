# Extracted from test42alldiffsasr.r:2337

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("prediction_alldiffs")
cat("#### Test for allDifferences.data.frame sort.alldiffs on Oats with asreml42\n")
cat("#### Test for findLSDminerrors with LSDby and singular predictions\n")
cat("#### Test for LSDs and halfLSIs on system data with asreml42\n")
cat("#### Test for LSD on Oats with asreml42\n")
cat("#### Test for sort.alldiffs on Smarthouse with asreml42\n")
cat("#### Test for LSD with sort.alldiffs on Smarthouse with asreml42\n")
cat("#### Test for LSDsupplied on Oats with asreml42\n")
cat("#### Test for single-prediction LSDs in 821 Barley with asreml42\n")
cat("#### Test for LSD on WaterRunoff with asreml42\n")
cat("#### Test for exploreLSDs on WaterRunoff with asreml42\n")
cat("#### Test for exploreLSDs on Oats with asreml42\n")
cat("#### Test for sort.alldiffs on WaterRunoff with asreml42\n")
cat("#### Test for sort.alldiffs on Oats with asreml42\n")
cat("#### Test for subset.alldiffs on Smarthouse with asreml42\n")
cat("#### Test for facCombine.alldiffs on Ladybird with asreml42\n")
cat("#### Test for facRecast.alldiffs and approx.se on Ladybird with asreml42\n")
cat("#### Test for linear.transformation on Oats with asreml42\n")
cat("#### Test for linear.transformation on dat699 with asreml42\n")
cat("#### Test for linear.transformation on WaterRunoff with asreml42\n")
cat("#### Test for addBacktransforms on WaterRunoff with asreml42\n")
cat("#### Test for ratioTansforms on system data with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asremlPlus)
library(dae)
load(system.file("extdata", "testDiffs.rda", package = "asremlPlus", mustWork = TRUE))
Preds.ratio.RGR <- ratioTransform.alldiffs(alldiffs.obj = diffs.RGR,
                                             ratio.factor = "Salinity", 
                                             numerator.levels = "Salt",
                                             denominator.levels = "Control")
