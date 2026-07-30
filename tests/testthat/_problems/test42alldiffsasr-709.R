# Extracted from test42alldiffsasr.r:709

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

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(dae)
data(Oats.dat)
LSD.hdr <- c("c", "minLSD", "meanLSD", "maxLSD", "assignedLSD", "accuracyLSD")
m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
testthat::expect_equal(length(m1.asr$vparameters),3)
current.asrt <- as.asrtests(m1.asr)
diffs <- predictPlus(m1.asr, classify = "Nitrogen:Variety", 
                       wald.tab = current.asrt$wald.tab, 
                       error.intervals = "half", 
                       tables = "none", Vmatrix = TRUE)
testthat::expect_true(validAlldiffs(diffs))
testthat::expect_true(all(abs(c(66, 15.47426, 18.54066, 19.56707, 18.54066, 0.1653879, 2, 4) - 
                                  diffs$LSD) < 1e-05))
