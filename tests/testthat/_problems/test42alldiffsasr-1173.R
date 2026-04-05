# Extracted from test42alldiffsasr.r:1173

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("prediction_alldiffs")
cat("#### Test for allDifferences.data.frame sort.alldiffs on Oats with asreml42\n")
cat("#### Test for LSDs and halfLSIs on system data with asreml42\n")
cat("#### Test for LSD on Oats with asreml42\n")
cat("#### Test for sort.alldiffs on Smarthouse with asreml42\n")
cat("#### Test for LSD with sort.alldiffs on Smarthouse with asreml42\n")
cat("#### Test for LSDsupplied on Oats with asreml42\n")
cat("#### Test for single-prediction LSDs in 821 Barley with asreml42\n")
cat("#### Test for LSD on WaterRunoff with asreml42\n")
cat("#### Test for exploreLSDs on WaterRunoff with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(dae)
data(WaterRunoff.dat)
m1.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                   random = ~ Benches:MainPlots,
                   keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- as.asrtests(m1.asr, NULL, NULL)
testthat::expect_equal(length(m1.asr$vparameters),2)
current.asrt <- as.asrtests(m1.asr)
current.asrt <- rmboundary(current.asrt)
current.asr <- current.asrt$asreml.obj
TS.diffs <- predictPlus(classify = "Sources:Type", 
                          asreml.obj = current.asr, 
                          wald.tab = current.asrt$wald.tab, 
                          present = c("Sources", "Type", "Species"),
                          tables = "none")
LSDstat <- exploreLSDs(TS.diffs, LSDtype = "factor.combinations", LSDby = "Sources")
