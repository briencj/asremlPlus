# Extracted from test42alldiffsasr.r:127

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("prediction_alldiffs")
cat("#### Test for allDifferences.data.frame sort.alldiffs on Oats with asreml42\n")
cat("#### Test for findLSDminerrors with LSDby and singular predictions\n")
cat("#### Test for LSDs and halfLSIs on system data with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asremlPlus)
library(dae)
load(system.file("extdata", "testDiffs.rda", package = "asremlPlus", mustWork = TRUE))
diffs.new <- redoErrorIntervals(diffs.ClUp, error.intervals = "half", 
                                 LSDtype = "factor", LSDby = attr(diffs.ClUp, which = "LSDby"))
