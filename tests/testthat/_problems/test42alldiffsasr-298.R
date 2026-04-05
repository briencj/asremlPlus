# Extracted from test42alldiffsasr.r:298

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("prediction_alldiffs")
cat("#### Test for allDifferences.data.frame sort.alldiffs on Oats with asreml42\n")
cat("#### Test for LSDs and halfLSIs on system data with asreml42\n")
cat("#### Test for LSD on Oats with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(dae)
data(Oats.dat)
m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat, pworkspace = "1Gb")
current.asrt <- as.asrtests(m1.asr)
wald.tab <-  current.asrt$wald.tab
den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", vcov=TRUE)
Var.diffs <- allDifferences(predictions = Var.pred$pvals,
                              classify = "Nitrogen:Variety", 
                              vcov = Var.pred$vcov, tdf = den.df)
testthat::expect_true(all("LSDtype" %in% names(attributes(Var.diffs))))
testthat::expect_true(all(attr(Var.diffs, which = "LSDtype") == "overall"))
testthat::expect_true(all(attr(Var.diffs, which = "LSDstatistic") == "mean"))
Int.pred <- predict(m1.asr, classify="Intercept", vcov=TRUE)
Int.diffs <- allDifferences(predictions = Int.pred$pvals,
                              classify = "Intercept", 
                              vcov = Int.pred$vcov, tdf = den.df)
testthat::expect_true(all("LSDtype" %in% names(attributes(Int.diffs))))
