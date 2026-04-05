# Extracted from test42alldiffsasr.r:338

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
testthat::expect_true(all(attr(Int.diffs, which = "LSDtype") == "overall"))
testthat::expect_true(all(attr(Int.diffs, which = "LSDstatistic") == "mean"))
Int.diffs <- linTransform(Var.diffs, linear.transformation = ~ 1,
                            error.intervals = "half", 
                            LSDtype = "overall", tables = "none")
testthat::expect_equal(names(Int.diffs$predictions), 
                         c("Nitrogen", "Variety", "predicted.value", "standard.error", 
                           "upper.halfLeastSignificant.limit", 
                           "lower.halfLeastSignificant.limit", "est.status"))
testthat::expect_true(all(abs(Int.diffs$predictions$predicted.value - mean(Oats.dat$Yield)) < 1e-04))
testthat::expect_true(all(abs(Int.diffs$predictions$upper.halfLeastSignificant.limit - 114.4348) < 1e-04))
testthat::expect_true(all(abs(c(0, 20.92506, 20.92506, 20.92506, 20.92506) - 
                                  Int.diffs$LSD[1:5]) < 1e-05))
testthat::expect_true(all(is.na(Int.diffs$LSD[6:8])))
lsd1 <- exploreLSDs(Int.diffs)
testthat::expect_true(all(abs(c(0, rep(20.92506, 8)) - 
                                  lsd1$statistics) < 1e-05))
testthat::expect_true(all(sapply(lsd1[c("accuracy","false.pos","false.neg")], 
                                   function(x) all(is.na(x[-1])) & x[1] == 0)))
testthat::expect_true(all(sapply(lsd1$per.pred.accuracy, function(x) all(is.na(x)))))
testthat::expect_true(all(abs(lsd1$LSD[upper.tri(lsd1$LSD)] - 20.92506 < 1e-05)))
Int.diffs <- linTransform(Var.diffs, linear.transformation = ~ 1,
                            error.intervals = "half", 
                            LSDtype = "factor", LSDby = "Nitrogen", 
                            tables = "none")
testthat::expect_equal(names(Int.diffs$predictions), 
                         c("Nitrogen", "Variety", "predicted.value", "standard.error", 
                           "upper.halfLeastSignificant.limit", 
                           "lower.halfLeastSignificant.limit", "est.status"))
testthat::expect_true(all(abs(Int.diffs$predictions$predicted.value - mean(Oats.dat$Yield)) < 1e-04))
testthat::expect_true(all(is.na(Int.diffs$predictions$upper.Confidence.limit)))
testthat::expect_true(all(abs(c(0, rep(20.92506, 4)) - 
                                  Int.diffs$LSD[1:5]) < 1e-05))
