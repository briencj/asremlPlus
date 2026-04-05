# Extracted from test42alldiffsasr.r:1381

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
cat("#### Test for exploreLSDs on Oats with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(dae)
data(Oats.dat)
m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
current.asrt <- as.asrtests(m1.asr)
wald.tab <-  current.asrt$wald.tab
den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", vcov=TRUE)
Var.diffs <- allDifferences(predictions = Var.pred$pvals,
                              classify = "Nitrogen:Variety", 
                              vcov = Var.pred$vcov, tdf = den.df)
kLSD <- Var.diffs$sed
kLSD <- qt(0.975, attr(Var.diffs, which = "tdf")) * kLSD[upper.tri(kLSD)]
kdif <- Var.diffs$differences
kdif <- abs(kdif[upper.tri(kdif)])
kLSD <- kdif >= kLSD
minLSD <- kdif >= Var.diffs$LSD$assignedLSD
fpos <- sum(minLSD & !kLSD)
fneg <- sum(!minLSD & kLSD)
testthat::expect_true(fpos == Var.diffs$LSD$falsePos)
testthat::expect_true(fneg == Var.diffs$LSD$falseNeg)
lsd <- exploreLSDs(Var.diffs)
testthat::expect_equal(names(lsd), c("frequencies", "distinct.vals", "statistics", "accuracy", 
                                       "false.pos", "false.neg","per.pred.accuracy", "LSD"))
testthat::expect_true(all(lapply(c("statistics", "accuracy"), function(k, lsd) nrow(lsd[[k]]), lsd = lsd) == 1))
testthat::expect_equal(names(lsd$frequencies), as.character(seq(17.25, 21.75, 0.5)))
testthat::expect_equal(lsd$distinct.vals, c(17.1, 21.6))
testthat::expect_true(all(abs(lsd$statistics[1,] - 
                                  c(66,17.11869,17.11869,17.11869,20.51095,
                                    21.64642,21.64642,21.64642,21.64642)) < 1e-05))
testthat::expect_true(all(abs(lsd$accuracy[1,] - 
                                  c(66,0.2644909,0.2644909,0.2644909,0.1653879,
                                    0.2091679,0.2091679,0.2091679,0.2091679)) < 1e-05))
testthat::expect_true(fpos == lsd$false.pos["mean"])
testthat::expect_true(fneg == lsd$false.neg["mean"])
testthat::expect_true(all(lsd$false.pos == c(66,3,3,3.0,0,0,0,0,0)))
testthat::expect_true(all(lsd$false.neg == c(66,0,0,0,3,4,4,4,4)))
testthat::expect_true(all(lapply(c("per.pred.accuracy", "LSD"), function(k, lsd) nrow(lsd[[k]]), lsd = lsd) == 12))
testthat::expect_equal(rownames(lsd$per.pred.accuracy), 
                         as.character(fac.combine(as.list(Var.diffs$predictions[c("Nitrogen","Variety")]), 
                                                  combine.levels = TRUE)))
testthat::expect_true(all(abs(lsd$per.pred.accuracy[1,] - 
                                  c(0.2644909,0.2644909,0.2644909,0.1653879,
                                    0.2091679,0.2091679,0.2091679,0.2091679)) < 1e-05))
lsd <- exploreLSDs(Var.diffs, LSDtype = "fact", LSDby = "Nitrogen")
