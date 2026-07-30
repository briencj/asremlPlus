# Extracted from test42alldiffsasr.r:1244

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
testthat::expect_equal(names(LSDstat), c("frequencies", "distinct.vals", "statistics", "accuracy", 
                                           "false.pos", "false.neg", "per.pred.accuracy", "LSD"))
testthat::expect_true(all(lapply(c("frequencies", "statistics", "accuracy"), 
                            function(k, LSDstat) nrow(LSDstat[[k]]), LSDstat = LSDstat) == 6))
testthat::expect_true(all(unlist(lapply(c("frequencies", "statistics", "accuracy"), function(k, LSDstat, dat) 
    all(rownames(LSDstat[[k]])== levels(WaterRunoff.dat$Sources)), LSDstat = LSDstat, dat = WaterRunoff.dat))))
testthat::expect_equal(names(LSDstat$frequencies), as.character(seq(0.17, 0.39, 0.02)))
testthat::expect_equal(LSDstat$distinct.vals$`Rain+Basalt`, c(0.197,0.294,0.307))
testthat::expect_true(all(abs(LSDstat$statistics[1,] - 
                                  c(3, 0.1982634,0.2175164, 0.246396, 
                                    0.2708314,0.2945287, 0.300556,0.3041724,0.3065833)) < 1e-05))
LSDst <- pickLSDstatistics(TS.diffs, LSDtype = "factor.combinations", LSDby = "Sources")
testthat::expect_true(all(LSDst == c("min","median","mean","q75","q75","mean")))
LSDstw1 <- pickLSDstatistics(TS.diffs, LSDtype = "factor.combinations", LSDby = "Sources", 
                               false.pos.wt = 1)
testthat::expect_true(all(LSDstw1 == c("min","median","mean","q75","q75","mean")))
TS.diffs.var <- recalcLSD(TS.diffs, LSDtype = "factor.combinations", LSDby = "Sources", 
                            LSDstatistic = c("q10","med","med","q75","q75", "med"))
testthat::expect_true(all(TS.diffs.var$LSD$falsePos) == 0)
testthat::expect_equal(sum(TS.diffs.var$LSD$falseNeg), 3)
LSDall <- findLSDminerrors(TS.diffs)
testthat::expect_equal(rownames(LSDall), "overall")
testthat::expect_true(all(abs(LSDall - c(0.3238331, 1, 24, 34)) < 0.00001))
