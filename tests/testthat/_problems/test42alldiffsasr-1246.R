# Extracted from test42alldiffsasr.r:1246

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
TS.diffs.all <- redoErrorIntervals(TS.diffs, LSDtype = "supplied", LSDsupplied = LSDall["LSD"])
TS.diffs.all$LSD["assignedLSD"]
testthat::expect_true(abs(TS.diffs.all$LSD["assignedLSD"] - 0.3238331) < 0.00001)
LSDallwt <- findLSDminerrors(TS.diffs, false.pos.wt = 30)
testthat::expect_equal(rownames(LSDallwt), "overall")
testthat::expect_true(all(abs(LSDallwt - c(0.3870779, 0, 36, 36)) < 0.00001))
LSDmin <- findLSDminerrors(TS.diffs, LSDtype = "factor.combinations", LSDby = "Sources")
testthat::expect_equal(rownames(LSDmin), c("Rainwater", "Recycled water", "Tap water", 
                                             "Rain+Basalt", "Rain+Dolomite", "Rain+Quartzite"))
testthat::expect_equal(LSDmin$false.criterion, c(0, 0, 1, 0, 1, 1))
TS.diffs.min <- redoErrorIntervals(TS.diffs, LSDtype = "supplied", LSDby = "Sources",
                                     LSDsupplied = LSDmin["LSD"])
testthat::expect_true(all(abs(TS.diffs.min$LSD$assignedLSD - 
                                  c(0.1982634, 0.2770138, 0.2680011, 0.2978220, 
                                    0.2949564, 0.2464056)) < 0.0001))
LSDminwt <- findLSDminerrors(TS.diffs, LSDtype = "factor.combinations", LSDby = "Sources",
                             false.pos.wt = c(5,5,0,5,1,1))
testthat::expect_equal(rownames(LSDminwt), c("Rainwater", "Recycled water", "Tap water", 
                                               "Rain+Basalt", "Rain+Dolomite", "Rain+Quartzite"))
testthat::expect_equal(LSDminwt$false.pos, c(0, 0, 1, 0, 0, 0))
testthat::expect_equal(LSDminwt$false.criterion, c(0, 0, 0, 0, 1, 1))
kLSD <- TS.diffs$sed[17:20,17:20]
kLSD <- qt(0.975, attr(TS.diffs, which = "tdf")) * kLSD[upper.tri(kLSD)]
kdif <- TS.diffs$differences[17:20,17:20]
kdif <- abs(kdif[upper.tri(kdif)])
kLSD <- kdif >= kLSD
minLSD <- kdif >= LSDstat$statistics$min[6]
fpos <- sum(minLSD & !kLSD)
fneg <- sum(!minLSD & kLSD)
testthat::expect_true(fpos == LSDstat$false.pos$min[6])
testthat::expect_true(fneg == LSDstat$false.neg$min[6])
testthat::expect_true(all(abs(LSDstat$accuracy[1,] - 
                                  c(3, 0.5463438,0.4094721, 0.2442706, 
                                    0.2679453,0.3268454,0.340344,0.3481875,0.3533133)) < 1e-05))
testthat::expect_true(all(lapply(c("per.pred.accuracy", "LSD"), 
                                   function(k, LSDstat) nrow(LSDstat[[k]]), LSDstat = LSDstat) == 20))
testthat::expect_equal(rownames(LSDstat$per.pred.accuracy), 
                         as.character(fac.combine(as.list(TS.diffs$predictions[c("Sources","Type")]), 
                                                  combine.levels = TRUE)))
