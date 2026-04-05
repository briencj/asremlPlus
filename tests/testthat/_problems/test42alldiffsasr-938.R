# Extracted from test42alldiffsasr.r:938

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

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(dae)
data(WaterRunoff.dat)
LSD.hdr <- c("c", "minLSD", "meanLSD", "maxLSD", "assignedLSD", "accuracyLSD", "falsePos", "falseNeg")
m1.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                   random = ~ Benches:MainPlots,
                   keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- as.asrtests(m1.asr, NULL, NULL)
testthat::expect_equal(length(m1.asr$vparameters),2)
current.asrt <- as.asrtests(m1.asr)
current.asrt <- rmboundary(current.asrt)
m1.asr <- current.asrt$asreml.obj
testthat::expect_equal(length(m1.asr$vparameters),1)
diffs.full.LSD <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", LSDtype = "factor", 
                                       LSDby = c("Type", "Species"),
                                       tables = "none", Vmatrix = TRUE)
tdf <- attr(diffs.full.LSD, which = "tdf")
t.val <- qt(0.975, tdf)
testthat::expect_true(setequal(names(diffs.full.LSD$predictions), 
                                 c("Sources", "Type", "Species", "predicted.value", 
                                   "standard.error", "upper.halfLeastSignificant.limit", 
                                   "lower.halfLeastSignificant.limit", "est.status")))
testthat::expect_true(!("LSDwarning" %in% names(diffs.full.LSD$predictions)))
testthat::expect_true(all(c("LSDtype", "LSDby", "LSDvalues", "avsed.tolerance", 
                              "accuracy.threshold") %in% names(attributes(diffs.full.LSD$predictions))))
testthat::expect_true((attr(diffs.full.LSD$predictions, which = "avsed.tolerance") == 0.25))
testthat::expect_true(is.na(attr(diffs.full.LSD$predictions, which = "accuracy.threshold")))
testthat::expect_true(all(abs(attr(diffs.full.LSD$predictions, which = "LSDvalues") - 
                                  c(0.3052944, 0.3052944, 0.3099100, 0.3076127, 
                                    0.3076127, 0.3739078, 0.3739078, 0.3783488)) < 1e-05))
testthat::expect_true(all(diffs.full.LSD$LSD$falsePos == c(rep(0,7),1)))
testthat::expect_true(all(diffs.full.LSD$LSD$falseNeg == c(rep(0,7),1)))
medianLSD.dat <-  data.frame(meanLSD = c(0.3052944, 0.3052944, 0.3121896, 0.3052944, 
                                           0.3052944, 0.3739078, 0.3739078, 0.3739078),
                               row.names = c("Landscape,S. iqscjbogxah", "Landscape,S. oymwjrcnepv", 
                                             "Landscape,S. ocphawvtlgi", "Medicinal,S. hbpgtylxqku", 
                                             "Medicinal,S. orcxszbujml", "Culinary,S. xeqackngdrt", 
                                             "Culinary,S. tkujbvipoyr",  "Control,Non-planted"))
levs <- strsplit(rownames(medianLSD.dat), ",", fixed = TRUE)
medianLSDfacs.dat <- data.frame(Type = factor(unlist(lapply(levs, function(lev) lev[1])), 
                                                levels = levels(WaterRunoff.dat$Type)), 
                                  Species = factor(unlist(lapply(levs, function(lev) lev[2])), 
                                                levels = levels(WaterRunoff.dat$Species)),
                                  assignedLSD = medianLSD.dat$meanLSD)
diffs.over.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", 
                                       avsed.tolerance = NA,
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
testthat::expect_true(abs(diffs.over.Acc$LSD$accuracyLSD - 0.1860607) < 01e-05)
testthat::expect_true(diffs.over.Acc$LSD$falsePos ==  15)
testthat::expect_true(diffs.over.Acc$LSD$falseNeg ==  10)
diffs.over.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", 
                                       LSDaccuracy = "maxDev", avsed.tolerance = NA,
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
testthat::expect_true(abs(diffs.over.Acc$LSD$accuracyLSD - 0.1860607) < 01e-05)
testthat::expect_true(diffs.over.Acc$LSD$falsePos ==  15)
testthat::expect_true(diffs.over.Acc$LSD$falseNeg ==  10)
diffs.over.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", 
                                       LSDaccuracy = "q90Dev", avsed.tolerance = NA,
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
testthat::expect_true(abs(diffs.over.Acc$LSD$accuracyLSD - 0.06892194) < 01e-05)
testthat::expect_true(diffs.over.Acc$LSD$falsePos ==  15)
testthat::expect_true(diffs.over.Acc$LSD$falseNeg ==  10)
diffs.over.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", 
                                       LSDaccuracy = "Root", avsed.tolerance = NA,
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
testthat::expect_true(abs(diffs.over.Acc$LSD$accuracyLSD - 0.06833254) < 01e-05)
testthat::expect_true(diffs.over.Acc$LSD$falsePos ==  15)
testthat::expect_true(diffs.over.Acc$LSD$falseNeg ==  10)
diffs.reLSD.q90 <- redoErrorIntervals(diffs.full.LSD, error.intervals = "half", 
                                        LSDtype = "factor", LSDby = c("Type", "Species"), 
                                        LSDtstatistic = "q90")
testthat::expect_true(all(diffs.reLSD.q90$LSD$assignedLSD <= diffs.reLSD.q90$LSD$maxLSD))
testthat::expect_true(all(abs(diffs.reLSD.q90$LSD$assignedLSD - c(0.3052944,0.3052944,0.3099100,0.3076127,
                                                                    0.3076127,0.3739078,0.3739078,0.3783488)) < 1e-05))
testthat::expect_true(all(diffs.reLSD.q90$LSD$falsePos == c(rep(0,7),1)))
testthat::expect_true(all(diffs.reLSD.q90$LSD$falseNeg == c(rep(0,7),1)))
diffs.reLSD <- redoErrorIntervals(diffs.full.LSD, error.intervals = "half", 
                                    LSDtype = "supplied", LSDby = c("Type", "Species"), 
                                    LSDsupplied = medianLSD.dat)
testthat::expect_is(diffs.reLSD, "alldiffs")
testthat::expect_true(validAlldiffs(diffs.reLSD))
testthat::expect_equal(nrow(diffs.reLSD$predictions),40)
testthat::expect_equal(ncol(diffs.reLSD$predictions),8)
testthat::expect_true(all(names(diffs.reLSD$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDby", "LSDstatistic") %in% 
                              names(attributes(diffs.reLSD))))
testthat::expect_true(all(c( "LSDtype", "LSDstatistic", "LSDby", "LSDvalues") %in% 
                              names(attributes(diffs.reLSD$predictions))))
testthat::expect_true(all(rownames(diffs.reLSD$LSD) == rownames(medianLSD.dat)))
testthat::expect_true(all(LSD.hdr %in% names(diffs.reLSD$LSD)))
testthat::expect_true(all(abs(medianLSD.dat$meanLSD - diffs.reLSD$LSD$assignedLSD) < 1e-05))
testthat::expect_true(all(diffs.reLSD$LSD$falsePos == c(rep(0,7),1)))
testthat::expect_true(all(diffs.reLSD$LSD$falseNeg == rep(0,8)))
tmp <- merge(diffs.reLSD$predictions, medianLSDfacs.dat)
testthat::expect_true(all(abs((tmp$upper.halfLeastSignificant.limit - tmp$lower.halfLeastSignificant.limit) 
                                - tmp$assignedLSD) < 1e-05))
exploreLSDs(diffs.reLSD.q90, LSDtype = "factor", LSDby = c("Type", "Species"))
