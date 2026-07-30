# Extracted from test42alldiffsasr.r:753

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
LSD.dat <- diffs$LSD
LSD.dat$assignedLSD <- qt(0.975, attr(diffs, which = "tdf"))*median(diffs$sed, na.rm = TRUE)
diffs.reLSD <- redoErrorIntervals(diffs, error.intervals = "half", 
                                    LSDtype = "supplied", LSDsupplied = LSD.dat["assignedLSD"])
testthat::expect_is(diffs.reLSD, "alldiffs")
testthat::expect_true(validAlldiffs(diffs.reLSD))
testthat::expect_equal(nrow(diffs.reLSD$predictions),12)
testthat::expect_equal(ncol(diffs.reLSD$predictions),7)
testthat::expect_true(all(names(diffs.reLSD$predictions)[5:6] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic") %in% names(attributes(diffs.reLSD))))
testthat::expect_true(is.null(attr(diffs.reLSD, which = "LSDby")))
testthat::expect_true(all(c( "LSDtype", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.reLSD$predictions))))
testthat::expect_true(rownames(diffs.reLSD$LSD) == "overall")
testthat::expect_true(all(LSD.hdr %in% names(diffs.reLSD$LSD)))
testthat::expect_true(all(abs(c(66, 15.47426, 18.54066, 19.56707, 19.56707, 0.2091679, 0, 4) - 
                                  diffs.reLSD$LSD) < 1e-05))
testthat::expect_true(all(abs((diffs.reLSD$predictions$upper.halfLeastSignificant.limit - 
                                   diffs.reLSD$predictions$lower.halfLeastSignificant.limit) 
                                - diffs.reLSD$LSD$assignedLSD) < 1e-05))
medianLSD <- LSD.dat[["meanLSD"]]
names(medianLSD) <- "overall"
diffs.med.reLSD <- redoErrorIntervals(diffs, error.intervals = "half", 
                                        LSDtype = "supplied", LSDsupplied = medianLSD)
testthat::expect_is(diffs.med.reLSD, "alldiffs")
testthat::expect_true(validAlldiffs(diffs.med.reLSD))
testthat::expect_equal(nrow(diffs.med.reLSD$predictions),12)
testthat::expect_equal(ncol(diffs.med.reLSD$predictions),7)
testthat::expect_true(all(names(diffs.med.reLSD$predictions)[5:6] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic") %in% names(attributes(diffs.med.reLSD))))
testthat::expect_true(is.null(attr(diffs.med.reLSD, which = "LSDby")))
testthat::expect_true(all(c( "LSDtype", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.med.reLSD$predictions))))
testthat::expect_true(rownames(diffs.med.reLSD$LSD) == "overall")
testthat::expect_true(all(LSD.hdr %in% names(diffs.med.reLSD$LSD)))
testthat::expect_true(all(abs(c(66, 15.47426, 18.54066, 19.56707, 18.54066, 0.1653879,2,4) - 
                                  diffs.med.reLSD$LSD) < 1e-05))
