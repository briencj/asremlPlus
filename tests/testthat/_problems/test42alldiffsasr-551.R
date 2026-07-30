# Extracted from test42alldiffsasr.r:551

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

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(dae)
data(Smarthouse.dat)
m1.asr <- asreml(y1 ~ Genotype*A*B, 
                   random=~Replicate/Mainplot/Subplot,
                   data=Smarthouse.dat)
testthat::expect_equal(length(m1.asr$vparameters),4)
current.asrt <- as.asrtests(m1.asr)
current.asrt <- rmboundary(current.asrt)
m <- current.asrt$asreml.obj
testthat::expect_equal(length(m$vparameters),3)
diffs <- predictPlus(m, classify = "Genotype:A:B", 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
testthat::expect_true(is.alldiffs(diffs))
testthat::expect_true(validAlldiffs(diffs))
testthat::expect_equal(nrow(diffs$predictions),120)
testthat::expect_equal(ncol(diffs$predictions),8)
testthat::expect_equal(as.character(diffs$predictions$Genotype[1]),"Axe")
testthat::expect_true(is.null(attr(diffs, which = "sortOrder")))
diffs.reord <- renewClassify(diffs, newclassify = "A:B:Genotype")
testthat::expect_equal(as.character(diffs.reord$predictions$Genotype[1]),"Axe")
testthat::expect_equal(as.character(diffs.reord$predictions$Genotype[2]),"Espada")
testthat::expect_true(abs(diffs.reord$predictions$predicted.value[2] - -0.2265723017) < 1e-06)
testthat::expect_true(all(names(diffs.reord$predictions)[6:7] == 
                              c("upper.StandardError.limit", "lower.StandardError.limit")))
testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic") %in% names(attributes(diffs.reord))))
testthat::expect_true(!("LSDby" %in% names(attributes(diffs.reord))))
testthat::expect_true(attr(diffs.reord, which = "LSDtype") == "overall")
testthat::expect_true(attr(diffs.reord, which = "LSDstatistic") == "mean")
testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.reord$predictions))))
testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.reord$backtransforms))))
testthat::expect_silent(plotPredictions(data = diffs$predictions, 
                                          classify = "Genotype:A:B", 
                                          y = "predicted.value", 
                                          error.intervals = "StandardError",  
                                          y.title = attr(diffs, 
                                                         which = "response.title")))
testthat::expect_warning(plotPredictions(data = diffs$predictions, 
                                          classify = "Genotype:A:B", 
                                          y = "predicted.value", 
                                          error.intervals = "StandardError",  
                                          y.title = attr(diffs, 
                                                         which = "response.title"),
                                          sortFactor = "Genotype"))
