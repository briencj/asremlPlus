# Extracted from test42Sortalldiffs.r:40

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("prediction_presentation")
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
testthat::expect_is(diffs, "alldiffs")
testthat::expect_true(validAlldiffs(diffs))
testthat::expect_equal(nrow(diffs$predictions),120)
testthat::expect_equal(ncol(diffs$predictions),8)
testthat::expect_equal(as.character(diffs$predictions$Genotype[1]),"Axe")
testthat::expect_equal(length(attributes(diffs)),11)
testthat::expect_true(is.null(attr(diffs, which = "sortOrder")))
testthat::expect_warning(plotPredictions(data = diffs$predictions, 
                                          classify = "Genotype:A:B", 
                                          y = "predicted.value", 
                                          error.intervals = "StandardError",  
                                          y.title = attr(diffs, 
                                                         which = "response.title"),
                                          sortFactor = "Genotype"))
