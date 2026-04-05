# Extracted from test42alldiffsasr.r:1954

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
cat("#### Test for sort.alldiffs on WaterRunoff with asreml42\n")
cat("#### Test for sort.alldiffs on Oats with asreml42\n")
cat("#### Test for subset.alldiffs on Smarthouse with asreml42\n")
cat("#### Test for facCombine.alldiffs on Ladybird with asreml42\n")
cat("#### Test for facRecast.alldiffs on Ladybird with asreml42\n")
cat("#### Test for linear.transformation on Oats with asreml42\n")

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
testthat::expect_equal(length(m1.asr$vparameters),3)
current.asrt <- as.asrtests(m1.asr)
diffs <- predictPlus(m1.asr, classify = "Nitrogen:Variety", Vmatrix = TRUE, 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
testthat::expect_is(diffs, "alldiffs")
testthat::expect_equal(length(attr(diffs$predictions, which = "heading")),1)
testthat::expect_true("asreml.predict" %in% class(diffs$predictions))
testthat::expect_equal(nrow(diffs$vcov),12)
testthat::expect_true(all(colnames(diffs$vcov)[1:2] %in% c("0,Victory", "0,Golden Rain")))
testthat::expect_true(all(abs((diffs$vcov[1,1:2] - c(82.93704, 35.74618))) < 1e-4))
nv <- length(levels(diffs$predictions$Variety))
L <- cbind(kronecker(matrix(rep(1, 3), nrow = 3), diag(1, nrow=nv)), 
             diag(-1, nrow=3*nv))
rownames(L) <- colnames(diffs$vcov)[4:12]
diffs.L <- predictPlus(m1.asr, classify = "Nitrogen:Variety", 
                         linear.transformation = L, Vmatrix = TRUE,  
                         wald.tab = current.asrt$wald.tab,
                         error.intervals = "Conf", tables = "none")
testthat::expect_is(diffs.L, "alldiffs")
testthat::expect_equal(length(attr(diffs.L$predictions, which = "heading")),2)
testthat::expect_true("asreml.predict" %in% class(diffs.L$predictions))
testthat::expect_true(abs((diffs$vcov[1,1] - diffs$vcov[4,1]) + 
                              (diffs$vcov[4,4] - diffs$vcov[1,4]) - diffs.L$vcov[1,1]) < 1e-04)
testthat::expect_equal(as.character(diffs.L$predictions$Combination[1]), "0.2,Victory")
testthat::expect_true(abs((diffs.L$predictions$predicted.value[1] - 
                               qt(0.975,attr(diffs.L, which = "tdf")) * 
                               diffs.L$predictions$standard.error[1]) - 
                              diffs.L$predictions$lower.Confidence.limit[1]) < 1e-04)
diffs.mod <- predictPlus(m1.asr, classify = "Nitrogen:Variety", 
                           linear.transformation = ~ Variety + Nitrogen, 
                           wald.tab = current.asrt$wald.tab,
                           error.intervals = "half", 
                           LSDtype = "factor.comb",
                           LSDby = "Nitrogen",
                           tables = "none")
testthat::expect_is(diffs.mod, "alldiffs")
testthat::expect_equal(length(attr(diffs.mod$predictions, which = "heading")),2)
testthat::expect_true("asreml.predict" %in% class(diffs.mod$predictions))
testthat::expect_true(is.null(diffs.mod$vcov))
