# Extracted from test42alldiffsasr.r:2087

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
cat("#### Test for linear.transformation on dat699 with asreml42\n")
cat("#### Test for linear.transformation on WaterRunoff with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(dae)
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE)
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        data= WaterRunoff.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
diffs <- predictPlus(classify = "Sources:Species", Vmatrix = TRUE, 
                       asreml.obj = current.asr, tables = "none", 
                       wald.tab = current.asrt$wald.tab, 
                       present = c("Type","Species","Sources"))
diffs.sub <- linTransform(diffs, classify = "Sources:Species", Vmatrix = TRUE,
                            linear.transformation = ~ Sources + Species,
                            tables = "none")
testthat::expect_equal(diffs.sub$predictions$predicted.value[1] - 
                           diffs.sub$predictions$predicted.value[6],
                         diffs.sub$predictions$predicted.value[7] - 
                           diffs.sub$predictions$predicted.value[12])
L <- kronecker(diag(1, nrow = 4), 
                 cbind(diag(1, nrow = 5), matrix(rep(-1, 5), ncol = 1)))
L <- mat.dirsum(list(L, 
                       kronecker(diag(1, nrow = 2), 
                                 cbind(diag(1, nrow = 7), 
                                       matrix(rep(-1, 7), ncol = 1)))))
testthat::expect_silent(diffs.L.EGLS <- linTransform(diffs.sub, 
                                                   classify = "Sources:Species",
                                                   linear.transformation = L,
                                                   tables = "none"))
testthat::expect_silent(diffs.L <- linTransform(diffs.sub, 
                                                  classify = "Sources:Species",
                                                  linear.transformation = L, 
                                                  EGLS.linTransform = FALSE,
                                                  tables = "none"))
