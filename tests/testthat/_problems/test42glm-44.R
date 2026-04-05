# Extracted from test42glm.r:44

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("prediction_alldiffs")
cat("#### Test for alldiffs with GLM on budworm using asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
df <- data.frame(ldose = rep(0:5, 2),
                   numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16),
                   sex = factor(rep(c("M", "F"), c(6, 6))),
                   N=rep(20,12))
df$numalive <- df$N-df$numdead
df$p <- df$numdead/df$N
asreml.options(pworkspace = "1gb")
as1 <- asreml(p ~ ldose + sex, family=asr_binomial(total=N), data=df)
testthat::expect_warning(
    diffs <- predictPlus(as1, classify = "sex:ldose", levels = list(ldose = 0:5),
                         transform.function = "logit", tables = "none"),
    regexp = "Denominator degrees of freedom obtained using dDF.fault method residual")
testthat::expect_equal(length(diffs), 7)
testthat::expect_true(!any(c("transformed.value", "approx.se") %in% names(diffs$predictions)))
testthat::expect_true(abs(diffs$predictions$predicted.value[1] - -3.473155) < 1e-05)
testthat::expect_true(abs(diffs$predictions$standard.error[1] - 0.4685204) < 1e-05)
testthat::expect_equal(nrow(diffs$backtransforms), 12)
testthat::expect_true(!any(c("transformed.value", "approx.se") %in% names(diffs$backtransforms)))
testthat::expect_true(abs(diffs$backtransforms$backtransformed.predictions[1] - 0.03008577) < 1e-05)
testthat::expect_true(abs(diffs$backtransforms$standard.error[1] - 0.01103991) < 1e-05)
testthat::expect_warning(
    diffs <- predictPlus(as1, classify = "sex:ldose", levels = list(ldose = 0:5),
                         tables = "none"),
    regexp = "Denominator degrees of freedom obtained using dDF.fault method residual")
testthat::expect_true(is.null(diffs$backtransforms))
testthat::expect_true(all(c("transformed.value", "approx.se") %in% names(diffs$predictions)))
testthat::expect_true(abs(diffs$predictions$predicted.value[1] - -3.473155) < 1e-05)
testthat::expect_true(abs(diffs$predictions$standard.error[1] - 0.4685204) < 1e-05)
testthat::expect_true(abs(diffs$predictions$transformed.value[1] - 0.03008577) < 1e-05)
testthat::expect_true(abs(diffs$predictions$approx.se[1] - 0.01103991) < 1e-05)
