# Extracted from test42Selection.r:583

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("model_selection")
cat("#### Test for chooseModel.data.frame with asreml42\n")
cat("#### Test for chooseModel.asrtests with asreml42\n")
cat("#### Test for testing fixed at terms with asreml42\n")
cat("#### Test for changeTerms with random at terms with asreml42\n")
cat("#### Test for testing MET at terms with asreml42\n")
cat("#### Test for at terms in testswapran with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
data(longit.dat)
asreml.options(fail = "soft", ai.sing = TRUE)
current.asr <- do.call(asreml, 
                         args=list(fixed = Area ~ Block + Treatments + Treatments:xDAP,
                                   random = ~ Block:Cart + at(Treatments):spl(xDAP, k = 10) + 
                                     Treatments:DAP + 
                                     Block:Cart:spl(xDAP) + Block:Cart:xDAP,
                                   residual = ~ Block:Cart:ar1h(DAP),
                                   keep.order=TRUE, data = longit.dat, maxit=100))
current.call <- current.asr$call
vpR <- grepl("Block:Cart:DAP!DAP", names(current.asr$vparameters.con))
vpR <- current.asr$vparameters.con[vpR]
(terms <- names(vpR[vpR == "B"]))
fixBoundResidualVariances <-function(current.asr)
  {
    repeat
    {
      asreml.options(fail = "soft", ai.sing = TRUE)
      current.call <- current.asr$call
      vpR <- grepl("Block:Cart:DAP!DAP", names(current.asr$vparameters.con))
      vpR <- current.asr$vparameters.con[vpR]
      (terms <- names(vpR[vpR == "B"]))
      if (length(terms) == 0 || length(sum(vpR == "F")) > 5) break
      current.asr <- setvarianceterms(call = current.call, terms = terms, 
                                      bounds = "F", initial.values = 0.0001,
                                      ignore.suffices = FALSE)
    }
    invisible(current.asr)
  }
current.asr <- fixBoundResidualVariances(current.asr)
testthat::expect_true(all(table(summary(current.asr)$varcomp$bound) ==  c(2,46,1)))
current.asrt <- as.asrtests(current.asr, NULL, NULL, label = "Selected variance model")
testthat::expect_true(current.asrt$asreml.obj$converge)
current.asrt <- testranfix(current.asrt, term = "Treatments:DAP",
                             positive.zero = TRUE)
testthat::expect_equal(current.asrt$test.summary$action[2], "Dropped")
testthat::expect_true(all(table(summary(current.asrt$asreml.obj)$varcomp$bound) ==  c(2,45,1)))
setvpars <- current.asrt$asreml.obj$call$setvparameters
testthat::expect_equal(nrow(setvpars), 1)
testthat::expect_equal(ncol(setvpars), 4)
testthat::expect_equal(setvpars$set.terms, "Block:Cart:DAP!DAP_17")
testthat::expect_equal(setvpars$bounds, "F")
testthat::expect_true(abs(setvpars$initial.values - 1e-04) < 1e-04)
current.asrt <- testswapran(current.asrt, oldterms = "at(Treatments):spl(xDAP, k = 10)",
                              newterms = "at(AMF):Zn:spl(xDAP, k = 10)",
                              simpler = TRUE,
                              label = "Heterogeneous Treatment splines")
testthat::expect_true(getTestPvalue(current.asrt, label = "Heterogeneous Treatment splines") < 0.05)
testthat::expect_true(names(current.asrt$asreml.obj$vparameters[1]) == 
                          "at(Treatments, '-,0'):spl(xDAP, k = 10)")
testthat::expect_true(names(current.asrt$asreml.obj$vparameters[5]) == 
                          "at(Treatments, '+,0'):spl(xDAP, k = 10)")
vpar.vals <- c(233.152932, 502.667930, 74.955973, 2.540186, 42.003197, 61.206138, 
                 32.367734, 36.902978)
names(vpar.vals) <- names(current.asrt$asreml.obj$vparameters[1:8])
testthat::expect_true(all.equal(current.asrt$asreml.obj$vparameters[1:8], 
                                  vpar.vals, tolerance = 1e-02))
current.asr <- do.call(asreml, 
                         args=list(fixed = Area ~ Block + Treatments + Treatments:xDAP,
                                   random = ~ Block:Cart + at(Treatments):spl(xDAP, k = 10) + 
                                     Treatments:DAP + 
                                     Block:Cart:spl(xDAP) + Block:Cart:xDAP + 
                                     idv(Block):Cart:ar1h(DAP),
                                   residual = ~ Block:Cart:DAP,
                                   keep.order=TRUE, data = longit.dat, maxit=100))
summary(current.asr)$varcomp
current.asr <- fixBoundResidualVariances(current.asr)
testthat::expect_true(all(table(summary(current.asr)$varcomp$bound) ==  c(9,40,1)))
current.asrt <- as.asrtests(current.asr, NULL, NULL, label = "Selected variance model")
testthat::expect_false(current.asrt$asreml.obj$converge)
setvpars <- current.asrt$asreml.obj$call$setvparameters
testthat::expect_equal(nrow(setvpars), 9)
testthat::expect_equal(ncol(setvpars), 4)
testthat::expect_equal(sum(grepl("Block:Cart:DAP!DAP", setvpars$set.terms)), 9)
