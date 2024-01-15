#devtools::test("asremlPlus")
context("model_selection")

cat("#### Test for REMLRT with asreml42\n")
test_that("REMLRT_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(Wheat.dat)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  # Fit initial model
  m1.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column + units,
                   residual = ~ ar1(Row):ar1(Column), 
                   maxit = 30, data=Wheat.dat)
  summary(m1.asr)$varcomp
  info <- infoCriteria(m1.asr)
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1346.76764), 1e-02)
  
  #Fit model without the units term
  m2.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column,
                   residual = ~ ar1(Row):ar1(Column), 
                   maxit = 30, data=Wheat.dat)
  summary(m2.asr)$varcomp
  info <- infoCriteria(m2.asr)
  testthat::expect_equal(info$varDF, 4)
  testthat::expect_lt(abs(info$AIC - 1352.941), 1e-03)
  testthat::expect_warning(
    test <- REMLRT(m2.asr, m1.asr), 
    regexp = "There were a total of 1 bound terms. These bound terms occur in both models")
  testthat::expect_lt(abs(test$p - 0.004232946), 1e-03)
  testthat::expect_equal(test$DF, 1)
  testthat::expect_warning(
    test <- REMLRT(m2.asr, m1.asr, DF = 1), 
    regexp = "There were a total of 1 bound terms. These bound terms occur in both models")
  testthat::expect_lt(abs(test$p - 0.004232946), 1e-03)
  testthat::expect_equal(test$DF, 1)

  m3.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column,
                   residual = ~ Row:Column, 
                   maxit = 30, data=Wheat.dat)
  summary(m3.asr)$varcomp
  test3 <- REMLRT(m3.asr, m1.asr)
  testthat::expect_lt(abs(test3$p - 2.596812e-13), 1e-03)
  testthat::expect_equal(test3$DF, 2)
  test3 <- REMLRT(m3.asr, m1.asr, DF = 3)
  testthat::expect_lt(abs(test3$p - 1.603828e-12), 1e-03)
  testthat::expect_equal(test3$DF, 3)
  
  info <- infoCriteria(m3.asr, IClikelihood = "full")
  testthat::expect_equal(info$fixedDF, 31)
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1720.888), 5e-03)
  testthat::expect_lt(abs(info$BIC - 1823.25), 5e-03)
  testthat::expect_lt(abs(info$loglik - m3.asr$loglik), 130)
})

cat("#### Test for wheat76 example with asreml42\n")
test_that("Wheat_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  ## Fit several models to the wheat data and caclulate their ICs
  data(Wheat.dat)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  # Fit initial model
  m.max <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                  random = ~ Row + Column + units,
                  residual = ~ ar1(Row):ar1(Column), 
                  maxit = 30, data=Wheat.dat)

  #Drop term for within Column pairs
  m1 <- asreml(yield ~ Rep + Variety, 
               random = ~ Row + Column + units,
               residual = ~ ar1(Row):ar1(Column), 
               maxit = 30, data=Wheat.dat)
  
  #Drop nugget term
  m2 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
               random = ~ Row + Column,
               residual = ~ ar1(Row):ar1(Column), 
               maxit = 30, data=Wheat.dat)

  #Drop Row autocorrelation
  m3 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                  random = ~ Row + Column + units,
                  residual = ~ Row:ar1(Column), 
                  data=Wheat.dat)

  #Drop Col autocorrelation
  m4 <- asreml(yield ~ Rep + WithinColPairs + Variety, 
               random = ~ Row + Column + units,
               residual = ~ ar1(Row):Column, 
               maxit = 30, data=Wheat.dat)

  mods.asr <- list(m.max, m1, m2, m3, m4)
  ic <- infoCriteria(mods.asr, IClikelihood = "full")
  testthat::expect_equal(nrow(ic), 5)
  testthat::expect_true(all(ic$fixedDF == c(31, 30, 31, 31, 31)))
  testthat::expect_true(all(ic$varDF == c(5, 5, 4, 4, 5)))
  testthat::expect_true(all(abs(ic$AIC - c(1653.100,1651.294,1654.613,1669.928,1708.997)) < 1e-01))
  testthat::expect_true(abs(ic$BIC[1] - 1761.483) < 1)
  testthat::expect_true(abs(ic$loglik[1] - (-790.5502)) < 1e-01)

})

cat("#### Test for IC with wheat94 using asreml42\n")
test_that("IC_wheat94_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(wheat94.dat)
  
  ### Start with a simple model
  fm0 <- asreml(yield ~ 1,
                random = ~ Variety + Block,
                data = wheat94.dat)
  
  fm0 <- update(fm0)
  
  current.asrt <- as.asrtests(fm0, NULL, NULL, 
                              label = "Simple model", IClikelihood = "full")
  testthat::expect_equal(nrow(current.asrt$wald.tab), 1)
  
  #Add autocorrelation
  current.asrt <- changeTerms(current.asrt, newResidual = "ar1(Col):ar1(Row)", 
                              label = "Add autocorrelation", IClikelihood = "full")
  testthat::expect_true(abs(diff(current.asrt$test.summary$AIC)) - 312.1018 < 1e-03)
  
  #Add units term
  current.asrt <- changeTerms(current.asrt, addRandom = "units", 
                              label = "Add units", IClikelihood = "full")
  
  vpar3 <- current.asrt$asreml.obj$vparameters[1:3]
  current.asrt <- iterate(current.asrt)
  testthat::expect_true(all(abs(current.asrt$asreml.obj$vparameters[1:3] - vpar3) < 1e-03))
  testthat::expect_true(abs(current.asrt$asreml.obj$loglik - -1563.459) < 1)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 1)
  
  #Add random Row and Col terms
  current.asrt <- changeTerms(current.asrt, addRandom = "Row + Col", 
                              label = "Add Row + Col", IClikelihood = "full")
  
  current.asrt <- iterate(current.asrt)
  #check that denDf for current model is the same as the number of variance parameters
  testthat::expect_true(
    length(current.asrt$asreml.obj$vparameters.con[!(current.asrt$asreml.obj$vparameters.con 
                                               %in% c("F","S","B"))])  == 
      current.asrt$test.summary$denDF[current.asrt$test.summary$terms == "Add Row + Col"])
  
  #Add fixed lin(Row) and lin(Col) terms
  current.asrt <- changeTerms(current.asrt, addFixed = "lin(Row) + lin(Col)", 
                              label = "Add lin(Row) + lin(Col)", IClikelihood = "full")
  #three fixed parameters?
  testthat::expect_equal(nrow(current.asrt$wald.tab), 
                         current.asrt$test.summary$DF[current.asrt$test.summary$terms == 
                                                        "Add lin(Row) + lin(Col)"])
  
  #Add random spl(Col) term
  current.asrt <- changeTerms(current.asrt, addRandom = "spl(Col)", 
                              label = "Add spl(Col)", 
                              IClikelihood = "full")
  testthat::expect_equal(
    length(current.asrt$asreml.obj$vparameters.con[!(current.asrt$asreml.obj$vparameters.con 
                                               %in% c("F","S","B"))]),  8)
  testthat::expect_equal(
    current.asrt$test.summary$denDF[current.asrt$test.summary$terms == "Add spl(Col)"][1], 7)
  
  #Restart with fixed Rowcode and Colcode covariates, units and autocorrelation
  fm6 <- asreml(yield ~ Rowcode + Colcode,
                random = ~ Variety + Block + units,
                residual = ~ ar1(Col):ar1(Row),
                data = wheat94.dat)
  
  fm6 <- update(fm6)
  
  current.asrt <- as.asrtests(fm6, wald.tab = NULL, 
                              test.summary = current.asrt$test.summary,
                              label = "Basic + Row/Col covariates", IClikelihood = "full")
  testthat::expect_true(tail(current.asrt$test.summary$action,1) == "Starting model")
  
  #Add random Row and Col terms
  current.asrt <- changeTerms(current.asrt, addRandom = "Row + Col", 
                              label = "Add Row + Col", IClikelihood = "full")
  
  current.asrt <- iterate(current.asrt)
  
  #Add fixed lin(Row) and lin(Col) terms
  current.asrt <- changeTerms(current.asrt, addFixed = "lin(Row) + lin(Col)", 
                              label = "Add lin(Row) + lin(Col)", IClikelihood = "full")
  
  #Add random spl(Col) term
  current.asrt <- changeTerms(current.asrt, addRandom = "spl(Col)", 
                              label = "Add spl(Col)", IClikelihood = "full")
  
  current.asrt <- iterate(current.asrt)
  testthat::expect_true(nrow(current.asrt$test.summary) %in% c(13,14))
  print(current.asrt$test.summary, omit.columns = "p")
  
  
  #Start with Maximal model
  fm.max <- asreml(yield ~ lin(Row) + lin(Col) + Rowcode + Colcode,
                   random = ~ Variety + Block + Row + spl(Col) + Col + units,
                   residual = ~ ar1(Col):ar1(Row),
                   data = wheat94.dat)
  
  current.asrt <- as.asrtests(fm.max, NULL, NULL, 
                              label = "Maximal model", IClikelihood = "full")
  current.asrt <- iterate(current.asrt)
  testthat::expect_true(tail(current.asrt$test.summary$action,1) == "Starting model")
  testthat::expect_equal(current.asrt$test.summary$DF, 7)
  testthat::expect_equal(current.asrt$test.summary$denDF, 8)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 9) #includes bound Block
  
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 
                         current.asrt$test.summary$denDF[1])
  
  #Drop random Row and Col terms
  current.asrt <- changeTerms(current.asrt, dropRandom = "Row + Col", 
                              label = "Drop Row + Col", IClikelihood = "full")
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 
                         current.asrt$test.summary$denDF[3])
  
  #Drop random spl(Col) term
  current.asrt <- changeTerms(current.asrt, dropRandom = "spl(Col)", 
                              label = "Drop spl(Col)", IClikelihood = "full")
  testthat::expect_equal(
    length(current.asrt$asreml.ob$vparameters.con[!(current.asrt$asreml.obj$vparameters.con 
                                               %in% c("F","S","B"))]), 
    current.asrt$test.summary$denDF[4])
  testthat::expect_true((abs(diff(current.asrt$test.summary$BIC[3:4])) - 4.062308) < 1e-05)
  
  #Use hypothesis testing with the maximal model
  current.asrt <- as.asrtests(fm.max, NULL, test.summary = current.asrt$test.summary, 
                              label = "Maximal model", IClikelihood = "full")
  current.asrt <- iterate(current.asrt)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(current.asrt$test.summary), 6)
  
  #Test random Row term
  current.asrt <- testranfix(current.asrt, term = "Row", alpha = 0.20)
  
  #Test random Col term
  current.asrt <- testranfix(current.asrt, term = "Col", alpha = 0.20)
  current.asrt <- iterate(current.asrt)
  
  #test random spl(Col) term
  if (getTestPvalue(current.asrt, label = "Col") > 0.05)
    current.asrt <- testranfix(current.asrt, term = "spl(Col)", alpha = 0.20)
  
  #tests for getTestPvalue
  testthat::expect_true(abs(getTestPvalue(current.asrt, label = "Col") - .5944761) < 1e-04)
  testthat::expect_error(getTestPvalue(current.asrt, label = "Co"))
  
  #Test units term
  current.asrt <- testranfix(current.asrt, term = "units", alpha = 0.20)
  testthat::expect_equal(nrow(current.asrt$test.summary), 10)
  
})

cat("#### Test for glmm ICs with budworm using asreml42\n")
test_that("GLMM_budworm_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  ##  1. the data - the MASS budworm data from function dose.p
  ##     in 'grouped' binomial format
  df <- data.frame(ldose = rep(0:5, 2),
                   numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16),
                   sex = factor(rep(c("M", "F"), c(6, 6))),
                   N=rep(20,12))
  df$numalive <- df$N-df$numdead
  df$p <- df$numdead/df$N
  
  
  as0 <- asreml(p ~ ldose, data=df, family=asr_binomial(total=N))
  as1 <- asreml(p ~ ldose + sex, data=df, family=asr_binomial(total=N))
  
  # asreml AIC agrees with glm
  info <- infoCriteria(list(as0, as1))
  testthat::expect_true(all(abs(info$AIC - c(20.98403, 12.75706)) < 1e-05))
  testthat::expect_true(all(abs(info[1,] - c(2, 0, 0, 20.98403, 21.95385, -8.492016)) < 1e-05))
  #test deviance & AIC diff
  testthat::expect_true(abs(with(info, loglik[1] - loglik[2])*(-2) - 10.22697) < 1e-05)
  testthat::expect_true(abs(with(info, AIC[1] - AIC[2]) - 8.226968) < 1e-05)
  
  ## 2. binary/bernoulli format: 
  # convert number alive and number dead to a set of 1/0 observations at each dose
  df.bin1 <- data.frame(ldose=rep(df$ldose, df$numdead), sex=rep(df$sex, df$numdead),
                        y=rep(rep(1, nrow(df)), df$numdead))
  df.bin2 <- data.frame(ldose=rep(df$ldose, df$numalive), sex=rep(df$sex, df$numalive),
                        y=rep(rep(0, nrow(df)), df$numalive))
  df.bin <- rbind(df.bin1, df.bin2)
  
  # asreml shows AIC lowest for dose model (compared to dose +sex)
  bin.as0 <- asreml(y ~ ldose, data=df.bin, family=asr_binomial())
  bin.as1 <- asreml(y ~ ldose + sex, data=df.bin, family=asr_binomial())
  info <- infoCriteria(list(bin.as0, bin.as1))
  #test deviance & AIC diff
  testthat::expect_true(abs(with(info, loglik[1] - loglik[2])*(-2) - 10.22697) < 1e-05)
  testthat::expect_true(abs(with(info, AIC[1] - AIC[2]) - 8.226968) < 1e-05)
})

cat("#### Test for getFormulae with wheat94 using asreml42\n")
test_that("Formulae_wheat94_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(wheat94.dat)

  fm.max <- asreml(yield ~ lin(Row) + lin(Col) + Rowcode + Colcode,
                   random = ~ Variety + Block + Row + spl(Col) + Col + units,
                   residual = ~ ar1(Col):ar1(Row),
                   data = wheat94.dat)
  
  mod <- getFormulae(fm.max, which = "all", envir = sys.frame(sys.nframe()))
  testthat::expect_true(all(unlist(lapply(mod, function(form) is.null(form) | 
                                            inherits(form, what = "formula")))))
  testthat::expect_true(all(names(mod) == c("fixed", "random", "residual", "sparse")))
  
  #Print fitted model
  testthat::expect_equal(length(mod), 4)
  testthat::expect_true(is.null(mod$sparse))
  
  p <- printFormulae(fm.max, which = "all")
  testthat::expect_true(all(nchar(p) > 11))
  testthat::expect_equal(length(p), 4)
  p <- printFormulae(fm.max, expanded = TRUE)
  testthat::expect_equal(length(p), 3)
  p <- printFormulae(fm.max, which = c("fixed", "random"))
  testthat::expect_equal(length(p), 2)
  p <- printFormulae(fm.max, which = "fixed")
  testthat::expect_equal(length(p), 1)
  
  #Test when have formulae in a character or list
  fix.mod <- mod$fixed
  fm.max <- asreml(fixed = fix.mod,
                   random = mod$random,
                   residual = mod$residual,
                   data = wheat94.dat)
  mod.fm <- getFormulae(fm.max, which = "all", envir = sys.frame(sys.nframe()))
  testthat::expect_true(all(unlist(lapply(mod.fm, function(form) is.null(form) | 
                                            inherits(form, what = "formula")))))
  testthat::expect_equal(length(mod), 4)
  testthat::expect_true(is.null(mod$sparse))
  p <- printFormulae(fm.max, which = "all")
  testthat::expect_true(all(nchar(p) > 11))
  
})

cat("#### Test for R2adj.asreml on Oats with asreml42\n")
test_that("R2adj_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  data(Oats.dat)
  
  #Test a model with no random terms
  m0.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   data=Oats.dat)
  R2.adj <- R2adj.asreml(m0.asr)
  wald.tab <- wald.asreml(m0.asr)
  summary(tab <- aov(Yield ~ Nitrogen*Variety, 
                     data=Oats.dat))
  treatSSqs <- wald.tab[grepl("Nitrogen" ,rownames(wald.tab)) | grepl("Variety" ,rownames(wald.tab)), "Sum of Sq"]
  treatSSq <- sum(treatSSqs)
  totalSSq <- 71 * var(Oats.dat$Yield)
  R2 <- treatSSq/totalSSq*100
  R2adj.calc <- (1 - ((1-(R2/100))*(71/(72-12))))*100
  R2.adj <- R2adj(m0.asr)
  testthat::expect_true(abs(R2.adj - R2adj.calc) < 1e-05)
  R2.adj.N <- R2adj(m0.asr, include.which.fixed = ~ Nitrogen)
  R2.adj.V <- R2adj(m0.asr, include.which.fixed = ~ Variety)
  R2.adj.NV <- R2adj(m0.asr, include.which.fixed = ~ Nitrogen:Variety)
  #Check that sum of R2s for the individual fixed terms equal the overall R2
  testthat::expect_true(abs(R2.adj - (R2.adj.N + R2.adj.V + R2.adj.NV)) < 1e-05)
  #Test two terms
  R2.adj.N_V <- R2adj(m0.asr, include.which.fixed = ~ Nitrogen + Variety)
  testthat::expect_true(abs(R2.adj.N_V - (R2.adj.N + R2.adj.V)) < 1e-05)
  
  #Trying to calculate the adjusted R2 manually - not far off
  R2adj.calc.N <- (1 - ((1-(treatSSqs[1]/(totalSSq)))*(71/(72-3))))*100
  R2adj.calc.V <- (1 - ((1-(treatSSqs[2]/(totalSSq)))*(71/(72-2))))*100
  R2adj.calc.NV <- (1 - ((1-(treatSSqs[3]/(totalSSq)))*(71/(72-6))))*100
  (R2adj.calc.N + R2adj.calc.V + R2adj.calc.NV)
  
  #Test a model with random terms
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  R2.adj.fix <- R2adj(m1.asr)
  testthat::expect_true(abs(R2.adj.fix - 37.18736) < 1e-02)
  R2.adj.ran <- R2adj(m1.asr, include.which.fixed = NULL, include.which.random = ~ .)
  testthat::expect_true(abs(R2.adj.ran - 38.62742) < 1e-02)
  R2.adj <- R2adj(m1.asr, include.which.random = ~ .)
  testthat::expect_true(abs(R2.adj - 75.81478) < 1e-03)
  testthat::expect_true(abs(R2.adj - (R2.adj.fix + R2.adj.ran)) < 1e-05)
})
