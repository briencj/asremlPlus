#devtools::test("asremlPlus")
context("model_selection")

cat("#### Test for REMLRT with asreml4\n")
test_that("REMLRT_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(Wheat.dat)
  
  # Fit initial model
  m1.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column + units,
                   residual = ~ ar1(Row):ar1(Column), 
                   data=Wheat.dat)
  summary(m1.asr)$varcomp
  info <- infoCriteria(m1.asr)
  testthat::expect_equal(info$DF, 5)
  testthat::expect_lt(abs(info$AIC - 1346.76), 1e-03)
  
  #Fit model without the units term
  m2.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column,
                   residual = ~ ar1(Row):ar1(Column), 
                   data=Wheat.dat)
  summary(m2.asr)$varcomp
  info <- infoCriteria(m2.asr)
  testthat::expect_equal(info$DF, 4)
  testthat::expect_lt(abs(info$AIC - 1352.941), 1e-03)
  test <- REMLRT(m2.asr, m1.asr)
  testthat::expect_lt(abs(test$p - 0.004232946), 1e-03)
  testthat::expect_equal(test$DF, 1)
  test <- REMLRT(m2.asr, m1.asr, DF = 1)
  testthat::expect_lt(abs(test$p - 0.004232946), 1e-03)
  testthat::expect_equal(test$DF, 1)

  m3.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                   random = ~ Row + Column,
                   residual = ~ Row:Column, 
                   data=Wheat.dat)
  summary(m3.asr)$varcomp
  test3 <- REMLRT(m3.asr, m1.asr)
  testthat::expect_lt(abs(test3$p - 2.596812e-13), 1e-03)
  testthat::expect_equal(test3$DF, 2)
  test3 <- REMLRT(m3.asr, m1.asr, DF = 3)
  testthat::expect_lt(abs(test3$p - 1.603828e-12), 1e-03)
  testthat::expect_equal(test3$DF, 3)
})
