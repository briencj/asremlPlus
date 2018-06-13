#devtools::test("asremlPlus")
context("model_selection3")
asr3.lib <- "D:\\Analyses\\R oldpkg" 

cat("#### Test estimateV not available for version 3\n")
test_that("chickpea_estimateV_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)

  data(chkpeadat)
  asreml.obj <- asreml(fixed = Biomass.plant ~ Lines * TRT + Smarthouse/(vLanes + vPos), 
                       random = ~Smarthouse:Zone + Smarthouse:spl(vLanes), 
                       residual = ~Smarthouse:ar1(Lane):ar1(Position), 
                       data = chkpeadat, trace = FALSE)
  
  #'## estimate with fixed spline - no G terms in V matrix
  testthat::expect_error(Vnospl <- estimateV(asreml.obj, 
                                             fixed.spline.terms = "Smarthouse:spl(vLanes)"))
  
})

