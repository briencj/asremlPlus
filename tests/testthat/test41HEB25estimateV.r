asr41.lib <- "D:\\Analyses\\R ASReml4.1" 

cat("#### Test estimateV specials for HEB25 with asreml41\n")
test_that("HEB25_estimateV_asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)

  data(cart.dat)
  
  #us(Treatment):Genotype
  asreml.options(keep.order = TRUE) #required for asreml4 only
  HEB25.asr <- asreml(fixed = Height ~ Smarthouse + Check + Treatment.1 + 
                        Smarthouse:xMainPosn + Smarthouse:xLane + Check:Treatment.1, 
                      random = ~ us(Treatment.1):Genotype.ID + Smarthouse:Zones:Mainplots, 
                      residual = ~idh(Treat.Smarthouse):Zones:Mainplots, 
                      data = cart.dat, na.action=na.method(y="include", x="include"), 
                      maxiter = 1000, trace = FALSE)
  summary(HEB25.asr)$varcomp
  Vus <- estimateV(HEB25.asr)
  Vus[1:10, 1:9]
  
  #corv(Treatment):Genotype
  HEB25.asr <- asreml(fixed = Height ~ Smarthouse + Check + Treatment.1 + 
                        Smarthouse:xMainPosn + Smarthouse:xLane + Check:Treatment.1, 
                      random = ~ corv(Treatment.1):Genotype.ID + Smarthouse:Zones:Mainplots, 
                      residual = ~idh(Treat.Smarthouse):Zones:Mainplots, 
                      data = cart.dat, na.action=na.method(y="include", x="include"), 
                      maxiter = 1000, trace = FALSE)
  summary(HEB25.asr)$varcomp
  V <- estimateV(HEB25.asr)
  V[1:10, 1:9]
  testthat::expect_false(abs(V[1,1] - 7524245) > 10)
  
  #corh(Treatment):Genotype
  HEB25.asr <- asreml(fixed = Height ~ Smarthouse + Check + Treatment.1 + 
                        Smarthouse:xMainPosn + Smarthouse:xLane + Check:Treatment.1, 
                      random = ~ corh(Treatment.1):Genotype.ID + Smarthouse:Zones:Mainplots, 
                      residual = ~idh(Treat.Smarthouse):Zones:Mainplots, 
                      data = cart.dat, na.action=na.method(y="include", x="include"), 
                      maxiter = 1000, trace = FALSE)
  summary(HEB25.asr)$varcomp
  V <- estimateV(HEB25.asr)
  V[1:10, 1:9]
  testthat::expect_false(abs(V[1,1] - 8234025) > 10)
  
  #corgh(Treatment):Genotype
  HEB25.asr <- asreml(fixed = Height ~ Smarthouse + Check + Treatment.1 + 
                        Smarthouse:xMainPosn + Smarthouse:xLane + Check:Treatment.1, 
                      random = ~ corgh(Treatment.1):Genotype.ID + Smarthouse:Zones:Mainplots, 
                      residual = ~idh(Treat.Smarthouse):Zones:Mainplots, 
                      data = cart.dat, na.action=na.method(y="include", x="include"), 
                      maxiter = 1000, trace = FALSE)
  summary(HEB25.asr)$varcomp
  V <- estimateV(HEB25.asr)
  V[1:10, 1:9]
  testthat::expect_false(any(abs(Vus - V) > 1e+03))
  
  #corgh(Treatment):Genotype + diag
  HEB25.asr <- asreml(fixed = Height ~ Smarthouse + Check + Treatment.1 + 
                        Smarthouse:xMainPosn + Smarthouse:xLane + Check:Treatment.1, 
                      random = ~ corgh(Treatment.1):Genotype.ID + Smarthouse:Zones:Mainplots, 
                      residual = ~diag(Treat.Smarthouse):Zones:Mainplots, 
                      data = cart.dat, na.action=na.method(y="include", x="include"), 
                      maxiter = 1000, trace = FALSE)
  summary(HEB25.asr)$varcomp
  Vdiag <- estimateV(HEB25.asr)
  Vdiag[1:10, 1:9]
  testthat::expect_false(any(abs(Vdiag - V) > 1e+03))
  
})