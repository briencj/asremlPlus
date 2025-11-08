#devtools::test("asremlPlus")
context("model_selection")

cat("#### Test estimateV for chick pea example with asreml42\n")
test_that("chickpea_estimateV_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)

  data(chkpeadat)
  asreml.options(design = TRUE, fail = "soft")
  asreml.obj <- asreml(fixed = Biomass.plant ~ Lines * TRT + Smarthouse/(vLanes + vPos), 
                       random = ~ Smarthouse:Zone + Smarthouse:spl(vLanes), 
                       residual = ~ idv(Smarthouse):ar1(Lane):ar1(Position), 
                       data = chkpeadat, maxit = 13)

  
  #'## estimate with fixed spline - no G terms in V matrix
  Vnospl <- estimateV(asreml.obj, fixed.spline.terms = "Smarthouse:spl(vLanes)")
  
  # Form variance matrix based on estimated variance parameters
  V <- kronecker(mat.ar1(asreml.obj$vparameters[5], length(levels(chkpeadat$Lane))),
                 mat.ar1(asreml.obj$vparameters[6], length(levels(chkpeadat$Position))))
  V <- asreml.obj$vparameters[4] * kronecker(diag(1,2), V)
  testthat::expect_false(any(abs(Vnospl - V) > 1e-08))
  
  #Estimate G and R separately
  Gnospl <- estimateV(asreml.obj, which.matrix = "G", 
                      fixed.spline.terms = "Smarthouse:spl(vLanes)")
  Rnospl <- estimateV(asreml.obj, which.matrix = "R", 
                      fixed.spline.terms = "Smarthouse:spl(vLanes)")
  testthat::expect_false(any(abs(Gnospl + Rnospl - V) > 1e-08))
  
  #'## estimate with random spline
  Vranspl <- estimateV(asreml.obj)
  Zspl <- as.matrix(asreml.obj$design[ , grepl("Smarthouse", colnames(asreml.obj$design)) &
                                         grepl("spl\\(vLanes\\)", colnames(asreml.obj$design))])
  V <- V + asreml.obj$sigma2 * asreml.obj$vparameters[2] * (Zspl %*% t(Zspl))
  testthat::expect_false(any(abs(Vranspl - V) > 1e-08))
  
  #Estimate G and R separately
  Granspl <- estimateV(asreml.obj, which.matrix = "G")
  Rranspl <- estimateV(asreml.obj, which.matrix = "R")
  testthat::expect_false(any(abs(Rranspl - Rnospl) > 1e-06))
  testthat::expect_false(any(abs(Granspl + Rranspl - Vranspl) > 1e-06))
  testthat::expect_false(any(abs(Granspl + Rranspl - V) > 1e-06))
  
  #Estimate with random spline and k = 6
  asreml.obj <- asreml(fixed = Biomass.plant ~ Lines * TRT + Smarthouse/(vLanes + vPos), 
                       random = ~ Smarthouse:spl(vLanes, k = 6), 
                       residual = ~Smarthouse:ar1(Lane):ar1(Position), 
                       data = chkpeadat, trace = FALSE)
  summary(asreml.obj)$varcomp
  V <- kronecker(mat.ar1(asreml.obj$vparameters[3], length(levels(chkpeadat$Lane))),
                 mat.ar1(asreml.obj$vparameters[4], length(levels(chkpeadat$Position))))
  V <- kronecker(diag(1,2), V)
  Vranspl6 <- estimateV(asreml.obj)
  Zspl6 <- as.matrix(asreml.obj$design[ , grepl("Smarthouse", colnames(asreml.obj$design)) &
                                          grepl("spl\\(vLanes\\, k = 6\\)", colnames(asreml.obj$design))])
  V <- asreml.obj$sigma2 * (V + asreml.obj$vparameters[1] * (Zspl6 %*% t(Zspl6)))
  testthat::expect_false(any(abs(Vranspl6 - V) > 1e-08))
  
  #Estimate G and R separately
  Granspl <- estimateV(asreml.obj, which.matrix = "G")
  Rranspl <- estimateV(asreml.obj, which.matrix = "R")
  testthat::expect_false(any(abs(Granspl + Rranspl - V) > 1e-06))
  
  asreml.options(design = FALSE)
  #no residual
  asreml.obj <- asreml(fixed = Biomass.plant ~ Lines * TRT + Smarthouse/(vLanes + vPos), 
                       random = ~ Smarthouse:spl(vLanes), 
                       data = chkpeadat, trace = FALSE)
  Vnospl <- estimateV(asreml.obj, fixed.spline.terms = "Smarthouse:spl(vLanes)")
  V <- diag(asreml.obj$sigma2, nrow = nrow(chkpeadat))
  testthat::expect_false(any(abs(Vnospl - V) > 1e-08))

  #Estimate G and R separately
  Gnospl <- estimateV(asreml.obj, which.matrix = "G", 
                      fixed.spline.terms = "Smarthouse:spl(vLanes)")
  Rnospl <- estimateV(asreml.obj, which.matrix = "R", 
                      fixed.spline.terms = "Smarthouse:spl(vLanes)")
  testthat::expect_false(any(abs(Gnospl + Rnospl - V) > 1e-08))

  asreml.obj <- asreml(fixed = Biomass.plant ~ Lines * TRT + Smarthouse/(vLanes + vPos), 
                       random = ~ Smarthouse:spl(vLanes), 
                       residual = ~ idv(units),
                       data = chkpeadat, trace = FALSE)
  Vnospl <- estimateV(asreml.obj, fixed.spline.terms = "Smarthouse:spl(vLanes)")
  V <- diag(asreml.obj$vparameters[2], nrow = nrow(chkpeadat))
  testthat::expect_false(any(abs(Vnospl - V) > 1e-08))
  
  #single dsum    
  asreml.obj <- asreml(fixed = Biomass.plant ~ Smarthouse + Lines * TRT, 
                       random = ~Smarthouse:Lane, 
                       residual = ~dsum(~ ar1(Lane):ar1(Position) + Lane:ar1(Position) | 
                                          Smarthouse, levels = list(c(1), c(2))),  
                       data = chkpeadat, trace = FALSE)
  Vdsum <- estimateV(asreml.obj)
  gamma.Lane <- asreml.obj$vparameters[1]
  s2.SW <- asreml.obj$vparameters[2]
  rho.lSW <- asreml.obj$vparameters[3]
  rho.pSW <- asreml.obj$vparameters[4]
  lane.ar1SW <- mat.ar1(order=24, rho=rho.lSW) 
  pos.ar1SW <- mat.ar1(order=22, rho=rho.pSW)
  RSW <- s2.SW * mat.dirprod(lane.ar1SW, pos.ar1SW)
  s2.SE <- asreml.obj$vparameters[5]
  rho.pSE <- asreml.obj$vparameters[6]
  pos.ar1SE <- mat.ar1(order=22, rho=rho.pSE)
  RSE <- s2.SE * mat.dirprod(mat.I(24), pos.ar1SE)
  V <- mat.dirsum(list(RSW, RSE))
  #2 Smarthouses = 1056
  V <- gamma.Lane * with(chkpeadat, fac.sumop(fac.combine(list(Smarthouse,Lane)))) + V 
  testthat::expect_false(any(abs(Vdsum - V) > 1e-08))

  #Estimate G and R separately
  Gdsum <- estimateV(asreml.obj, which.matrix = "G")
  Rdsum <- estimateV(asreml.obj, which.matrix = "R")
  testthat::expect_false(any(abs(Gdsum + Rdsum - V) > 1e-08))
  
  #multiple dsum - same as single dsum
  asreml.obj <- asreml(fixed = Biomass.plant ~ Smarthouse + Lines * TRT, 
                       random = ~Smarthouse:Lane, 
                       residual = ~ dsum(~ ar1(Lane):ar1(Position) | Smarthouse, levels = c(1)) +
                         + dsum(~ Lane:ar1(Position) | Smarthouse, levels = c(2)),  
                       data = chkpeadat, trace = FALSE)
  Vdsum <- estimateV(asreml.obj)
  testthat::expect_false(any(abs(Vdsum - V) > 1e-08))
  
  #Estimate G and R separately
  Gdsum <- estimateV(asreml.obj, which.matrix = "G")
  Rdsum <- estimateV(asreml.obj, which.matrix = "R")
  testthat::expect_false(any(abs(Gdsum + Rdsum - V) > 1e-08))

  #Spatial correlation in random model and nuggets in residual model 
  asreml.obj <- asreml(fixed = Biomass.plant ~ Lines * TRT + Smarthouse/(vLanes + vPos), 
                       random = ~ 
                         at(Smarthouse, 'SW'):spl(vLanes) + at(Smarthouse, 'SE'):spl(vLanes) +
                         at(Smarthouse, 'SW'):ar1(Lane):ar1(Position) + 
                         at(Smarthouse, 'SE'):ar1(Lane):ar1(Position), 
                       residual =  ~ dsum(~ units | Smarthouse, levels = list(c('SW'))) +
                         dsum(~ units | Smarthouse, levels = list(c('SE'))), 
                       data = chkpeadat, maxit = 30)
  
  #'## estimate with fixed spline - no G terms in V matrix
  Vnospl <- estimateV(asreml.obj, 
                      fixed.spline.terms = c("at(Smarthouse, 'SW'):spl(vLanes)",  
                                             "at(Smarthouse, 'SE'):spl(vLanes)"))
  
  # Form variance matrix based on estimated variance parameters
  V <- mat.dirsum(list(asreml.obj$vparameters[3] * 
                         kronecker(mat.ar1(asreml.obj$vparameters[4], length(levels(chkpeadat$Lane))),
                                   mat.ar1(asreml.obj$vparameters[5], length(levels(chkpeadat$Position)))),
                       asreml.obj$vparameters[6] * 
                         kronecker(mat.ar1(asreml.obj$vparameters[7], length(levels(chkpeadat$Lane))),
                                   mat.ar1(asreml.obj$vparameters[8], length(levels(chkpeadat$Position))))))
  V <- diag(rep(asreml.obj$vparameters[9:10], each = 528),1056) %*% V
  testthat::expect_true(any(abs(Vnospl - V) > 1e-08))
})


cat("#### Test estimateV for Wheat example with asreml42\n")
test_that("Wheat_estimateV_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)

  data(Wheat.dat)

  #No residual
  asreml.obj <- asreml(fixed = yield ~ Rep + Variety, 
                       random = ~Row, 
                       data = Wheat.dat, trace = FALSE)
  VWheat <- estimateV(asreml.obj)
  s2 <- asreml.obj$sigma2
  gamma.Row <- asreml.obj$vparameters[1]
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) +
    diag(1, nrow=150, ncol=150)
  V <- s2*V
  testthat::expect_false(any(abs(VWheat - V) > 1e-08))
  
  #units residual
  asreml.obj <- asreml(fixed = yield ~ Rep + Variety, 
                       random = ~Row, 
                       residual = ~ units, 
                       data = Wheat.dat, trace = FALSE)
  VWheat <- estimateV(asreml.obj)
  s2 <- asreml.obj$sigma2
  gamma.Row <- asreml.obj$vparameters[1]
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) +
    diag(1, nrow=150, ncol=150)
  V <- s2*V
  testthat::expect_false(any(abs(VWheat - V) > 1e-08))
  
  #units residual
  asreml.obj <- asreml(fixed = yield ~ Rep + Variety, 
                       random = ~Row, 
                       residual = ~ idv(units), 
                       data = Wheat.dat, trace = FALSE)
  VWheat <- estimateV(asreml.obj)
  s2 <- asreml.obj$sigma2
  gamma.Row <- asreml.obj$vparameters[1]
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) +
    diag(asreml.obj$vparameters[2], nrow=150, ncol=150)
  testthat::expect_false(any(abs(VWheat - V) > 1e-08))
  
  #named residual
  asreml.obj <- asreml(fixed = yield ~ Rep + Variety, 
                       random = ~Row, 
                       residual = ~ Row:Column, 
                       data = Wheat.dat, trace = FALSE)
  VWheat <- estimateV(asreml.obj)
  s2 <- asreml.obj$sigma2
  gamma.Row <- asreml.obj$vparameters[1]
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) +
    diag(1, nrow=150, ncol=150)
  V <- s2*V
  testthat::expect_false(any(abs(VWheat - V) > 1e-08))

  #residual with specials  
  asreml.obj <- asreml(fixed = yield ~ Rep + Variety, 
                       random = ~Row + units, 
                       residual = ~ar1(Row):ar1(Column), 
                       data = Wheat.dat, maxit = 25, trace = FALSE)
  VWheat <- estimateV(asreml.obj)
  s2 <- asreml.obj$sigma2
  gamma.Row <- asreml.obj$vparameters[1]
  gamma.unit <- asreml.obj$vparameters[2]
  rho.r <- asreml.obj$vparameters[4]
  rho.c <- asreml.obj$vparameters[5]
  row.ar1 <- mat.ar1(order=10, rho=rho.r)
  col.ar1 <- mat.ar1(order=15, rho=rho.c)
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) +
    gamma.unit * diag(1, nrow=150, ncol=150) +
    mat.dirprod(row.ar1, col.ar1)
  V <- s2*V
  testthat::expect_true(all(abs(VWheat - V) < 1e-08))

  #residual with Row corb  
  asreml.obj <- asreml(fixed = yield ~ Rep + Variety, 
                       random = ~Row, 
                       residual = ~corb(Row, b = 1):ar1(Column), 
                       data = Wheat.dat, maxit = 25, trace = FALSE)
  VWheat <- estimateV(asreml.obj)
  s2 <- asreml.obj$sigma2
  gamma.Row <- asreml.obj$vparameters[1]
  rho.r <- asreml.obj$vparameters[3]
  rho.c <- c(asreml.obj$vparameters[4])
  row.corb <- mat.banded(x = c(1,rho.r), nrow=10, ncol=10)
  col.ar1 <- mat.ar1(order=15, rho=rho.c)
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) +
    mat.dirprod(row.corb, col.ar1)
  V <- s2*V
  testthat::expect_true(all.equal(VWheat, V, tolerance =  1e-05))

  #residual with Col corb  
  asreml.obj <- asreml(fixed = yield ~ Rep + Variety, 
                       random = ~Row, 
                       residual = ~ar1(Row):corb(Column, b = 4), 
                       data = Wheat.dat, maxit = 25, trace = FALSE)
  VWheat <- estimateV(asreml.obj)
  s2 <- asreml.obj$sigma2
  gamma.Row <- asreml.obj$vparameters[1]
  rho.r <- asreml.obj$vparameters[3]
  rho.c <- c(asreml.obj$vparameters[4:7])
  row.ar1 <- mat.ar1(order=10, rho=rho.r)
  col.corb <- mat.banded(x = c(1,rho.c), nrow=15, ncol=15)
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) +
    mat.dirprod(row.ar1, col.corb)
  V <- s2*V
  testthat::expect_false(any(abs(VWheat - V) > 1e-08))

  #residual with two corb  
  asreml.obj <- asreml(fixed = yield ~ Rep + Variety, 
                       random = ~Row, 
                       residual = ~corb(Row, b = 1):corb(Column, b = 4), 
                       data = Wheat.dat, maxit = 25, trace = FALSE)
  VWheat <- estimateV(asreml.obj)
  s2 <- asreml.obj$sigma2
  gamma.Row <- asreml.obj$vparameters[1]
  rho.r <- asreml.obj$vparameters[3]
  rho.c <- c(asreml.obj$vparameters[c(4:7)])
  row.corb <- mat.banded(x=c(1,rho.r), nrow=10, ncol=10)
  col.corb <- mat.banded(x = c(1,rho.c), nrow=15, ncol=15)
  V <- fac.vcmat(Wheat.dat$Row, gamma.Row) +
    mat.dirprod(row.corb, col.corb)
  V <- s2*V
  testthat::expect_false(any(abs(VWheat - V) > 1e-08))
  RWheat <- estimateV(asreml.obj, which.matrix = "R")
  R <- s2*mat.dirprod(row.corb, col.corb)
  testthat::expect_false(any(abs(RWheat - R) > 1e-08))
  
})

cat("#### Test estimateV str for orange with asreml42\n")
test_that("Orange_estimateV_str_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  # Orange tree data from asreml examples
  data(orange)
  
  orange <- within(orange, 
                   {
                     x2 <- x^2
                     cx <- as.numfac(x, center = TRUE)
                     cx2 <- cx^2
                   })
  asreml.options(design = TRUE)
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  ##Independent quadratic terms - with diag
  asreml.obj <- asreml(circ ~x, 
                       random= ~ str( ~Tree/(x + x2), ~diag(3):id(5)),
                       data = orange, maxit = 30)
  (asreml.obj$vparameters)
  
  G.g <- kronecker(diag(asreml.obj$vparameters[1:3]), mat.I(5))
  Z <- as.matrix(asreml.obj$design[,3:17] )
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  G <- estimateV(asreml.obj, which.matrix = "G")
  testthat::expect_true(all(abs(G - V.g) < 1e-06))
  
  #idh
  asreml.obj <- asreml(circ ~x, 
                       random= ~ str( ~Tree/(x + x2), ~idh(3):id(Tree)),
                       data = orange, maxit = 30)
  (asreml.obj$vparameters)
  G.g <- kronecker(diag(asreml.obj$vparameters[1:3]), mat.I(5))
  Z <- as.matrix(asreml.obj$design[,3:17] )
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  G <- estimateV(asreml.obj, which.matrix = "G")
  testthat::expect_true(all(abs(G - V.g) < 1e-06))
  
  ##Correlated slope and intercept + fixed Season
  asreml.obj <- asreml(circ ~ x,
                       random = ~ str( ~Tree/(x + x2), ~us(3, init = c(1.2, 0.8, 1.2, 0.8,0.8,1.2)):id(5)),
                       data = orange, maxit=500)     
  (asreml.obj$vparameters)
  G.g <- kronecker(matrix(asreml.obj$vparameters[c(1,2,4, 2,3,5, 4:6)], nrow = 3), mat.I(5))
  Z.3_17 <- as.matrix(asreml.obj$design[,3:17] )
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  G <- estimateV(asreml.obj, which.matrix = "G", bound.exclusions = NULL)
  testthat::expect_true(all.equal(G, V.g))
  
  ##Correlated slope and intercept - corg
  asreml.obj <- asreml(circ ~ x,
                       random = ~ str( ~Tree/(x + x2), ~corg(3):id(5)),
                       data = orange, maxit=20)   
  (vpar <- asreml.obj$vparameters)
  g <- diag(1, nrow = 3)
  g[1,2] <- g[2,1] <- vpar[2]
  g[1,3] <- g[3,1] <- vpar[3]
  g[3,2] <- g[2,3] <- vpar[4]
  g <- vpar[1]*g
  G.g <- kronecker(matrix(g, nrow = 3), mat.I(5))
  Z <- as.matrix(asreml.obj$design[,3:17] )
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  G <- estimateV(asreml.obj, which.matrix = "G", bound.exclusions = NULL)
  testthat::expect_true(all.equal(G, V.g, tol = 1e-03))
  
  ##Correlated slope and intercept - cor
  asreml.obj <- asreml(circ ~ cx,
                       random = ~ str( ~Tree/(cx + cx2), ~cor(3):id(5)),
                       data = orange, maxit=50)   
  (vpar <- asreml.obj$vparameters)
  g <- matrix(vpar[2], nrow = 3, ncol = 3)
  diag(g) <- 1
  G.g <- vpar[1] * kronecker(matrix(g, nrow = 3), mat.I(5))
  Z <- as.matrix(asreml.obj$design[,3:17] )
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  G <- estimateV(asreml.obj, which.matrix = "G", bound.exclusions = NULL)
  #mulltiply G by 10 because G.param has the initial value 0.1 stored instead of 1, the correlations being bound
  testthat::expect_true(all.equal(G,V.g, tol = 1e-03)) #problem with bound variance 08/2025
  
  ##Correlated slope and intercept - ar1
  asreml.obj <- asreml(circ ~ cx,
                       random = ~ str( ~Tree/(cx + cx2), ~ar1(3):id(5)),
                       data = orange, maxit=20)   
  (vpar <- asreml.obj$vparameters)
  g <- mat.ar1(vpar[2], 3)
  G.g <-vpar[1] * kronecker(matrix(g, nrow = 3), mat.I(5))
  Z <- as.matrix(asreml.obj$design[,3:17] )
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  G <- estimateV(asreml.obj, which.matrix = "G", bound.exclusions = NULL)
  testthat::expect_true(all.equal(G, V.g, tol = 1e-03))
  
  ##Correlated slope and intercept - ar2
  asreml.obj <- asreml(circ ~ x,
                       random = ~ str( ~Tree/(x + x2), ~ar2(3):id(5)),
                       data = orange, maxit=20)   
  (vpar <- asreml.obj$vparameters)
  g <- mat.ar2(vpar[2:3], 3)
  G.g <- vpar[1] * kronecker(matrix(g, nrow = 3), mat.I(5))
  Z <- as.matrix(asreml.obj$design[,3:17] )
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  G <- estimateV(asreml.obj, which.matrix = "G", bound.exclusions = NULL)
  testthat::expect_true(all.equal(G, V.g, tol = 1e-03))
  
  ##Correlated slope and intercept - ma1
  asreml.obj <- asreml(circ ~ x,
                       random = ~ str( ~Tree/(x + x2), ~ma1(3):id(5)),
                       data = orange, maxit=20)   
  (vpar <- asreml.obj$vparameters)
  g <- mat.ma1(vpar[2], 3)
  G.g <- vpar[1] * kronecker(matrix(g, nrow = 3), mat.I(5))
  Z <- as.matrix(asreml.obj$design[,3:17] )
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  G <- estimateV(asreml.obj, which.matrix = "G", bound.exclusions = NULL)
  testthat::expect_true(all.equal(G, V.g, tol = 1e-03))
  
  ##Correlated slope and intercept - sar
  asreml.obj <- asreml(circ ~ x,
                       random = ~ str( ~Tree/(x + x2), ~sar(3):id(5)),
                       data = orange, maxit=20)   
  (vpar <- asreml.obj$vparameters)
  g <- mat.sar(vpar[2], 3)
  G.g <- vpar[1] * kronecker(matrix(g, nrow = 3), mat.I(5))
  Z <- as.matrix(asreml.obj$design[,3:17] )
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  G <- estimateV(asreml.obj, which.matrix = "G", bound.exclusions = NULL)
  testthat::expect_true(all.equal(G, V.g, tol = 1e-03))
  
})


cat("#### Test for corr random models with asreml42\n")
test_that("corr_random_barley_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  asreml.options(design = TRUE)
  
  #This data is for the Durban 2003 barley data cited in Piepho, Boer and Williams (2022)
  data("barley.dat")
  
  #Fit initial model - Row and column random
  current.asr <- do.call(asreml, 
                         list(yield ~ rep + gen, 
                              random = ~ row + col,
                              data=barley.dat, maxit = 50))
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - -484.1135), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random row and col effects")
  init.asrt <- rmboundary(init.asrt)
  
  corr.asrt <- addSpatialModel(init.asrt, spatial.model = "corr", 
                               row.covar = "crow", col.covar = "ccol",
                               row.factor = "row", col.factor = "col")
  asreml.obj <- corr.asrt$asreml.obj
  (asreml.obj$vparameters)
  (asreml.obj$vparameters.con)

  #Test of estimate V
  G <- estimateV(asreml.obj, which.matrix = "G")
  G.g <- asreml.obj$vparameters[1]*kronecker(mat.ar1(asreml.obj$vparameters[2], 16), 
                                             mat.ar1(asreml.obj$vparameters[3], 34))
  Z <- as.matrix(asreml.obj$design[,276:length(colnames(asreml.obj$design))])
  V.g <- asreml.obj$sigma2 *(Z %*% G.g %*% t(Z))
  testthat::expect_true(all.equal(V.g, G, check.attributes = FALSE))
})

