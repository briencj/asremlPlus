
cat("#### Test estimateV str, spl & dev for orange with asreml42\n")
test_that("Orange_estimateV_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  # Orange tree data from asreml examples
  data(orange)
  
  ##Indepedent slope and intercept - with str
  asreml.options(design = TRUE)
  asreml.obj <- asreml(circ ~x, 
                       random= ~ str( ~Tree/x, ~diag(2):id(5)) + spl(x) + spl(x):Tree,
                       knot.points = list(x = c(118,484,664,1004,1231,1372,1582)), 
                       data = orange, maxiter = 30)
  summary(asreml.obj)$varcomp
  G.g <- kronecker(diag(asreml.obj$vparameters[2:3]), mat.I(5))
  V.g <- asreml.obj$design[,8:17] %*% G.g %*% t(as.matrix(asreml.obj$design[,8:17]))
  V.g <- V.g + asreml.obj$vparameters["spl(x)"] * asreml.obj$design[,3:7] %*%
    t(as.matrix(asreml.obj$design[,3:7]))
  V.g <- V.g + asreml.obj$vparameters["spl(x):Tree"] * asreml.obj$design[,18:42] %*%
    t(as.matrix(asreml.obj$design[,18:42]))
  V.g <- asreml.obj$sigma2 * (V.g + mat.I(nrow(orange)))
  V <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(V - V.g) < 1e-06))
  
  asreml.obj <- asreml(circ ~x, 
                       random= ~ str( ~Tree/x, ~idh(2):id(Tree)) + spl(x) + spl(x):Tree,
                       knot.points = list(x = c(118,484,664,1004,1231,1372,1582)), 
                       data = orange, maxiter = 30)
  summary(asreml.obj)$varcomp
  V <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(V - V.g) < 1e-06))
  
  #Add dev
  asreml.obj <- asreml(circ ~ x,
                       random = ~ str( ~Tree/x, ~diag(2):id(5)) + spl(x) + spl(x):Tree + dev(x),
                       knot.points = list(x = c(118,484,664,1004,1231,1372,1582)),
                       data = orange, maxiter=20)
  summary(asreml.obj)$varcomp
  G.g <- kronecker(diag(asreml.obj$vparameters[c("Tree+Tree:x!diag(2)_1","Tree+Tree:x!diag(2)_2")]), 
                   mat.I(5))
  colnos <- match(c(paste0("Tree_", 1:5), paste0("Tree_", 1:5, ":x")), 
                  dimnames(asreml.obj$design)[[2]])
  V.g <- asreml.obj$design[,colnos] %*% G.g %*% t(as.matrix(asreml.obj$design[,colnos]))
  colnos <- sort(match(outer(paste0("spl(x)_", 1:5), paste0(":Tree_", 1:5), paste0), 
                       dimnames(asreml.obj$design)[[2]]))
  V.g <- V.g + asreml.obj$vparameters["spl(x):Tree"] * asreml.obj$design[,colnos] %*%
    t(as.matrix(asreml.obj$design[,colnos]))
  colnos <- grep("dev\\(x\\)", dimnames(asreml.obj$design)[[2]])
  V.g <- V.g + asreml.obj$vparameters["dev(x)"] * asreml.obj$design[,colnos] %*%
    t(as.matrix(asreml.obj$design[,colnos]))
  V.g <- asreml.obj$sigma2 * (V.g + mat.I(nrow(orange)))
  V <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(V - V.g) < 1e-06))
  
  ##Correlated slope and intercept + fixed Season
  asreml.obj <- asreml(circ ~ x + Season,
                       random = ~ str( ~Tree/x, ~us(2,init=c(5.0,-0.01,0.0001)):id(5)) + 
                         spl(x) + spl(x):Tree + dev(x),
                       knot.points = list(x = c(118,484,664,1004,1231,1372,1582)),
                       data = orange, maxiter=20)     
  summary(asreml.obj)$varcomp
  colnos <- grep("Tree\\+Tree", names(asreml.obj$vparameters))
  G.g <- kronecker(matrix(asreml.obj$vparameters[colnos[c(1,2,2,3)]], nrow = 2), mat.I(5))
  colnos <- match(c(paste0("Tree_", 1:5), paste0("Tree_", 1:5, ":x")), 
                  dimnames(asreml.obj$design)[[2]])
  V.g <- asreml.obj$design[,colnos] %*% G.g %*% t(as.matrix(asreml.obj$design[,colnos]))
  colnos <- match(paste0("spl(x)_", 1:5), dimnames(asreml.obj$design)[[2]])
  V.g <- V.g + asreml.obj$vparameters["spl(x)"] * asreml.obj$design[,colnos] %*%
    t(as.matrix(asreml.obj$design[,colnos]))
  colnos <- sort(match(outer(paste0("spl(x)_", 1:5), paste0(":Tree_", 1:5), paste0), 
                       dimnames(asreml.obj$design)[[2]]))
  V.g <- V.g + asreml.obj$vparameters["spl(x):Tree"] * asreml.obj$design[,colnos] %*%
    t(as.matrix(asreml.obj$design[,colnos]))
  V.g <- asreml.obj$sigma2 * (V.g + mat.I(nrow(orange)))
  V <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(V - V.g) < 1e-06))
  
  #random slope
  asreml.obj <- asreml(circ ~x, 
                       random= ~ Tree + Tree:x + spl(x) + dev(x),
                       knot.points = list(x = c(118,484,664,1004,1231,1372,1582)), 
                       data = orange, maxiter=30)
  summary(asreml.obj)$varcomp
  V.g <- asreml.obj$vparameters["Tree"] * asreml.obj$design[,3:7] %*%
    t(as.matrix(asreml.obj$design[,3:7]))
  V.g <- V.g + asreml.obj$vparameters["Tree:x"] * asreml.obj$design[,8:12] %*%
    t(as.matrix(asreml.obj$design[,8:12]))
  V.g <- V.g + asreml.obj$vparameters["dev(x)"] * asreml.obj$design[,18:24] %*%
    t(as.matrix(asreml.obj$design[,18:24]))
  V.g <- asreml.obj$sigma2 * (V.g + mat.I(nrow(orange)))
  testthat::expect_warning(V <- estimateV(asreml.obj))
  testthat::expect_true(all(abs(V - V.g) < 1e-06))
  

  #Overall spline and deviations based on factor - fails because cannot have a factor in dev
  testthat::expect_error(asreml.obj <- 
                           asreml(circ ~x , 
                                  random= ~ str( ~Tree/x, ~diag(2):id(5)) + spl(X) + dev(X),
                                  knot.points = list(x = c(118,484,664,1004,1231,1372,1582)), 
                                  data = orange))

  asreml.options(design = FALSE) 
  
})
  