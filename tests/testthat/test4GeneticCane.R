
cat("#### Test estimateV specials for Cane with asreml4\n")
test_that("HEB25_estimateV_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
### Parana Local 851
### 200 test lines, 40 check plots
### 10 Columns x 24 rows
  data(local851)
  test.specials <- c( "ar1", "ar2", "ar3", "sar","sar2",
                      "ma1", "ma2", "arma", "exp", "gau", 
                      "cor")
  vpar.vals <- c(-0.1274424, 0.0669068, -0.2265849,  -0.1238088, -0.157972,  
                 0.1200424, -0.04090029, 0.1, 3.116814e-08, 2.721824e-08, 
                 0.04134244)
  names(vpar.vals) <- test.specials
  V.el <- c(-14.69756, -14.82033, -12.26416, -13.75973, -13.78844, 
            -12.3181, -11.08994, 42.3259, 10.94163, 10.94113,
            10.41591)
  names(V.el) <- test.specials
  
  ### Model with genetic variance only
  for (func in test.specials[test.specials != "arma"])
  {
    ranform <- as.formula(paste("~ ar1(Col):", func, "(Row)", sep = ""))
    asreml.obj <- asreml(tch ~ Control/Check, 
                         random = ~ Col + Row + New,
                         residual = ranform, 
                         data=site2, 
                         na.action=na.method(y="include", x="include"))
    print(asreml.obj$vparameters[length(asreml.obj$vparameters)])
    testthat::expect_true(abs(asreml.obj$vparameters[length(asreml.obj$vparameters)] - 
                                vpar.vals[func]) < 1e-06)
    V <- estimateV(asreml.obj)
    print(V[2, 1])
    testthat::expect_true(abs(V[2, 1] - V.el[func]) < 1e-04)
  }
  
  ### Because of a bug, cannot yest test corb
  testthat::expect_error(asreml.obj <- asreml(tch ~ Control/Check, 
                                              random = ~ Col + Row + New,
                                              residual = ~ ar1(Col):corb(Row, b = 3), 
                                              data=site2, 
                                              na.action=na.method(y="include", x="include")))
  #summary(asreml.obj)$varcomp 
  
  
  
  ### Model with genetic competition without pedigree
  asreml.options(design = TRUE)
  asreml.obj <- asreml(tch ~ Control/Check, 
                       random=~ Col + Row + str(~New+ grp(neighbour),~us(2):id(201)),
                       residual =~ ar1(Col):ar1(Row),
                       group=list(neighbour=c(11:211)),data=site2,
                       na.action=na.method(y="include", x="include"))
  summary.asreml(asreml.obj)$varcomp
  asreml.obj$vparameters
  V.g <- with(site2, fac.vcmat(Col, asreml.obj$vparameters["Col"]) + 
                     fac.vcmat(Row, asreml.obj$vparameters["Row"]))
  G.g <- kronecker(matrix(asreml.obj$vparameters[c(3,4,4,5)], nrow = 2, ncol = 2), mat.I(201))
  cols <- c(grep("New", colnames(asreml.obj$design), fixed = TRUE), 
            grep("grp", colnames(asreml.obj$design), fixed = TRUE))
  V.g <- V.g + (asreml.obj$design[, cols] %*% G.g %*% t(as.matrix(asreml.obj$design[, cols])))
  ar1C <- fac.ar1mat(site2$Col, asreml.obj$vparameters["Col:Row!Col!cor"])
  ar1R <- fac.ar1mat(site2$Row, asreml.obj$vparameters["Col:Row!Row!cor"])
  V.g <- asreml.obj$sigma2 * (V.g + ar1C * ar1R)
  V <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(V - V.g) < 1e-06))
  
  #Test for grp function  
  asreml.obj <- asreml(tch ~ Control/Check, 
                       random=~ Col + Row + New+ grp(neighbour),
                       residual =~ ar1(Col):ar1(Row),
                       group=list(neighbour=c(11:211)),data=site2,
                       na.action=na.method(y="include", x="include"))
  summary.asreml(asreml.obj)$varcomp
  asreml.obj$vparameters
  V.g <- with(site2, fac.vcmat(Col, asreml.obj$vparameters["Col"]) + 
                fac.vcmat(Row, asreml.obj$vparameters["Row"]) + 
                fac.vcmat(New, asreml.obj$vparameters["New"]))
  G.g <- asreml.obj$vparameters["grp(neighbour)"] * mat.I(201)
  cols <- grep("grp", colnames(asreml.obj$design), fixed = TRUE)
  V.g <- V.g + (asreml.obj$design[, cols] %*% G.g %*% t(as.matrix(asreml.obj$design[, cols])))
  ar1C <- fac.ar1mat(site2$Col, asreml.obj$vparameters["Col:Row!Col!cor"])
  ar1R <- fac.ar1mat(site2$Row, asreml.obj$vparameters["Col:Row!Row!cor"])
  V.g <- asreml.obj$sigma2 * (V.g + ar1C * ar1R)
  V <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(V - V.g) < 1e-03))

  asreml.options(design = FALSE) 
 
})
