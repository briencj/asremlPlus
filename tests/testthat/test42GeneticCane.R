
cat("#### Test estimateV specials for Cane with asreml42\n")
test_that("HEB25_estimateV_asreml42", {
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
  vpar.vals <- c(-0.1295442, 0.0669068, -2.265849e-01, -1.238088e-01, -1.57972e-01, 
                 1.200424e-01, -4.0854110e-02, 0, 3.116814e-08,  2.721824e-08, 4.134244e-02)
  # vpar.vals <- c(-1.278067e-01, 6.676302e-02, -2.265831e-01, -1.242148e-01, -1.580255e-01, 
  #                1.186700e-01, -4.063130e-02, 0, 1.937006e-07, 1.549264e-07, 4.164659e-02 )
  names(vpar.vals) <- test.specials
  V.el <- c(-14.6976, -14.82033, -12.26416, -13.75973, -13.78844, 
            -12.3181, -11.08808, 0, 10.94163, 10.94113, 10.41591)
  # V.el <- c(-14.72583, -14.82721, -12.17283, -13.78792, -13.80899, 
  #           -12.23504, -11.12781, 0, 10.91292, 10.91283,  10.48848)
  names(V.el) <- test.specials
  
  ### Model with genetic variance only
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft", design = TRUE)
  for (func in test.specials[test.specials != "arma"])
  {
    ranform <- as.formula(paste("~ ar1(Col):", func, "(Row)", sep = ""))
    asreml.obj <- asreml(tch ~ Control/Check, 
                         random = ~ Col + Row + New,
                         residual = ranform, 
                         data=site2, maxit = 30,
                         na.action=na.method(y="include", x="include"))
    print(asreml.obj$vparameters[length(asreml.obj$vparameters)])
    cat("\n",func,": ", asreml.obj$vparameters[length(asreml.obj$vparameters)], " and ",vpar.vals[func],"\n\n")
   if (abs(asreml.obj$vparameters[length(asreml.obj$vparameters)] - 
      vpar.vals[func]) > 1e-03)
     cat("Difference is ", abs(asreml.obj$vparameters[length(asreml.obj$vparameters)] - 
                                 vpar.vals[func]), "\n\n" )
    testthat::expect_true(abs(asreml.obj$vparameters[length(asreml.obj$vparameters)] - 
                                vpar.vals[func]) < 1e-02)
    V <- estimateV(asreml.obj)
    cat("\n",func,": ", V[2, 1], " and ",V.el[func],"\n\n")
    testthat::expect_true(abs(V[2, 1] - V.el[func]) < 0.5)
  }
  
  ### Model with genetic variance only - test R estimated using estimateV
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  models <- list(ar1 = mat.ar1, ar2 = mat.ar2, ar3 = mat.ar3, cor = mat.cor,
                 sar = mat.sar, sar2 = mat.sar2, ma1 = mat.ma1, ma2 = mat.ma2, 
                 exp = mat.exp, gau = mat.gau)
  for (func in names(models))
  {
    ranform <- as.formula(paste("~ ar1(Col):", func, "(Row)", sep = ""))
    asreml.obj <- asreml(tch ~ Control/Check, 
                         random = ~ Col + Row + New,
                         residual = ranform, 
                         data=site2, maxit = 30,
                         na.action=na.method(y="include", x="include"))
    colcorr <- asreml.obj$vparameters[grepl("Col!cor", names(asreml.obj$vparameters))]
    
    if (func %in% c("exp", "gau"))
    { 
      rowcorrs <- asreml.obj$vparameters[grepl("Row!pow", names(asreml.obj$vparameters))]
      R.calc <- asreml.obj$sigma2 * kronecker(mat.ar1(colcorr, 10), models[[func]](rowcorrs, c(1:24)))
    }
    else
    { 
      rowcorrs <- asreml.obj$vparameters[grepl("Row!cor", names(asreml.obj$vparameters))]
      R.calc <- asreml.obj$sigma2 * kronecker(mat.ar1(colcorr, 10), models[[func]](rowcorrs, 24))
    }
    V <- estimateV(asreml.obj, which.matrix = "R")
    testthat::expect_true(all.equal(R.calc, V))
  }

  ### This had a bug, but seems to be working - it is now converging
  testthat::expect_warning(
    asreml.obj <- asreml(tch ~ Control/Check, 
                         random = ~ Col + Row + New,
                         residual = ~ ar1(Col):corb(Row, b = 3), 
                         data=site2, 
                         na.action=na.method(y="include", x="include")),
    regexp = "Some components changed by more than 1% on the last iteration")
  testthat::expect_true(asreml.obj$converge)
  testthat::expect_equal(nrow(summary(asreml.obj)$varcomp), 8)
  colcorr <- asreml.obj$vparameters[grepl("Col!cor", names(asreml.obj$vparameters))]
  rowcorrs <- asreml.obj$vparameters[grepl("Row!cor", names(asreml.obj$vparameters))]
  R.calc <- asreml.obj$sigma2 * kronecker(mat.ar1(colcorr, 10), mat.banded(c(1, rowcorrs), ncol = 24, nrow = 24))
  V <- estimateV(asreml.obj, which.matrix = "R")
  testthat::expect_true(all.equal(R.calc, V))
  
  ### Random model with New + cor - test G estimated using estimateV
  for (func in names(models))
  {
    tmp <- site2 #needed to make site2 local to the for loop
    ranform <- as.formula(paste0("~ New + ar1(Col):", func, "(Row)"))
    asreml.obj <- asreml(tch ~ Control/Check, 
                         random = ranform, 
                         data=tmp, maxit = 30,
                         na.action=na.method(y="include", x="include"))
    colcorr <- asreml.obj$vparameters[grepl("Col!cor", names(asreml.obj$vparameters))]
    
    if (func %in% c("exp", "gau"))
    { 
      rowcorrs <- asreml.obj$vparameters[grepl("Row!pow", names(asreml.obj$vparameters))]
      G.calc <- asreml.obj$vparameters["Col:Row"] * kronecker(mat.ar1(colcorr, 10), models[[func]](rowcorrs, c(1:24)))
    } else
    { 
      rowcorrs <- asreml.obj$vparameters[grepl("Row!cor", names(asreml.obj$vparameters))]
      G.calc <- asreml.obj$vparameters["Col:Row"] * kronecker(mat.ar1(colcorr, 10), models[[func]](rowcorrs, 24))
    }
    Z.new <- as.matrix(asreml.obj$design[,grepl("New", colnames(asreml.obj$design))])
    D <- diag(sqrt(asreml.obj$vparameters["New"]), nrow = 240, ncol = 240)
    G.calc <- G.calc + D %*% (Z.new %*% t(Z.new)) %*% D
    G.calc <- asreml.obj$sigma2 * G.calc
    V <- estimateV(asreml.obj, which.matrix = "G")
    dimnames(V) <- NULL
    testthat::expect_true(all.equal(G.calc, V))
  }
  
  ### Model with genetic competition without pedigree
  asreml.options(design = TRUE)
  asreml.obj <- asreml(tch ~ Control/Check, 
                       random=~ Col + Row + str(~New + grp(neighbour),~us(2):id(201)),
                       residual =~ ar1(Col):ar1(Row),
                       group=list(neighbour=c(11:211)),data=site2,
                       na.action=na.method(y="include", x="include"))
  summary(asreml.obj)$varcomp
  asreml.obj$vparameters
  V.g <- with(site2, fac.vcmat(Col, asreml.obj$vparameters["Col"]) + 
                     fac.vcmat(Row, asreml.obj$vparameters["Row"]))
  G.g <- kronecker(matrix(asreml.obj$vparameters[c(3,4,4,5)], nrow = 2, ncol = 2), mat.I(201))
  cols <- c(grep("^New", colnames(asreml.obj$design)), 
            grep("^grp", colnames(asreml.obj$design)))
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
  summary(asreml.obj)$varcomp
  asreml.obj$vparameters
  V.g <- with(site2, fac.vcmat(Col, asreml.obj$vparameters["Col"]) + 
                fac.vcmat(Row, asreml.obj$vparameters["Row"]) + 
                fac.vcmat(New, asreml.obj$vparameters["New"]))
  G.g <- asreml.obj$vparameters["grp(neighbour)"] * mat.I(201)
  cols <- grep("^grp", colnames(asreml.obj$design))
  V.g <- V.g + (asreml.obj$design[, cols] %*% G.g %*% t(as.matrix(asreml.obj$design[, cols])))
  ar1C <- fac.ar1mat(site2$Col, asreml.obj$vparameters["Col:Row!Col!cor"])
  ar1R <- fac.ar1mat(site2$Row, asreml.obj$vparameters["Col:Row!Row!cor"])
  V.g <- asreml.obj$sigma2 * (V.g + ar1C * ar1R)
  V <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(V - V.g) < 1e-03))

  asreml.options(design = FALSE) 
 
})

