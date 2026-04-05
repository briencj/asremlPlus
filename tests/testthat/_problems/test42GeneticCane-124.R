# Extracted from test42GeneticCane.R:124

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
cat("#### Test estimateV specials for Cane with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(dae)
library(asreml)
library(asremlPlus)
data(local851)
test.specials <- c( "ar1", "ar2", "ar3", "sar","sar2",
                      "ma1", "ma2", "arma", "exp", "gau", 
                      "cor")
vpar.vals <- c(-0.1295442, 0.0669068, -2.265849e-01, -1.238088e-01, -1.57972e-01, 
                 1.200424e-01, -4.0854110e-02, 0, 3.116814e-08,  2.721824e-08, 4.134244e-02)
names(vpar.vals) <- test.specials
V.el <- c(-14.6976, -14.82033, -12.26416, -13.75973, -13.78844, 
            -12.3181, -11.08808, 0, 10.94163, 10.94113, 10.41591)
names(V.el) <- test.specials
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
      R.calc <- asreml.obj$sigma2 * kronecker(mat.ar1(colcorr, 10), 
                                              models[[func]](rowcorrs, c(1:24)))
    }
    else
    { 
      rowcorrs <- asreml.obj$vparameters[grepl("Row!cor", names(asreml.obj$vparameters))]
      R.calc <- asreml.obj$sigma2 * kronecker(mat.ar1(colcorr, 10), 
                                              models[[func]](rowcorrs, 24))
    }
    V <- estimateV(asreml.obj, which.matrix = "R")
    testthat::expect_true(all.equal(R.calc, V))
  }
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
      G.calc <- asreml.obj$vparameters["Col:Row"] * kronecker(mat.ar1(colcorr, 10),
                                                              models[[func]](rowcorrs, c(1:24)))
    } else
    { 
      rowcorrs <- asreml.obj$vparameters[grepl("Row!cor", names(asreml.obj$vparameters))]
      G.calc <- asreml.obj$vparameters["Col:Row"] * kronecker(mat.ar1(colcorr, 10), 
                                                              models[[func]](rowcorrs, 24))
    }
    Z.new <- as.matrix(asreml.obj$design[,grepl("New", colnames(asreml.obj$design))])
    D <- diag(sqrt(asreml.obj$vparameters["New"]), nrow = 240, ncol = 240)
    G.calc <- G.calc + D %*% (Z.new %*% t(Z.new)) %*% D
    G.calc <- asreml.obj$sigma2 * G.calc
    V <- estimateV(asreml.obj, which.matrix = "G")
    dimnames(V) <- NULL
    testthat::expect_true(all.equal(G.calc, V))
  }
