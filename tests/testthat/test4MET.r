
cat("#### Test estimateV at, fa, rr functions for MET data with asreml4\n")
test_that("MET_estimateV_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(MET)
  asreml.options(design = TRUE, keep.order=TRUE)
  #Test at
  asreml.obj <-asreml(fixed = GY.tha ~  + at(expt, c(1:5)):rep + at(expt, c(1)):vrow + 
                        at(expt, c(2,3,6,7)):colblocks + 
                        at(expt, c(1:5,7)):vcol + Genotype*Condition*expt,
                      random = ~  at(expt, c(1)):dev(vrow) + at(expt, c(2)):spl(vcol) +  
                        at(expt, c(3,5,7)):dev(vcol) + at(expt, c(7)):units,
                      data=comb.dat, maxiter = 100, workspace = "1Gb")
  
  summary(asreml.obj)$varcomp
  ranterms <- names(asreml.obj$G.param)
  n <- nrow(comb.dat)
  V.g <- matrix(0, nrow = n, ncol = n)
  for (term in ranterms)
  {
    cols <- grep(term, colnames(asreml.obj$design), fixed = TRUE)
    V.g <- V.g + asreml.obj$vparameters[term] * (asreml.obj$design[, cols] %*% 
                                                   t(as.matrix(asreml.obj$design[, cols])))
  }
  V.g <- asreml.obj$sigma2 * (V.g + mat.I(n))
  Vat <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(Vat - V.g) < 1e-06))
  
  #Test fa
  asreml.obj <-asreml(fixed = GY.tha ~  + at(expt, c(1:5)):rep + at(expt, c(1)):vrow + 
                        at(expt, c(2,3,6,7)):colblocks + 
                        at(expt, c(1:5,7)):vcol + Condition*expt,
                      random = ~  fa(exptCond, k = 2):Genotype + 
                        at(expt, c(1)):dev(vrow) + at(expt, c(2)):spl(vcol) +  
                        at(expt, c(3,5,7)):dev(vcol) + at(expt, c(7)):units,
                      data=comb.dat, maxiter = 100, workspace = "1Gb")
  
  summary(asreml.obj)$varcomp
  ranterms <- names(asreml.obj$G.param)
  n <- nrow(comb.dat)
  V.g <- matrix(0, nrow = n, ncol = n)
  for (term in ranterms[2:7])
  {
    cols <- grep(term, colnames(asreml.obj$design), fixed = TRUE)
    V.g <- V.g + asreml.obj$vparameters[term] * (asreml.obj$design[, cols] %*% 
                                                   t(as.matrix(asreml.obj$design[, cols])))
  }
  term <- ranterms[1]
  cols <- grep(term, colnames(asreml.obj$design), fixed = TRUE)[1:1364]
  vp <- asreml.obj$vparameters[names(asreml.obj$vparameters)[grep(term, 
                                           names(asreml.obj$vparameters), fixed = TRUE)]]
  spec.var <- diag(vp[grepl("!var", names(vp), fixed = TRUE)])
  loads <- matrix(vp[grepl("!fa", names(vp), fixed = TRUE)], ncol = 2)
  Gfa <- loads %*% t(loads) + spec.var
  V.g <- V.g + (asreml.obj$design[, cols] %*% kronecker(Gfa, mat.I(62)) %*%
                  t(as.matrix(asreml.obj$design[, cols])))
  V.g <- asreml.obj$sigma2 * (V.g + mat.I(n))
  Vfa <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(Vfa - V.g) < 1e-06))

  #Test rr
  asreml.obj <-asreml(fixed = GY.tha ~  + at(expt, c(1:5)):rep + at(expt, c(1)):vrow + 
                        at(expt, c(2,3,6,7)):colblocks + 
                        at(expt, c(1:5,7)):vcol + Condition*expt,
                      random = ~  rr(exptCond, k = 2):Genotype + 
                        at(expt, c(1)):dev(vrow) + at(expt, c(2)):spl(vcol) +  
                        at(expt, c(3,5,7)):dev(vcol) + at(expt, c(7)):units,
                      data=comb.dat, maxiter = 100, workspace = "1Gb")
  
  summary(asreml.obj)$varcomp
  ranterms <- names(asreml.obj$G.param)
  n <- nrow(comb.dat)
  V.g <- matrix(0, nrow = n, ncol = n)
  for (term in ranterms[2:7])
  {
    cols <- grep(term, colnames(asreml.obj$design), fixed = TRUE)
    V.g <- V.g + asreml.obj$vparameters[term] * (asreml.obj$design[, cols] %*% 
                                                   t(as.matrix(asreml.obj$design[, cols])))
  }
  term <- ranterms[1]
  cols <- grep(term, colnames(asreml.obj$design), fixed = TRUE)[1:1364]
  vp <- asreml.obj$vparameters[names(asreml.obj$vparameters)[grep(term, 
                                                                  names(asreml.obj$vparameters), fixed = TRUE)]]
  loads <- matrix(vp[grepl("!fa", names(vp), fixed = TRUE)], ncol = 2)
  Gfa <- loads %*% t(loads)
  V.g <- V.g + (asreml.obj$design[, cols] %*% kronecker(Gfa, mat.I(62)) %*%
                  t(as.matrix(asreml.obj$design[, cols])))
  V.g <- asreml.obj$sigma2 * (V.g + mat.I(n))
  Vrr <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(Vrr - V.g) < 1e-06))
  
})