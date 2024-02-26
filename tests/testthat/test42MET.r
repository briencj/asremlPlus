
cat("#### Test estimateV at, fa, rr functions for MET data with asreml42\n")
test_that("MET_estimateV_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  #Function too replace design matrix colnames when computing V
  replaceColnames <- function(des)
  { 
    colnam <- colnames(des)
    colnam <- strsplit(colnam, split = ":")
    colnam <- lapply(colnam, strsplit, split = "_")
    colnam <- lapply(colnam, 
                     function(nam) 
                       paste(sapply(nam, function(fac) fac[1]), collapse = ":"))
    colnames(des) <- colnam
    return(des)
  }
  
  data(MET)
  asreml.options(design = TRUE, keep.order=TRUE)
  #Test at
  asreml.obj <-asreml(fixed = GY.tha ~  + at(expt, c(1:5)):rep + at(expt, c(1)):vrow + 
                        at(expt, c(2,3,6,7)):colblocks + 
                        at(expt, c(1:5,7)):vcol + Genotype*Condition*expt,
                      random = ~  at(expt, c(1)):dev(vrow) + at(expt, c(2)):spl(vcol) +  
                        at(expt, c(3,5,7)):dev(vcol) + at(expt, c(7)):units,
                      data=comb.dat, maxit = 100, workspace = "1Gb")
  
  summary(asreml.obj)$varcomp
  ranterms <- names(asreml.obj$G.param)
  n <- nrow(comb.dat)
  design <- asreml.obj$design
  design <- replaceColnames(design)
  V.g <- matrix(0, nrow = n, ncol = n)
  for (term in ranterms)
  {
    cols <- grep(term, colnames(design), fixed = TRUE)
    V.g <- V.g + asreml.obj$vparameters[term] * (design[, cols] %*% 
                                                   t(as.matrix(design[, cols])))
  }
  V.g <- as.matrix(asreml.obj$sigma2 * (V.g + mat.I(n)))
  Vat <- estimateV(asreml.obj)
  testthat::expect_true(all.equal(Vat, V.g))
  set.daeTolerance(1e-04, 1e-04)
  #All defaults
  R2.adj <- R2adj(asreml.obj)
  testthat::expect_true(all(abs(R2.adj - 86.26415) < 0.01))
  testthat::expect_equal(attr(R2.adj, which = "fixed"), ~.)
  #Specify a set of fixed terms
  set.daeTolerance(1e-04, 1e-04)
  fix.mod <- as.formula(paste("~", paste0("at(expt, '", levels(comb.dat$expt)[1:5], 
                                          "'):rep", collapse = " + ")))
  R2.adj <- R2adj(asreml.obj, orthogonalize = "eigenmethods", 
                  include.which.fixed = fix.mod)
  testthat::expect_true(all(abs(R2.adj - 50.50893) < 0.01))
  fix.mod <- as.formula(paste("~", paste0(c(paste0("at(expt, '", levels(comb.dat$expt)[1], 
                                                   "'):rep"), 
                                            paste0("rep:at(expt, '", levels(comb.dat$expt)[2:5], 
                                                   "')")), collapse = " + ")))
  testthat::expect_equal(attr(R2.adj, which = "fixed"), fix.mod)
  
  ## Testing random term removal
  #No fixed terms and a specified random term
  R2.adj <- R2adj(asreml.obj, include.which.fixed = NULL, 
                  include.which.random = ~ at(expt, 'tcnue10'):dev(vcol))
  testthat::expect_true(all(abs(R2.adj - 0.7028481) < 1e-02))
  testthat::expect_equal(attr(R2.adj, which = "random"), ~at(expt, 'tcnue10'):dev(vcol))
  #No fixed terms and three specified random terms
  R2.adj <- R2adj(asreml.obj, include.which.fixed = NULL, 
                  include.which.random = ~ at(expt, 'tcnue10'):dev(vcol) + 
                    at(expt, 'rsnue11'):dev(vcol) + at(expt, 'tarlee13'):dev(vcol))
  testthat::expect_true(all(abs(R2.adj - 1.019012) < 1e-04))
  testthat::expect_equal(attr(R2.adj, which = "random"), 
                         ~at(expt, 'tcnue10'):dev(vcol)+dev(vcol):at(expt, 'rsnue11')+dev(vcol):at(expt, 'tarlee13'))
   #No fixed terms and minus a single bound random at term using actual level as in coefficients$random
  R2.adj <- R2adj(asreml.obj, include.which.fixed = NULL, 
                  include.which.random = ~ . - at(expt, 'tarlee13'):units)
  testthat::expect_true(all(abs(R2.adj - 2.741433) < 1e-02))
  testthat::expect_true(attr(R2.adj, which = "random") == 
                           ~at(expt, "mtnue10"):dev(vrow) + at(expt, "pnnue10"):spl(vcol) + 
                          at(expt, "tcnue10"):dev(vcol) + 
                          dev(vcol):at(expt, "rsnue11") + dev(vcol):at(expt, "tarlee13"))
  #No fixed terms and minus a single random at term using ordinal level
  R2.adj <- R2adj(asreml.obj, include.which.fixed = NULL, 
                  include.which.random = ~ . - at(expt, c(7)):units) #not removed because ordinal levels
  testthat::expect_true(grepl('at\\(expt, \\"tarlee13\\"\\):units', 
                              as.character(attr(R2.adj, which = "random"))[2])) #term still in formula
  #No fixed terms and all random terms - provides benchmark
  R2.adj <- R2adj(asreml.obj, include.which.fixed = NULL, include.which.random = ~ .) 
  testthat::expect_true(all(abs(R2.adj - 2.741433) < 1e-02))
  testthat::expect_true(attr(R2.adj, which = "random") == ~.)
  #No fixed terms and minus one of a compound random at term using actual level
  R2.adj <- R2adj(asreml.obj, include.which.fixed = NULL, 
                  include.which.random = ~ . - at(expt, 'tarlee13'):dev(vcol))
  testthat::expect_true(all(abs(R2.adj - 2.577726) < 1e-02))
  testthat::expect_true(attr(R2.adj, which = "random") == 
                           ~at(expt, "mtnue10"):dev(vrow) + at(expt, "pnnue10"):spl(vcol) + 
                           at(expt, "tcnue10"):dev(vcol) + dev(vcol):at(expt, "rsnue11") + 
                           at(expt, "tarlee13"):units)

  
  ## Test fixed model term removal
  #No random terms and minus a single fixed term using actual level
  R2.adj <- R2adj(asreml.obj, include.which.fixed = ~ . - at(expt, 'mtnue10'):vrow, 
                  orthogonalize = "eigen")
  testthat::expect_true(all(abs(R2.adj - 86.27269) < 1e-02))
  testthat::expect_false(grepl('at\\(expt, \\"mtnue10\\"\\):vrow', 
                              as.character(attr(R2.adj, which = "fixed"))[2]))
  #No random terms and minus two single fixed terms using actual level from a compound term
  R2.adj <- R2adj(asreml.obj, 
                  include.which.fixed = ~ . - (at(expt, 'pnnue10'):colblocks + 
                                                 at(expt, 'geranium13'):colblocks), 
                  orthogonalize = "eigen")
  testthat::expect_true(all(abs(R2.adj - 63.41818) < 1e-02))
  testthat::expect_false(grepl('at\\(expt, \\"pnnue10\\"\\):colblocks', 
                               as.character(attr(R2.adj, which = "fixed"))[2]))
  testthat::expect_false(grepl('at\\(expt, \\"geranium13\\"\\):colblocks', 
                               as.character(attr(R2.adj, which = "fixed"))[2]))
  testthat::expect_true(grepl('at\\(expt, \\"tcnue10\\"\\):colblocks', 
                               as.character(attr(R2.adj, which = "fixed"))[2]))
  testthat::expect_true(grepl('colblocks:at\\(expt, \\"tarlee13\\"\\)', 
                              as.character(attr(R2.adj, which = "fixed"))[2]))
  
  ## Test fixed model term removal
  fix.mod <- c("at(expt, 'mtnue10'):rep", "at(expt, 'pnnue10'):rep", 
               "at(expt, 'tcnue10'):rep", "at(expt, 'csnue11'):rep", 
               "at(expt, 'rsnue11'):rep", "at(expt, 'mtnue10'):vrow", 
               "at(expt, 'pnnue10'):colblocks", "at(expt, 'tcnue10'):colblocks", 
               "at(expt, 'geranium13'):colblocks", "at(expt, 'tarlee13'):colblocks", 
               "at(expt, 'mtnue10'):vcol","at(expt, 'pnnue10'):vcol", 
               "at(expt, 'tcnue10'):vcol", "at(expt, 'csnue11'):vcol", 
               "at(expt, 'rsnue11'):vcol", "at(expt, 'tarlee13'):vcol")
  
  #Test at with actual levels and single terms in the model
  asreml.obj <-asreml(fixed = GY.tha ~  at(expt, c(1:5)):rep + at(expt, c(1)):vrow + 
                        at(expt, c('pnnue10', "tcnue10", "geranium13", "tarlee13")):colblocks + 
                        at(expt, c(1:5,7)):vcol + Genotype*Condition*expt,
                      random = ~  at(expt, 'mtnue10'):dev(vrow) + at(expt, 'pnnue10'):spl(vcol) + 
                        at(expt, 'tcnue10'):dev(vcol) + at(expt, 'rsnue11'):dev(vcol) +  
                        at(expt, 'tarlee13'):dev(vcol) + at(expt, 'tarlee13'):units,
                      data=comb.dat, maxit = 100, workspace = "1Gb")
  
  summary(asreml.obj)$varcomp
  #No fixed terms and minus a single bound random at term using actual level as in coefficients$random
  R2.adj <- R2adj(asreml.obj, include.which.fixed = NULL, 
                  include.which.random = ~ . - at(expt, 'tarlee13'):units)
  testthat::expect_true(all(abs(R2.adj - 2.741433) < 1e-02))
  testthat::expect_false(grepl('at\\(expt, \\"tarlee13\\"\\):units', 
                              as.character(attr(R2.adj, which = "random"))[2])) #term still in formula
  #No random terms and minus two single fixed terms using actual level from a compound term
  R2.adj <- R2adj(asreml.obj, 
                  include.which.fixed = ~ . - (at(expt, 'pnnue10'):colblocks + 
                                                 at(expt, 'geranium13'):colblocks), 
                  orthogonalize = "eigen")
  testthat::expect_true(all(abs(R2.adj - 63.41818) < 1e-02))
  testthat::expect_false(grepl('at\\(expt, \\"pnnue10\\"\\):colblocks', 
                               as.character(attr(R2.adj, which = "fixed"))[2]))
  testthat::expect_false(grepl('at\\(expt, \\"geranium13\\"\\):colblocks', 
                               as.character(attr(R2.adj, which = "fixed"))[2]))
  testthat::expect_true(grepl('at\\(expt, \\"tcnue10\\"\\):colblocks', 
                              as.character(attr(R2.adj, which = "fixed"))[2]))
  testthat::expect_true(grepl('colblocks:at\\(expt, \\"tarlee13\\"\\)', 
                              as.character(attr(R2.adj, which = "fixed"))[2]))
  
 
  #Test fa
  asreml.obj <-asreml(fixed = GY.tha ~  + at(expt, c(1:5)):rep + at(expt, c(1)):vrow + 
                        at(expt, c(2,3,6,7)):colblocks + 
                        at(expt, c(1:5,7)):vcol + Condition*expt,
                      random = ~  fa(exptCond, k = 2):Genotype + 
                        at(expt, c(1)):dev(vrow) + at(expt, c(2)):spl(vcol) +  
                        at(expt, c(3,5,7)):dev(vcol) + at(expt, c(7)):units,
                      data=comb.dat, maxit = 100, workspace = "1Gb")
  
  summary(asreml.obj)$varcomp
  ranterms <- names(asreml.obj$G.param)
  n <- nrow(comb.dat)
  V.g <- matrix(0, nrow = n, ncol = n)
  # for (term in ranterms[2:7])
  # {
  #   cols <- grep(term, colnames(design), fixed = TRUE)
  #   print(length(cols))
  #   V.g <- V.g + asreml.obj$vparameters[term] * (asreml.obj$design[, cols] %*% 
  #                                                  t(as.matrix(asreml.obj$design[, cols])))
  # }
  # term <- ranterms[1]
  # term <- "fa(exptCond, k = 2)"
  # cols <- grep(term, colnames(asreml.obj$design), fixed = TRUE)[1:1364]
  # vp <- asreml.obj$vparameters[names(asreml.obj$vparameters)[grep(term, 
  #                                                                 names(asreml.obj$vparameters), fixed = TRUE)]]
  # spec.var <- diag(vp[grepl("!var", names(vp), fixed = TRUE)])
  # loads <- matrix(vp[grepl("!fa", names(vp), fixed = TRUE)], ncol = 2)
  # Gfa <- loads %*% t(loads) + spec.var
  # V.g <- V.g + (asreml.obj$design[, cols] %*% kronecker(Gfa, mat.I(62)) %*%
  #                 t(as.matrix(asreml.obj$design[, cols])))
  # V.g <- asreml.obj$sigma2 * (V.g + mat.I(n))
  
  design <- asreml.obj$design
  design <- replaceColnames(design)
  for (term in ranterms[2:7])
  {
    cols <- grep(term, colnames(design), fixed = TRUE)
    V.g <- V.g + asreml.obj$vparameters[term] * (design[, cols] %*% t(as.matrix(design[, cols])))
  }
  term <- ranterms[1]
  cols <- grep(term, colnames(design), fixed = TRUE)[1:1364]
  vp <- asreml.obj$vparameters[names(asreml.obj$vparameters)[grep(term, 
                                           names(asreml.obj$vparameters), fixed = TRUE)]]
  spec.var <- diag(vp[grepl("!var", names(vp), fixed = TRUE)])
  loads <- matrix(vp[grepl("!fa", names(vp), fixed = TRUE)], ncol = 2)
  Gfa <- loads %*% t(loads) + spec.var
  V.g <- V.g + (design[, cols] %*% kronecker(Gfa, mat.I(62)) %*%
                  t(as.matrix(design[, cols])))
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
                      data=comb.dat, maxit = 100, workspace = "1Gb")
  
  summary(asreml.obj)$varcomp
  ranterms <- names(asreml.obj$G.param)
  n <- nrow(comb.dat)
  V.g <- matrix(0, nrow = n, ncol = n)
  design <- asreml.obj$design
  design <- replaceColnames(design)
  for (term in ranterms[2:7])
  {
    cols <- grep(term, colnames(design), fixed = TRUE)
    V.g <- V.g + asreml.obj$vparameters[term] * (design[, cols] %*% t(as.matrix(design[, cols])))
  }
  term <- ranterms[1]
  cols <- grep(term, colnames(design), fixed = TRUE)[1:1364]
  vp <- asreml.obj$vparameters[names(asreml.obj$vparameters)[grep(term, 
                                                                  names(asreml.obj$vparameters), fixed = TRUE)]]
  loads <- matrix(vp[grepl("!fa", names(vp), fixed = TRUE)], ncol = 2)
  Gfa <- loads %*% t(loads)
  V.g <- V.g + (design[, cols] %*% kronecker(Gfa, mat.I(62)) %*%
                  t(as.matrix(design[, cols])))
  V.g <- asreml.obj$sigma2 * (V.g + mat.I(n))
  Vrr <- estimateV(asreml.obj)
  testthat::expect_true(all(abs(Vrr - V.g) < 1e-06))
  
  asreml.options(design = FALSE) 

})





