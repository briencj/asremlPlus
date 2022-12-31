#devtools::test("asremlPlus")
context("model_selection")

cat("#### Test for chooseModel.data.frame with asreml4\n")
test_that("choose.model.data.frame_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  print(packageVersion("asreml"))
  #'## Load data
  data("Ladybird.dat")
  
  #'## ANOVA of logits
  Ladybird.aov <- aov(logitP ~ Host*Cadavers*Ladybird + Error(Run/Plant), 
                      data=Ladybird.dat)
  summary(Ladybird.aov)
  
  #'## Mixed model analysis of logits 
  m <- asreml(logitP ~ Host*Cadavers*Ladybird, 
              random = ~ Run,
              residual = ~ Run:Plant,
              data = Ladybird.dat)
  testthat::expect_true(all(summary(m)$varcomp$bound == c("B", "P"))) #shows bound Run component
  
  #'### Unconstrain Reps to make the analysis equivalent to ANOVA
  m <- setvarianceterms(m$call, terms = "Run", bounds = "U")
  summary(m)$varcomp #shows negative Run component
  testthat::expect_true(m$vparameters["Run"] < 0)
  
  #'### Use chooseModel.data.frame
  wald.tab <- wald.asreml(m, denDF = "numeric")$Wald
  testthat::expect_equal(nrow(wald.tab), 8)
  
  #'### Choose marginality-compliant model from wald.tab, obtaining marginality using pstructure
  Ladybird.pstr <- pstructure(formula = ~ Host*Cadavers*Ladybird, 
                              data = Ladybird.dat)
  HCL.marg <- marginality(Ladybird.pstr)
  testthat::expect_equal(nrow(HCL.marg), 7)
  sigmod <- chooseModel(wald.tab, terms.marginality = HCL.marg)
  testthat::expect_true(all(unlist(lapply(sigmod$sig.terms, function(term) term)) == 
                              c("Cadavers:Ladybird", "Host")))
  testthat::expect_equal(nrow(sigmod$choose.summary), 5)
  testthat::expect_true(all(names(sigmod$choose.summary) == c("terms", "DF", "denDF", "p", "action")))
  
  #'### Rechoose marginality-compliant model from wald.tab, but omit DF and denDF
  sigmod <- chooseModel(wald.tab, omit.DF = TRUE, terms.marginality = HCL.marg)
  testthat::expect_equal(nrow(sigmod$choose.summary), 5)
  testthat::expect_true(all(names(sigmod$choose.summary) == c("terms", "p", "action")))
  sigmod <- chooseModel(wald.tab, denDF = NA, terms.marginality = HCL.marg)
  testthat::expect_equal(nrow(sigmod$choose.summary), 5)
  testthat::expect_true(all(names(sigmod$choose.summary) == c("terms", "DF", "denDF", "p", "action")))

  #'### Specify the denDF argument in various ways, all of which result in the overwriting of denDF
  sigmod <- chooseModel(wald.tab, denDF = wald.tab$denDF, terms.marginality = HCL.marg)
  testthat::expect_equal(nrow(sigmod$choose.summary), 5)
  testthat::expect_true(all(names(sigmod$choose.summary) == c("terms", "DF", "denDF", "p", "action")))
  
  den.df <- wald.tab$denDF
  sigmod <- chooseModel(wald.tab, denDF = den.df, terms.marginality = HCL.marg)
  testthat::expect_equal(nrow(sigmod$choose.summary), 5)
  testthat::expect_true(all(names(sigmod$choose.summary) == c("terms", "DF", "denDF", "p", "action")))
  
  sigmod <- chooseModel(wald.tab, denDF = 59, terms.marginality = HCL.marg)
  testthat::expect_equal(nrow(sigmod$choose.summary), 5)
  testthat::expect_true(all(names(sigmod$choose.summary) == c("terms", "DF", "denDF", "p", "action")))
  
  new.tab <- wald.tab
  names(new.tab)[match("denDF", names(new.tab))] <- "den.df"
  sigmod <- chooseModel(wald.tab, denDF = 59, terms.marginality = HCL.marg)
  testthat::expect_equal(nrow(sigmod$choose.summary), 5)
  testthat::expect_true(all(names(sigmod$choose.summary) == c("terms", "DF", "denDF", "p", "action")))
  
})

cat("#### Test for chooseModel.asrtests with asreml4\n")
test_that("choose.model.asrtests_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  print(packageVersion("asreml"))
  data(WaterRunoff.dat)
  asreml::asreml.options(keep.order = TRUE)
  current.asr <- do.call("asreml",
                         args = list(fixed = log.Turbidity ~ Benches +
                                       (Sources * (Type + Species)) * Date,
                                     random = ~Benches:MainPlots:SubPlots:spl(xDay),
                                     data = quote(WaterRunoff.dat)))
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  
  #some tests for validWaldTab
  testthat::expect_error(test.wald <- as.asrtests(current.asr, 
                                                  wald.tab = WaterRunoff.dat))
  asrt.wald <- testranfix(current.asrt, term = "Sources:Species", ssType = "conditional")
  testthat::expect_equal(ncol(asrt.wald$wald.tab), 6)
  testthat::expect_true("F.con" %in% colnames(asrt.wald$wald.tab))
  
  terms.treat <- c("Sources", "Type", "Species", 
                   "Sources:Type", "Sources:Species")
  terms <- sapply(terms.treat, 
                  FUN=function(term){paste("Date:",term,sep="")}, 
                  simplify=TRUE)
  terms <- c("Date", terms)
  terms <- unname(terms)
  marginality <-  matrix(c(1,0,0,0,0,0, 1,1,0,0,0,0,  1,0,1,0,0,0, 
                           1,0,1,1,0,0, 1,1,1,0,1,0, 1,1,1,1,1,1), nrow=6)
  rownames(marginality) <- terms
  colnames(marginality) <- terms
  choose <- chooseModel(current.asrt, marginality, denDF="algebraic")
  current.asrt <- choose$asrtests.obj
  sig.terms <- choose$sig.terms
  testthat::expect_equal(length(sig.terms), 2)
  testthat::expect_equal(sig.terms[[1]], "Date:Species")
  testthat::expect_equal(sig.terms[[2]], "Date:Sources")
})

cat("#### Test for testing at terms with asreml4\n")
test_that("at_testing_testranfix_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  print(packageVersion("asreml"))
  #'## Load data
  data("Exp355.Control.dat")
  
  dat <- with(dat, dat[order(Genotype, Salt, NP_AMF, InTreat), ])
  
  #Fit model with quotes around AMF_plus 
  # NB must have quotes for character levels, 
  # but cannot have in testranfix term because not in wald.tab or varcomp rownames
  current.asr <- asreml(fixed = TSP ~ Lane + xPosn + AMF*Genotype*NP + 
                          at(AMF, "AMF_plus"):per.col + (Genotype*NP):at(AMF, "AMF_plus"):per.col,
                        random = ~ spl(xPosn) + Position ,
                        residual = ~ Genotype:idh(NP_AMF):InTreat,
                        keep.order=TRUE, data = dat, 
                        maxiter=50, na.action = na.method(x="include"))
  
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(current.asrt$wald.tab),14)
  testthat::expect_true(all(c("AMF:Genotype:NP", "at(AMF, AMF_plus):per.col", "Genotype:at(AMF, AMF_plus):per.col", 
                              "NP:at(AMF, AMF_plus):per.col", "Genotype:NP:at(AMF, AMF_plus):per.col") %in% 
                              rownames(current.asrt$wald.tab)))
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp),7)
  
  t.asrt <- testranfix(current.asrt, term = "Genotype:NP:at(AMF, AMF_plus):per.col", 
                       drop.fix.ns = TRUE)
  t.asrt$wald.tab
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13) 
  testthat::expect_true(all(c("AMF:Genotype:NP", "at(AMF, AMF_plus):per.col", "Genotype:at(AMF, AMF_plus):per.col", 
                              "NP:at(AMF, AMF_plus):per.col") %in% rownames(current.asrt$wald.tab)))
  #Change position of at term in testranfix
  t.asrt <- testranfix(current.asrt, term = "at(AMF, AMF_plus):per.col:Genotype:NP", 
                       drop.fix.ns = TRUE)
  t.asrt$wald.tab
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  testthat::expect_true(all(c("AMF:Genotype:NP", "at(AMF, AMF_plus):per.col", "Genotype:at(AMF, AMF_plus):per.col", 
                              "NP:at(AMF, AMF_plus):per.col") %in% rownames(current.asrt$wald.tab)))
  
  #Fit model with level index of 2
  current.asr <- asreml(fixed = TSP ~ Lane + xPosn + AMF*Genotype*NP + 
                          at(AMF, "AMF_plus"):per.col + (Genotype*NP):at(AMF, 2):per.col,
                        random = ~ spl(xPosn) + Position ,
                        residual = ~ Genotype:idh(NP_AMF):InTreat,
                        keep.order=TRUE, data = dat, 
                        maxiter=50, na.action = na.method(x="include"))
  
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(current.asrt$wald.tab),14)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp),7)
  
  t.asrt <- testranfix(current.asrt, term = "at(AMF, AMF_plus):per.col:Genotype:NP", 
                       drop.fix.ns = TRUE)
  t.asrt$wald.tab
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  testthat::expect_true(all(c("AMF:Genotype:NP", "at(AMF, AMF_plus):per.col", "Genotype:at(AMF, AMF_plus):per.col", 
                              "NP:at(AMF, AMF_plus):per.col") %in% rownames(current.asrt$wald.tab)))
  
  #Test for a numeric level that is not the same as the levels index (1:no.levels)
  current.asr <- asreml(fixed = TSP ~ at(Lane, 4) + xPosn + AMF*Genotype*NP + 
                          at(AMF, "AMF_plus"):per.col + (Genotype*NP):at(AMF, c(2)):per.col,
                        random = ~ spl(xPosn) + Position ,
                        residual = ~ Genotype:idh(NP_AMF):InTreat,
                        keep.order=TRUE, data = dat, 
                        maxiter=50, na.action = na.method(x="include"))
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- rmboundary(current.asrt)
  current.asrt$wald.tab
  testthat::expect_error(
    t.asrt <- testranfix(current.asrt, term = "at(Lane, 8)", 
                         drop.fix.ns = TRUE), 
    regexp = "at has numeric values that are more than the number of levels")

  #Test adding multiple terms
  current.asr <- asreml(fixed = TSP ~ Lane + xPosn,
                        random = ~ spl(xPosn) + Position ,
                        residual = ~ Genotype:idh(NP_AMF):InTreat,
                        keep.order=TRUE, data = dat, 
                        maxiter=50, na.action = na.method(x="include"))
  
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(current.asrt$wald.tab),3)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp),7)
  
  full.asrt <- changeTerms(current.asrt, 
                           addFixed = 'AMF*Genotype*NP + at(AMF, "AMF_plus"):per.col + (Genotype*NP):at(AMF, 2):per.col')
  testthat::expect_equal(nrow(full.asrt$wald.tab), 14)
  testthat::expect_true(all(c("AMF:Genotype:NP", "at(AMF, AMF_plus):per.col", "Genotype:per.col:at(AMF, AMF_plus)", 
                              "NP:per.col:at(AMF, AMF_plus)", "Genotype:NP:per.col:at(AMF, AMF_plus)") %in% 
                              rownames(full.asrt$wald.tab)))
  
  full.asrt <- changeTerms(current.asrt, 
                              addFixed = 'AMF*Genotype*NP + at(AMF, "AMF_plus"):per.col + (Genotype*NP):at(AMF, 2):per.col')
  testthat::expect_equal(nrow(full.asrt$wald.tab), 14)
  testthat::expect_true(all(c("AMF:Genotype:NP", "at(AMF, AMF_plus):per.col", "Genotype:per.col:at(AMF, AMF_plus)", 
                              "NP:per.col:at(AMF, AMF_plus)", "Genotype:NP:per.col:at(AMF, AMF_plus)") %in% 
                              rownames(full.asrt$wald.tab)))

  #Try different specification of the at level for removing a fixed term that had a level when added
  t.asrt <- changeTerms(full.asrt, dropFixed = 'at(AMF, "AMF_plus"):per.col')
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  testthat::expect_true(!("at(AMF, AMF_plus):per.col)" %in% rownames(t.asrt$wald.tab)))

  t.asrt <- changeTerms(full.asrt, dropFixed = 'at(AMF, AMF_plus):per.col')
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  testthat::expect_true(!("at(AMF, AMF_plus):per.col)" %in% rownames(t.asrt$wald.tab)))
  
  t.asrt <- changeTerms(full.asrt, dropFixed = 'at(AMF, 2):per.col')
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  testthat::expect_true(!("at(AMF, AMF_plus):per.col)" %in% rownames(t.asrt$wald.tab)))
  
  #Try different specification of the at level for removing a fixed term that had a level index when added
  t.asrt <- changeTerms(full.asrt, dropFixed = 'Genotype:at(AMF, "AMF_plus"):per.col')
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  testthat::expect_true(!("Genotype:per.col:at(AMF, AMF_plus)" %in% rownames(t.asrt$wald.tab)))
  
  t.asrt <- changeTerms(full.asrt, dropFixed = 'Genotype:at(AMF, AMF_plus):per.col')
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  testthat::expect_true(!("Genotype:per.col:at(AMF, AMF_plus)" %in% rownames(t.asrt$wald.tab)))
  
  t.asrt <- changeTerms(full.asrt, dropFixed = 'Genotype:at(AMF, 2):per.col')
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  testthat::expect_true(!("Genotype:per.col:at(AMF, AMF_plus)" %in% rownames(t.asrt$wald.tab)))
  
  #Remove parentheses from (Genotype*NP):at(...)
  t.asrt <- changeTerms(current.asrt, 
                        addFixed = 'AMF*Genotype*NP + at(AMF, "AMF_plus"):per.col + Genotype*NP:at(AMF, 2):per.col')
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  testthat::expect_true(all(c("AMF:Genotype:NP", "at(AMF, AMF_plus):per.col", 
                              "NP:per.col:at(AMF, AMF_plus)", "Genotype:NP:per.col:at(AMF, AMF_plus)") %in% 
                              rownames(t.asrt$wald.tab)))
  
})



cat("#### Test for changeTerms with at functions with asreml4\n")
test_that("at_testing_changeTerms_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  print(packageVersion("asreml"))
  #'## Load data
  data(sPSA.DASTs.dat)


  indv.dat <- with(sPSA.DASTs.dat, sPSA.DASTs.dat[order(SHZone, ZMainunit, Salinity), ])
  
  
  #Set up fixed model and fit initial model
  mod.ch <- "sPSA.2 ~ Smarthouse + at(Smarthouse):cMainPosn + Genotype*Salinity"
  mod <- as.formula(mod.ch)
  cat("\n\n#### ",mod.ch,"\n\n")
  asreml.options(keep.order = TRUE)

  #Fit model with levels and levels indices
  current.asr <- do.call(asreml, 
                         args=list(fixed = mod,
                                   random = ~ at(Smarthouse):Zone/Lane + 
                                     at(Smarthouse, "NE"):spl(cMainPosn) + at(Smarthouse, "NE"):dev(cMainPosn) + 
                                     at(Smarthouse, 2):spl(cMainPosn) + at(Smarthouse, 2):dev(cMainPosn) + 
                                     SHZone:ZMainunit,
                                   residual = ~ SHZone:ZMainunit:Salinity,
                                   data = indv.dat, maxiter=50))
  current.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                              label = "Starting with homogeneous variances model")
  testthat::expect_equal(length(current.asrt$asreml.obj$vparameters), 10)
  testthat::expect_true(all(names(current.asrt$asreml.obj$vparameters) %in% 
                          c("at(Smarthouse, NE):Zone", "at(Smarthouse, NW):Zone", "at(Smarthouse, NE):spl(cMainPosn)", 
                            "at(Smarthouse, NW):spl(cMainPosn)", "at(Smarthouse, NE):dev(cMainPosn)", 
                            "at(Smarthouse, NW):dev(cMainPosn)", "at(Smarthouse, NE):Zone:Lane", 
                            "at(Smarthouse, NW):Zone:Lane", "SHZone:ZMainunit", "SHZone:ZMainunit:Salinity!R")))
  
  #remove boundaries, if any
  t.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 9)
  testthat::expect_false("at(Smarthouse, NE):dev(cMainPosn)" %in% names(t.asrt$asreml.obj$vparameters))
  
  #Test removing the spline term for NE
  t.asrt <- changeTerms(current.asrt, dropRandom = 'at(Smarthouse, "NE"):spl(cMainPosn)', label = "Drop unbound spl")
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 8)
  testthat::expect_true(!any(c("at(Smarthouse, NE):dev(cMainPosn)", "at(Smarthouse, NE):spl(cMainPosn)") %in% 
                          names(t.asrt$asreml.obj$vparameters)))
  
  t.asrt <- changeTerms(current.asrt, dropRandom = 'at(Smarthouse, NE):spl(cMainPosn)', label = "Drop unbound spl")
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 8)
  testthat::expect_true(!any(c("at(Smarthouse, NE):dev(cMainPosn)", "at(Smarthouse, NE):spl(cMainPosn)") %in% 
                               names(t.asrt$asreml.obj$vparameters)))
  
  t.asrt <- changeTerms(current.asrt, dropRandom = 'at(Smarthouse, 1):spl(cMainPosn)', label = "Drop unbound spl")
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 8)
  testthat::expect_true(!any(c("at(Smarthouse, NE):dev(cMainPosn)", "at(Smarthouse, NE):spl(cMainPosn)") %in% 
                               names(t.asrt$asreml.obj$vparameters)))

  #Test removing the devn term for NW
  t.asrt <- changeTerms(current.asrt, dropRandom = 'at(Smarthouse, "NW"):dev(cMainPosn)', label = "Drop unbound dev")
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 8)
  testthat::expect_true(!any(c("at(Smarthouse, NE):dev(cMainPosn)", "at(Smarthouse, NW):dev(cMainPosn)") %in% 
                               names(t.asrt$asreml.obj$vparameters)))
  
  t.asrt <- changeTerms(current.asrt, dropRandom = 'at(Smarthouse, NW):dev(cMainPosn)', label = "Drop unbound dev")
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 8)
  testthat::expect_true(!any(c("at(Smarthouse, NE):dev(cMainPosn)", "at(Smarthouse, NW):dev(cMainPosn)") %in% 
                               names(t.asrt$asreml.obj$vparameters)))
  
  t.asrt <- changeTerms(current.asrt, dropRandom = 'at(Smarthouse, 2):dev(cMainPosn)', label = "Drop unbound dev")
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 8)
  testthat::expect_true(!any(c("at(Smarthouse, NE):dev(cMainPosn)", "at(Smarthouse, NW):dev(cMainPosn)") %in% 
                               names(t.asrt$asreml.obj$vparameters)))
  
  #The next two examples show that the term corresponding to a single level cannot be removed unless it was fitted as a single term
  current.asr <- do.call(asreml, 
                         args=list(fixed = mod,
                                   random = ~ at(Smarthouse):Zone/Lane + 
                                     at(Smarthouse):spl(cMainPosn) + at(Smarthouse):dev(cMainPosn) + 
                                     SHZone:ZMainunit,
                                   residual = ~ SHZone:ZMainunit:Salinity,
                                   data = indv.dat, maxiter=50))
  current.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                              label = "Starting with homogeneous variances model")
  testthat::expect_equal(length(current.asrt$asreml.obj$vparameters), 10)
  testthat::expect_true(all(names(current.asrt$asreml.obj$vparameters) %in% 
                              c("at(Smarthouse, NE):Zone", "at(Smarthouse, NW):Zone", "at(Smarthouse, NE):spl(cMainPosn)", 
                                "at(Smarthouse, NW):spl(cMainPosn)", "at(Smarthouse, NE):dev(cMainPosn)", 
                                "at(Smarthouse, NW):dev(cMainPosn)", "at(Smarthouse, NE):Zone:Lane", 
                                "at(Smarthouse, NW):Zone:Lane", "SHZone:ZMainunit", "SHZone:ZMainunit:Salinity!R")))
  
  #boundary term cannot be removed
  t.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 10)
  testthat::expect_true("at(Smarthouse, NE):dev(cMainPosn)" %in% names(t.asrt$asreml.obj$vparameters))
  
  #Remove whole term
  t.asrt <- changeTerms(current.asrt, dropRandom = "at(Smarthouse):spl(cMainPosn)")
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 8)
  testthat::expect_true(!any(c("at(Smarthouse, NE):spl(cMainPosn)", "at(Smarthouse, NW):spl(cMainPosn)") %in% 
                               names(t.asrt$asreml.obj$vparameters)))

  #Fit with multiple levels indices
  current.asr <- do.call(asreml, 
                         args=list(fixed = mod,
                                   random = ~ at(Smarthouse):Zone/Lane + 
                                     at(Smarthouse, 1:2):spl(cMainPosn) + at(Smarthouse, 1:2):dev(cMainPosn) + 
                                     SHZone:ZMainunit,
                                   residual = ~ SHZone:ZMainunit:Salinity,
                                   data = indv.dat, maxiter=50))
  current.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                              label = "Starting with homogeneous variances model")
  testthat::expect_equal(length(current.asrt$asreml.obj$vparameters), 10)
  testthat::expect_true(all(names(current.asrt$asreml.obj$vparameters) %in% 
                              c("at(Smarthouse, NE):Zone", "at(Smarthouse, NW):Zone", "at(Smarthouse, NE):spl(cMainPosn)", 
                                "at(Smarthouse, NW):spl(cMainPosn)", "at(Smarthouse, NE):dev(cMainPosn)", 
                                "at(Smarthouse, NW):dev(cMainPosn)", "at(Smarthouse, NE):Zone:Lane", 
                                "at(Smarthouse, NW):Zone:Lane", "SHZone:ZMainunit", "SHZone:ZMainunit:Salinity!R")))
  
  #boundary term cannobt br removed
  t.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 10)
  testthat::expect_true("at(Smarthouse, NE):dev(cMainPosn)" %in% names(t.asrt$asreml.obj$vparameters))
  
  #Remove whole term
  t.asrt <- changeTerms(current.asrt, dropRandom = "at(Smarthouse, 1:2):spl(cMainPosn)")
  testthat::expect_equal(length(t.asrt$asreml.obj$vparameters), 8)
  testthat::expect_true(!any(c("at(Smarthouse, NE):spl(cMainPosn)", "at(Smarthouse, NW):spl(cMainPosn)") %in% 
                               names(t.asrt$asreml.obj$vparameters)))
})

cat("#### Test for testing MET at terms with asreml4\n")
test_that("at_multilevel_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  print(packageVersion("asreml"))

  #'## Load data
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
  current.asrt <- as.asrtests(asreml.obj, NULL, NULL)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 24)
  
  asreml.options(step.size = 0.0001)
  
  #Single term in at expresion with the level and drop.fix.ns = TRUE -- works
  t.asrt <- testranfix(current.asrt, 
                       term = "at(expt, mtnue10):vrow", 
                       drop.fix.ns = TRUE,
                       dDF.na = "residual", update = FALSE)
  testthat::expect_equal(nrow(t.asrt$wald.tab), 23)
  testthat::expect_false("at(expt, mtnue10):vrow" %in% rownames(t.asrt$wald.tab))
  
  
  #Multiple fixed terms in an at expresion generates an error
  testthat::expect_error(t.asrt <- testranfix(current.asrt, 
                                              term = "at(expt, c(1:5)):rep", 
                                              drop.fix.ns = TRUE,
                                              dDF.na = "residual", update = FALSE))
  
  #Multiple random terms in an at expression  generates an error
  testthat::expect_error(t.asrt <- testranfix(current.asrt, 
                                              term = "at(expt, c(3,5,7)):dev(vcol)", 
                                              drop.ran.ns = TRUE,
                                              dDF.na = "residual", update = FALSE))

  #Single random term in an at expression - thinks absent
  t.asrt <- testranfix(current.asrt, 
                       term = "at(expt, tarlee13):dev(vcol)", 
                       drop.ran.ns = TRUE,
                       dDF.na = "residual", update = FALSE)
  testthat::expect_equal(t.asrt$test.summary$action[1], "Absent")
  
  #Test multiple at term with changeTerms
  t.asrt <- changeTerms(current.asrt, dropFixed = "at(expt, c(1:5)):rep", update = FALSE)
  testthat::expect_equal(nrow(t.asrt$wald.tab), 19)
  testthat::expect_false(all(grepl("\\:rep", rownames(t.asrt$wald.tab))))
})


cat("#### Test for at terms in testswapran with asreml4\n")
test_that("at_testswapran_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  print(packageVersion("asreml"))
  #'## Load data
  data(longit.dat)
  
  asreml.options(fail = "soft", upsd = TRUE, pxem = 1, step.size = 0.1, ai.sing = TRUE)
  current.asr <- do.call(asreml, 
                         args=list(fixed = Area ~ Block + Treatments + Treatments:xDAP,
                                   random = ~ Block:Cart + at(Treatments):spl(xDAP, k = 10) + Treatments:DAP + 
                                     Block:Cart:spl(xDAP) + Block:Cart:xDAP,
                                   residual = ~ Block:Cart:ar1h(DAP),
                                   keep.order=TRUE, data = longit.dat, maxiter=100))
  
  #'## Function to deal with bound variances - set to 1e-04
  fixBoundResidualVariances <-function(current.asr)
  {
    repeat
    {
      current.call <- current.asr$call
      vpR <- grepl("Block:Cart:DAP!DAP", 
                   names(current.asr$vparameters.con), fixed = TRUE)
      vpR <- current.asr$vparameters.con[vpR]
      (terms <- names(vpR[vpR == 7]))
      if (length(terms) == 0 || length(sum(vpR == 4)) > 5) break
      current.asr <- setvarianceterms(current.call, terms = terms, 
                                      bounds = "F", initial.values = 0.0001,
                                      ignore.suffices = FALSE)
    }
    invisible(current.asr)
  }
  current.asr <- fixBoundResidualVariances(current.asr)
  testthat::expect_true(all(table(summary(current.asr)$varcomp$bound) ==  c(2,46,1)))
  
  #'## Load starting model into an asrtests object
  current.asrt <- as.asrtests(current.asr, NULL, NULL, label = "Selected variance model")
  testthat::expect_true(current.asrt$asreml.obj$converge)
  
  #'### Test for Treatments:DAP deviations terms
  current.asrt <- testranfix(current.asrt, term = "Treatments:DAP",
                             positive.zero = TRUE)
  testthat::expect_equal(current.asrt$test.summary$action[2], "Dropped")
  testthat::expect_true(all(table(summary(current.asrt$asreml.obj)$varcomp$bound) ==  c(2,45,1)))
  
  #'### Test for different curvatures in splines
  current.asrt <- testswapran(current.asrt, oldterms = "at(Treatments):spl(xDAP, k = 10)",
                              newterms = "at(AMF):Zn:spl(xDAP, k = 10)",
                              simpler = TRUE,
                              label = "Heterogeneous Treatment splines")
  testthat::expect_true(getTestPvalue(current.asrt, label = "Heterogeneous Treatment splines") < 0.05)
  testthat::expect_true(names(current.asrt$asreml.obj$vparameters[1]) == 
                          "at(Treatments, -,0):spl(xDAP, k = 10)")
  testthat::expect_true(names(current.asrt$asreml.obj$vparameters[5]) == 
                          "at(Treatments, +,0):spl(xDAP, k = 10)")
  testthat::expect_true(all(abs(current.asrt$asreml.obj$vparameters[1:8] - 
                                  c(233.152932, 502.667930, 74.955973, 
                                    2.540186, 42.003197, 61.206138, 
                                    32.367734, 36.902978)) < 1e-02))
})


cat("#### Test for spline testing with asreml4\n")
test_that("spl.asrtests_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  print(packageVersion("asreml"))
  data(WaterRunoff.dat)
  asreml::asreml.options(keep.order = TRUE)
  current.asr <- do.call("asreml", 
                         args = list(fixed = log.Turbidity ~ Benches + 
                                       (Sources * (Type + Species)) * Date, 
                                     random = ~Benches:MainPlots:SubPlots:spl(xDay, k = 6), 
                                     data = WaterRunoff.dat))
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  
  #Test random splines
  current.asrt <- testranfix(current.asrt, term = "Benches:MainPlots:SubPlots:spl(xDay, k = 6)")
  current.asrt$test.summary
  
  testthat::expect_equal(nrow(current.asrt$test.summary), 1)
  testthat::expect_true(abs(current.asrt$test.summary$p - 0.08013755) < 1e-06)
  
  data(Wheat.dat)
  #'## Add cubic trend to Row so that spline is not bound
  Wheat.dat <- within(Wheat.dat, 
                      {
                        vRow <- as.numeric(Row)
                        vRow <- vRow - mean(unique(vRow))
                        yield <- yield + 10*vRow + 5 * (vRow^2) + 5 * (vRow^3)
                      })
  
  #'## Fit model using asreml4
  asreml.obj <- asreml(fixed = yield ~ Rep + vRow + Variety, 
                       random = ~spl(vRow, k=6) + units, 
                       residual = ~ar1(Row):ar1(Column), 
                       data = Wheat.dat, trace = FALSE)
  testthat::expect_true(summary(asreml.obj)$varcomp$bound[1] == "B")
  asreml.obj <- asreml(fixed = yield ~ Rep + vRow + Variety, 
                       random = ~spl(vRow) + units, 
                       residual = ~ar1(Row):ar1(Column), 
                       data = Wheat.dat, trace = FALSE)
  testthat::expect_true(summary(asreml.obj)$varcomp$bound[1] == "B")

})


cat("#### Test for reparamSigDevn.asrtests with asreml4\n")
test_that("reparamSigDevn.asrtests_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  data(WaterRunoff.dat)
  asreml::asreml.options(keep.order = TRUE)
  current.asr <- asreml(fixed = log.Turbidity ~ Benches + Sources + Type + Species + 
                          Sources:Type + Sources:Species + Sources:Species:xDay + 
                          Sources:Species:Date, 
                        data = WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  
  #Examine terms that describe just the interactions of Date and the treatment factors
  terms.treat <- c("Sources", "Type", "Species", "Sources:Type", "Sources:Species")
  date.terms <- sapply(terms.treat, 
                       FUN=function(term){paste("Date:",term,sep="")}, 
                       simplify=TRUE)
  date.terms <- c("Date", date.terms)
  date.terms <- unname(date.terms)
  treat.marginality <-  matrix(c(1,0,0,0,0,0, 1,1,0,0,0,0,  1,0,1,0,0,0, 
                                 1,0,1,1,0,0, 1,1,1,0,1,0, 1,1,1,1,1,1), nrow=6)
  rownames(treat.marginality) <- date.terms
  colnames(treat.marginality) <- date.terms
  choose <- chooseModel(current.asrt, treat.marginality, denDF="algebraic")
  current.asrt <- choose$asrtests.obj
  current.asr <- current.asrt$asreml.obj
  sig.date.terms <- choose$sig.terms
  
  #Remove all Date terms left in the fixed model
  terms <- "(Date/(Sources * (Type + Species)))"
  current.asrt <- changeTerms(current.asrt, dropFixed = terms)
  #if there are significant date terms, reparameterize to xDays + spl(xDays) + Date
  if (length(sig.date.terms) != 0)
  { #add lin + spl + devn for each to fixed and random models
    trend.date.terms <- sapply(sig.date.terms, 
                               FUN=function(term){sub("Date","xDay",term)}, 
                               simplify=TRUE)
    trend.date.terms <- paste(trend.date.terms,  collapse=" + ")
    current.asrt <- changeTerms(current.asrt, addFixed=trend.date.terms)
    trend.date.terms <- sapply(sig.date.terms, 
                               FUN=function(term){sub("Date","spl(xDay)",term)}, 
                               simplify=TRUE)
    trend.date.terms <- c(trend.date.terms, sig.date.terms)
    trend.date.terms <- paste(trend.date.terms,  collapse=" + ")
    current.asrt <- changeTerms(current.asrt, addRandom = trend.date.terms)
    current.asrt <- rmboundary.asrtests(current.asrt)
  }
  
  #Now test terms for sig date terms
  spl.terms <- sapply(terms.treat, 
                      FUN=function(term){paste("spl(xDay):",term,sep="")}, 
                      simplify=TRUE)
  spl.terms <- c("spl(xDay)",spl.terms)
  lin.terms <- sapply(terms.treat, 
                      FUN=function(term){paste(term,":xDay",sep="")}, 
                      simplify=TRUE)
  lin.terms <- c("xDay",lin.terms)
  systematic.terms <- c(terms.treat, lin.terms, spl.terms, date.terms)
  systematic.terms <- unname(systematic.terms)
  treat.marginality <-  matrix(c(1,0,0,0,0,0, 1,1,0,0,0,0,  1,0,1,0,0,0, 
                                 1,0,1,1,0,0, 1,1,1,1,1,0, 1,1,1,1,1,1), nrow=6)
  systematic.marginality <- kronecker(matrix(c(1,0,0,0, 1,1,0,0, 
                                               1,1,1,0, 1,1,1,1), nrow=4), 
                                      treat.marginality)
  systematic.marginality <- systematic.marginality[-1, -1]
  rownames(systematic.marginality) <- systematic.terms
  colnames(systematic.marginality) <- systematic.terms
  choose <- chooseModel(current.asrt, systematic.marginality, 
                        denDF="algebraic", pos=TRUE)
  current.asrt <- choose$asrtests.obj
  
  #Check if any deviations are significant and, for those that are, go back to 
  #fixed dates
  current.asrt <- reparamSigDevn(current.asrt, choose$sig.terms, 
                                 trend.num = "xDay", devn.fac = "Date", 
                                 denDF = "algebraic")
  k <- match("Sources:Species:Date",rownames(current.asrt$wald.tab))
  testthat::expect_equal(nrow(current.asrt$wald.tab), 9)
  testthat::expect_equal(nrow(current.asrt$test.summary), 6)
  testthat::expect_true(!is.na(k))
})


cat("#### Test for changeModelOnIC with wheat94 using asreml4\n")
test_that("changeModelOnIC_wheat94_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the 1994 wheat example from Gilmour et al. (1995)
  data(wheat94.dat)
  
  
  #Start with Maximal model
  fm.max <- asreml(yield ~ lin(Row) + lin(Col) + Rowcode + Colcode,
                   random = ~ Variety + Block + Row + spl(Col) + Col + units,
                   residual = ~ ar1(Col):ar1(Row),
                   data = wheat94.dat)
  
  current.asrt <- as.asrtests(fm.max, NULL, NULL, 
                              label = "Maximal model", IClikelihood = "full")
  current.asrt <- iterate(current.asrt)
  testthat::expect_true(tail(current.asrt$test.summary$action,1) == "Starting model")
  testthat::expect_equal(current.asrt$test.summary$DF, 7)
  testthat::expect_equal(current.asrt$test.summary$denDF, 8)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 9) #includes bound Block
  
  
  current.asrt <- changeModelOnIC(current.asrt, dropRandom = "Block",
                                  IClikelihood = "full", checkboundaryonly = TRUE, which.IC="AIC")
  
  
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 
                         current.asrt$test.summary$denDF[1])
  
  #Drop random Row and Col terms
  current.asrt <- changeModelOnIC(current.asrt, dropRandom = "Row + Col", 
                                  label = "Drop Row + Col", 
                                  which.IC = "AIC", IClikelihood = "full")
  testthat::expect_equal(getTestEntry(current.asrt, label = "Drop Row + Col")[["denDF"]], -2)
  testthat::expect_equal(getTestEntry(current.asrt, label = "Drop Row + Col")[["action"]], "Unswapped")

  #Drop random spl(Col) term
  current.asrt <- changeModelOnIC(current.asrt, dropRandom = "spl(Col)", 
                                  label = "Drop spl(Col)", IClikelihood = "full")
  testthat::expect_true(getTestEntry(current.asrt, label = "Drop spl(Col)")[["denDF"]] %in% c(-2,-3))
  testthat::expect_equal(getTestEntry(current.asrt, label = "Drop spl(Col)")[["action"]], "Unswapped")
  testthat::expect_true(abs(getTestEntry(current.asrt, label = "Drop spl(Col)")[["AIC"]] - 6.981351) < 1e-05)

  #Drop random units term
  current.asrt <- changeModelOnIC(current.asrt, dropRandom = "units", 
                                  label = "Drop units", IClikelihood = "full")
  testthat::expect_equal(getTestEntry(current.asrt, label = "Drop units")[["denDF"]], -1)
  testthat::expect_equal(getTestEntry(current.asrt, label = "Drop units")[["action"]], "Unswapped")
  testthat::expect_true(abs(getTestEntry(current.asrt, label = "Drop units")[["AIC"]] - 9.511413) < 1e-05)
  
  mod <- printFormulae(current.asrt$asreml.obj)
  testthat::expect_equal(length(mod), 3)
  
  
  #Use REML likelihood and BIC
  current.asrt <- as.asrtests(fm.max, NULL, label = "Maximal model", 
                              IClikelihood = "REML")
  current.asrt <- iterate(current.asrt)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(current.asrt$test.summary), 2)
  
  #Drop random Row and Col terms
  current.asrt <- changeModelOnIC(current.asrt, dropRandom = "Row + Col", 
                                  label = "Drop Row + Col", 
                                  which.IC = "BIC", IClikelihood = "REML")
  testthat::expect_equal(current.asrt$test.summary$denDF[3], -2)
  testthat::expect_equal(current.asrt$test.summary$action[current.asrt$test.summary$terms == 
                                                            "Drop Row + Col"], "Swapped")
  testthat::expect_true((abs(current.asrt$test.summary$BIC[3]) - 8.598262) < 1e-02)
  
  #Drop random spl(Col) term
  current.asrt <- changeModelOnIC(current.asrt, dropRandom = "spl(Col)", 
                                  label = "Drop spl(Col)", 
                                  which.IC = "BIC", IClikelihood = "REML")
  testthat::expect_true(current.asrt$test.summary$denDF[4] %in% c(-1, -2))
  testthat::expect_equal(current.asrt$test.summary$action[current.asrt$test.summary$terms == 
                                                            "Drop spl(Col)"], "Swapped")

  #Drop random units term
  current.asrt <- changeModelOnIC(current.asrt, dropRandom = "units", 
                                  label = "Drop units", 
                                  which.IC = "BIC", IClikelihood = "REML")
  testthat::expect_equal(current.asrt$test.summary$denDF[5], -1)
  testthat::expect_equal(current.asrt$test.summary$action[current.asrt$test.summary$terms == 
                                                            "Drop units"], "Unswapped")
  
  mod <- printFormulae(current.asrt$asreml.obj)
  testthat::expect_equal(length(mod), 3)
  
})


cat("#### Test for changeModelOnIC example using asreml4\n")
test_that("changeModelOnIC_Example_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(Wheat.dat)
  
  #'## Fit maximal model
  current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                        random = ~ Row + Column + units,
                        residual = ~ ar1(Row):ar1(Column), 
                        data=Wheat.dat)
  current.asr <- update(current.asr)
  current.asrt <- as.asrtests(current.asr, NULL, NULL, 
                              label = "Maximal model", IClikelihood = "full")
  #current.asrt <- rmboundary(current.asrt)
  #testthat::expect_true(current.asrt$asreml.obj$converge)
  testthat::expect_true(current.asrt$test.summary$action[1] == "Starting model")
  testthat::expect_equal(current.asrt$test.summary$DF[1], 31)
  testthat::expect_equal(current.asrt$test.summary$denDF[1], 5)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp), 6)
  
  # Drop both Row and Column
  current.asrt <- changeModelOnIC(current.asrt, 
                                  dropRandom = "Row + Column", label = "Drop Row + Column",
                                  checkboundaryonly = TRUE,
                                  which.IC = "AIC", IClikelihood = "full")
  testthat::expect_true(current.asrt$asreml.obj$converge)
  testthat::expect_equal(current.asrt$test.summary$denDF[2], -1)
  
  # Replace residual with model without Row autocorrelation
  current.asrt <- changeModelOnIC(current.asrt, 
                                  newResidual = "Row:ar1(Column)", 
                                  label="Row autocorrelation",
                                  IClikelihood = "full")
  testthat::expect_true(current.asrt$asreml.obj$converge)
  testthat::expect_equal(current.asrt$test.summary$denDF[3], -2)
  testthat::expect_true((abs(current.asrt$test.summary$AIC[3]) - 21.709898) < 1e-03)
  
  mod <- printFormulae(current.asrt$asreml.obj)
  testthat::expect_equal(length(mod), 3)
  testthat::expect_true(grepl("units", mod[2], fixed = TRUE))
  
})



cat("#### Test for fixedcorrelations using asreml4\n")
test_that("Fixedcorrelations_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  ## use asremlPlus to analyse the wheat (barley) example from section 8.6 of the asreml manual (Butler et al. 2010)
  data(PSA.27.dat)
  
  asreml.options(ai.sing = TRUE)
  m.asr <- do.call(asreml, 
                   args=list(fixed = PSA.27 ~ Lane + Position,
                             residual = ~ ar1(Lane):Position,
                             data = PSA.27.dat, maxiter=50))
  m.asrt <- as.asrtests(m.asr, NULL, NULL, label = "Start with Lane autocorrelation",
                        IClikelihood = "full")
  m.asrt <- rmboundary(m.asrt)
  testthat::expect_true(m.asrt$asreml.obj$converge)
  
  testthat::expect_silent(
    m1.asrt <- changeModelOnIC(m.asrt, addRandom = "units", label = "units", allow.fixedcorrelation = FALSE,
                               IClikelihood = "full"))
  tests<- m1.asrt$test.summary
  testthat::expect_equal(m1.asrt$test.summary$action[2], "Unchanged - fixed correlation")
  testthat::expect_true(is.null(getFormulae(m1.asrt$asreml.obj)$random))
  
  m2.asrt <- changeModelOnIC(m.asrt, addRandom = "units", label = "units", allow.fixedcorrelation = TRUE,
                             IClikelihood = "full")
  testthat::expect_equal(m2.asrt$test.summary$action[2], "Swapped")
  testthat::expect_true(grepl("units", as.character(getFormulae(m2.asrt$asreml.obj)$random)[2], fixed = TRUE))
  summary(m2.asrt$asreml.obj)$varcomp
  testthat::expect_equal(unname(vpc.char(m2.asrt$asreml.obj)["Lane:Position!Lane!cor"]), "F")
  
  m3.asrt <- changeTerms(m.asrt, addRandom = "units", label = "Add units", allow.fixedcorrelation = FALSE)
  testthat::expect_equal(m3.asrt$test.summary$action[2], "Unchanged - fixed correlation")
  testthat::expect_true(is.null(getFormulae(m3.asrt$asreml.obj)$random))
  
  m4.asrt <- changeTerms(m.asrt, addRandom = "units", label = "Add units", allow.fixedcorrelation = TRUE)
  testthat::expect_equal(m4.asrt$test.summary$action[2], "Changed random")
  testthat::expect_true(grepl("units", as.character(getFormulae(m4.asrt$asreml.obj)$random)[2], fixed = TRUE))
  
  m4.asrt <- testranfix(m4.asrt, term = "units", positive.zero = TRUE, allow.fixedcorrelation = TRUE)
  testthat::expect_equal(m4.asrt$test.summary$action[3], "Retained")
  testthat::expect_true(grepl("units", as.character(getFormulae(m4.asrt$asreml.obj)$random)[2], fixed = TRUE))
  testthat::expect_equal(unname(vpc.char(m4.asrt$asreml.obj)["Lane:Position!Lane!cor"]), "F")

  m5.asrt <- testranfix(m4.asrt, term = "units", positive.zero = TRUE, allow.fixedcorrelation = TRUE,
                        IClikelihood = "REML")
  testthat::expect_equal(m5.asrt$test.summary$action[4], "Retained")
  testthat::expect_true(grepl("units", as.character(getFormulae(m5.asrt$asreml.obj)$random)[2], fixed = TRUE))
  testthat::expect_equal(unname(vpc.char(m5.asrt$asreml.obj)["Lane:Position!Lane!cor"]), "F")
  testthat::expect_true(all(abs(c(m5.asrt$test.summary$AIC[4],m5.asrt$test.summary$BIC[4]) - 
                                  c(2352.823, 2361.365)) < 1e-04))
  
  m6.asrt <- testranfix(m4.asrt, term = "Lane", allow.fixedcorrelation = TRUE,
                        IClikelihood = "REML")
  testthat::expect_equal(m6.asrt$test.summary$action[4], "Significant")
  
  testthat::expect_warning(
    m6.asrt <- testranfix(m4.asrt, term = "Lane", allow.fixedcorrelation = FALSE,
                          IClikelihood = "REML"), 
    regexp = "The estimated value of one or more correlations in the supplied asreml fit for PSA.27 is fixed")
  testthat::expect_equal(m6.asrt$test.summary$action[4], "Significant")

  #The fixed correlation is in m4.asrt and do not know how to remove it.
  testthat::expect_warning(
    m6.asrt <- testranfix(m4.asrt, term = "units", positive.zero = TRUE, allow.fixedcorrelation = FALSE), 
    regexp = "The estimated value of one or more correlations in the supplied asreml fit for PSA.27 is fixed")
  testthat::expect_equal(m6.asrt$test.summary$action[4], "Unchanged - fixed correlation")
  testthat::expect_true(grepl("units", as.character(getFormulae(m6.asrt$asreml.obj)$random)[2], fixed = TRUE))

  #Start with both Lane and Position autocorrelation
  m.asr <- do.call(asreml, 
                   args=list(fixed = PSA.27 ~ Lane + Position,
                             random = ~ units,
                             residual = ~ ar1(Lane):ar1(Position),
                             data = PSA.27.dat, maxiter=75))
  m.asrt <- as.asrtests(m.asr, NULL, NULL, label = "Start with all autocorrelation",
                        IClikelihood = "full")
  m.asrt <- iterate(m.asrt)
  m.asrt <- rmboundary(m.asrt)
  testthat::expect_true(m.asrt$asreml.obj$converge)
  
  m1.asrt <- changeModelOnIC(m.asrt, newResidual = "ar1(Lane):Position", label = "Lane autocorrelation", 
                             allow.fixedcorrelation = FALSE, update = FALSE,
                             IClikelihood = "full")
  testthat::expect_equal(m1.asrt$test.summary$action[2], "Unchanged - fixed correlation")
  testthat::expect_true(grepl("ar1(Lane):ar1(Position)", 
                              as.character(getFormulae(m1.asrt$asreml.obj)$residual)[2], fixed = TRUE))

  m2.asrt <- changeModelOnIC(m.asrt, newResidual = "ar1(Lane):Position", label = "Lane autocorrelation", 
                             allow.fixedcorrelation = TRUE, update = FALSE, 
                             IClikelihood = "full")
  testthat::expect_equal(m2.asrt$test.summary$action[2], "Swapped")
  testthat::expect_true(grepl("ar1(Lane):Position", 
                              as.character(getFormulae(m2.asrt$asreml.obj)$residual)[2], fixed = TRUE))
  
  m3.asrt <- testresidual(m.asrt, terms = "ar1(Lane):Position", label = "Lane autocorrelation", 
                          simpler = TRUE, allow.fixedcorrelation = FALSE, update = FALSE)
  testthat::expect_equal(m3.asrt$test.summary$action[2], "Unchanged - fixed correlation")
  testthat::expect_true(grepl("ar1(Lane):ar1(Position)", 
                              as.character(getFormulae(m1.asrt$asreml.obj)$residual)[2], fixed = TRUE))

  m4.asrt <- testresidual(m.asrt, terms = "ar1(Lane):Position", label = "Lane autocorrelation", 
                          simpler = TRUE, allow.fixedcorrelation = TRUE, update = FALSE)
  testthat::expect_equal(m4.asrt$test.summary$action[2], "Swapped")
  testthat::expect_true(grepl("ar1(Lane):Position", 
                              as.character(getFormulae(m4.asrt$asreml.obj)$residual)[2], fixed = TRUE))
  
  #Check warning message when supplied asreml.obj has a fixed correlation
  testthat::expect_output(testthat::expect_warning(
    m.asr <- do.call(asreml, 
                     args=list(fixed = PSA.27 ~ 1,
                               random = ~ Lane + Position + units,
                               residual = ~ ar1(Lane):Position,
                               data = PSA.27.dat, maxiter=50))))
  testthat::expect_warning(
    m.asrt <- as.asrtests(m.asr, NULL, NULL, label = "Start with all autocorrelation",
                          IClikelihood = "full"))
  m.asrt <- rmboundary(m.asrt)
  testthat::expect_true(m.asrt$asreml.obj$converge)
  
  testthat::expect_warning(
    m1.asrt <- changeModelOnIC(m.asrt, dropRandom = "units", allow.fixedcorrelation = FALSE),
    regexp = "The estimated value of one or more correlations in the supplied asreml fit for PSA.27 is fixed")
  testthat::expect_warning(
    m2.asrt <- testresidual(m.asrt, terms = "ar1(Lane):ar1(Position)", allow.fixedcorrelation = FALSE),
    regexp = "The estimated value of one or more correlations in the supplied asreml fit for PSA.27 is fixed")
  m1.asr <- newfit(m.asr, random. = ~ . - units, allow.fixedcorrelation = TRUE)
  testthat::expect_false(any("units" == rownames(attr(m1.asr$formulae$random, which = "factors"))))
  
  testthat::expect_warning(
    m2.asr <- newfit(m.asr, random. = ~ . - units, allow.fixedcorrelation = FALSE),
    regexp = "The estimated value of one or more correlations in the supplied asreml fit for PSA.27 is fixed")
  testthat::expect_true(any("units" == rownames(attr(m2.asr$formulae$random, which = "factors"))))

  #Test repararmSigDevn
  PSA.27.dat <- within(PSA.27.dat, 
                       {
                         xPosn <- dae::as.numfac(Position)
                         xPosn <- xPosn - mean(unique(xPosn))
                       })
  m.asr <- do.call(asreml, 
                   args=list(fixed = PSA.27 ~ Lane  + xPosn,
                             random = ~ spl(xPosn) + Position + units,
                             residual = ~ ar1(Lane):Position,
                             data = PSA.27.dat, maxiter=75))
  m.asrt <- as.asrtests(m.asr, NULL, NULL, label = "Start with all autocorrelation",
                        IClikelihood = "full")
  asreml.options(ai.sing = TRUE)
  m1.asrt <- reparamSigDevn(m.asrt, terms = "Position", trend.num = "xPosn", devn.fac = "Position", 
                            allow.fixedcorrelation = TRUE, update = FALSE)
  m1.asrt <- iterate(m1.asrt)
  testthat::expect_equal(m1.asrt$test.summary$action[2], "Changed fixed, random")
  testthat::expect_true(unname(vpc.char(m1.asrt$asreml.obj)["Lane:Position!Lane!cor"]) %in% c("F", "B"))
  
  m2.asrt <- reparamSigDevn(m.asrt, terms = "Position", trend.num = "xPosn", devn.fac = "Position", 
                            allow.fixedcorrelation = FALSE, update = FALSE)
  testthat::expect_equal(m2.asrt$test.summary$action[2], "Unchanged - fixed correlation")
  
  #Test testswapran
  m.asr <- do.call(asreml, 
                   args=list(fixed = PSA.27 ~ 1,
                             random = ~ Lane + units,
                             residual = ~ ar1(Lane):Position,
                             data = PSA.27.dat, maxiter=100))
  m.asrt <- as.asrtests(m.asr, NULL, NULL, label = "Start with all autocorrelation",
                        IClikelihood = "full")
  m.asrt <- rmboundary(m.asrt)
  testthat::expect_true(m.asrt$asreml.obj$converge)
  
  m1.asrt <- testswapran(m.asrt, oldterms = "Lane", newterms = "Position", allow.fixedcorrelation = TRUE)
  testthat::expect_equal(m1.asrt$test.summary$action[2], "Rejected")
  testthat::expect_true(unname(vpc.char(m1.asrt$asreml.obj)["Lane:Position!Lane!cor"]) %in% c("F", "B"))
  
  testthat::expect_warning(
    m2.asrt <- testswapran(m.asrt, oldterms = "Lane", newterms = "Position", allow.fixedcorrelation = FALSE),
    regexp = "The estimated value of one or more correlations in the supplied asreml fit for PSA.27 is fixed")

})
