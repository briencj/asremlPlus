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
  wald.tab <- wald(m, denDF = "numeric")$Wald
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
test_that("at_testing_asreml4", {
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
  # but cannot have in testranfix term because not in wald.tab rownames
  current.asr <- asreml(fixed = TSP ~ Lane + xPosn + AMF*Genotype*NP + 
                          at(AMF, "AMF_plus"):per.col + (Genotype*NP):at(AMF, "AMF_plus"):per.col,
                        random = ~ spl(xPosn) + Position ,
                        residual = ~ Genotype:idh(NP_AMF):InTreat,
                        keep.order=TRUE, data = dat, 
                        maxiter=50, na.action = na.method(x="include"))
  
  current.asrt <- asrtests(current.asr, NULL, NULL)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(current.asrt$wald.tab),14)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp),7)
  
  t.asrt <- testranfix(current.asrt, term = "Genotype:NP:at(AMF, AMF_plus):per.col", 
                       drop.fix.ns = TRUE)
  t.asrt$wald.tab
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  #Change position of at term in testranfix
  t.asrt <- testranfix(current.asrt, term = "at(AMF, AMF_plus):per.col:Genotype:NP", 
                       drop.fix.ns = TRUE)
  t.asrt$wald.tab
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  
  #Fit model with level index of 2
  current.asr <- asreml(fixed = TSP ~ Lane + xPosn + AMF*Genotype*NP + 
                          at(AMF, "AMF_plus"):per.col + (Genotype*NP):at(AMF, 2):per.col,
                        random = ~ spl(xPosn) + Position ,
                        residual = ~ Genotype:idh(NP_AMF):InTreat,
                        keep.order=TRUE, data = dat, 
                        maxiter=50, na.action = na.method(x="include"))
  
  current.asrt <- asrtests(current.asr, NULL, NULL)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_equal(nrow(current.asrt$wald.tab),14)
  testthat::expect_equal(nrow(summary(current.asrt$asreml.obj)$varcomp),7)
  
  t.asrt <- testranfix(current.asrt, term = "at(AMF, AMF_plus):per.col:Genotype:NP", 
                       drop.fix.ns = TRUE)
  t.asrt$wald.tab
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
  
  #Test for a numeric level that is not the same as the levels index (1:no.levels)
  current.asr <- asreml(fixed = TSP ~ at(Lane, 4) + xPosn + AMF*Genotype*NP + 
                          at(AMF, "AMF_plus"):per.col + (Genotype*NP):at(AMF, c(2)):per.col,
                        random = ~ spl(xPosn) + Position ,
                        residual = ~ Genotype:idh(NP_AMF):InTreat,
                        keep.order=TRUE, data = dat, 
                        maxiter=50, na.action = na.method(x="include"))
  current.asrt <- asrtests(current.asr, NULL, NULL)
  current.asrt <- rmboundary(current.asrt)
  current.asrt$wald.tab
  t.asrt <- testranfix(current.asrt, term = "at(Lane, 8)", 
                       drop.fix.ns = TRUE)
  t.asrt$wald.tab
  testthat::expect_equal(nrow(t.asrt$wald.tab), 13)
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
  current.asrt <- asrtests(asreml.obj, NULL, NULL)
  testthat::expect_equal(nrow(current.asrt$wald.tab), 24)
  
  asreml.options(step.size = 0.0001)
  
  #Single term in at expresion with the level and drop.fix.ns = TRUE -- works
  t.asrt <- testranfix(current.asrt, 
                       term = "at(expt, mtnue10):vrow", 
                       drop.fix.ns = TRUE,
                       dDF.na = "residual", update = FALSE)
  testthat::expect_equal(nrow(t.asrt$wald.tab), 23)
  
  
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
  testthat::expect_false(abs(summary(asreml.obj)$varcomp$component[1] - 4.495472e-02) > 1e-06)
  testthat::expect_true(summary(asreml.obj)$varcomp$bound[1] == "B")
  asreml.obj <- asreml(fixed = yield ~ Rep + vRow + Variety, 
                       random = ~spl(vRow) + units, 
                       residual = ~ar1(Row):ar1(Column), 
                       data = Wheat.dat, trace = FALSE)
  testthat::expect_false(abs(summary(asreml.obj)$varcomp$component[1] - 8.676175e-03) > 1e-06)
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
