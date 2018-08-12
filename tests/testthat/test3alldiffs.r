#devtools::test("asremlPlus")
context("prediction_presentation3")
asr3.lib <- "D:\\Analyses\\R oldpkg" 

cat("#### Test for allDifferences.data.frame sort.alldiffs on Oats with asreml3\n")
test_that("allDifferences_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  testthat::expect_equal(length(m1.asr$gammas),3)
  current.asrt <- asrtests(m1.asr)
  
  #Test for as.alldiffs
  Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", 
                              sed=TRUE)$predictions
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  Var.diffs <- as.alldiffs(predictions = Var.pred$pvals, 
                           sed = Var.pred$sed, 
                           tdf = den.df)
  testthat::expect_equal(length(attributes(Var.diffs)),3)
  testthat::expect_null(attr(Var.diffs, which = "sortOrder"))
  
  #Test for allDifferences
  Var.sort.diffs <- allDifferences(predictions = Var.pred$pvals,
                                   classify = "Nitrogen:Variety", 
                                   sed = Var.pred$sed, tdf = den.df, 
                                   sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(as.character(Var.sort.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.sort.diffs$predictions$Variety[2]) == "Golden Rain")
  testthat::expect_equal(length(attr(Var.sort.diffs, which = "sortOrder")),3)
  
  #Test for re-order factors
  Var.reord.diffs <- allDifferences(predictions = Var.pred$pvals,
                              classify = "Variety:Nitrogen", 
                              sed = Var.pred$sed, tdf = den.df)
  testthat::expect_true(as.character(Var.reord.diffs$predictions$Variety[1]) == "Victory" &
                          as.character(Var.reord.diffs$predictions$Variety[2]) == "Victory")
  
  #Test for re-order factors with reorderClassify
  Var.reord.diffs <- reorderClassify(Var.diffs, newclassify = "Variety:Nitrogen")
  testthat::expect_true(as.character(Var.reord.diffs$predictions$Variety[1]) == "Victory" &
                          as.character(Var.reord.diffs$predictions$Variety[2]) == "Victory")
  
  #Test for re-order factors and sort
  Var.both.diffs <- allDifferences(predictions = Var.pred$pvals,
                                    classify = "Variety:Nitrogen", 
                                    sed = Var.pred$sed, tdf = den.df, 
                                    sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(as.character(Var.both.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.both.diffs$predictions$Variety[2]) == "Marvellous")
  Var.both.diffs <- reorderClassify(Var.diffs, newclassify = "Variety:Nitrogen", 
                                    sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(as.character(Var.both.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.both.diffs$predictions$Variety[2]) == "Marvellous")
})


cat("#### Test for sort.alldiffs on Smarthouse with asreml3\n")
test_that("sort.alldiffs_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  library(dae)
  data(Smarthouse.dat)
  
  #Set up without any sorting
  m1.asr <- asreml(y1 ~ Genotype*A*B, 
                   random=~Replicate/Mainplot/Subplot,
                   data=Smarthouse.dat)
  testthat::expect_equal(length(m1.asr$gammas),4)
  current.asrt <- asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m <- current.asrt$asreml.obj
  testthat::expect_equal(length(m$gammas),3)
  
  diffs <- predictPlus(m, classify = "Genotype:A:B", 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_equal(nrow(diffs$predictions),120)
  testthat::expect_equal(ncol(diffs$predictions),8)
  testthat::expect_equal(as.character(diffs$predictions$Genotype[1]),"Axe")
  testthat::expect_equal(length(attributes(diffs)),8)
  testthat::expect_null(attr(diffs, which = "sortOrder"))
  
  #Test reodering of the classify
  diffs.reord <- reorderClassify(diffs, newclassify = "A:B:Genotype")
  testthat::expect_equal(as.character(diffs.reord$predictions$Genotype[1]),"Axe")
  testthat::expect_equal(as.character(diffs.reord$predictions$Genotype[2]),"Espada")
  testthat::expect_true(abs(diffs.reord$predictions$predicted.value[2] - -0.2265723017) < 1e-06)
  
  testthat::expect_silent(plotPredictions(data = diffs$predictions, 
                                          classify = "Genotype:A:B", 
                                          y = "predicted.value", 
                                          error.intervals = "StandardError",  
                                          y.title = attr(diffs, 
                                                         which = "response.title")))
  
  #Test sort in plotPredictions
  testthat::expect_silent(plotPredictions(data = diffs$predictions, 
                                          classify = "Genotype:A:B", 
                                          y = "predicted.value", 
                                          error.intervals = "StandardError",  
                                          y.title = attr(diffs, 
                                                         which = "response.title"),
                                          sortFactor = "Genotype"))
  
  #Test sort.alldiffs and save order for use with other response variables
  diffs.sort <- sort(diffs, sortFactor = "Genotype")
  sort.order <- attr(diffs.sort, which = "sortOrder")
  testthat::expect_is(diffs.sort, "alldiffs")
  testthat::expect_equal(nrow(diffs.sort$predictions),120)
  testthat::expect_equal(ncol(diffs.sort$predictions),8)
  testthat::expect_equal(as.character(diffs.sort$predictions$Genotype[1]),"Gladius")
  testthat::expect_equal(length(attributes(diffs.sort)),10)
  testthat::expect_equal(length(attr(diffs.sort, which = "sortOrder")),10)
  
  #Test sort.alldiffs with supplied sortOrder
  m2.asr <- asreml(y2 ~ Genotype*A*B, 
                   random=~Replicate/Mainplot/Subplot,
                   data=Smarthouse.dat)
  testthat::expect_equal(length(m1.asr$gammas),4)
  current.asrt <- asrtests(m2.asr)
  diffs2.sort <- predictPlus(m2.asr, classify = "Genotype:A:B", 
                             wald.tab = current.asrt$wald.tab,
                             error.intervals = "Stand", tables = "none",
                             sortFactor = "Genotype", 
                             sortOrder = sort.order)
  testthat::expect_equal(as.character(diffs.sort$predictions$Genotype[1]),
                         as.character(diffs2.sort$predictions$Genotype[1]))
  testthat::expect_equal(attr(diffs.sort, which = "sortOrder"),
                         attr(diffs2.sort, which = "sortOrder"))
  
  #Test sort.alldiffs with sortWithinVals and increasing order
  diffs1.sort <- sort(diffs, sortFactor = "Genotype", 
                      sortWithinVals = list(A = "N3", B = "D4"),
                      decreasing = TRUE)
  testthat::expect_is(diffs1.sort, "alldiffs")
  testthat::expect_equal(as.character(diffs1.sort$predictions$Genotype[2]),"Wyalkatchem")
  testthat::expect_equal(length(attr(diffs1.sort, which = "sortOrder")), 10)
  #Check plot
  testthat::expect_silent(plotPredictions(data = diffs1.sort$predictions, 
                                          classify = "Genotype:A:B", 
                                          y = "predicted.value", 
                                          error.intervals = "StandardError",  
                                          y.title = attr(diffs, 
                                                         which = "response.title"),
                                          sortFactor = "Genotype"))
})

cat("#### Test for sort.alldiffs on Oats with asreml3\n")
test_that("sort.alldiffs_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  testthat::expect_equal(length(m1.asr$gammas),3)
  current.asrt <- asrtests(m1.asr)
  
  #Test for as.alldiffs
  Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", 
                      sed=TRUE)$predictions
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  Var.diffs <- as.alldiffs(predictions = Var.pred$pvals, 
                           sed = Var.pred$sed, 
                           tdf = den.df)
  testthat::expect_equal(length(attributes(Var.diffs)),3)
  testthat::expect_null(attr(Var.diffs, which = "sortOrder"))
  
  #Test for allDifferences
  Var.diffs <- allDifferences(predictions = Var.pred$pvals,
                              classify = "Nitrogen:Variety", 
                              sed = Var.pred$sed, tdf = den.df, 
                              sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_equal(as.character(Var.diffs$predictions$Variety[1]),"Marvellous")
  testthat::expect_equal(length(attr(Var.diffs, which = "sortOrder")),3)
  
  
  
  #Test for predictPlus no sorting
  diffs <- predictPlus(m1.asr, classify = "Nitrogen:Variety", 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_equal(nrow(diffs$predictions),12)
  testthat::expect_equal(ncol(diffs$predictions),7)
  testthat::expect_equal(as.character(diffs$predictions$Variety[1]),"Victory")
  testthat::expect_equal(length(attributes(diffs)),8)
  testthat::expect_equal(attr(diffs, which = "tdf"),45)
  testthat::expect_null(attr(diffs, which = "sortOrder"))
  
  testthat::expect_silent(plotPredictions(data = diffs$predictions, 
                                          classify ="Nitrogen:Variety", 
                                          y = "predicted.value", 
                                          error.intervals = "StandardError",  
                                          y.title = attr(diffs, 
                                                         which = "response.title")))
  #TestPlot with sort
  testthat::expect_silent(plotPvalues(diffs, gridspacing = 3, 
                                      sortFactor = "Variety", decreasing = TRUE))
  
  #Test for sort.alldiffs
  diffs.sort <- sort(diffs, sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_equal(as.character(diffs.sort$predictions$Variety[1]),"Marvellous")
  testthat::expect_silent(plotPvalues(diffs.sort, gridspacing = 3))
  
  #Test for predictPresent
  mx.asr <- asreml(Yield ~ xNitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  
  testthat::expect_equal(length(mx.asr$gammas),3)
  current.asrt <- asrtests(mx.asr)
  print(current.asrt)
  
  diffs <- predictPresent(mx.asr, terms = "xNitrogen:Variety", 
                          x.num = "xNitrogen", x.fac = "Nitrogen", 
                          x.pred.values = sort(unique(Oats.dat$xNitrogen)),
                          wald.tab = current.asrt$wald.tab,
                          error.intervals = "Stand", tables = "none",
                          sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_is(diffs[[1]], "alldiffs")
  testthat::expect_equal(nrow(diffs[[1]]$predictions),12)
  testthat::expect_equal(ncol(diffs[[1]]$predictions),7)
  testthat::expect_equal(as.character(diffs[[1]]$predictions$Variety[[1]]),"Marvellous")
  testthat::expect_equal(length(attributes(diffs$xNitrogen.Variety)),10)
  testthat::expect_equal(length(attr(diffs[[1]], which = "sortOrder")),3)
  
  #Test for predictPlus with sortFactor
  diffs <- predictPlus(mx.asr, classify = "xNitrogen:Variety", 
                       x.num = "xNitrogen", x.fac = "Nitrogen", 
                       x.pred.values = sort(unique(Oats.dat$xNitrogen)),
                       wald.tab = current.asrt$wald.tab, 
                       error.intervals = "Stand", tables = "none",
                       sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_equal(nrow(diffs$predictions),12)
  testthat::expect_equal(ncol(diffs$predictions),7)
  testthat::expect_equal(as.character(diffs$predictions$Variety[1]),"Marvellous")
  testthat::expect_equal(length(attributes(diffs)),10)
  testthat::expect_true(all(attr(diffs, which = "sortOrder") == 
                              levels(diffs$predictions$Variety)))
  testthat::expect_true(all(attr(diffs, which = "sortOrder") == 
                              c("Marvellous","Golden Rain", "Victory")))
  testthat::expect_silent(plotPredictions(data = diffs$predictions, 
                                          classify ="xNitrogen:Variety", 
                                          x.num = "xNitrogen", x.fac = "Nitrogen", 
                                          y = "predicted.value", 
                                          error.intervals = "StandardError",  
                                          y.title = attr(diffs, 
                                                         which = "response.title")))
  testthat::expect_silent(plotPvalues(diffs, gridspacing = 3))
  
})

cat("#### Test for subset.alldiffs on Smarthouse with asreml3\n")
test_that("subset.alldiffs_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  library(dae)
  data(Smarthouse.dat)
  
  #Run analysis and produce alldiffs object
  m1.asr <- asreml(y1 ~ Genotype*A*B, 
                   random=~Replicate/Mainplot/Subplot,
                   data=Smarthouse.dat)
  testthat::expect_equal(length(m1.asr$gammas),4)
  current.asrt <- asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m <- current.asrt$asreml.obj
  testthat::expect_equal(length(m$gammas),3)
  
  diffs <- predictPlus(m, classify = "Genotype:A:B", 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_equal(nrow(diffs$predictions),120)
  testthat::expect_equal(ncol(diffs$predictions),8)
  testthat::expect_equal(as.character(diffs$predictions$Genotype[1]),"Axe")
  testthat::expect_equal(length(attributes(diffs)),8)
  testthat::expect_null(attr(diffs, which = "sortOrder"))
  
  #Form subset
  diffs.subs <- subset(diffs, 
                       subset = !grepl("E",Genotype, fixed = TRUE) & 
                         B %in% c("D1","D2"))
  testthat::expect_is(diffs.subs, "alldiffs")
  testthat::expect_equal(nrow(diffs.subs$predictions),48)
  testthat::expect_equal(ncol(diffs.subs$predictions),8)
  testthat::expect_false(any(diffs.subs$predictions$Genotype %in% c("Excalibur","Espada")))
  testthat::expect_false(any(diffs.subs$predictions$B %in% c("D3","D4")))
  testthat::expect_equal(length(attributes(diffs.subs)),8)
  
  #Test subset with removal of vars
  diffs.subs <- subset(diffs, subset = A == "N1" & B == "D2", rmClassifyVars = c("A","B"))
  testthat::expect_is(diffs.subs, "alldiffs")
  testthat::expect_equal(nrow(diffs.subs$predictions),10)
  testthat::expect_equal(ncol(diffs.subs$predictions),6)
  testthat::expect_false(any(c("A","B") %in% names(diffs.subs$predictions)))

  data(WaterRunoff.dat)
  #Run analysis and produce alldiffs object
  current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- asrtests(current.asr, NULL, NULL)
  diffs <- predictPlus.asreml(classify = "Sources:Type", 
                              asreml.obj = current.asr, tables = "none", 
                              wald.tab = current.asrt$wald.tab, 
                              present = c("Type","Species","Sources"))

  #Use subset.alldiffs to select a subset of the alldiffs object
  diffs.subs <- subset(diffs, 
                       subset = grepl("R", Sources, fixed = TRUE) & 
                         Type %in% c("Control","Medicinal"))
  testthat::expect_is(diffs.subs, "alldiffs")
  testthat::expect_equal(nrow(diffs.subs$predictions),10)
  testthat::expect_equal(ncol(diffs.subs$predictions),7)
  testthat::expect_false(any(diffs.subs$predictions$Sources == "Tap water"))
  testthat::expect_false(any(diffs.subs$predictions$B %in% c("Landscape","Culinary")))
  testthat::expect_equal(length(attributes(diffs.subs)),8)
})

cat("#### Test for linear.transformation Oats with asreml3\n")
test_that("linear.transformation_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  testthat::expect_equal(length(m1.asr$gammas),3)
  current.asrt <- asrtests(m1.asr)
  
  #Test store of vcov by predictPlus
  diffs <- predictPlus(m1.asr, classify = "Nitrogen:Variety", Vmatrix = TRUE, 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_equal(nrow(diffs$vcov),12)
  testthat::expect_true(all(colnames(diffs$vcov)[1:2] %in% c("0,Victory", "0,Golden Rain")))
  testthat::expect_true(all(abs((diffs$vcov[1,1:2] - c(82.93704, 35.74618))) < 1e-4))

  #Test linear transformation that compares with other N levels
  nv <- length(levels(diffs$predictions$Variety))
  L <- cbind(kronecker(matrix(rep(1, 3), nrow = 3), diag(1, nrow=nv)), 
             diag(-1, nrow=3*nv))
  rownames(L) <- colnames(diffs$vcov)[4:12]
  diffs.L <- predictPlus(m1.asr, classify = "Nitrogen:Variety", 
                         linear.transformation = L, Vmatrix = TRUE, 
                         wald.tab = current.asrt$wald.tab,
                         error.intervals = "Conf", tables = "none")
  testthat::expect_is(diffs.L, "alldiffs")
  testthat::expect_true(abs((diffs$vcov[1,1] - diffs$vcov[4,1]) + 
                              (diffs$vcov[4,4] - diffs$vcov[1,4]) - diffs.L$vcov[1,1]) < 1e-04)
  testthat::expect_equal(as.character(diffs.L$predictions$Combination[1]), "0.2,Victory")
  testthat::expect_true(abs((diffs.L$predictions$predicted.value[1] - 
                               qt(0.975,attr(diffs.L, which = "tdf")) * 
                               diffs.L$predictions$standard.error[1]) - 
                              diffs.L$predictions$lower.Confidence.limit[1]) < 1e-04)
  
  #Test model for linear transformation
  diffs.mod <- predictPlus(m1.asr, classify = "Nitrogen:Variety", 
                           linear.transformation = ~ Variety + Nitrogen, 
                           wald.tab = current.asrt$wald.tab,
                           error.intervals = "half", 
                           meanLSD.type = "factor.comb",
                           LSDby = "Nitrogen",
                           tables = "none")
  testthat::expect_is(diffs.mod, "alldiffs")
  testthat::expect_true(is.null(diffs.mod$vcov))

  m2.asr <- asreml(Yield ~ Nitrogen+Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  preds <- predict(m2.asr, classify = "Nitrogen:Variety")$predictions$pvals 
  testthat::expect_true(all(abs(diffs.mod$predictions$predicted.value - 
                                  preds$predicted.value) < 1e-04))
  testthat::expect_true(abs(diffs.mod$predictions$standard.error[1] - 
                              preds$standard.error[1]) > 0.01)
  
  #Test example in manual
  data(WaterRunoff.dat)
  #Run analysis and produce alldiffs object
  current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- asrtests(current.asr, NULL, NULL)
  diffs <- predictPlus.asreml(classify = "Sources:Species", Vmatrix = TRUE, 
                              asreml.obj = current.asr, tables = "none", 
                              wald.tab = current.asrt$wald.tab, 
                              present = c("Type","Species","Sources"))
  #Obtain predictions additive for Source and Species
  diffs.sub <- linTransform(diffs, classify = "Sources:Species", Vmatrix = TRUE,
                            linear.transformation = ~ Sources + Species,
                            tables = "none")
  testthat::expect_equal(diffs.sub$predictions$predicted.value[1] - 
                           diffs.sub$predictions$predicted.value[6],
                         diffs.sub$predictions$predicted.value[7] - 
                           diffs.sub$predictions$predicted.value[12])

  #Form contrast matrix for all sources
  L <- kronecker(diag(1, nrow = 4), 
                 cbind(diag(1, nrow = 5), matrix(rep(-1, 5), ncol = 1)))
  L <- mat.dirsum(list(L, 
                       kronecker(diag(1, nrow = 2), 
                                 cbind(diag(1, nrow = 7), 
                                       matrix(rep(-1, 7), ncol = 1)))))
  #Will get NaNs because differences between every 6th contrast are zero 
  #because the predictions are additive
  testthat::expect_warning(diffs.L <- linTransform(diffs.sub, 
                                                   classify = "Sources:Species",
                                                   linear.transformation = L,
                                                   tables = "none"))
  testthat::expect_true(all(abs(diffs.L$predictions$predicted.value[c(1,6,11,16,21,28)] - 
  (diffs.sub$predictions$predicted.value[1] - diffs.sub$predictions$predicted.value[6])) < 1e-06))
  
  #More efficient version for manual
  data(WaterRunoff.dat)
  current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- asrtests(current.asr, NULL, NULL)
  #Get additive predictions directly using predictPlus
  diffs.sub <- predictPlus.asreml(classify = "Sources:Species", Vmatrix = TRUE, 
                                  linear.transformation = ~ Sources + Species,
                                  asreml.obj = current.asr, tables = "none", 
                                  wald.tab = current.asrt$wald.tab, 
                                  present = c("Type","Species","Sources"))
  #Contrast matrix for differences between each species and non-planted for the last source
  L <- cbind(matrix(rep(0,7*32), nrow = 7, ncol = 32),
             diag(1, nrow = 7), 
             matrix(rep(-1, 7), ncol = 1))
  rownames(L) <- as.character(diffs.sub$predictions$Species[33:39])
  diffs.L <- linTransform(diffs.sub, 
                          classify = "Sources:Species",
                          linear.transformation = L,
                          tables = "predictions")
  testthat::expect_true(abs(diffs.L$predictions$predicted.value[1] + 0.0406963) < 1e-04)
  testthat::expect_true(diffs.L$predictions$Combination[1] == "S. iqscjbogxah")
  
  #Test a single contrast
  L1 <- matrix(L[1,1:40], nrow = 1)
  diffs.L1 <- linTransform(diffs.sub, 
                          classify = "Sources:Species",
                          linear.transformation = L1,
                          tables = "predictions")
  
  #Test for unbalanced two-way factorial 
  #- demonstrates will not make predictions for missing cells because NA in two-way table
  twoway <- data.frame(A = factor(rep(1:2, c(4,5))),
                       B = factor(rep(c(1:3,1:3), c(2,1,1,2,1,2))),
                       y = c(6,4,3,3,3,5,5,5,7))
  twoway <- twoway[-4,]
  mod <- asreml(fixed = y ~ A*B, data = twoway)
  pred.full <- do.call("predict",
                       args = list(object = mod, classify = "A:B", data = quote(twoway)))
  testthat::expect_true(any(is.na(pred.full$predictions$pvals$predicted.value)))
  
  #Make predictions additive
  diffs.sub <- do.call("predictPlus",
                       args = list(asreml.obj = mod, classify = "A:B", 
                                   linear.transformation = ~ A + B,
                                   tables = "predictions", 
                                   error.intervals = "Stand", 
                                   data = quote(twoway)))
  testthat::expect_equal(nrow(diffs.sub$predictions), 5)
  #Fit additive mode
  mod.add <- asreml(y ~ A+B, data = twoway)
  pred.add <- do.call("predict",
                      args = list(object = mod.add, classify = "A:B", data = quote(twoway)))
  testthat::expect_true(!any(is.na(pred.add$predictions$pvals$predicted.value)))
  #Compare projected and fitted additive predictions
  testthat::expect_true(abs(diffs.sub$predictions$predicted.value[1] - 
                              pred.add$predictions$pvals$predicted.value[1]) > 0.01)

  #Test backtransforms
  data(cart.dat)
  suffices <- list("32_42","42_50","50_59")
  responses.RGR <- paste("RGR_sm", suffices, sep = "_")
  responses.lRGR <- gsub("RGR", "lRGR", responses.RGR, fixed = TRUE)
  lRGR <- as.data.frame(lapply(responses.RGR, 
                               function(resp, dat)
                               {
                                 x <- log(dat[[resp]])
                               }, dat = cart.dat))
  names(lRGR) <- responses.lRGR
  cart.dat <- cbind(cart.dat, lRGR)

  
  nresp <- length(responses.lRGR)
  save <- vector(mode = "list", length = nresp)
  names(save) <- responses.lRGR
  nresp <- 1
  for (k in 1:nresp)
  {
    fix <- paste(responses.lRGR[k], " ~ Smarthouse/(xMainPosn + xLane) + Genotype.ID*Treatment.1",
                 sep= "")
    HEB25.asr <- do.call("asreml", 
                         list(fixed = as.formula(fix), 
                              random = ~ Smarthouse:Zones:Mainplots, 
                              rcov = ~ idh(Treat.Smarthouse):Zones:Mainplots, 
                              data = cart.dat, workspace=32e+06, 
                              keep.order=TRUE, na.method.X="include", maxiter=50))
    summary(HEB25.asr)$varcomp
    current.asrt <- do.call("asrtests", 
                            args = list(asreml.obj = HEB25.asr, 
                                        wald.tab = NULL, test.summary = NULL))
    current.asrt <- rmboundary(current.asrt)
    current.asr <- current.asrt$asreml.obj
    wald.tab <- recalcWaldTab(asrtests.obj=current.asrt, dDF.na="residual")
    HEB25.diffs <- do.call("predictPlus",
                           list(asreml.obj = current.asr, 
                                classify="Treatment.1:Genotype.ID",
                                wald.tab = wald.tab, Vmatrix = TRUE, 
                                transform.power = 0,
                                tables= "none", pworkspace=32E+06))
    #deal with missing value for HEB-20-125,W
    genos <- HEB25.diffs$predictions$Genotype.ID[HEB25.diffs$predictions$Treatment.1=="D"]
    genos <- as.character(genos)
    ng <- length(genos)
    L <- diag(1, nrow=(ng))
    rownames(L) <- colnames(L) <- genos
    L <- L[-match("HEB-20-125", genos),] #remove row
    L[,match("HEB-20-125", genos)] <- 0  #zero column
    L <- cbind(L, diag(-1, nrow=(ng-1)))
    testthat::expect_true(is.na(match("HEB-20-125", rownames(L))))
    save[[responses.lRGR[k]]] <- linTransform(HEB25.diffs, 
                                             classify = "Treatment.1:Genotype.ID", 
                                             linear.transformation = L, 
                                             transform.power = 0, 
                                             error.intervals = "Conf", 
                                             tables = "none")
  }
  testthat::expect_equal(sum(unlist(lapply(save$lRGR_sm_32_42, is.null))),1)
  testthat::expect_equal(dim(save$lRGR_sm_32_42$predictions), c(446,6))
  testthat::expect_equal(as.character(save$lRGR_sm_32_42$predictions$Combination[1]), 
                         "Barke")
  testthat::expect_true(all(abs(save$lRGR_sm_32_42$predictions[1, 2:5] - 
                                  c(-0.4131313, 0.04177124, -0.3306084, -0.4956541)) < 1e-04))
  testthat::expect_true(all(is.na(save$lRGR_sm_32_42$backtransforms[, "standard.error"])))
  testthat::expect_true(all(abs(exp(save$lRGR_sm_32_42$predictions[1, c(2,4:5)]) - 
                                  save$lRGR_sm_32_42$backtransforms[1, c(2,4:5)]) < 1e-04))
  })

