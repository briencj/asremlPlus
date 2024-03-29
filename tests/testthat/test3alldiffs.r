#devtools::test("asremlPlus")
context("prediction_alldiffs3")
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
  current.asrt <- as.asrtests(m1.asr)
  
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]

  #Test for as.alldiffs
  Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", 
                              sed=TRUE)$predictions
  Var.pred$pvals$Nitrogen <- factor(Var.pred$pvals$Nitrogen)
  Var.diffs <- as.alldiffs(predictions = Var.pred$pvals, 
                           sed = Var.pred$sed, 
                           classify = "Nitrogen:Variety", response = "Yield", tdf = den.df)
  testthat::expect_true(is.alldiffs(Var.diffs))
  testthat::expect_equal(nrow(Var.diffs$predictions),12)
  testthat::expect_true(validAlldiffs(Var.diffs))
  testthat::expect_true(is.null(attr(Var.diffs, which = "sortOrder")))
  
  #Test for allDifferences
  Var.sort.diffs <- allDifferences(predictions = Var.pred$pvals,
                                   classify = "Nitrogen:Variety", 
                                   sed = Var.pred$sed, tdf = den.df, 
                                   sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(is.alldiffs(Var.sort.diffs))
  testthat::expect_true(validAlldiffs(Var.sort.diffs))
  testthat::expect_equal(length(attr(Var.sort.diffs$predictions, which = "heading")),1)
  testthat::expect_true("asreml.predict" %in% class(Var.sort.diffs$predictions))
  testthat::expect_equal(length(attr(Var.sort.diffs, which = "sortOrder")),3)
  testthat::expect_true(as.character(Var.sort.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.sort.diffs$predictions$Variety[2]) == "Golden Rain")

  #Test for re-order factors with allDifferences
  Var.reord.diffs <- allDifferences(predictions = Var.pred$pvals,
                                    classify = "Variety:Nitrogen", 
                                    sed = Var.pred$sed, tdf = den.df)
  testthat::expect_true(as.character(Var.reord.diffs$predictions$Variety[1]) == "Victory" &
                          as.character(Var.reord.diffs$predictions$Variety[2]) == "Victory")

  #Test for re-order factors with renewClassify
  Var.reord.diffs <- renewClassify(Var.diffs, newclassify = "Variety:Nitrogen")
  testthat::expect_true(as.character(Var.reord.diffs$predictions$Variety[1]) == "Victory" &
                          as.character(Var.reord.diffs$predictions$Variety[2]) == "Victory")
  testthat::expect_equal(length(attr(Var.reord.diffs$predictions, which = "heading")),1)
  testthat::expect_true("asreml.predict" %in% class(Var.reord.diffs$predictions))
 
  #Test for re-order factors and sort
  Var.both.diffs <- allDifferences(predictions = Var.pred$pvals,
                                   classify = "Variety:Nitrogen", 
                                   sed = Var.pred$sed, tdf = den.df, 
                                   sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(as.character(Var.both.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.both.diffs$predictions$Variety[2]) == "Marvellous")
  Var.both.diffs <- renewClassify(Var.diffs, newclassify = "Variety:Nitrogen", 
                                     sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(as.character(Var.both.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.both.diffs$predictions$Variety[2]) == "Marvellous")
  
  #Test a single factor prediction
  diffsN <- predictPlus(m1.asr, classify = "Nitrogen", tables = "none")
  testthat::expect_true(validAlldiffs(diffsN))
  
  #Test single factor linear.transform
  Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", vcov=TRUE)$predictions
  Var.pred$pvals$Nitrogen <- factor(Var.pred$pvals$Nitrogen)
  Var.diffs <- allDifferences(predictions = Var.pred$pvals,
                              classify = "Nitrogen:Variety", 
                              vcov = Var.pred$vcov, tdf = den.df)
  Var.diffs.one <- linTransform(Var.diffs, linear.transformation = ~Nitrogen,
                                error.intervals = "half", tables = "none")
  testthat::expect_true(all(abs(Var.diffs.one$LSD - 9.883479) < 1e-06))
  testthat::expect_true(all(abs(Var.diffs.one$LSD - 
                              attr(Var.diffs.one$predictions, which = "meanLSD")) < 1E-06))
  #Test LSDby not in linear.transformation
  testthat::expect_warning(Var.diffs.by <- linTransform(Var.diffs, 
                                                        linear.transformation = ~Nitrogen,
                                                        error.intervals = "half", 
                                                        LSDtype = "factor", 
                                                        LSDby = "Variety", 
                                                        tables = "none"))
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
  current.asrt <- as.asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m <- current.asrt$asreml.obj
  testthat::expect_equal(length(m$gammas),3)
  
  diffs <- predictPlus(m, classify = "Genotype:A:B", 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_true(is.alldiffs(diffs))
  testthat::expect_true(validAlldiffs(diffs))
  testthat::expect_equal(nrow(diffs$predictions),120)
  testthat::expect_equal(ncol(diffs$predictions),8)
  testthat::expect_equal(as.character(diffs$predictions$Genotype[1]),"Axe")
  testthat::expect_true(is.null(attr(diffs, which = "sortOrder")))
  
  #Test reodering of the classify
  diffs.reord <- renewClassify(diffs, newclassify = "A:B:Genotype")
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
  testthat::expect_true(validAlldiffs(diffs.sort))
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
  current.asrt <- as.asrtests(m2.asr)
  diffs2.sort <- predictPlus(m2.asr, classify = "Genotype:A:B", 
                             wald.tab = current.asrt$wald.tab,
                             error.intervals = "Stand", tables = "none",
                             sortFactor = "Genotype", 
                             sortOrder = sort.order)
  testthat::expect_equal(as.character(diffs.sort$predictions$Genotype[1]),
                         as.character(diffs2.sort$predictions$Genotype[1]))
  testthat::expect_equal(attr(diffs.sort, which = "sortOrder"),
                         attr(diffs2.sort, which = "sortOrder"))
  
  #Test sort.alldiffs with sortParallelToCombo and increasing order
  diffs1.sort <- sort(diffs, sortFactor = "Genotype", 
                      sortParallelToCombo = list(A = "N3", B = "D4"),
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

cat("#### Test for sort.alldiffs on WaterRunoff with asreml3\n")
test_that("sort.alldiffsWater3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)
  
  #Analyse pH  
  m1.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- as.asrtests(m1.asr, NULL, NULL)
  testthat::expect_equal(length(m1.asr$gammas),2)
  current.asrt <- as.asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m1.asr <- current.asrt$asreml.obj
  testthat::expect_equal(length(m1.asr$gammas),1)
  
  TS.diffs <- predictPlus.asreml(classify = "Sources:Type", 
                                 asreml.obj = m1.asr, tables = "none", 
                                 wald.tab = current.asrt$wald.tab, 
                                 present = c("Type","Species","Sources"))
  testthat::expect_true(is.alldiffs(TS.diffs))
  testthat::expect_true(validAlldiffs(TS.diffs))
  testthat::expect_equal(nrow(TS.diffs$predictions),20)
  testthat::expect_equal(ncol(TS.diffs$predictions),7)
  testthat::expect_equal(as.character(TS.diffs$predictions$Type[1]),"Landscape")
  testthat::expect_true(is.null(attr(TS.diffs, which = "sortOrder")))
  
  #Test reodering of the classify
  TS.diffs.reord <- renewClassify(TS.diffs, newclassify = "Type:Sources")
  testthat::expect_equal(as.character(TS.diffs.reord$predictions$Sources[1]),"Rainwater")
  testthat::expect_equal(as.character(TS.diffs.reord$predictions$Sources[2]),"Recycled water")
  testthat::expect_true(abs(TS.diffs.reord$predictions$predicted.value[2] - 7.646389) < 1e-06)
  
  #Test sort.alldiffs and save order for use with other response variables
  TS.diffs.sort <- sort(TS.diffs, sortFactor = "Sources", sortParallelToCombo = list(Type = "Control"))
  sort.order <- attr(TS.diffs.sort, which = "sortOrder")
  testthat::expect_is(TS.diffs.sort, "alldiffs")
  testthat::expect_true(validAlldiffs(TS.diffs.sort))
  testthat::expect_equal(nrow(TS.diffs.sort$predictions),20)
  testthat::expect_equal(ncol(TS.diffs.sort$predictions),7)
  testthat::expect_equal(as.character(TS.diffs.sort$predictions$Sources[1]),"Recycled water")
  testthat::expect_equal(length(attributes(TS.diffs.sort)),10)
  testthat::expect_equal(length(attr(TS.diffs.sort, which = "sortOrder")),6)
  
  #Test sort.alldiffs with supplied sortOrder
  m2.asr <- asreml(fixed = Turbidity ~ Benches + (Sources * (Type + Species)), 
                   random = ~ Benches:MainPlots,
                   keep.order=TRUE, data= WaterRunoff.dat)
  testthat::expect_equal(length(m2.asr$gammas),2)
  current.asrt <- as.asrtests(m2.asr)
  diffs2.sort <- predictPlus(m2.asr, classify = "Sources:Type", 
                             pairwise = FALSE, error.intervals = "Stand", 
                             tables = "none", present = c("Type","Species","Sources"),
                             sortFactor = "Sources", 
                             sortOrder = sort.order)
  testthat::expect_equal(as.character(TS.diffs.sort$predictions$Sources[1]),
                         as.character(diffs2.sort$predictions$Sources[1]))
  testthat::expect_equal(attr(TS.diffs.sort, which = "sortOrder"),
                         attr(diffs2.sort, which = "sortOrder"))
  
  #Test removing a multilevel classifying factor that is marginal to another factor using subset
  diffs.full <- predictPlus.asreml(asreml.obj = m1.asr, 
                                   classify = "Sources:Type:Species", 
                                   wald.tab = current.asrt$wald.tab, 
                                   present = c("Type","Species","Sources"),
                                   tables = "none", Vmatrix = TRUE)
  testthat::expect_true(setequal(names(diffs.full$predictions), 
                                 c("Sources", "Type", "Species", "predicted.value", 
                                   "standard.error", "upper.Confidence.limit", 
                                   "lower.Confidence.limit", "est.status")))
  diffs.fit <- linTransform(diffs.full, classify = "Sources:Type:Species",
                            linear.transformation = ~ Sources:Type, 
                            error.intervals="half", 
                            LSDtype="factor", LSDby="Type", 
                            tables = "none")
  testthat::expect_true(setequal(names(diffs.fit$predictions), 
                                 c("Sources", "Type", "Species", "predicted.value", 
                                   "standard.error", "upper.halfLeastSignificant.limit", 
                                   "lower.halfLeastSignificant.limit", "est.status")))
  diffs.fit <- subset(diffs.fit, rmClassifyVars = "Type")
  testthat::expect_true(setequal(names(diffs.fit$predictions), 
                                 c("Sources", "Species", "predicted.value", 
                                   "standard.error", "upper.halfLeastSignificant.limit", 
                                   "lower.halfLeastSignificant.limit", "est.status")))
  
  #Check that renewClassify also works when the full classify is not supplied
  diffs.red <- diffs.full
  diffs.red$predictions <- diffs.red$predictions[,
                                            -match("Type", names(diffs.red$predictions))]
  diffs.red <- renewClassify(diffs.red, newclassify = "Sources:Species")
  testthat::expect_true(setequal(names(diffs.red$predictions), 
                                 c("Sources", "Species", "predicted.value", 
                                   "standard.error", "upper.Confidence.limit", 
                                   "lower.Confidence.limit", "est.status")))

  #Check that renewClassify fails when newclassify does not uniquely index the predictions
  diffs.red <- diffs.full
  diffs.red$predictions <- diffs.red$predictions[,
                                          -match("Species", names(diffs.red$predictions))]
  testthat::expect_error(diffs.red <- renewClassify(diffs.red, 
                                                    newclassify = "Sources:Type"))
  
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
  current.asrt <- as.asrtests(m1.asr)
  
  #Test for as.alldiffs
  Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", 
                      sed=TRUE)$predictions
  Var.pred$pvals$Nitrogen <- factor(Var.pred$pvals$Nitrogen)
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  Var.diffs <- as.alldiffs(predictions = Var.pred$pvals, 
                           sed = Var.pred$sed, 
                           tdf = den.df, classify = "Nitrogen:Variety")
  testthat::expect_true(validAlldiffs(Var.diffs))
  testthat::expect_equal(length(attributes(Var.diffs)),4)
  testthat::expect_true(is.null(attr(Var.diffs, which = "sortOrder")))
  
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
  testthat::expect_true(validAlldiffs(diffs))
  testthat::expect_equal(nrow(diffs$predictions),12)
  testthat::expect_equal(ncol(diffs$predictions),7)
  testthat::expect_equal(as.character(diffs$predictions$Variety[1]),"Victory")
  testthat::expect_equal(length(attributes(diffs)),8)
  testthat::expect_true(is.null(attr(diffs, which = "sortOrder")))
  
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
  current.asrt <- as.asrtests(mx.asr)
  print(current.asrt)
  
  diffs <- predictPresent(mx.asr, terms = "xNitrogen:Variety", 
                          x.num = "xNitrogen", x.fac = "Nitrogen", 
                          x.pred.values = sort(unique(Oats.dat$xNitrogen)),
                          wald.tab = current.asrt$wald.tab,
                          error.intervals = "Stand", tables = "none",
                          sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_is(diffs[[1]], "alldiffs")
  testthat::expect_true(validAlldiffs(diffs$xNitrogen.Variety))
  testthat::expect_equal(nrow(diffs$xNitrogen.Variety$predictions),12)
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
  testthat::expect_true(validAlldiffs(diffs))
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
  current.asrt <- as.asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m <- current.asrt$asreml.obj
  testthat::expect_equal(length(m$gammas),3)
  
  diffs <- predictPlus(m, classify = "Genotype:A:B", 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs))
  testthat::expect_equal(nrow(diffs$predictions),120)
  testthat::expect_equal(ncol(diffs$predictions),8)
  testthat::expect_equal(as.character(diffs$predictions$Genotype[1]),"Axe")
  testthat::expect_equal(length(attributes(diffs)),8)
  testthat::expect_true(is.null(attr(diffs, which = "sortOrder")))
  
  #Form subset
  diffs.subs <- subset(diffs, 
                       subset = !grepl("E",Genotype, fixed = TRUE) & 
                         B %in% c("D1","D2"))
  testthat::expect_is(diffs.subs, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs.subs))
  testthat::expect_equal(nrow(diffs.subs$predictions),48)
  testthat::expect_equal(ncol(diffs.subs$predictions),8)
  testthat::expect_false(any(diffs.subs$predictions$Genotype %in% c("Excalibur","Espada")))
  testthat::expect_false(any(diffs.subs$predictions$B %in% c("D3","D4")))
  testthat::expect_equal(length(attributes(diffs.subs)),8)
  
  #Test subset with removal of vars
  diffs.subs <- subset(diffs, subset = A == "N1" & B == "D2", rmClassifyVars = c("A","B"))
  testthat::expect_is(diffs.subs, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs.subs))
  testthat::expect_equal(nrow(diffs.subs$predictions),10)
  testthat::expect_equal(ncol(diffs.subs$predictions),6)
  testthat::expect_false(any(c("A","B") %in% names(diffs.subs$predictions)))
  
  data(WaterRunoff.dat)
  #Run analysis and produce alldiffs object
  current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  diffs <- predictPlus.asreml(classify = "Sources:Type", 
                              asreml.obj = current.asr, tables = "none", 
                              wald.tab = current.asrt$wald.tab, 
                              present = c("Type","Species","Sources"))

  #Use subset.alldiffs to select a subset of the alldiffs object
  diffs.subs <- subset(diffs, 
                       subset = grepl("R", Sources, fixed = TRUE) & 
                         Type %in% c("Control","Medicinal"))
  testthat::expect_is(diffs.subs, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs.subs))
  testthat::expect_equal(nrow(diffs.subs$predictions),10)
  testthat::expect_equal(ncol(diffs.subs$predictions),7)
  testthat::expect_false(any(diffs.subs$predictions$Sources == "Tap water"))
  testthat::expect_false(any(diffs.subs$predictions$B %in% c("Landscape","Culinary")))
  testthat::expect_equal(length(attributes(diffs.subs)),8)
})

cat("#### Test for facCombine.alldiffs on Ladybird with asreml3\n")
test_that("facCombine.alldiffs3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  library(dae)
  data("Ladybird.dat")
  
  #Mixed model analysis of logits 
  m1.asr <- asreml(logitP ~ Host*Cadavers*Ladybird, 
                   random = ~ Run,
                   data = Ladybird.dat)
  
  testthat::expect_equal(length(m1.asr$gammas),2)
  current.asrt <- as.asrtests(m1.asr)
  testthat::expect_true(validAsrtests(current.asrt))
  
  HCL.pred <- predict(m1.asr, classify="Host:Cadavers:Ladybird", 
                                     sed=TRUE)
  if (getASRemlVersionLoaded(nchar = 1) == "3")
    HCL.pred <- HCL.pred$predictions
  HCL.preds <- HCL.pred$pvals
  HCL.sed <- HCL.pred$sed
  HCL.vcov <- NULL
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Host:Cadavers:Ladybird", rownames(wald.tab)), "denDF"]
  testthat::expect_equal(den.df,60)
  
  ## Use lme4 and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(logitP ~ Host*Cadavers*Ladybird + (1|Run),
                              data=Ladybird.dat)
    HCL.emm <- emmeans::emmeans(m1.lmer, specs = ~ Host:Cadavers:Ladybird)
    HCL.preds <- summary(HCL.emm)
    den.df <- min(HCL.preds$df)
    ## Modify HCL.preds to be compatible with a predictions.frame
    names(HCL.preds)[match(c("emmean", "SE", "lower.CL", "upper.CL"), 
                           names(HCL.preds))] <- c("predicted.value", 
                                                   "standard.error", 
                                                   "lower.Confidence.limit",
                                                   "upper.Confidence.limit")
    HCL.preds$est.status <- "Estimable"
    HCL.preds$est.status[is.na(HCL.preds$predicted.value)] <- "Aliased"
    HCL.vcov <- vcov(HCL.emm)
    HCL.sed <- NULL
  }
  
  ## Form an all.diffs object with predictions obtained with either asreml or lmerTest
  HCL.diffs <- as.alldiffs(predictions = HCL.preds, classify = "Host:Cadavers:Ladybird", 
                           sed = HCL.sed, vcov = HCL.vcov, tdf = den.df)
    
  ## check the class and validity of the alldiffs object
  is.alldiffs(HCL.diffs)
  validAlldiffs(HCL.diffs)
  testthat::expect_equal(nrow(HCL.diffs$predictions),12)
  testthat::expect_equal(ncol(HCL.diffs$predictions),9)
  
  ## Combine Cadavers and Ladybird
  Comb.diffs <- facCombine(HCL.diffs, factors = c("Cadavers","Ladybird"))
  
  ## check the validity of Comb.diffs
  validAlldiffs(Comb.diffs)
  testthat::expect_equal(nrow(Comb.diffs$predictions),12)
  testthat::expect_equal(ncol(Comb.diffs$predictions),8)
  testthat::expect_true(all(c("Host", "Cadavers_Ladybird", "predicted.value") %in% 
                              names(Comb.diffs$predictions)))
    
  ## Rename Cadavers
  HCL.rename.diffs <- facRename(HCL.diffs, factor.names = "Cadavers", newnames = "Cadaver.nos")
  testthat::expect_true(validAlldiffs(HCL.rename.diffs))
  testthat::expect_true("Cadaver.nos" %in% names(HCL.rename.diffs$predictions))
  
  ## Recast Ladybird
  HCL.recast.diffs <- facRecast(HCL.diffs, factor = "Ladybird", newlabels = c("none", "present"))
  testthat::expect_true(validAlldiffs(HCL.recast.diffs))
  testthat::expect_true(all(levels(HCL.recast.diffs$predictions$Ladybird) == c("none", "present")))
  HCL.recast.diffs <- facRecast.alldiffs(HCL.recast.diffs, factor = "Host", levels.order = c("trefoil", "bean"))
  testthat::expect_true(validAlldiffs(HCL.recast.diffs))
  testthat::expect_true(all(levels(HCL.recast.diffs$predictions$Host) == c("trefoil", "bean")))
  HCL.recast.diffs <- facRecast(HCL.recast.diffs, factor = "Ladybird", levels.order = c("present", "none"), 
                                newlabels = c("yes","no"))
  testthat::expect_true(validAlldiffs(HCL.recast.diffs))
  testthat::expect_true(all(levels(HCL.recast.diffs$predictions$Ladybird) == c("yes", "no")))
  
})

cat("#### Test for linear.transformation on Oats with asreml3\n")
test_that("linear.transform_Oats_asreml3", {
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
  current.asrt <- as.asrtests(m1.asr)
  
  #Test store of vcov by predictPlus
  diffs <- predictPlus(m1.asr, classify = "Nitrogen:Variety", Vmatrix = TRUE, 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_equal(length(attr(diffs$predictions, which = "heading")),1)
  testthat::expect_true("asreml.predict" %in% class(diffs$predictions))
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
  testthat::expect_equal(length(attr(diffs.L$predictions, which = "heading")),2)
  testthat::expect_true("asreml.predict" %in% class(diffs.L$predictions))
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
                           LSDtype = "factor.comb",
                           LSDby = "Nitrogen",
                           tables = "none")
  testthat::expect_is(diffs.mod, "alldiffs")
  testthat::expect_equal(length(attr(diffs.mod$predictions, which = "heading")),2)
  testthat::expect_true("asreml.predict" %in% class(diffs.mod$predictions))
  testthat::expect_true(is.null(diffs.mod$vcov))

  m2.asr <- asreml(Yield ~ Nitrogen+Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  preds <- predict(m2.asr, classify = "Nitrogen:Variety")$predictions$pvals 
  testthat::expect_true(all(abs(diffs.mod$predictions$predicted.value - 
                                  preds$predicted.value) < 1e-04))
  testthat::expect_true(abs(diffs.mod$predictions$standard.error[1] - 
                              preds$standard.error[1]) > 0.01)
})

cat("#### Test for linear.transformation on WaterRunoff with asreml3\n")
test_that("linear.transform_WaterRunoff_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  library(dae)
  #Test example in manual
  data(WaterRunoff.dat)
  #Run analysis and produce alldiffs object
  current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
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
  #check for zero seds and their removal
  ksed <- na.omit(as.vector(diffs.L$sed))
  testthat::expect_true(abs(diffs.L$LSD["minLSD"] - 0.1246359) < 1e-06)
  testthat::expect_true(all(abs(diffs.L$predictions$predicted.value[c(1,6,11,16,21,28)] - 
                                  (diffs.sub$predictions$predicted.value[1] - 
                                     diffs.sub$predictions$predicted.value[6])) < 1e-06))
  
  #More efficient version for manual
  data(WaterRunoff.dat)
  current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
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
                              keep.order=TRUE, na.method.X="include", maxit=50))
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

test_that("addBacktransforms_WaterRunoff_asreml3", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr3.lib)
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)
  
  ##Use asreml to get predictions and associated statistics
  
  current.asr <- asreml(fixed = log.Turbidity ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  TS.diffs <- predictPlus(classify = "Sources:Type", 
                          asreml.obj = current.asr, 
                          wald.tab = current.asrt$wald.tab, 
                          present = c("Sources", "Type", "Species"))
  
  ## Plot p-values for predictions obtained using asreml3
  if (exists("TS.diffs"))
  {
    ##Add the backtransforms component for predictions obtained using asreml or lmerTest  
    TS.diffs <- addBacktransforms.alldiffs(TS.diffs, transform.power = 0)
    testthat::expect_false(is.null(TS.diffs$backtransforms))
    testthat::expect_true(all(abs(exp(TS.diffs$predictions$predicted.value)-
                                    TS.diffs$backtransforms$backtransformed.predictions) < 1e-06))
    testthat::expect_true(all(abs(exp(TS.diffs$predictions$upper.Confidence.limit)-
                                    TS.diffs$backtransforms$upper.Confidence.limit) < 1e-06))
  }
})
