#devtools::test("asremlPlus")
context("prediction_presentation3")
asr3.lib <- "D:\\Analyses\\R oldpkg" 

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
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs))
  testthat::expect_equal(nrow(diffs$predictions),120)
  testthat::expect_equal(ncol(diffs$predictions),8)
  testthat::expect_equal(as.character(diffs$predictions$Genotype[1]),"Axe")
  testthat::expect_equal(length(attributes(diffs)),8)
  testthat::expect_true(is.null(attr(diffs, which = "sortOrder")))
  
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

cat("#### Test for sort.alldiffs on Oats with asremls\n")
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
  Var.pred <- asreml:::predict.asreml(m1.asr, classify="Nitrogen:Variety", 
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

cat("#### Test for sort in standard order for classify on Oats with asreml3\n")
test_that("classify.sort_asreml3", {
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
  Var.pred <- asreml:::predict.asreml(m1.asr, classify="Variety:Nitrogen", 
                                      sed=TRUE)$predictions
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  Var.diffs <- as.alldiffs(predictions = Var.pred$pvals, 
                           sed = Var.pred$sed, 
                           tdf = den.df)
  testthat::expect_equal(length(attributes(Var.diffs)),3)
  testthat::expect_true(is.null(attr(Var.diffs, which = "sortOrder")))
  #asreml-R3 does not ensure that factors are in the same order as the classify - 
  #they are in the same order as the term in asreml call
  testthat::expect_true(all(names(Var.diffs$predictions)[1:2] ==  c("Nitrogen","Variety")))
  
  #Test for allDifferences without sortFactor
  Var.diffs <- allDifferences(predictions = Var.pred$pvals,
                              classify = "Variety:Nitrogen", 
                              sed = Var.pred$sed, tdf = den.df)
  testthat::expect_true(all(names(Var.diffs$predictions)[1:2] ==  c("Variety","Nitrogen")))
  testthat::expect_equal(as.character(Var.diffs$predictions$Variety[1]),"Victory")
  
  #TestPlot with sort
  testthat::expect_silent(plotPvalues(Var.diffs, gridspacing = 4))
  
  #Test for allDifferences with sortFactor
  Var.diffs <- allDifferences(predictions = Var.pred$pvals,
                              classify = "Variety:Nitrogen", 
                              sed = Var.pred$sed, tdf = den.df, 
                              sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(all(names(Var.diffs$predictions)[1:2] ==  c("Variety","Nitrogen")))
  testthat::expect_equal(as.character(Var.diffs$predictions$Variety[1]),"Marvellous")
  testthat::expect_equal(length(attr(Var.diffs, which = "sortOrder")),3)
  
  #TestPlot with sort
  testthat::expect_silent(plotPvalues(Var.diffs, gridspacing = 3, 
                                      sortFactor = "Variety", decreasing = TRUE))
  
  #Test for sort.alldiffs
  diffs <- allDifferences(predictions = Var.pred$pvals,
                          classify = "Variety:Nitrogen", 
                          sed = Var.pred$sed, tdf = den.df)
  diffs.sort <- sort(diffs, sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_equal(as.character(diffs.sort$predictions$Variety[1]),"Marvellous")
  testthat::expect_silent(plotPvalues(diffs.sort, gridspacing = 3))
  
})

