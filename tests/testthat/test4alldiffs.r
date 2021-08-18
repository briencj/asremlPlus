#devtools::test("asremlPlus")
context("prediction_alldiffs")

cat("#### Test for allDifferences.data.frame sort.alldiffs on Oats with asreml4\n")
test_that("allDifferences_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  testthat::expect_equal(length(m1.asr$vparameters),3)
  current.asrt <- as.asrtests(m1.asr)
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  
  #Test for as.alldiffs - only loads into alldiffs object
  Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", sed=TRUE)
  Var.diffs <- as.alldiffs(predictions = Var.pred$pvals, 
                           sed = Var.pred$sed, 
                           classify = "Nitrogen:Variety", response = "Yield", tdf = den.df)
  testthat::expect_true(is.alldiffs(Var.diffs))
  testthat::expect_equal(nrow(Var.diffs$predictions),12)
  testthat::expect_true(validAlldiffs(Var.diffs))
  testthat::expect_true(is.null(attr(Var.diffs, which = "sortOrder")))
  testthat::expect_true(setequal(names(Var.diffs$predictions), 
                                 c("Nitrogen", "Variety", "predicted.value", "standard.error", "est.status")))
  testthat::expect_true(attr(Var.diffs, which = "alpha") == 0.05)
  testthat::expect_true(all(!(c( "LSDtype", "LSDstatistic") %in% names(attributes(Var.diffs)))))
  testthat::expect_true(is.null(Var.diffs$LSD))
  
  #Test for allDifferences without the sort
  Var.diffs <- allDifferences(predictions = Var.pred$pvals, 
                           sed = Var.pred$sed, 
                           classify = "Nitrogen:Variety", response = "Yield", tdf = den.df)
  testthat::expect_true(is.alldiffs(Var.diffs))
  testthat::expect_equal(nrow(Var.diffs$predictions),12)
  testthat::expect_true(validAlldiffs(Var.diffs))
  testthat::expect_true(is.null(attr(Var.diffs, which = "sortOrder")))
  testthat::expect_true(setequal(names(Var.diffs$predictions), 
                                 c("Nitrogen", "Variety", "predicted.value", "standard.error", "est.status")))
  testthat::expect_true(attr(Var.diffs, which = "alpha") == 0.05)
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic") %in% names(attributes(Var.diffs))))
  testthat::expect_true(!is.null(Var.diffs$LSD))
  
  #Test for allDifferences
  Var.sort.diffs <- allDifferences(predictions = Var.pred$pvals,
                                   classify = "Nitrogen:Variety", 
                                   sed = Var.pred$sed, tdf = den.df, 
                                   sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(is.alldiffs(Var.sort.diffs))
  testthat::expect_true(validAlldiffs(Var.sort.diffs))
  testthat::expect_equal(length(attr(Var.sort.diffs$predictions, which = "heading")),4)
  testthat::expect_true("asreml.predict" %in% class(Var.sort.diffs$predictions))
  testthat::expect_equal(length(attr(Var.sort.diffs, which = "sortOrder")),3)
  testthat::expect_true(as.character(Var.sort.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.sort.diffs$predictions$Variety[2]) == "Golden Rain")
  testthat::expect_true(setequal(names(Var.sort.diffs$predictions), 
                                 c("Nitrogen", "Variety", "predicted.value", "standard.error", "est.status")))
  testthat::expect_true(all(c("alpha", "LSDtype", "LSDstatistic") %in% names(attributes(Var.sort.diffs))))

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
  testthat::expect_equal(length(attr(Var.reord.diffs$predictions, which = "heading")),4)
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
  
 })


cat("#### Test for LSDs and halfLSIs on system data with asreml4\n")
test_that("LSD_LSI_SystemData_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asremlPlus)
  library(dae)
  load(system.file("extdata", "testDiffs.rda", package = "asremlPlus", mustWork = TRUE))
  
  #Rectify diff.Clup
  diffs.new <- redoErrorIntervals(diffs.ClUp, error.intervals = "half", 
                                 LSDtype = "factor", LSDby = attr(diffs.ClUp, which = "LSDby"))
  testthat::expect_true(all(names(diffs.new$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDby", "LSDstatistic") %in% 
                              names(attributes(diffs.new))))
  testthat::expect_true(attr(diffs.new, which = "LSDtype") == "factor.combinations")
  testthat::expect_true(all(attr(diffs.new, which = "LSDby") == attr(diffs.ClUp, which = "LSDby")))
  testthat::expect_true(attr(diffs.new, which = "LSDstatistic") == "mean")
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.new$predictions))))
  testthat::expect_true(length(attr(diffs.new$predictions, which = "LSDvalues")) ==  20)
  testthat::expect_true(all(abs(attr(diffs.new$predictions, which = "LSDvalues")[1:3] - 
                                  c(0.6946730,0.6969075,0.6965694)) < 1e-05))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.new$backtransforms))))
  testthat::expect_true(is.null(attr(diffs.new$backtransforms, which = "LSDvalues")))
  

  #Change half-LSIs to CIs and LSD component from overall to factor.combinations
  diffs.CI <- redoErrorIntervals(diffs.ClUp, LSDtype = "factor", LSDby = attr(diffs.ClUp, which = "LSDby"))
  testthat::expect_true(all(names(diffs.CI$predictions)[6:7] == 
                              c("upper.Confidence.limit", "lower.Confidence.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.CI))))
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.CI$predictions))))
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.CI$backtransforms))))
  testthat::expect_true(attr(diffs.CI, which = "LSDtype") == "factor.combinations")
  testthat::expect_true(all(attr(diffs.CI, which = "LSDby") == attr(diffs.ClUp, which = "LSDby")))
  testthat::expect_true(attr(diffs.CI, which = "LSDstatistic") == "mean")
  
  #Just recalc the LSDs component without updating the intervals
  diffs.LSD <- recalcLSD(diffs.ClUp, LSDtype = "factor", LSDby = attr(diffs.ClUp, which = "LSDby"))
  testthat::expect_true(all(names(diffs.LSD$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD))))
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD$predictions))))
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD$backtransforms))))
  testthat::expect_true(all(attr(diffs.LSD, which = "LSDby") == attr(diffs.ClUp, which = "LSDby")))
  testthat::expect_true(all(c( "transform.power", "offset", "scale") %in% names(attributes(diffs.LSD$backtransforms))))
  testthat::expect_true(attr(diffs.LSD$backtransforms, which = "transform.power") == 0)
  
  #facRename
  diffs.LSD.rename <- facRename(diffs.new, factor.names = "Temperature", newnames = "Climate")
  testthat::expect_true(all(names(diffs.LSD.rename$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.rename))))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.LSD.rename$predictions))))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% 
                              names(attributes(diffs.LSD.rename$backtransforms))))
  testthat::expect_true(all(attr(diffs.LSD.rename, which = "LSDby") == c("Climate", "Genotype")))
  testthat::expect_true(all(attr(diffs.LSD.rename$predictions, which = "LSDby") == c("Climate", "Genotype")))
  testthat::expect_true(all(attr(diffs.LSD.rename$backtransforms, which = "LSDby") == c("Climate", "Genotype")))
  testthat::expect_true(all(c( "transform.power", "offset", "scale") %in% 
                              names(attributes(diffs.LSD.rename$backtransforms))))
  testthat::expect_true(attr(diffs.LSD.rename$backtransforms, which = "transform.power") == 0)
  
  #facRecast
  diffs.LSD.recast <- facRecast(diffs.LSD.rename, factor = "Climate", newlabels = c("Cool","Warm"))
  testthat::expect_true(all(names(diffs.LSD.recast$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.recast))))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.LSD.recast$predictions))))
  testthat::expect_true(any(grepl("Warm", names(attr(diffs.LSD.recast$predictions, which = "LSDvalues")))))
  testthat::expect_true(any(grepl("Warm", rownames(diffs.LSD.recast$LSD), fixed = TRUE)))
  testthat::expect_true(!any(grepl("Hot", rownames(diffs.LSD.recast$LSD), fixed = TRUE)))
  testthat::expect_true(any(grepl("Warm", 
                                  names(attr(diffs.LSD.recast$predictions, which = "LSDvalues")), 
                                  fixed = TRUE)))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% 
                              names(attributes(diffs.LSD.recast$backtransforms))))
  testthat::expect_true(all(attr(diffs.LSD.recast, which = "LSDby") == c("Climate", "Genotype")))
  testthat::expect_true(all(levels(diffs.LSD.recast$predictions$Climate) == c("Cool", "Warm")))
  testthat::expect_true(all(attr(diffs.LSD.recast$predictions, which = "LSDby") == c("Climate", "Genotype")))
  testthat::expect_true(all(attr(diffs.LSD.recast$backtransforms, which = "LSDby") == c("Climate", "Genotype")))
  testthat::expect_true(all(c( "transform.power", "offset", "scale") %in% 
                              names(attributes(diffs.LSD.recast$backtransforms))))
  testthat::expect_true(attr(diffs.LSD.recast$backtransforms, which = "transform.power") == 0)
  
  #facCombine
  #Combine factors in the LSDBY
  diffs.LSD.comb <- facCombine(diffs.new, factors = c("Temperature", "Genotype"))
  testthat::expect_true(all(names(diffs.LSD.comb$predictions)[c(1,5:6)] == 
                              c("Temperature_Genotype", 
                                "upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.comb))))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.LSD.comb$predictions))))
  testthat::expect_true(any(grepl("Cool_1", names(attr(diffs.LSD.comb$predictions, which = "LSDvalues")))))
  testthat::expect_true(any(grepl("Cool_9", rownames(diffs.LSD.comb$LSD), fixed = TRUE)))
  testthat::expect_true(any(grepl("Hot_7", 
                                  names(attr(diffs.LSD.comb$predictions, which = "LSDvalues")), 
                                  fixed = TRUE)))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% 
                              names(attributes(diffs.LSD.comb$backtransforms))))
  testthat::expect_true(all(attr(diffs.LSD.comb, which = "LSDby") == c("Temperature_Genotype")))
  testthat::expect_true(all(attr(diffs.LSD.comb$predictions, which = "LSDby") == c("Temperature_Genotype")))
  testthat::expect_true(all(attr(diffs.LSD.comb$backtransforms, which = "LSDby") == c("Temperature_Genotype")))
  testthat::expect_true(all(c( "transform.power", "offset", "scale") %in% 
                              names(attributes(diffs.LSD.comb$backtransforms))))
  testthat::expect_true(attr(diffs.LSD.comb$backtransforms, which = "transform.power") == 0)
  
  #Test combining a factor in and a factor outside the LSDby
  #Removes Temperature from the LSDby, leaving Genotype only, and results in avsed.tolerance being exceeded, revert to CIs
  diffs.LSD.comb <- facCombine(diffs.new, factors = c("Temperature", "Salinity"))
  testthat::expect_true(all(names(diffs.LSD.comb$predictions)[c(1)] == "Temperature_Salinity"))
  testthat::expect_true(all(names(diffs.LSD.comb$predictions)[c(5:6)] != 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.comb))))
  testthat::expect_true(attr(diffs.LSD.comb, which = "LSDby") == "Genotype")
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.comb$predictions))))
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.comb$backtransforms))))
  testthat::expect_true(all(as.character(1:10) == rownames(diffs.LSD.comb$LSD)))
  
  #Test combining a factor in and a factor outside the LSDby, with avsed.tolerance set to NA
  testthat::expect_error(diffs.LSD.comb <- facCombine(diffs.new, factors = c("Temperature", "Salinity"), 
                                                      avsed.tolerance = NA))
  diffs.LSD.force <- redoErrorIntervals(diffs.new, error.intervals = "half", 
                                        LSDtype = "factor", LSDby = attr(diffs.ClUp, which = "LSDby"),
                                        avsed.tolerance = NA)
  diffs.LSD.comb <- facCombine(diffs.LSD.force, factors = c("Temperature", "Salinity"))
  testthat::expect_true(all(names(diffs.LSD.comb$predictions)[c(1,5:6)] == 
                              c("Temperature_Salinity", 
                                "upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.comb))))
  testthat::expect_true(attr(diffs.LSD.comb, which = "LSDby") == "Genotype")
  testthat::expect_true(all(as.character(1:10) == rownames(diffs.LSD.comb$LSD)))
  testthat::expect_true(attr(diffs.LSD.comb, which = "LSDby") == "Genotype")
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.comb$predictions))))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.comb$backtransforms))))
  testthat::expect_true(all(as.character(1:10) == names(attr(diffs.LSD.comb$predictions, which = "LSDby"))))
  

  #Combine all factors, reducing the LSDby to NULL, reverts to CIs
  diffs.LSD.comb <- facCombine(diffs.new, factors = c("Temperature", "Salinity", "Genotype"))
  testthat::expect_true(all(names(diffs.LSD.comb$predictions)[c(1)] == "Temperature_Salinity_Genotype"))
  testthat::expect_true(all(names(diffs.LSD.comb$predictions)[c(5:6)] != 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic") %in% names(attributes(diffs.LSD.comb))))
  testthat::expect_true(is.null(attr(diffs.LSD.comb, which = "LSDby")))
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.comb$predictions))))
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD.comb$backtransforms))))
  testthat::expect_true(nrow(diffs.LSD.comb$LSD) == 1)

  #Combine all factors, reducing the LSDby to NULL, set avsed.tolerance to NA to avoid reversion to CIs
  diffs.LSD.comb <- facCombine(diffs.LSD.force, factors = c("Temperature", "Salinity", "Genotype"))
  testthat::expect_true(all(names(diffs.LSD.comb$predictions)[c(1,4:5)] == 
                              c("Temperature_Salinity_Genotype", 
                                "upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic") %in% names(attributes(diffs.LSD.comb))))
  testthat::expect_true(is.null(attr(diffs.LSD.comb, which = "LSDby")))
  testthat::expect_true((attr(diffs.LSD.comb, which = "LSDtype") == "overall"))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic") %in% names(attributes(diffs.LSD.comb$predictions))))
  testthat::expect_true(is.null(attr(diffs.LSD.comb$predictions, which = "LSDby")))
  testthat::expect_true((attr(diffs.LSD.comb$predictions, which = "LSDtype") == "overall"))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic") %in% names(attributes(diffs.LSD.comb$backtransforms))))
  testthat::expect_true(is.null(attr(diffs.LSD.comb$backtransforms, which = "LSDby")))
  testthat::expect_true((attr(diffs.LSD.comb$backtransforms, which = "LSDtype") == "overall"))
})


cat("#### Test for LSD on Oats with asreml4\n")
test_that("LSD_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  current.asrt <- as.asrtests(m1.asr)
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  
  #Test single factor linear.transform
  Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", vcov=TRUE)
  Var.diffs <- allDifferences(predictions = Var.pred$pvals,
                              classify = "Nitrogen:Variety", 
                              vcov = Var.pred$vcov, tdf = den.df)
  testthat::expect_true(all("LSDtype" %in% names(attributes(Var.diffs))))
  testthat::expect_true(all(attr(Var.diffs, which = "LSDtype") == "overall"))
  testthat::expect_true(all(attr(Var.diffs, which = "LSDstatistic") == "mean"))
  Var.diffs.one <- linTransform(Var.diffs, linear.transformation = ~Nitrogen,
                                error.intervals = "half", tables = "none")
  testthat::expect_true(all("LSDtype" %in% names(attributes(Var.diffs.one))))
  testthat::expect_true(all(abs(Var.diffs.one$LSD[-match(c("n","accuracyLSD"), names(Var.diffs.one$LSD))] - 
                                  9.883479) < 1e-06))
  testthat::expect_true(all(abs(Var.diffs.one$LSD$assignedLSD - 
                                  attr(Var.diffs.one$predictions, which = "LSDvalues")) < 1E-06))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic", "LSDaccuracy") %in% 
                              names(attributes(Var.diffs.one))))
  testthat::expect_true(all(c("LSDtype", "LSDstatistic", "LSDvalues", "LSDaccuracy", "avsed.tolerance", 
                              "accuracy.threshold") %in% names(attributes(Var.diffs.one$predictions))))
  testthat::expect_true(all(abs(attr(Var.diffs.one$predictions, which = "LSDvalues") - 9.883479) < 1e-05))  
  #Test LSDstatistic = "min"
  Var.diffs.one.min <- redoErrorIntervals(Var.diffs.one, error.intervals = "half", 
                                          LSDtype = "overall", LSDstatistic = "min")
  testthat::expect_true(all(c("LSDtype", "LSDstatistic", "LSDvalues", "LSDaccuracy", "avsed.tolerance", 
                              "accuracy.threshold") %in% names(attributes(Var.diffs.one.min$predictions))))
  testthat::expect_true(attr(Var.diffs.one.min$predictions, which = "LSDstatistic") == "minimum")
  testthat::expect_true(attr(Var.diffs.one.min$predictions, which = "LSDvalues") == Var.diffs.one.min$LSD["minLSD"])
  
  #Test LSDstatistic = "q90"
  Var.diffs.one.q90 <- allDifferences(predictions = Var.pred$pvals,
                                      classify = "Nitrogen:Variety", LSDstatistic = "q90", 
                                      vcov = Var.pred$vcov, tdf = den.df)
  testthat::expect_true(all("LSDtype" %in% names(attributes(Var.diffs.one.q90))))
  testthat::expect_true(all(attr(Var.diffs.one.q90, which = "LSDtype") == "overall"))
  testthat::expect_true(all(attr(Var.diffs.one.q90, which = "LSDstatistic") == "q90"))

  #Test LSDby not in linear.transformation
  testthat::expect_warning(Var.diffs.by <- linTransform(Var.diffs, 
                                                        linear.transformation = ~Nitrogen,
                                                        error.intervals = "half", 
                                                        LSDtype = "factor", 
                                                        LSDby = "Variety", 
                                                        tables = "none"))
  
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic", "LSDaccuracy") %in% 
                              names(attributes(Var.diffs.by))))
  testthat::expect_true(all(nrow(Var.diffs.by$LSD) == 3))
  testthat::expect_true(all(c("LSDtype", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(Var.diffs.by$predictions))))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic", "LSDaccuracy", "avsed.tolerance", 
                               "accuracy.threshold") %in% names(attributes(Var.diffs.by$predictions))))
  testthat::expect_true(all(abs(Var.diffs.by$LSD["meanLSD"] - 
                                  attr(Var.diffs.by$predictions, which = "LSDvalues")) < 1E-06))
  testthat::expect_true(all(abs(attr(Var.diffs.one$predictions, which = "LSDvalues") - 9.883479) < 1e-05))  

  #Test for predictPlus with numeric in classify
  mx.asr <- asreml(Yield ~ xNitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  current.asrt <- as.asrtests(mx.asr)
  
  testthat::expect_equal(length(mx.asr$vparameters),3)
  
  diffs <- predictPlus(mx.asr, classify = "xNitrogen:Variety", 
                       x.num = "xNitrogen", 
                       x.pred.values = sort(unique(Oats.dat$xNitrogen)),
                       wald.tab = current.asrt$wald.tab, 
                       error.intervals = "half", LSDtype = "factor", 
                       LSDby = "xNitrogen", 
                       tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs))
  testthat::expect_equal(nrow(diffs$predictions),12)
  testthat::expect_equal(ncol(diffs$predictions),7)
  testthat::expect_equal(length(attributes(diffs)),12)
  testthat::expect_true(all(nrow(diffs$LSD) == 4))
  testthat::expect_true(all(abs(diffs$LSD["meanLSD"] - diffs$LSD["assignedLSD"]) < 1E-06))
  testthat::expect_true(all(diffs$LSD["accuracyLSD"] < 1E-10))

})

cat("#### Test for sort.alldiffs on Smarthouse with asreml4\n")
test_that("sort.alldiffs4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(Smarthouse.dat)
  
  #Set up without any sorting
  m1.asr <- asreml(y1 ~ Genotype*A*B, 
                   random=~Replicate/Mainplot/Subplot,
                   data=Smarthouse.dat)
  testthat::expect_equal(length(m1.asr$vparameters),4)
  current.asrt <- as.asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m <- current.asrt$asreml.obj
  testthat::expect_equal(length(m$vparameters),3)

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
  testthat::expect_true(all(names(diffs.reord$predictions)[6:7] == 
                              c("upper.StandardError.limit", "lower.StandardError.limit")))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic") %in% names(attributes(diffs.reord))))
  testthat::expect_true(!("LSDby" %in% names(attributes(diffs.reord))))
  testthat::expect_true(attr(diffs.reord, which = "LSDtype") == "overall")
  testthat::expect_true(attr(diffs.reord, which = "LSDstatistic") == "mean")
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.reord$predictions))))
  testthat::expect_true(!any(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.reord$backtransforms))))

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
  testthat::expect_equal(length(attributes(diffs.sort)),13)
  testthat::expect_equal(length(attr(diffs.sort, which = "sortOrder")),10)
  
  #Test sort.alldiffs with supplied sortOrder
  m2.asr <- asreml(y2 ~ Genotype*A*B, 
                   random=~Replicate/Mainplot/Subplot,
                   data=Smarthouse.dat)
  testthat::expect_equal(length(m1.asr$vparameters),4)
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
  testthat::expect_equal(as.character(attr(diffs1.sort, which = "sortOrder")[2]),"Wyalkatchem")
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

cat("#### Test for LSD with sort.alldiffs on Smarthouse with asreml4\n")
test_that("LSDsort.alldiffs4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(Smarthouse.dat)
  
  #Set up without any sorting
  m1.asr <- asreml(y1 ~ Genotype*A*B, 
                   random=~Replicate/Mainplot/Subplot,
                   data=Smarthouse.dat)
  testthat::expect_equal(length(m1.asr$vparameters),4)
  current.asrt <- as.asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m <- current.asrt$asreml.obj
  testthat::expect_equal(length(m$vparameters),3)
  
  #Include LSIs
  diffs <- predictPlus(m, classify = "Genotype:A:B", 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "half", 
                       LSDtype = "factor", LSDby = c("Genotype", "A"), 
                       LSDstatistic = "max", 
                       tables = "none")
  testthat::expect_true(is.alldiffs(diffs))
  testthat::expect_true(validAlldiffs(diffs))
  testthat::expect_equal(nrow(diffs$predictions),120)
  testthat::expect_equal(ncol(diffs$predictions),8)
  testthat::expect_equal(as.character(diffs$predictions$Genotype[1]),"Axe")
  testthat::expect_true(is.null(attr(diffs, which = "sortOrder")))
  testthat::expect_true(all(names(diffs$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs))))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs$predictions))))
  testthat::expect_true(any(grepl("Axe,N1", names(attr(diffs$predictions, which = "LSDvalues")))))
  testthat::expect_true(any(grepl("Axe,N1", rownames(diffs$LSD), fixed = TRUE)))
  testthat::expect_true(any(grepl("Axe,N1", 
                                  names(attr(diffs$predictions, which = "LSDvalues")), 
                                  fixed = TRUE)))
  testthat::expect_true(is.null(diffs$backtransforms))
  testthat::expect_true(all(attr(diffs, which = "LSDby") == c("Genotype", "A")))
  testthat::expect_true(all(attr(diffs$predictions, which = "LSDby") == c("Genotype", "A")))

  #Test re-odering of the classify - no reordering of the LSDs
  diffs.reord <- renewClassify(diffs, newclassify = "A:B:Genotype")
  testthat::expect_true(all(names(diffs.reord$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.reord))))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.reord$predictions))))
  testthat::expect_true(any(grepl("Axe,N1", names(attr(diffs.reord$predictions, which = "LSDvalues")))))
  testthat::expect_true(any(grepl("Axe,N1", rownames(diffs.reord$LSD), fixed = TRUE)))
  testthat::expect_true(all(abs(diffs.reord$LSD$maxLSD[1:2] - c(0.07689865, 0.08286324)) < 1e-05))
  testthat::expect_true(any(grepl("Axe,N1", 
                                  names(attr(diffs.reord$predictions, which = "LSDvalues")), 
                                  fixed = TRUE)))
  testthat::expect_true(all(attr(diffs.reord, which = "LSDby") == c("Genotype", "A")))
  testthat::expect_true(all(attr(diffs.reord$predictions, which = "LSDby") == c("Genotype", "A")))

  diffs.reLSD <- recalcLSD(diffs.reord, LSDtype = "factor", LSDby = c("A", "Genotype"), 
                           LSDstatistic = "min")
  testthat::expect_true(any(grepl("N1,Axe", rownames(diffs.reLSD$LSD), fixed = TRUE)))
  testthat::expect_true(all(abs(diffs.reLSD$LSD$maxLSD[1:2] - c(0.07689865, 0.07689865)) < 1e-05))
  testthat::expect_true(("LSDby" %in% names(attributes(diffs.reLSD))))
  testthat::expect_true(attr(diffs.reLSD, which = "LSDtype") == "factor.combinations")
  testthat::expect_true(attr(diffs.reLSD, which = "LSDstatistic") == "minimum")
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                               names(attributes(diffs.reLSD$predictions))))

  #Test sort.alldiffs and save order for use with other response variables
  diffs.sort <- sort(diffs, sortFactor = "Genotype")
  testthat::expect_true(all(names(diffs.sort$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.sort))))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.sort$predictions))))
  testthat::expect_true(any(grepl("Axe,N1", names(attr(diffs.sort$predictions, which = "LSDvalues")))))
  testthat::expect_true(any(grepl("Axe,N1", rownames(diffs.sort$LSD), fixed = TRUE)))
  testthat::expect_true(all(abs(diffs.reord$LSD$maxLSD[1:2] - c(0.07689865, 0.08286324)) < 1e-05))
  testthat::expect_true(any(grepl("Axe,N1", 
                                  names(attr(diffs.sort$predictions, which = "LSDvalues")), 
                                  fixed = TRUE)))
  testthat::expect_true(all(attr(diffs.sort, which = "LSDby") == c("Genotype", "A")))
  testthat::expect_true(all(attr(diffs.sort$predictions, which = "LSDby") == c("Genotype", "A")))
  
})

cat("#### Test for LSDsupplied on Oats with asreml4\n")
test_that("sort.alldiffs4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  LSD.hdr <- c("n", "minLSD", "meanLSD", "maxLSD", "assignedLSD", "accuracyLSD")
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  testthat::expect_equal(length(m1.asr$vparameters),3)
  current.asrt <- as.asrtests(m1.asr)
  
  #test supplying a median LSD value
  diffs <- predictPlus(m1.asr, classify = "Nitrogen:Variety", 
                       wald.tab = current.asrt$wald.tab, 
                       error.intervals = "half", tables = "none")
  testthat::expect_true(validAlldiffs(diffs))
  LSD.dat <- diffs$LSD
  LSD.dat$assignedLSD <- qt(0.975, attr(diffs, which = "tdf"))*median(diffs$sed, na.rm = TRUE)
  
  diffs.reLSD <- redoErrorIntervals(diffs, error.intervals = "half", 
                                    LSDtype = "supplied", LSDsupplied = LSD.dat["assignedLSD"])
  testthat::expect_is(diffs.reLSD, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs.reLSD))
  testthat::expect_equal(nrow(diffs.reLSD$predictions),12)
  testthat::expect_equal(ncol(diffs.reLSD$predictions),7)
  testthat::expect_true(all(names(diffs.reLSD$predictions)[5:6] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic") %in% names(attributes(diffs.reLSD))))
  testthat::expect_true(is.null(attr(diffs.reLSD, which = "LSDby")))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.reLSD$predictions))))
  testthat::expect_true(rownames(diffs.reLSD$LSD) == "overall")
  testthat::expect_true(all(LSD.hdr %in% names(diffs.reLSD$LSD)))
  testthat::expect_true(all(abs(c(132, 15.47426, 18.54066, 19.56707, 19.56707, 0.2091679) - diffs.reLSD$LSD) < 1e-05))
  #Check limit difference equals the LSD$meanLSD
  testthat::expect_true(all(abs((diffs.reLSD$predictions$upper.halfLeastSignificant.limit - 
                                   diffs.reLSD$predictions$lower.halfLeastSignificant.limit) 
                                - diffs.reLSD$LSD$assignedLSD) < 1e-05))
  
  medianLSD <- LSD.dat[["meanLSD"]]
  names(medianLSD) <- "overall"
  diffs.med.reLSD <- redoErrorIntervals(diffs, error.intervals = "half", 
                                        LSDtype = "supplied", LSDsupplied = medianLSD)
  testthat::expect_is(diffs.med.reLSD, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs.med.reLSD))
  testthat::expect_equal(nrow(diffs.med.reLSD$predictions),12)
  testthat::expect_equal(ncol(diffs.med.reLSD$predictions),7)
  testthat::expect_true(all(names(diffs.med.reLSD$predictions)[5:6] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic") %in% names(attributes(diffs.med.reLSD))))
  testthat::expect_true(is.null(attr(diffs.med.reLSD, which = "LSDby")))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.med.reLSD$predictions))))
  testthat::expect_true(rownames(diffs.med.reLSD$LSD) == "overall")
  testthat::expect_true(all(LSD.hdr %in% names(diffs.med.reLSD$LSD)))
  testthat::expect_true(all(abs(c(132, 15.47426, 18.54066, 19.56707, 18.54066, 0.1653879) - 
                                  diffs.med.reLSD$LSD) < 1e-05))
  #Check limit difference equals the LSD$meanLSD
  testthat::expect_true(all(abs((diffs.med.reLSD$predictions$upper.halfLeastSignificant.limit - 
                                   diffs.med.reLSD$predictions$lower.halfLeastSignificant.limit) 
                                - diffs.med.reLSD$LSD$assignedLSD) < 1e-05))
  
  })

cat("#### Test for LSD on WaterRunoff with asreml4\n")
test_that("LSDWater4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)
  LSD.hdr <- c("n", "minLSD", "meanLSD", "maxLSD", "assignedLSD", "accuracyLSD")
  
  #Analyse pH  
  m1.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                   random = ~ Benches:MainPlots,
                   keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- as.asrtests(m1.asr, NULL, NULL)
  testthat::expect_equal(length(m1.asr$vparameters),2)
  current.asrt <- as.asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m1.asr <- current.asrt$asreml.obj
  testthat::expect_equal(length(m1.asr$vparameters),1)
  
  diffs.full.LSD <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", LSDtype = "factor", 
                                       LSDby = c("Type", "Species"),
                                       tables = "none", Vmatrix = TRUE)
  tdf <- attr(diffs.full.LSD, which = "tdf")
  t.val <- qt(0.975, tdf)
  testthat::expect_true(setequal(names(diffs.full.LSD$predictions), 
                                 c("Sources", "Type", "Species", "predicted.value", 
                                   "standard.error", "upper.halfLeastSignificant.limit", 
                                   "lower.halfLeastSignificant.limit", "est.status")))
  testthat::expect_true(!("LSDwarning" %in% names(diffs.full.LSD$predictions)))
  testthat::expect_true(all(c("LSDtype", "LSDby", "LSDvalues", "avsed.tolerance", 
                              "accuracy.threshold") %in% names(attributes(diffs.full.LSD$predictions))))
  testthat::expect_true((attr(diffs.full.LSD$predictions, which = "avsed.tolerance") == 0.25))
  testthat::expect_true(is.na(attr(diffs.full.LSD$predictions, which = "accuracy.threshold")))
  testthat::expect_true(all(abs(attr(diffs.full.LSD$predictions, which = "LSDvalues") - 
                                  c(0.3052944, 0.3052944, 0.3099100, 0.3076127, 
                                    0.3076127, 0.3739078, 0.3739078, 0.3783488)) < 1e-05))
  
  
  medianLSD.dat <-  data.frame(meanLSD = c(0.3052944, 0.3052944, 0.3121896, 0.3052944, 
                                           0.3052944, 0.3739078, 0.3739078, 0.3739078),
                               row.names = c("Landscape,S. iqscjbogxah", "Landscape,S. oymwjrcnepv", 
                                             "Landscape,S. ocphawvtlgi", "Medicinal,S. hbpgtylxqku", 
                                             "Medicinal,S. orcxszbujml", "Culinary,S. xeqackngdrt", 
                                             "Culinary,S. tkujbvipoyr",  "Control,Non-planted"))
  levs <- strsplit(rownames(medianLSD.dat), ",", fixed = TRUE)
  medianLSDfacs.dat <- data.frame(Type = factor(unlist(lapply(levs, function(lev) lev[1])), 
                                                levels = levels(WaterRunoff.dat$Type)), 
                                  Species = factor(unlist(lapply(levs, function(lev) lev[2])), 
                                                levels = levels(WaterRunoff.dat$Species)),
                                  assignedLSD = medianLSD.dat$meanLSD)
  
  #### Test different LSDaccuracy with overall

  #Test of LSDaccuracy = "maxAbs"
  diffs.over.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", 
                                       avsed.tolerance = NA,
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
  testthat::expect_true(abs(diffs.over.Acc$LSD$accuracyLSD - 0.1860607) < 01e-05)
  
  #Test of LSDaccuracy = "maxDev"
  diffs.over.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", 
                                       LSDaccuracy = "maxDev", avsed.tolerance = NA,
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
  testthat::expect_true(abs(diffs.over.Acc$LSD$accuracyLSD - 0.1860607) < 01e-05)
  
  #Test of LSDaccuracy = "q90Dev"
  diffs.over.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", 
                                       LSDaccuracy = "q90Dev", avsed.tolerance = NA,
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
  testthat::expect_true(abs(diffs.over.Acc$LSD$accuracyLSD - 0.06892194) < 01e-05)
  
  #Test of LSDaccuracy = "Root"
  diffs.over.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", 
                                       LSDaccuracy = "Root", avsed.tolerance = NA,
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
  testthat::expect_true(abs(diffs.over.Acc$LSD$accuracyLSD - 0.06833254) < 01e-05)
  

  #Test LSDstatistic = "q90"
  diffs.reLSD.q90 <- redoErrorIntervals(diffs.full.LSD, error.intervals = "half", 
                                        LSDtype = "factor", LSDby = c("Type", "Species"), 
                                        LSDtstatistic = "q90")
  testthat::expect_true(all(diffs.reLSD.q90$LSD$assignedLSD <= diffs.reLSD.q90$LSD$maxLSD))
  
  testthat::expect_true(all(abs(diffs.reLSD.q90$LSD$assignedLSD - c(0.3052944,0.3052944,0.3099100,0.3076127,
                                                                    0.3076127,0.3739078,0.3739078,0.3783488)) < 1e-05))
    
  #supplied is a data.frame with rownames and LSD suppled 
  diffs.reLSD <- redoErrorIntervals(diffs.full.LSD, error.intervals = "half", 
                                    LSDtype = "supplied", LSDby = c("Type", "Species"), 
                                    LSDsupplied = medianLSD.dat)
  testthat::expect_is(diffs.reLSD, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs.reLSD))
  testthat::expect_equal(nrow(diffs.reLSD$predictions),40)
  testthat::expect_equal(ncol(diffs.reLSD$predictions),8)
  testthat::expect_true(all(names(diffs.reLSD$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDby", "LSDstatistic") %in% 
                              names(attributes(diffs.reLSD))))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic", "LSDby", "LSDvalues") %in% 
                              names(attributes(diffs.reLSD$predictions))))
  testthat::expect_true(all(rownames(diffs.reLSD$LSD) == rownames(medianLSD.dat)))
  testthat::expect_true(all(LSD.hdr %in% names(diffs.reLSD$LSD)))
  testthat::expect_true(all(abs(medianLSD.dat$meanLSD - diffs.reLSD$LSD$assignedLSD) < 1e-05))
  #Check limit difference equals the median LSDs
  tmp <- merge(diffs.reLSD$predictions, medianLSDfacs.dat)
  testthat::expect_true(all(abs((tmp$upper.halfLeastSignificant.limit - tmp$lower.halfLeastSignificant.limit) 
                                - tmp$assignedLSD) < 1e-05))
  
  #supplied is a data.frame with factors and LSDsupplied
  diffs.reLSD <- redoErrorIntervals(diffs.full.LSD, error.intervals = "half", 
                                    LSDtype = "supplied", LSDby = c("Type", "Species"), 
                                    LSDsupplied = medianLSDfacs.dat)
  testthat::expect_is(diffs.reLSD, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs.reLSD))
  testthat::expect_equal(nrow(diffs.reLSD$predictions),40)
  testthat::expect_equal(ncol(diffs.reLSD$predictions),8)
  testthat::expect_true(all(names(diffs.reLSD$predictions)[6:7] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDby", "LSDstatistic") %in% 
                              names(attributes(diffs.reLSD))))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic", "LSDby", "LSDvalues") %in% 
                              names(attributes(diffs.reLSD$predictions))))
  testthat::expect_true(all(rownames(diffs.reLSD$LSD) == rownames(medianLSD.dat)))
  testthat::expect_true(all(LSD.hdr %in% names(diffs.reLSD$LSD)))
  testthat::expect_true(all(abs(medianLSD.dat$meanLSD - diffs.reLSD$LSD$assignedLSD) < 1e-05))
  #Check limit difference equals the median LSDs
  tmp <- merge(diffs.reLSD$predictions, medianLSDfacs.dat)
  testthat::expect_true(all(abs((tmp$upper.halfLeastSignificant.limit - tmp$lower.halfLeastSignificant.limit) 
                                - tmp$assignedLSD) < 1e-05))
  

  ##### Test of new LSD component and accuracy.threshold
  diffs.full.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", 
                                       LSDtype = "factor", 
                                       LSDby = c("Type", "Species"),
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
  testthat::expect_true(setequal(c("names", "classify", "tdf", "alpha", "class", "LSDtype", "LSDby", 
                                   "LSDstatistic", "response", "response.title", "term", 
                                   "LSDaccuracy"), names(attributes(diffs.full.Acc))))
  testthat::expect_true(all(names(diffs.full.Acc$LSD) == LSD.hdr))
  testthat::expect_true(all(diffs.full.Acc$LSD$accuracyLSD < 0.03))
  testthat::expect_true(all(diffs.full.Acc$LSD$meanLSD == diffs.full.Acc$LSD$assignedLSD))
  testthat::expect_true(sum(diffs.full.Acc$predictions$LSDwarning) == 0)
  testthat::expect_true(setequal(c("names", "class", "row.names", "heading", "LSDtype", "LSDby", 
                              "LSDstatistic", "LSDaccuracy", "avsed.tolerance", "accuracy.threshold", 
                              "LSDvalues"), names(attributes(diffs.full.Acc$predictions))))
  testthat::expect_true((attr(diffs.full.Acc$predictions, which = "LSDaccuracy") == "maxAbsDeviation"))
  testthat::expect_true(all(unlist(attributes(diffs.full.Acc)[c("LSDtype", "LSDby", "LSDstatistic", "LSDaccuracy")]) 
                            == unlist(attributes(diffs.full.Acc$predictions)[c("LSDtype", "LSDby", "LSDstatistic", 
                                                                               "LSDaccuracy")])))

  #Test of LSDaccuracy = "Root"
  diffs.full.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                       classify = "Sources:Type:Species", 
                                       wald.tab = current.asrt$wald.tab, 
                                       present = c("Type","Species","Sources"),
                                       error.intervals = "halfLeast", LSDtype = "factor", 
                                       LSDby = c("Type", "Species"), LSDaccuracy = "Root",
                                       accuracy.threshold = 0.10,
                                       tables = "none", Vmatrix = TRUE)
  testthat::expect_true(setequal(c("names", "classify", "tdf", "alpha", "class", "LSDtype", "LSDby", 
                                   "LSDstatistic", "response", "response.title", "term", 
                                   "LSDaccuracy"), names(attributes(diffs.full.Acc))))
  testthat::expect_true(all(names(diffs.full.Acc$LSD) == LSD.hdr))
  testthat::expect_true(all(diffs.full.Acc$LSD$accuracyLSD < 0.02))
  testthat::expect_true(all(diffs.full.Acc$LSD$meanLSD == diffs.full.Acc$LSD$assignedLSD))
  testthat::expect_true(sum(diffs.full.Acc$predictions$LSDwarning) == 0)
  testthat::expect_true(setequal(c("names", "class", "row.names", "heading", "LSDtype", "LSDby", 
                                   "LSDstatistic", "LSDaccuracy", "avsed.tolerance", "accuracy.threshold", 
                                   "LSDvalues"), names(attributes(diffs.full.Acc$predictions))))
  testthat::expect_true((attr(diffs.full.Acc$predictions, which = "LSDaccuracy") == "RootMeanSqDeviation"))
  testthat::expect_true(all(unlist(attributes(diffs.full.Acc)[c("LSDtype", "LSDby", "LSDstatistic", "LSDaccuracy")]) 
                            == unlist(attributes(diffs.full.Acc$predictions)[c("LSDtype", "LSDby", "LSDstatistic", 
                                                                               "LSDaccuracy")])))
  
  #### Testing of different LSDby values with accuracy.threshold
  
  #Test of overall and accuracy.threshold
  diffs.full.all.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                           classify = "Sources:Type:Species", 
                                           wald.tab = current.asrt$wald.tab, 
                                           present = c("Type","Species","Sources"),
                                           error.intervals = "halfLeast", 
                                           LSDtype = "overall", avsed.tolerance = NA,
                                           accuracy.threshold = 0.10,
                                           tables = "none", Vmatrix = TRUE)
  testthat::expect_true(setequal(c("names", "classify", "tdf", "alpha", "class", "LSDtype",  
                                   "LSDstatistic", "response", "response.title", "term", 
                                   "LSDaccuracy"), names(attributes(diffs.full.all.Acc))))
  testthat::expect_true(all(names(diffs.full.all.Acc$LSD) == LSD.hdr))
  testthat::expect_true(all(diffs.full.all.Acc$LSD$accuracyLSD < 0.19))
  testthat::expect_true(all(diffs.full.all.Acc$LSD$meanLSD == diffs.full.all.Acc$LSD$assignedLSD))
  ksed <- na.omit(as.vector(diffs.full.all.Acc$sed))
  testthat::expect_true(abs(diffs.full.all.Acc$LSD$accuracyLSD - 
      max(abs(t.val *ksed - diffs.full.all.Acc$LSD$assignedLSD))/diffs.full.all.Acc$LSD$assignedLSD) < 1e-05)
  testthat::expect_true(sum(diffs.full.all.Acc$predictions$LSDwarning) == 14)
  testthat::expect_true(setequal(apply(diffs.full.all.Acc$sed, 
                                       FUN = function(ksed, t.val, assLSD)
                                         (max(abs(t.val*ksed - assLSD), na.rm = TRUE) / assLSD > 0.10), 
                                       t.val = t.val, assLSD = diffs.full.all.Acc$LSD$assignedLSD, MARGIN = 1), 
                                 diffs.full.all.Acc$predictions$LSDwarning)) #check LSDwarning
  
  testthat::expect_true(setequal(c("names", "class", "row.names", "heading", "LSDtype", 
                                   "LSDstatistic", "LSDaccuracy", "avsed.tolerance", "accuracy.threshold", 
                                   "LSDvalues"), names(attributes(diffs.full.all.Acc$predictions))))
  testthat::expect_true((attr(diffs.full.all.Acc$predictions, which = "LSDaccuracy") == "maxAbsDeviation"))
  testthat::expect_true(all(unlist(attributes(diffs.full.all.Acc)[c("LSDtype", "LSDstatistic", "LSDaccuracy")]) 
                            == unlist(attributes(diffs.full.all.Acc$predictions)[c("LSDtype", "LSDstatistic", 
                                                                               "LSDaccuracy")])))
  

  #Test of factor.comb and accuracy.threshold, when there are some single prediction combinations
  diffs.full.comb.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                            classify = "Sources:Type:Species", 
                                            wald.tab = current.asrt$wald.tab, 
                                            present = c("Type","Species","Sources"),
                                            error.intervals = "halfLeast", 
                                            LSDtype = "factor", 
                                            LSDby = c("Type", "Sources"),
                                            accuracy.threshold = 0.10,
                                            tables = "none", Vmatrix = TRUE)
  testthat::expect_true(setequal(c("names", "classify", "tdf", "alpha", "class", "LSDtype", "LSDby", 
                                   "LSDstatistic", "response", "response.title", "term", 
                                   "LSDaccuracy"), names(attributes(diffs.full.comb.Acc))))
  testthat::expect_true(all(names(diffs.full.comb.Acc$LSD) == LSD.hdr))
  testthat::expect_true(all(is.na(diffs.full.comb.Acc$LSD$accuracyLSD) | diffs.full.comb.Acc$LSD$accuracyLSD) < 0.015)
  testthat::expect_true(all(diffs.full.comb.Acc$LSD$meanLSD == diffs.full.comb.Acc$LSD$assignedLSD))
  testthat::expect_true(sum(diffs.full.comb.Acc$predictions$LSDwarning, na.rm = TRUE) == 0)
  assRL <-  unlist(lapply(1:3, function(k, diffs, t.value) 
    max(abs(t.value * diffs$sed[k,1:3] - diffs$LSD$assignedLSD[1]), na.rm = TRUE)/diffs$LSD$assignedLSD[1],
    diffs = diffs.full.comb.Acc, t.value = t.val))
  testthat::expect_true(setequal(assRL > 0.10, diffs.full.comb.Acc$predictions$LSDwarning[1:3]))
  testthat::expect_true(setequal(c("names", "class", "row.names", "heading", "LSDtype", "LSDby", 
                                   "LSDstatistic", "LSDaccuracy", "avsed.tolerance", "accuracy.threshold", 
                                   "LSDvalues"), names(attributes(diffs.full.comb.Acc$predictions))))
  testthat::expect_true((attr(diffs.full.comb.Acc$predictions, which = "LSDaccuracy") == "maxAbsDeviation"))
  testthat::expect_true(all(unlist(attributes(diffs.full.comb.Acc)[c("LSDtype", "LSDby", "LSDstatistic", "LSDaccuracy")]) 
                            == unlist(attributes(diffs.full.comb.Acc$predictions)[c("LSDtype", "LSDby", 
                                                                                    "LSDstatistic", "LSDaccuracy")])))

  #Test of supplied and accuracy.threshold, when there are some single prediction combinations
  LSDsupp <- rep(0.3, 20)
  names(LSDsupp) <- levels(with(WaterRunoff.dat, fac.combine(list(Type, Sources), combine.levels = TRUE)))
  diffs.full.supp.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                            classify = "Sources:Type:Species", 
                                            wald.tab = current.asrt$wald.tab, 
                                            present = c("Type","Species","Sources"),
                                            error.intervals = "halfLeast", 
                                            LSDtype = "supplied", 
                                            LSDsupplied = LSDsupp,
                                            LSDby = c("Type", "Sources"),
                                            accuracy.threshold = 0.10,
                                            tables = "none", Vmatrix = TRUE)
  testthat::expect_true(setequal(c("names", "classify", "tdf", "alpha", "class", "LSDtype", "LSDby", 
                                   "LSDstatistic", "response", "response.title", "term", 
                                   "LSDaccuracy"), names(attributes(diffs.full.supp.Acc))))
  testthat::expect_true(all(names(diffs.full.supp.Acc$LSD) == LSD.hdr))
  testthat::expect_true(all(is.na(diffs.full.supp.Acc$LSD$accuracyLSD) | diffs.full.supp.Acc$LSD$accuracyLSD < 0.25))
  testthat::expect_true(all(diffs.full.supp.Acc$LSD$assignedLSD == 0.3))
  testthat::expect_true(sum(diffs.full.supp.Acc$predictions$LSDwarning, na.rm = TRUE) == 4)
  assRL <- unlist(lapply(1:3, function(k, diffs, t.value) 
    max(abs(t.value * diffs$sed[k,1:3] - diffs$LSD$assignedLSD[1]), na.rm = TRUE)/diffs$LSD$assignedLSD[1],
    diffs = diffs.full.supp.Acc, t.value = t.val))
  testthat::expect_true(setequal(assRL > 0.10, diffs.full.supp.Acc$predictions$LSDwarning[1:3]))
  testthat::expect_true(setequal(c("names", "class", "row.names", "heading", "LSDtype", "LSDby", 
                                   "LSDstatistic", "LSDaccuracy", "avsed.tolerance", "accuracy.threshold", 
                                   "LSDvalues"), names(attributes(diffs.full.supp.Acc$predictions))))
  testthat::expect_true((attr(diffs.full.supp.Acc$predictions, which = "LSDaccuracy") == "maxAbsDeviation"))
  testthat::expect_true(all(unlist(attributes(diffs.full.supp.Acc)[c("LSDtype", "LSDby", "LSDstatistic", "LSDaccuracy")]) 
                            == unlist(attributes(diffs.full.supp.Acc$predictions)[c("LSDtype", "LSDby", 
                                                                                    "LSDstatistic", "LSDaccuracy")])))
  
  #Test of per.prediction and accuracy.threshold
  diffs.full.ppred.Acc <- predictPlus.asreml(asreml.obj = m1.asr, 
                                             classify = "Sources:Type:Species", 
                                             wald.tab = current.asrt$wald.tab, 
                                             present = c("Type","Species","Sources"),
                                             error.intervals = "halfLeast", 
                                             LSDtype = "per.pred", 
                                             accuracy.threshold = 0.10,
                                             tables = "none", Vmatrix = TRUE)
  testthat::expect_true(setequal(c("names", "classify", "tdf", "alpha", "class", "LSDtype",  
                                   "LSDstatistic", "response", "response.title", "term", 
                                   "LSDaccuracy"), names(attributes(diffs.full.ppred.Acc))))
  testthat::expect_true(all(names(diffs.full.ppred.Acc$LSD) == LSD.hdr))
  testthat::expect_true(all(unlist(attributes(diffs.full.ppred.Acc)[c("LSDtype", "LSDby", "LSDstatistic", "LSDaccuracy")]) 
                            == unlist(attributes(diffs.full.ppred.Acc$predictions)[c("LSDtype", "LSDstatistic", 
                                                                                     "LSDaccuracy")])))
  testthat::expect_equal(nrow(diffs.full.ppred.Acc$LSD),40)
  testthat::expect_true(abs(t.val*sqrt(mean(diffs.full.ppred.Acc$sed[1,]^2, na.rm = TRUE)) - 
                              diffs.full.ppred.Acc$LSD$assignedLSD[1]) < 1e-05) #check assigned LSD
  testthat::expect_true(abs(max(abs(t.val*diffs.full.ppred.Acc$sed[1,] - 
                                      diffs.full.ppred.Acc$LSD$assignedLSD[1]), na.rm = TRUE) / 
                              diffs.full.ppred.Acc$LSD$assignedLSD[1] -
                              diffs.full.ppred.Acc$LSD$accuracyLSD[1]) < 1e-05) #Test accuracyLSD
  testthat::expect_true(all(diffs.full.ppred.Acc$LSD$accuracyLSD < 0.15))
  testthat::expect_equal(sum(diffs.full.ppred.Acc$predictions$LSDwarning), 39)
  
})

cat("#### Test for exploreLSDs on WaterRunoff with asreml4\n")
test_that("exploreLSDWater4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)
  
  #Analyse pH  
  m1.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                   random = ~ Benches:MainPlots,
                   keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- as.asrtests(m1.asr, NULL, NULL)
  testthat::expect_equal(length(m1.asr$vparameters),2)
  current.asrt <- as.asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  current.asr <- current.asrt$asreml.obj
  
  TS.diffs <- predictPlus(classify = "Sources:Type", 
                          asreml.obj = current.asr, 
                          wald.tab = current.asrt$wald.tab, 
                          present = c("Sources", "Type", "Species"),
                          tables = "none")
  
  ##Explore the LSD values for predictions obtained using asreml or lmerTest  
  LSDstat <- exploreLSDs(TS.diffs, LSDtype = "factor.combinations", LSDby = "Sources")
  testthat::expect_equal(names(LSDstat), c("frequencies", "distinct.vals", "statistics", "accuracy", 
                                           "per.pred.accuracy", "LSD"))
  testthat::expect_true(all(lapply(c("frequencies", "statistics", "accuracy"), 
                            function(k, LSDstat) nrow(LSDstat[[k]]), LSDstat = LSDstat) == 6))
  testthat::expect_true(all(unlist(lapply(c("frequencies", "statistics", "accuracy"), function(k, LSDstat, dat) 
    all(rownames(LSDstat[[k]])== levels(WaterRunoff.dat$Sources)), LSDstat = LSDstat, dat = WaterRunoff.dat))))
  testthat::expect_equal(names(LSDstat$frequencies), as.character(seq(0.17, 0.39, 0.02)))
  testthat::expect_equal(LSDstat$distinct.vals$`Rain+Basalt`, c(0.197,0.294,0.307))
  testthat::expect_true(all(abs(LSDstat$statistics[1,] - 
                                  c(6, 0.1982634,0.1982634,0.2708314,0.2945287,0.3065833,0.3065833)) < 1e-05))
  testthat::expect_true(all(abs(LSDstat$accuracy[1,] - 
                                  c(6, 0.5463438,0.5463438,0.2679453,0.3268454,0.3533133,0.3533133)) < 1e-05))
  testthat::expect_true(all(lapply(c("per.pred.accuracy", "LSD"), 
                                   function(k, LSDstat) nrow(LSDstat[[k]]), LSDstat = LSDstat) == 20))
  testthat::expect_equal(rownames(LSDstat$per.pred.accuracy), 
                         as.character(fac.combine(as.list(TS.diffs$predictions[c("Sources","Type")]), 
                                                  combine.levels = TRUE)))
  testthat::expect_true(all(abs(LSDstat$per.pred.accuracy[1,] - 
                                  c(0.4855427,0.4855427,0.2679453,0.3268454,0.3533133,0.3533133)) < 1e-05))
  
  LSDstat <- exploreLSDs(TS.diffs, LSDtype = "factor.combinations", LSDby = "Sources", LSDaccuracy = "maxDev")
  testthat::expect_equal(names(LSDstat), c("frequencies", "distinct.vals", "statistics", "accuracy", 
                                           "per.pred.accuracy", "LSD"))
  testthat::expect_true(all(lapply(c("frequencies", "statistics", "accuracy"), 
                                   function(k, LSDstat) nrow(LSDstat[[k]]), LSDstat = LSDstat) == 6))
  testthat::expect_true(all(unlist(lapply(c("frequencies", "statistics", "accuracy"), function(k, LSDstat, dat) 
    all(rownames(LSDstat[[k]])== levels(WaterRunoff.dat$Sources)), LSDstat = LSDstat, dat = WaterRunoff.dat))))
  testthat::expect_equal(names(LSDstat$frequencies), as.character(seq(0.17, 0.39, 0.02)))
  testthat::expect_equal(LSDstat$distinct.vals$`Rain+Basalt`, c(0.197,0.294,0.307))
  testthat::expect_true(all(abs(LSDstat$statistics[1,] - 
                                  c(6,0.1982634,0.1982634,0.2708314,0.2945287,0.3065833,0.3065833)) < 1e-05))
  testthat::expect_true(all(abs(LSDstat$accuracy[1,] - 
                                  c(6,0.5463438,0.5463438,0.1320083,0.04092854,0,0)) < 1e-05))
  testthat::expect_true(all(lapply(c("per.pred.accuracy", "LSD"), 
                                   function(k, LSDstat) nrow(LSDstat[[k]]), LSDstat = LSDstat) == 20))
  testthat::expect_equal(rownames(LSDstat$per.pred.accuracy), 
                         as.character(fac.combine(as.list(TS.diffs$predictions[c("Sources","Type")]), 
                                                  combine.levels = TRUE)))
  testthat::expect_true(all(abs(LSDstat$per.pred.accuracy[1,] - 
                                  c(0.4855427,0.4855427,0.08749856,0,-0.03931926,-0.03931926) < 1e-05)))
  
  LSD.dat <- plotLSDs(TS.diffs, factors.per.grid = 1)
  testthat::expect_equal(nrow(LSD.dat),400)
  testthat::expect_equal(levels(LSD.dat$X1),rownames(TS.diffs$sed))
  testthat::expect_equal(length(unique(signif(LSD.dat$LSD, digits = 4))), 28)
})

cat("#### Test for exploreLSDs on Oats with asreml4\n")
test_that("exploreLSDOatsr4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  current.asrt <- as.asrtests(m1.asr)
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  
  #Test single factor linear.transform
  Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", vcov=TRUE)
  Var.diffs <- allDifferences(predictions = Var.pred$pvals,
                              classify = "Nitrogen:Variety", 
                              vcov = Var.pred$vcov, tdf = den.df)
  
  lsd <- exploreLSDs(Var.diffs)  
  testthat::expect_equal(names(lsd), c("frequencies", "distinct.vals", "statistics", "accuracy", 
                                       "per.pred.accuracy", "LSD"))
  testthat::expect_true(all(lapply(c("statistics", "accuracy"), function(k, lsd) nrow(lsd[[k]]), lsd = lsd) == 1))
  testthat::expect_equal(names(lsd$frequencies), as.character(seq(17.25, 21.75, 0.5)))
  testthat::expect_equal(lsd$distinct.vals, c(17.1, 21.6))
  testthat::expect_true(all(abs(lsd$statistics[1,] - 
                                  c(132,17.11869,17.11869,20.51095,21.64642,21.64642,21.64642)) < 1e-05))
  testthat::expect_true(all(abs(lsd$accuracy[1,] - 
                                  c(132,0.2644909,0.2644909,0.1653879,0.2091679,0.2091679,0.2091679)) < 1e-05))
  testthat::expect_true(all(lapply(c("per.pred.accuracy", "LSD"), function(k, lsd) nrow(lsd[[k]]), lsd = lsd) == 12))
  testthat::expect_equal(rownames(lsd$per.pred.accuracy), 
                         as.character(fac.combine(as.list(Var.diffs$predictions[c("Nitrogen","Variety")]), 
                                                  combine.levels = TRUE)))
  testthat::expect_true(all(abs(lsd$per.pred.accuracy[1,] - 
                                  c(0.2644909,0.2644909,0.1653879,0.2091679,0.2091679,0.2091679)) < 1e-05))
  
  lsd <- exploreLSDs(Var.diffs, LSDtype = "fact", LSDby = "Nitrogen")  
  testthat::expect_equal(names(lsd), c("frequencies", "distinct.vals", "statistics", "accuracy", 
                                       "per.pred.accuracy", "LSD"))
  testthat::expect_true(all(lapply(c("frequencies", "statistics", "accuracy"), 
                                   function(k, lsd) nrow(lsd[[k]]), lsd = lsd) == 4))
  testthat::expect_equal(names(lsd$frequencies), as.character(seq(17.25, 21.75, 0.5)))
  testthat::expect_true(all(unlist(lapply(lsd$distinct.vals, function(val) val == 21.6))))
  testthat::expect_true(all(abs(lsd$statistics[1,-1] - 21.64642) < 1e-05))
  testthat::expect_true(all(lsd$accuracy[1,-1] < 1e-08))
  testthat::expect_true(all(lapply(c("per.pred.accuracy", "LSD"), function(k, lsd) nrow(lsd[[k]]), lsd = lsd) == 12))
  testthat::expect_equal(rownames(lsd$per.pred.accuracy), 
                         as.character(fac.combine(as.list(Var.diffs$predictions[c("Nitrogen","Variety")]), 
                                                  combine.levels = TRUE)))
  testthat::expect_true(all(lsd$per.pred.accuracy[1,] < 1e-08))

  LSD.dat <- plotLSDs((Var.diffs))
  testthat::expect_equal(nrow(LSD.dat),144)
  testthat::expect_equal(levels(LSD.dat$X1),rownames(Var.diffs$sed))
  testthat::expect_true(all(abs(na.omit(LSD.dat$LSD) - 21.64642) < 1e-05 | 
                              abs(na.omit(LSD.dat$LSD) - 17.11869) < 1e-05))
})  

cat("#### Test for sort.alldiffs on WaterRunoff with asreml4\n")
test_that("sort.alldiffsWater4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)
  
  #Analyse pH  
  m1.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- as.asrtests(m1.asr, NULL, NULL)
  testthat::expect_equal(length(m1.asr$vparameters),2)
  current.asrt <- as.asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m1.asr <- current.asrt$asreml.obj
  testthat::expect_equal(length(m1.asr$vparameters),1)
  
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
  testthat::expect_equal(length(attributes(TS.diffs.sort)),13)
  testthat::expect_equal(length(attr(TS.diffs.sort, which = "sortOrder")),6)
  
  #Test sort.alldiffs with supplied sortOrder
  m2.asr <- asreml(fixed = Turbidity ~ Benches + (Sources * (Type + Species)), 
                   random = ~ Benches:MainPlots,
                   keep.order=TRUE, data= WaterRunoff.dat)
  testthat::expect_equal(length(m2.asr$vparameters),2)
  current.asrt <- as.asrtests(m2.asr)
  diffs2.sort <- predictPlus(m2.asr, classify = "Sources:Type", 
                             pairwise = FALSE, Vmatrix = TRUE, error.intervals = "Stand", 
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
  testthat::expect_true(all(c("LSDtype", "LSDby") %in% names(attributes(diffs.fit))))
  testthat::expect_true(attr(diffs.fit, which = "LSDtype") == "factor.combinations")
  testthat::expect_true(attr(diffs.fit, which = "LSDby") == "Type")
  testthat::expect_true(all(c("LSDtype", "LSDby", "LSDvalues") %in% names(attributes(diffs.fit$predictions))))
  testthat::expect_true(all(abs(attr(diffs.fit$predictions, which = "LSDvalues") - 
                                  c(0.1771545, 0.2175130, 0.2643927, 0.3783488)) < 1e-05))
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

cat("#### Test for sort.alldiffs on Oats with asreml4\n")
test_that("sort.alldiffs4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  testthat::expect_equal(length(m1.asr$vparameters),3)
  current.asrt <- as.asrtests(m1.asr)
  
  #Test for as.alldiffs
  Var.pred <- predict(m1.asr, classify="Nitrogen:Variety", sed=TRUE)
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  Var.diffs <- as.alldiffs(predictions = Var.pred$pvals, 
                           sed = Var.pred$sed, 
                           tdf = den.df, classify = "Nitrogen:Variety")
  testthat::expect_true(validAlldiffs(Var.diffs))
  testthat::expect_equal(length(attributes(Var.diffs)),5)
  testthat::expect_true(is.null(attr(Var.diffs, which = "sortOrder")))
  
  #Test for predictPlus no sorting
  diffs <- predictPlus(m1.asr, classify = "Nitrogen:Variety", 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs))
  testthat::expect_equal(nrow(diffs$predictions),12)
  testthat::expect_equal(ncol(diffs$predictions),7)
  testthat::expect_equal(as.character(diffs$predictions$Variety[1]),"Victory")
  testthat::expect_equal(length(attributes(diffs)),11)
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
  
  testthat::expect_equal(length(mx.asr$vparameters),3)
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
  testthat::expect_equal(length(attributes(diffs$xNitrogen.Variety)),13)
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
  testthat::expect_equal(length(attributes(diffs)),13)
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

cat("#### Test for subset.alldiffs on Smarthouse with asreml4\n")
test_that("subset.alldiffs4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(Smarthouse.dat)
  
  #Run analysis and produce diffs object
  m1.asr <- asreml(y1 ~ Genotype*A*B, 
                   random=~Replicate/Mainplot/Subplot,
                   data=Smarthouse.dat)
  testthat::expect_equal(length(m1.asr$vparameters),4)
  current.asrt <- as.asrtests(m1.asr)
  current.asrt <- rmboundary(current.asrt)
  m <- current.asrt$asreml.obj
  testthat::expect_equal(length(m$vparameters),3)
  
  diffs <- predictPlus(m, classify = "Genotype:A:B", 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs))
  testthat::expect_equal(nrow(diffs$predictions),120)
  testthat::expect_equal(ncol(diffs$predictions),8)
  testthat::expect_equal(as.character(diffs$predictions$Genotype[1]),"Axe")
  testthat::expect_equal(length(attributes(diffs)),11)
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
  testthat::expect_equal(length(attributes(diffs.subs)),11)
  
  #Test subset with removal of vars
  diffs.subs <- subset(diffs, subset = A == "N1" & B == "D2", rmClassifyVars = c("A","B"))
  testthat::expect_is(diffs.subs, "alldiffs")
  testthat::expect_true(validAlldiffs(diffs.subs))
  testthat::expect_equal(nrow(diffs.subs$predictions),10)
  testthat::expect_equal(ncol(diffs.subs$predictions),6)
  testthat::expect_false(any(c("A","B") %in% names(diffs.subs$predictions)))
  
  data(WaterRunoff.dat)
  #Run analysis and produce alldiffs object
  asreml.options(keep.order = TRUE) #required for asreml4 only
  current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        data= WaterRunoff.dat)
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
  testthat::expect_equal(length(attributes(diffs.subs)),11)
})

cat("#### Test for facCombine.alldiffs on Ladybird with asreml4\n")
test_that("facCombine.alldiffs4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data("Ladybird.dat")
  
  #Mixed model analysis of logits 
  m1.asr <- asreml(logitP ~ Host*Cadavers*Ladybird, 
                   random = ~ Run,
                   data = Ladybird.dat)

  testthat::expect_equal(length(m1.asr$vparameters),2)
  current.asrt <- as.asrtests(m1.asr)
  testthat::expect_true(validAsrtests(current.asrt))

  HCL.pred <- asreml::predict.asreml(m1.asr, classify="Host:Cadavers:Ladybird", 
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
  HCL.diffs <- allDifferences(predictions = HCL.preds, classify = "Host:Cadavers:Ladybird", 
                              sed = HCL.sed, vcov = HCL.vcov, tdf = den.df)
  
  ## check the class and validity of the alldiffs object
  is.alldiffs(HCL.diffs)
  validAlldiffs(HCL.diffs)
  testthat::expect_equal(nrow(HCL.diffs$predictions),12)
  testthat::expect_equal(ncol(HCL.diffs$predictions),9)
  testthat::expect_true(all(c("LSDtype", "LSDstatistic") %in% names(attributes(HCL.diffs))))
  testthat::expect_equal(attr(HCL.diffs, which = "LSDtype"), "overall")
  testthat::expect_equal(attr(HCL.diffs, which = "LSDstatistic"), "mean")
  testthat::expect_true(!is.null(HCL.diffs$LSD))
  testthat::expect_equal(nrow(HCL.diffs$LSD), 1)
  testthat::expect_true(all(abs(HCL.diffs$LSD[c("minLSD", "meanLSD", "maxLSD", "assignedLSD")] - 0.5534673) < 1e-05))
  testthat::expect_true(!any(c("LSDtype", "LSDstatistic", "LSDvalues") %in% 
                               names(attributes(HCL.diffs$predictions))))
  
  ## Combine Cadavers and Ladybird
  Comb.diffs <- facCombine(HCL.diffs, factors = c("Cadavers","Ladybird"))
  testthat::expect_true(all(c("LSDtype", "LSDstatistic") %in% names(attributes(Comb.diffs))))
  testthat::expect_equal(attr(Comb.diffs, which = "LSDtype"), "overall")
  testthat::expect_equal(attr(Comb.diffs, which = "LSDstatistic"), "mean")
  testthat::expect_true(!is.null(Comb.diffs$LSD))
  testthat::expect_equal(nrow(Comb.diffs$LSD), 1)
  testthat::expect_true(all(abs(Comb.diffs$LSD[c("minLSD", "meanLSD", "maxLSD", "assignedLSD")] - 0.5534673) < 1e-05))
  
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
  testthat::expect_true(all(c("LSDtype", "LSDstatistic") %in% names(attributes(HCL.rename.diffs))))
  testthat::expect_equal(attr(HCL.rename.diffs, which = "LSDtype"), "overall")
  testthat::expect_equal(attr(HCL.rename.diffs, which = "LSDstatistic"), "mean")
  testthat::expect_true(!is.null(HCL.rename.diffs$LSD))
  testthat::expect_equal(nrow(HCL.rename.diffs$LSD), 1)
  testthat::expect_true(all(abs(HCL.rename.diffs$LSD[c("minLSD", "meanLSD", "maxLSD", "assignedLSD")] - 0.5534673) < 1e-05))

  #Change to half-LSI intervals
  HCL.diffs.LSI <- redoErrorIntervals(HCL.diffs, error.intervals = "half")
  
  testthat::expect_true(is.alldiffs(HCL.diffs.LSI))
  testthat::expect_true(validAlldiffs(HCL.diffs.LSI))
  testthat::expect_equal(nrow(HCL.diffs.LSI$predictions),12)
  testthat::expect_equal(ncol(HCL.diffs.LSI$predictions),9)
  testthat::expect_true(all(c("LSDtype", "LSDstatistic") %in% names(attributes(HCL.diffs.LSI))))
  testthat::expect_equal(attr(HCL.diffs.LSI, which = "LSDtype"), "overall")
  testthat::expect_equal(attr(HCL.diffs.LSI, which = "LSDstatistic"), "mean")
  testthat::expect_true(!is.null(HCL.diffs.LSI$LSD))
  testthat::expect_equal(nrow(HCL.diffs.LSI$LSD), 1)
  testthat::expect_true(all(abs(HCL.diffs.LSI$LSD[c("minLSD", "meanLSD", "maxLSD", "assignedLSD")] - 0.5534673) < 1e-05))
  testthat::expect_true(all(c("LSDtype", "LSDstatistic", "LSDvalues") %in% 
                               names(attributes(HCL.diffs.LSI$predictions))))
  testthat::expect_equal(attr(HCL.diffs.LSI$predictions, which = "LSDtype"), "overall")
  testthat::expect_equal(attr(HCL.diffs.LSI$predictions, which = "LSDstatistic"), "mean")
  testthat::expect_true(all(abs(attr(HCL.diffs.LSI$predictions, which = "LSDvalues") - 0.5534673) < 1e-05))
  
  ## Combine Cadavers and Ladybird
  Comb.diffs.LSI <- facCombine(HCL.diffs.LSI, factors = c("Cadavers","Ladybird"))
  testthat::expect_true(validAlldiffs(Comb.diffs.LSI))
  testthat::expect_equal(nrow(Comb.diffs.LSI$predictions),12)
  testthat::expect_equal(ncol(Comb.diffs.LSI$predictions),8)
  testthat::expect_true(all(c("LSDtype", "LSDstatistic") %in% names(attributes(Comb.diffs.LSI))))
  testthat::expect_equal(attr(Comb.diffs.LSI, which = "LSDtype"), "overall")
  testthat::expect_equal(attr(Comb.diffs.LSI, which = "LSDstatistic"), "mean")
  testthat::expect_true(!is.null(Comb.diffs.LSI$LSD))
  testthat::expect_equal(nrow(Comb.diffs.LSI$LSD), 1)
  testthat::expect_true(all(abs(Comb.diffs.LSI$LSD[c("minLSD", "meanLSD", "maxLSD", "assignedLSD")] - 0.5534673) < 1e-05))
  
  testthat::expect_true(all(c("Host", "Cadavers_Ladybird", "predicted.value") %in% 
                              names(Comb.diffs.LSI$predictions)))
  testthat::expect_true(all(c("LSDtype", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(HCL.diffs.LSI$predictions))))
  testthat::expect_equal(attr(Comb.diffs.LSI$predictions, which = "LSDtype"), "overall")
  testthat::expect_equal(attr(Comb.diffs.LSI$predictions, which = "LSDstatistic"), "mean")
  testthat::expect_true(all(abs(attr(Comb.diffs.LSI$predictions, which = "LSDvalues") - 0.5534673) < 1e-05))
  
  ## Rename Cadavers
  HCL.rename.diffs.LSI <- facRename(HCL.diffs.LSI, factor.names = "Cadavers", newnames = "Cadaver.nos")
  testthat::expect_true(validAlldiffs(HCL.rename.diffs.LSI))
  testthat::expect_true("Cadaver.nos" %in% names(HCL.rename.diffs.LSI$predictions))
  testthat::expect_true(all(c("LSDtype", "LSDstatistic") %in% names(attributes(HCL.rename.diffs.LSI))))
  testthat::expect_equal(attr(HCL.rename.diffs.LSI, which = "LSDtype"), "overall")
  testthat::expect_equal(attr(HCL.rename.diffs.LSI, which = "LSDstatistic"), "mean")
  testthat::expect_true(!is.null(HCL.rename.diffs.LSI$LSD))
  testthat::expect_equal(nrow(HCL.rename.diffs.LSI$LSD), 1)
  testthat::expect_true(all(abs(HCL.rename.diffs.LSI$LSD[c("minLSD", "meanLSD", "maxLSD", "assignedLSD")] - 0.5534673) < 1e-05))
  
})

cat("#### Test for facRecast.alldiffs on Ladybird with asreml4\n")
test_that("facRecast.alldiffs4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data("Ladybird.dat")
  
  #Mixed model analysis of logs 
  Ladybird.dat$log.P <- log(Ladybird.dat$Prop*100 + 1)
  m1.asr <- do.call(asreml, 
                    args = list(fixed = log.P ~ Host*Cadavers*Ladybird, 
                                random = ~ Run,
                                data = Ladybird.dat))
  testthat::expect_equal(length(m1.asr$vparameters),2)
  current.asrt <- as.asrtests(m1.asr)
  testthat::expect_true(validAsrtests(current.asrt))
  
  HCL.diffs <- predictPlus(m1.asr, classify = "Host:Cadavers:Ladybird", tables = "none", 
                           wald.tab = current.asrt$wald.tab, transform.power = 0, offset = 1)
  
  ## Recast Ladybird
  HCL.recast.diffs <- facRecast(HCL.diffs, factor = "Ladybird", newlabels = c("none", "present"))
  testthat::expect_true(validAlldiffs(HCL.recast.diffs))
  testthat::expect_true(all(levels(HCL.recast.diffs$predictions$Ladybird) == c("none", "present")))
  testthat::expect_true(all(levels(HCL.recast.diffs$backtransforms$Ladybird) == c("none", "present")))
  testthat::expect_true(all(abs((exp(HCL.recast.diffs$predictions$predicted.value) - 1) -
                                  HCL.recast.diffs$backtransforms$backtransformed.predictions) < 1e-05))
  testthat::expect_true(all(rownames(HCL.recast.diffs$differences)[1:2] == c("bean,5,none", "bean,5,present")))
  HCL.recast.diffs <- facRecast.alldiffs(HCL.recast.diffs, factor = "Host", 
                                         levels.order = c("trefoil", "bean"))
  testthat::expect_true(validAlldiffs(HCL.recast.diffs))
  testthat::expect_true(all(levels(HCL.recast.diffs$predictions$Host) == c("trefoil", "bean")))
  testthat::expect_true(all(levels(HCL.recast.diffs$backtransforms$Host) == c("trefoil", "bean")))
  testthat::expect_true(all(rownames(HCL.recast.diffs$differences)[1:2] == c("trefoil,5,none", "trefoil,5,present")))
  testthat::expect_true(all(abs((exp(HCL.recast.diffs$predictions$predicted.value) - 1) -
                                  HCL.recast.diffs$backtransforms$backtransformed.predictions) < 1e-05))
  HCL.recast.diffs <- facRecast(HCL.recast.diffs, factor = "Ladybird", 
                                levels.order = c("present", "none"), 
                                newlabels = c("yes","no"))
  testthat::expect_true(validAlldiffs(HCL.recast.diffs))
  testthat::expect_true(all(levels(HCL.recast.diffs$predictions$Ladybird) == c("yes", "no")))
  testthat::expect_true(all(levels(HCL.recast.diffs$backtransforms$Ladybird) == c("yes", "no")))
  testthat::expect_true(all(rownames(HCL.recast.diffs$differences)[1:2] == c("trefoil,5,yes", "trefoil,5,no")))
  testthat::expect_true(all(abs((exp(HCL.recast.diffs$predictions$predicted.value) - 1) -
                                  HCL.recast.diffs$backtransforms$backtransformed.predictions) < 1e-05))
  
})

cat("#### Test for linear.transformation on Oats with asreml4\n")
test_that("linear.transform_Oats_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  testthat::expect_equal(length(m1.asr$vparameters),3)
  current.asrt <- as.asrtests(m1.asr)
  
  #Test store of vcov by predictPlus
  diffs <- predictPlus(m1.asr, classify = "Nitrogen:Variety", Vmatrix = TRUE, 
                       wald.tab = current.asrt$wald.tab,
                       error.intervals = "Stand", tables = "none")
  testthat::expect_is(diffs, "alldiffs")
  testthat::expect_equal(length(attr(diffs$predictions, which = "heading")),4)
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
  testthat::expect_equal(length(attr(diffs.L$predictions, which = "heading")),5)
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
  testthat::expect_equal(length(attr(diffs.mod$predictions, which = "heading")),5)
  testthat::expect_true("asreml.predict" %in% class(diffs.mod$predictions))
  testthat::expect_true(is.null(diffs.mod$vcov))
  
  m2.asr <- asreml(Yield ~ Nitrogen+Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  preds <- predict(m2.asr, classify = "Nitrogen:Variety")$pvals 
  testthat::expect_true(all(abs(diffs.mod$predictions$predicted.value - 
                                  preds$predicted.value) < 1e-04))
  testthat::expect_true(abs(diffs.mod$predictions$standard.error[1] - 
                              preds$std.error[1]) > 0.01)
})
  
cat("#### Test for linear.transformation on WaterRunoff with asreml4\n")
test_that("linear.transform_WaterRunoff_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  #Test example in manual
  data(WaterRunoff.dat)
  #Run analysis and produce alldiffs object
  asreml.options(keep.order = TRUE) #required for asreml4 only
  current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        data= WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  diffs <- predictPlus(classify = "Sources:Species", Vmatrix = TRUE, 
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
  ksed <- ksed[!(ksed < 1e-08)]
  testthat::expect_true(length(ksed) == diffs.L$LSD["n"] && diffs.L$LSD["n"] == 968)
  testthat::expect_true(all(abs(diffs.L$predictions$predicted.value[c(1,6,11,16,21,28)] - 
                                  (diffs.sub$predictions$predicted.value[1] - 
                                     diffs.sub$predictions$predicted.value[6])) < 1e-06))
  
  #More efficient version for manual
  data(WaterRunoff.dat)
    asreml.options(keep.order = TRUE) #required for asreml4 only
  current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        data= WaterRunoff.dat)
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
  mod <- do.call("asreml", 
                 list(y ~ A*B, data = twoway))
  pred.full <- predict(mod, classify = "A:B")
  testthat::expect_true(any(is.na(pred.full$pvals$predicted.value)))
  
  #Make predictions additive
  diffs.sub <- predictPlus.asreml(asreml.obj = mod, classify = "A:B", 
                                  linear.transformation = ~ A + B,
                                  tables = "predictions", 
                                  error.intervals = "Stand")
  testthat::expect_equal(nrow(diffs.sub$predictions), 5)
  #Fit additive mode
  mod.add <- asreml(y ~ A+B, data = twoway)
  pred.add <- predict(mod.add, classify = "A:B")
  testthat::expect_true(!any(is.na(pred.add$pvals$predicted.value)))
  #Compare projected and fitted additive predictions
  testthat::expect_true(abs(diffs.sub$predictions$predicted.value[1] - 
                              pred.add$pvals$predicted.value[1]) > 0.01)
  
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
  asreml.options(keep.order = TRUE) #required for asreml4 only
  for (k in 1:nresp)
  {
    fix <- paste(responses.lRGR[k], 
                 " ~ Smarthouse/(xMainPosn + xLane) + Genotype.ID*Treatment.1",
                 sep= "")
    HEB25.asr <- do.call("asreml",
                         args = list(fixed = as.formula(fix), 
                                     random = ~ Smarthouse:Zones:Mainplots, 
                                     residual = ~ idh(Treat.Smarthouse):Zones:Mainplots, 
                                     data = cart.dat, workspace="500mb", 
                                     na.action=na.method(y="include", x="include"),
                                     maxiter=50))
    summary(HEB25.asr)$varcomp
    current.asrt <- as.asrtests(HEB25.asr)
    current.asrt <- rmboundary(current.asrt)
    current.asr <- current.asrt$asreml.obj
    wald.tab <- recalcWaldTab(asrtests.obj=current.asrt, dDF.na="residual")
    HEB25.diffs <- predictPlus(current.asr, classify="Treatment.1:Genotype.ID",
                               wald.tab = wald.tab, Vmatrix = TRUE, 
                               transform.power = 0,
                               tables= "none", pworkspace=32E+06)
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

cat("#### Test for addBacktransforms on WaterRunoff with asreml4\n")
test_that("addBacktransforms_WaterRunoff_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)

  ##Use asreml to get predictions and associated statistics
  
  asreml.options(keep.order = TRUE) #required for asreml-R4 only
  current.asr <- asreml(fixed = log.Turbidity ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        keep.order=TRUE, data= WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  TS.diffs <- predictPlus(classify = "Sources:Type", 
                          asreml.obj = current.asr, 
                          wald.tab = current.asrt$wald.tab, 
                          present = c("Sources", "Type", "Species"))
  
  ## Plot p-values for predictions obtained using asreml4
  if (exists("TS.diffs"))
  {
    ##Add the backtransforms component for predictions obtained using asreml or lmerTest  
    TS.diffs <- addBacktransforms.alldiffs(TS.diffs, transform.power = 0)
    testthat::expect_false(is.null(TS.diffs$backtransforms))
    testthat::expect_true(all(abs(exp(TS.diffs$predictions$predicted.value)-
                                    TS.diffs$backtransforms$backtransformed.predictions) < 1e-06))
    testthat::expect_true(all(abs(exp(TS.diffs$predictions$upper.Confidence.limit)-
                                    TS.diffs$backtransforms$upper.Confidence.limit) < 1e-06))
    
    TS.diffs.LSD <- redoErrorIntervals(TS.diffs, error.intervals = "half", 
                                       LSDtype = "factor", LSDby = "Type")
    testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic") %in% names(attributes(TS.diffs.LSD))))
    testthat::expect_true(all(attr(TS.diffs.LSD, which = "LSDtype") == "factor.combinations"))
    testthat::expect_true(all(attr(TS.diffs.LSD, which = "LSDstatistic") == "mean"))
    testthat::expect_true(all(c("LSDtype", "LSDstatistic", "LSDvalues") %in% 
                                names(attributes(TS.diffs.LSD$predictions))))
    testthat::expect_true(all(c("LSDtype", "LSDstatistic") %in% 
                                names(attributes(TS.diffs.LSD$backtransforms))))
    testthat::expect_true(is.null(attr(TS.diffs.LSD$backtransforms, which = "LSDvalues")))
  }
})


cat("#### Test for ratioTansforms on system data with asreml4\n")
test_that("ratioTransforms_SystemData_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asremlPlus)
  library(dae)
  load(system.file("extdata", "testDiffs.rda", package = "asremlPlus", mustWork = TRUE))
  
  Preds.ratio.RGR <- ratioTransform.alldiffs(alldiffs.obj = diffs.RGR,
                                             ratio.factor = "Salinity", 
                                             numerator.levels = "Salt",
                                             denominator.levels = "Control")
  testthat::expect_true(all(abs(Preds.ratio.RGR$`Salt,Control`$predicted.value[1:3] - 
                                  c(0.7330241,0.9098487,0.9513331)) < 1e-05))
  testthat::expect_true(all(Preds.ratio.RGR$`Salt,Control`$Temperature[1:3] == "Cool"))
  testthat::expect_true(all(Preds.ratio.RGR$`Salt,Control`$Genotype[1:3] == as.character(1:3)))
  testthat::expect_true(all(names(Preds.ratio.RGR$`Salt,Control`)[5:6] == c("upper.Confidence.limit",
                                                                            "lower.Confidence.limit")))

  #Because diffs.ClUp was built using asremplus v4.2-xx, it needs to be renewed for the new attributes
  #For testing see 
  diffs.new <- redoErrorIntervals(diffs.ClUp, error.intervals = "half", 
                                  LSDtype = "factor", LSDby = c("Temperature", "Genotype"))

  Preds.ratio.ClUp <- pairdiffsTransform(diffs.new, method = "log",
                                         pairs.factor = "Temperature", 
                                         first.levels = "Hot",
                                         second.levels = "Cool",
                                         error.intervals = "halfLeast",
                                         LSDtype = "factor", LSDby = "Genotype",
                                         tables = "backtrans")
  testthat::expect_true(all(abs(Preds.ratio.ClUp$`Hot,Cool`$predictions$predicted.value[1:3] - 
                                  c(-0.05752483,0.12987766,0.06916038)) < 1e-05))
  testthat::expect_true(all(Preds.ratio.ClUp$`Hot,Cool`$predictions$Salinity[1:3] == "Control"))
  testthat::expect_true(all(Preds.ratio.ClUp$`Hot,Cool`$predictions$Genotype[1:3] == as.character(1:3)))
  testthat::expect_true(all(names(Preds.ratio.ClUp$`Hot,Cool`$predictions)[5:6] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(c("tdf", "alpha", "LSDtype", "LSDstatistic") %in% 
                              names(attributes(Preds.ratio.ClUp$`Hot,Cool`))))
  testthat::expect_true(all(attr(Preds.ratio.ClUp$`Hot,Cool`, which = "LSDtype") == "factor.combinations"))
  testthat::expect_true(all(attr(Preds.ratio.ClUp$`Hot,Cool`, which = "LSDstatistic") == "mean"))
  testthat::expect_true(all(c("LSDtype", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(Preds.ratio.ClUp$`Hot,Cool`$predictions))))
  testthat::expect_true(all(c("LSDtype", "LSDstatistic") %in% 
                              names(attributes(Preds.ratio.ClUp$`Hot,Cool`$backtransforms))))
  testthat::expect_true(is.null(attr(Preds.ratio.ClUp$`Hot,Cool`$backtransforms, which = "LSDvalues")))
})

cat("#### Test for ratioTansforms on the Oats data with asreml4\n")
test_that("ratioTransforms_SystemData_asreml4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  data("Oats.dat")
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  current.asrt <- as.asrtests(m1.asr)
  wald.tab <-  current.asrt$wald.tab
  Var.diffs <- predictPlus(m1.asr, classify="Nitrogen:Variety", pairwise = TRUE,
                          Vmatrix = TRUE, error.intervals = "halfLeast",
                          LSDtype = "factor", LSDby = "Variety",
                          wald.tab = wald.tab)

  #Test ratioTransform
  Preds.ratio.OatsN <- ratioTransform(alldiffs.obj = Var.diffs,
                                      ratio.factor = "Nitrogen", 
                                      numerator.levels = "0.6",
                                      denominator.levels = "0.2")
  testthat::expect_true(names(Preds.ratio.OatsN) == "0.6,0.2")
  testthat::expect_true(all(abs(Preds.ratio.OatsN$`0.6,0.2`$predicted.value - 
                                  c(1.321561,1.267343,1.168971)) < 1e-05))
  testthat::expect_true(all(Preds.ratio.OatsN$`0.6,0.2`$Variety == c("Victory", "Golden Rain","Marvellous")))
  testthat::expect_true(all(names(Preds.ratio.OatsN$`0.2,0`)[4:5] == c("upper.Confidence.limit",
                                                                       "lower.Confidence.limit")))
  
  #Test for ordering of the ratioTransforms
  diffs.sort <- sort(Var.diffs, sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_equal(as.character(diffs.sort$predictions$Variety[1]),"Marvellous")
  testthat::expect_equal(rownames(diffs.sort$differences)[1],"0,Marvellous")
  testthat::expect_equal(colnames(diffs.sort$p.differences)[1],"0,Marvellous")
  testthat::expect_silent(plotPvalues(diffs.sort, gridspacing = 3))
  
  Preds.ratio.sort.OatsN <- ratioTransform(alldiffs.obj = diffs.sort, 
                                           ratio.factor = "Nitrogen", 
                                           numerator.levels = "0.6",
                                           denominator.levels = "0.2")
  testthat::expect_true(names(Preds.ratio.sort.OatsN) == "0.6,0.2")
  testthat::expect_true(all(abs(Preds.ratio.sort.OatsN$`0.6,0.2`$predicted.value - 
                                  c(1.168971,1.267343,1.321561)) < 1e-05))
  testthat::expect_true(all(Preds.ratio.sort.OatsN$`0.6,0.2`$Variety == c("Marvellous", "Golden Rain","Victory")))
  testthat::expect_true(all(names(Preds.ratio.sort.OatsN$`0.6,0.2`)[4:5] == c("upper.Confidence.limit",
                                                                       "lower.Confidence.limit")))

  #Test pairdiffsTransform
  Preds.diffs.OatsN <- pairdiffsTransform(alldiffs.obj = Var.diffs,
                                          pairs.factor = "Nitrogen", 
                                          first.levels = "0.6",
                                          second.levels = "0.2", error.intervals = "halfLeast",
                                          tables = "none")
  testthat::expect_true(names(Preds.diffs.OatsN) == "0.6,0.2")
  
  testthat::expect_true(all(abs(Preds.diffs.OatsN$`0.6,0.2`$predictions$predicted.value - 
                                  (Var.diffs$predictions$predicted.value[10:12] - 
                                     Var.diffs$predictions$predicted.value[4:6])) < 1e-05))
  testthat::expect_true(all(Preds.diffs.OatsN$`0.6,0.2`$predictions$Variety == c("Victory", "Golden Rain","Marvellous")))
  testthat::expect_true(all(names(Preds.diffs.OatsN$`0.6,0.2`$predictions)[4:5] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))

  #Remove the Victory, 0.2 combination to test what happens when not all numerator combinations are present
  Var.red.diffs <- subset(Var.diffs, subset = !(Variety == "Victory" & Nitrogen == "0.2"))
  testthat::expect_equal(nrow(Var.red.diffs$predictions),  11)

  Preds.red.ratio.OatsN <- ratioTransform(alldiffs.obj = Var.red.diffs,
                                          ratio.factor = "Nitrogen", 
                                          numerator.levels = c("0.2","0.4","0.6"),
                                          denominator.levels = "0")
  testthat::expect_true(all(names(Preds.red.ratio.OatsN) == c("0.2,0", "0.4,0", "0.6,0")))
  testthat::expect_true(is.na(Preds.red.ratio.OatsN$`0.2,0`$predicted.value[1])) 
  testthat::expect_true(all(abs(Preds.red.ratio.OatsN$`0.2,0`$predicted.value[2:3] - 
                                  c(1.231250,1.251923)) < 1e-05))
  testthat::expect_true(all(Preds.red.ratio.OatsN$`0.2,0`$Variety == c("Victory", "Golden Rain","Marvellous")))
  testthat::expect_true(all(names(Preds.red.ratio.OatsN$`0.2,0`)[4:5] == c("upper.Confidence.limit",
                                                                       "lower.Confidence.limit")))
  
  Preds.red.diffs.OatsN <- pairdiffsTransform(alldiffs.obj = Var.red.diffs,
                                              pairs.factor = "Nitrogen", 
                                              first.levels = c("0.2","0.4","0.6"),
                                              second.levels = "0", error.intervals = "halfLeast",
                                              tables = "none")
  testthat::expect_true(all(names(Preds.red.diffs.OatsN) == c("0.2,0", "0.4,0", "0.6,0")))
  testthat::expect_true(all(abs(Preds.red.diffs.OatsN$`0.2,0`$predictions$predicted.value - 
                                  c(18.50000, 21.83333)) < 1e-05))
  testthat::expect_true(all(Preds.red.diffs.OatsN$`0.2,0`$predictions$Variety == c("Golden Rain","Marvellous")))
  testthat::expect_true(all(names(Preds.red.diffs.OatsN$`0.2,0`$predictions)[4:5] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  testthat::expect_true(all(abs(Preds.red.diffs.OatsN$`0.4,0`$predictions$predicted.value - 
                                  c(39.33333,34.66667,30.50000)) < 1e-05))
  testthat::expect_true(all(Preds.red.diffs.OatsN$`0.4,0`$predictions$Variety == c("Victory","Golden Rain","Marvellous")))
  
  #Remove the Victory, 0.4 combination to test what happens when not all denominator combinations are present
  Var.red.diffs <- subset(Var.diffs, subset = !(Variety == "Victory" & Nitrogen == "0"))
  testthat::expect_equal(nrow(Var.red.diffs$predictions),  11)
  
  Preds.red.ratio.OatsN <- ratioTransform(alldiffs.obj = Var.red.diffs,
                                          ratio.factor = "Nitrogen", 
                                          numerator.levels = c("0.2","0.4","0.6"),
                                          denominator.levels = "0")
  testthat::expect_true(all(names(Preds.red.ratio.OatsN) == c("0.2,0", "0.4,0", "0.6,0")))
  testthat::expect_true(all(unlist(lapply(Preds.red.ratio.OatsN, function(pd) is.na(pd$predicted.value[1])))))
  testthat::expect_true(all(abs(Preds.red.ratio.OatsN$`0.2,0`$predicted.value[2:3] - 
                                  (Var.red.diffs$predictions$predicted.value[4:5]/
                                     Var.red.diffs$predictions$predicted.value[1:2])) < 1e-05))
  testthat::expect_true(all(Preds.red.ratio.OatsN$`0.2,0`$Variety == c("Victory", "Golden Rain","Marvellous")))
  testthat::expect_true(all(names(Preds.red.ratio.OatsN$`0.2,0`)[4:5] == c("upper.Confidence.limit",
                                                                       "lower.Confidence.limit")))
  
  Preds.red.diffs.OatsN <- pairdiffsTransform(alldiffs.obj = Var.red.diffs,
                                              pairs.factor = "Nitrogen", 
                                              first.levels = c("0.2","0.4","0.6"),
                                              second.levels = "0", error.intervals = "halfLeast",
                                              tables = "none")
  testthat::expect_true(all(names(Preds.red.diffs.OatsN) == c("0.2,0", "0.4,0", "0.6,0")))
  testthat::expect_true(all(abs(Preds.red.diffs.OatsN$`0.2,0`$predictions$predicted.value - 
                                  c(18.50000, 21.83333)) < 1e-05))
  testthat::expect_true(all(Preds.red.diffs.OatsN$`0.2,0`$predictions$Variety == c("Golden Rain","Marvellous")))
  testthat::expect_true(all(names(Preds.red.diffs.OatsN$`0.2,0`$predictions)[4:5] == 
                              c("upper.halfLeastSignificant.limit", "lower.halfLeastSignificant.limit")))
  
  #Test ratioTransform when a a level occurs in both numerator.levels and denominator.levels
  Preds.ratio.OatsN <- ratioTransform(alldiffs.obj = Var.diffs,
                                      ratio.factor = "Nitrogen", 
                                      numerator.levels = c("0.2","0.6"),
                                      denominator.levels = c("0.2","0.6"))
  testthat::expect_true(all(names(Preds.ratio.OatsN) == c("0.2,0.2","0.6,0.2","0.2,0.6","0.6,0.6")))
  
  testthat::expect_true(all(unlist(lapply(Preds.ratio.OatsN, is.null)) == c(TRUE,FALSE,FALSE,TRUE)))
  testthat::expect_true(all(abs(Preds.ratio.OatsN$`0.6,0.2`$predicted.value - 
                                  c(1.321561,1.267343,1.168971)) < 1e-05))
  testthat::expect_true(all(Preds.ratio.OatsN$`0.6,0.2`$Variety == c("Victory", "Golden Rain","Marvellous")))
  testthat::expect_true(all(names(Preds.ratio.OatsN$`0.2,0`)[4:5] == c("upper.Confidence.limit",
                                                                       "lower.Confidence.limit")))
  
})

