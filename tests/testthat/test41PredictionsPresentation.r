#devtools::test("asremlPlus")
context("prediction_presentation")
asr41.lib <- "D:\\Analyses\\R ASReml4.1" 

cat("#### Test for Intercept prediction on Oats with asreml41\n")
test_that("predict_Intercept4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  testthat::expect_equal(length(m1.asr$vparameters),3)
  current.asrt <- as.asrtests(m1.asr)
  
  #Test for Intercept predict
  Int.pred <- predict(m1.asr, classify="(Intercept)")$pvals
  testthat::expect_equal(nrow(Int.pred), 1)
  testthat::expect_true(abs( Int.pred$predicted.value - 103.9722) < 1e-04)
  Int.diffs <- predictPlus(m1.asr, classify="(Intercept)")
  testthat::expect_equal(length(Int.diffs),7)
  testthat::expect_equal(nrow(Int.diffs$predictions), 1)
  testthat::expect_true(abs( Int.diffs$predictions$predicted.value - 103.9722) < 1e-04)
  
  xtitl <- "Overall mean"
  names(xtitl) <- "Intercept"
  testthat::expect_silent(plotPredictions(classify="(Intercept)", y = "predicted.value", 
                  data = Int.diffs$predictions, 
                  y.title = "Yield", titles = xtitl,
                  error.intervals = "Conf"))
})


cat("#### Test for predictPlus.asreml41\n")
test_that("predictPlus.asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)
  asreml.options(keep.order = TRUE) #required for asreml4 only
  testthat::expect_warning(current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                                                 random = ~ Benches:MainPlots,
                                                 data= WaterRunoff.dat))
  testthat::expect_output(current.asrt <- as.asrtests(current.asr, NULL, NULL),
                          regexp = "Calculating denominator DF")
  testthat::expect_silent(diffs <- predictPlus(classify = "Sources:Type", 
                                               asreml.obj = current.asr, tables = "none", 
                                               wald.tab = current.asrt$wald.tab, 
                                               present = c("Type","Species","Sources")))
  testthat::expect_is(diffs, "alldiffs")
  
  #### Get the observed combinations of the factors and variables in classify
  class.facs <- c("Species","Date","xDay")
  levs <- as.data.frame(table(WaterRunoff.dat[class.facs]))
  levs <- levs[do.call(order, levs), ]
  levs <- as.list(levs[levs$Freq != 0, class.facs])
  levs$xDay <- as.numfac(levs$xDay)
  
  current.asr <- asreml(fixed = log.Turbidity ~ Benches + Sources + Type + Species +
                          Sources:Type + Sources:Species + 
                          Sources:xDay + Species:xDay + Species:Date,
                        data = WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  
  diffs.p <- predictPlus(asreml.obj = current.asr, 
                         classify="Species:Date:xDay", 
                         term = "Species:Date", 
                         parallel = TRUE, levels=levs, 
                         present=c("Type","Species","Sources"), 
                         x.num = "xDay", x.fac = "Date", 
                         x.plot.values=c(0,28,56,84), tables = "none",
                         wald.tab = current.asrt$wald.tab)
  
  testthat::expect_is(diffs.p, "alldiffs")
  
})

cat("#### Test for plotPredictions.asreml41\n")
test_that("plotPredictions.asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(ggplot2)
  library(dae)
  data(WaterRunoff.dat)
  #### Get the observed combinations of the factors and variables in classify
  class.facs <- c("Species","Date","xDay")
  levs <- as.data.frame(table(WaterRunoff.dat[class.facs]))
  levs <- levs[do.call(order, levs), ]
  levs <- as.list(levs[levs$Freq != 0, class.facs])
  levs$xDay <- as.numfac(levs$xDay)

  asreml.options(keep.order = TRUE) #required for asreml4 only
  current.asr <- asreml(fixed = log.Turbidity ~ Benches + Sources + Type + Species +
                          Sources:Type + Sources:Species + 
                          Sources:xDay + Species:xDay + Species:Date,
                        data = WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  predictions <- predict(current.asr, class="Species:Date:xDay", 
                         parallel = TRUE, levels = levs, 
                         present = c("Type","Species","Sources"))$pvals
  predictions <- predictions[predictions$status == "Estimable",]
  
  x.title <- "Days since first observation"
  names(x.title) <- "xDay"
  #Get predictions without specifying levels
  plotPredictions(classify="Species:Date:xDay", y = "predicted.value", 
                  data = predictions, wald.tab = current.asrt$wald.tab, 
                  x.num = "xDay", x.fac = "Date", 
                  titles = x.title,
                  y.title = "Predicted log(Turbidity)",
                  present = c("Type","Species","Sources"),
                  error.intervals = "none", 
                  ggplotFuncs = list(ggtitle("Transformed turbidity over time")))
  
  #Specify the levs and parallel = TRUE
  diffs <- predictPlus(asreml.obj = current.asr, 
                       classify="Species:Date:xDay", 
                       term = "Species:Date", 
                       present=c("Type","Species","Sources"), 
                       x.num = "xDay", x.fac = "Date", 
                       parallel = TRUE, levels = levs,
                       x.plot.values=c(0,28,56,84),
                       wald.tab = current.asrt$wald.tab)
  plotPredictions(classify="Species:Date:xDay", y = "predicted.value", 
                  data = diffs$predictions, wald.tab = current.asrt$wald.tab, 
                  x.num = "xDay", x.fac = "Date", 
                  titles = x.title,
                  y.title = "Predicted log(Turbidity)")
  testthat::expect_silent("dummy")
})

cat("#### Test for predictPresent.asreml41\n")
test_that("predictPresent.asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  data(WaterRunoff.dat)
  #### Get the observed combinations of the factors and variables in classify
  class.facs <- c("Species","Date","xDay")
  levs <- as.data.frame(table(WaterRunoff.dat[class.facs]))
  levs <- levs[do.call(order, levs), ]
  levs <- as.list(levs[levs$Freq != 0, class.facs])
  levs$xDay <- as.numfac(levs$xDay)

  titles <- list("Days since first observation", "Days since first observation", "pH", "Turbidity (NTU)")
  names(titles) <- names(WaterRunoff.dat)[c(5,7,11:12)]
  asreml.options(keep.order = TRUE) #required for asreml4 only
  current.asr <- asreml(fixed = log.Turbidity ~ Benches + Sources + Type + Species + 
                          Sources:Type + Sources:Species + Sources:Species:xDay + 
                          Sources:Species:Date, 
                        data = WaterRunoff.dat)
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  #Example that fails because Date has levels that are not numeric in nature
  testthat::expect_error(diff.list <- predictPresent(terms = "Date:Sources:Species", 
                                                     asreml.obj = current.asrt$asreml.obj, 
                                                     wald.tab = current.asrt$wald.tab, 
                                                     x.fac = "Date", 
                                                     plots = "predictions", 
                                                     error.intervals = "StandardError", 
                                                     titles = titles, 
                                                     transform.power = 0, 
                                                     present = c("Type","Species","Sources"), 
                                                     tables = "differences", 
                                                     level.length = 6))
  #Example that does not produce predictions because has Date but not xDay
  testthat::expect_error(diff.list <- predictPresent(terms = "Date:Sources:Species", 
                                                     asreml.obj = current.asrt$asreml.obj, 
                                                     wald.tab = current.asrt$wald.tab, 
                                                     plots = "predictions", 
                                                     error.intervals = "StandardError", 
                                                     titles = titles, 
                                                     transform.power = 0, 
                                                     present = c("Type","Species","Sources","Date"), 
                                                     tables = "differences", 
                                                     level.length = 6))
  
  #### Get the observed combinations of the factors and variables in classify
  class.facs <- c("Sources","Species","Date","xDay")
  levs <- as.data.frame(table(WaterRunoff.dat[class.facs]))
  levs <- levs[do.call(order, levs), ]
  levs <- as.list(levs[levs$Freq != 0, class.facs])
  levs$xDay <- as.numfac(levs$xDay)
  
  # parallel and levels are arguments from predict.asreml
  diff.list <- predictPresent.asreml(asreml.obj = current.asrt$asreml.obj, 
                                     terms = "Date:Sources:Species:xDay",
                                     x.num = "xDay", x.fac = "Date", 
                                     parallel = TRUE, levels = levs, 
                                     wald.tab = current.asrt$wald.tab, 
                                     plots = "predictions", 
                                     error.intervals = "StandardError", 
                                     titles = titles, 
                                     transform.power = 0, 
                                     present = c("Type","Species","Sources"), 
                                     tables = "none", 
                                     level.length = 6)
  testthat::expect_equal(length(diff.list), 1)
  testthat::expect_match(names(diff.list), "Date.Sources.Species.xDay")

  # test that backtransforms have halfLSD intervals
  diff.list <- predictPresent.asreml(asreml.obj = current.asrt$asreml.obj, 
                                     terms = "Date:Sources:Species:xDay",
                                     x.num = "xDay", x.fac = "Date", 
                                     parallel = TRUE, levels = levs, 
                                     wald.tab = current.asrt$wald.tab, 
                                     plots = "backtransforms", 
                                     error.intervals = "halfLeast", 
                                     avsed.tolerance = 1,
                                     titles = titles, 
                                     transform.power = 0, 
                                     present = c("Type","Species","Sources"), 
                                     tables = "none", 
                                     level.length = 6)
  testthat::expect_equal(length(diff.list), 1)
  testthat::expect_match(names(diff.list), "Date.Sources.Species.xDay")
  testthat::expect_true(all(c("upper.halfLeastSignificant.limit", 
                              "lower.halfLeastSignificant.limit") %in% 
                              names(diff.list$Date.Sources.Species.xDay$backtransforms)))
  
  })


#### This test is not relevant to asreml3 because its saving of sed and vcov are different
cat("#### Test for error when no predictions.asreml41\n")
test_that("noPredictions.asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  data(gw.dat)
  current.asr <- do.call(asreml, 
                         args=list(fixed = y ~ Species*Substrate*Irrigation,
                                   random = ~ Row + Column,
                                   keep.order=TRUE, data = gw.dat, 
                                   maxit=50, workspace = 1e08, stepsize = 0.0001))
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  current.asrt <- rmboundary(current.asrt)
  testthat::expect_error(diffs <- predictPresent(current.asrt$asreml.obj,
                                                 terms = "Irrigation",
                                                 error.intervals = "Conf", 
                                                 wald.tab = current.asrt$wald.tab,
                                                 tables = "none")[[1]], 
                         regexp = "predict.asreml has not returned the sed component for the predictions as requested",
                         fixed = TRUE)
  
  testthat::expect_error(diffs <- predictPresent(current.asrt$asreml.obj,
                                                 terms = "Irrigation",
                                                 linear.transformation = ~ Irrigation,
                                                 error.intervals = "Conf", 
                                                 wald.tab = current.asrt$wald.tab,
                                                 tables = "none")[[1]], 
                         regexp = "predict.asreml has not returned the variance matrix of the predictions as requested",
                         fixed = TRUE)
})

cat("#### Test for plotPvalues.asreml41\n")
test_that("plotPvalues.asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(dae)
  library(reshape2)
  data(WaterRunoff.dat)
  asreml.options(keep.order = TRUE) #required for asreml4 only
  testthat::expect_output(current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                                                random = ~ Benches:MainPlots,
                                                data= WaterRunoff.dat))
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  diffs <- predictPlus.asreml(classify = "Sources:Type", 
                              asreml.obj = current.asr, tables = "none", 
                              wald.tab = current.asrt$wald.tab, 
                              present = c("Type","Species","Sources"))
  testthat::expect_is(diffs, "alldiffs")
  
  p <- diffs$p.differences
  p <- within(reshape2::melt(p), 
              { 
                Var1 <- factor(Var1, levels=dimnames(diffs$p.differences)[[1]])
                Var2 <- factor(Var2, levels=levels(Var1))
              })
  names(p) <- c("Rows","Columns","p")
  testthat::expect_silent(plotPvalues(p, x = "Rows", y = "Columns", 
                                      gridspacing = rep(c(3,4), c(4,2)), 
                                      show.sig = TRUE))

  #Test different size, face and colour
  testthat::expect_silent(plotPvalues(p, x = "Rows", y = "Columns", 
                                      gridspacing = rep(c(3,4), c(4,2)), 
                                      show.sig = TRUE, sig.size = 5, sig.colour = "blue"))
  testthat::expect_silent(plotPvalues(p, x = "Rows", y = "Columns", 
                                      gridspacing = rep(c(3,4), c(4,2)), 
                                      show.sig = TRUE, sig.size = 5, sig.face = "bold", 
                                      sig.family = "serif"))
  
  #Plot with sections
  pdata <- plotPvalues(diffs, sections = "Sources", show.sig = TRUE)
  testthat::expect_equal(nrow(pdata$pvalues), 400)
  testthat::expect_equal(ncol(pdata$pvalues), 5)
  testthat::expect_true(all(c("Rows","Columns","p","sections1","sections2") %in% names(pdata$pvalues)))
  testthat::expect_equal(length(pdata$plots), 6)
  testthat::expect_equal(names(pdata$plots), c("Rainwater","Recycled water","Tap water",
                                               "Rain+Basalt","Rain+Dolomite","Rain+Quartzite"))
  
  #Plot without sections, but automatic gridspacing
  pupdata <- plotPvalues(diffs, show.sig = TRUE, factors.per.grid = 1)
  testthat::expect_equal(nrow(pupdata$pvalues), 400)
  testthat::expect_equal(ncol(pupdata$pvalues), 3)
  testthat::expect_true(all(c("Rows","Columns","p") %in% names(pupdata$pvalues)))
  testthat::expect_equal(sum(!is.na(pupdata$pvalues$p)), 380)
  testthat::expect_equal(length(pupdata$plots), 1)
  
  #Plot without sections, but automatic gridspacing and upper triangle
  pupdata <- plotPvalues(diffs, show.sig = TRUE, factors.per.grid = 1, 
                         triangles = "upper")
  testthat::expect_equal(nrow(pupdata$pvalues), 400)
  testthat::expect_equal(ncol(pupdata$pvalues), 3)
  testthat::expect_true(all(c("Rows","Columns","p") %in% names(pupdata$pvalues)))
  testthat::expect_equal(sum(!is.na(pupdata$pvalues$p)), 190)
  
  
  #Plot without sections, but manual gridspacing and upper triangle
  pupdata <- plotPvalues(diffs, show.sig = TRUE, gridspacing = rep(c(3,4), c(4,2)), 
                         triangles = "upper")
  testthat::expect_equal(nrow(pupdata$pvalues), 400)
  testthat::expect_equal(ncol(pupdata$pvalues), 3)
  testthat::expect_true(all(c("Rows","Columns","p") %in% names(pupdata$pvalues)))
  testthat::expect_equal(sum(!is.na(pupdata$pvalues$p)), 190)
  
  #Plot without sections, but manual gridspacing and lower triangle
  pupdata <- plotPvalues(diffs, sections = "Sources", show.sig = TRUE, triangles = "upper")
  pupdata$pvalues <- na.omit(pupdata$pvalues)
  testthat::expect_equal(nrow(pupdata$pvalues), 190)
  testthat::expect_equal(ncol(pupdata$pvalues), 5)
  testthat::expect_true(all(c("Rows","Columns","p","sections1","sections2") %in% 
                                    names(pupdata$pvalues)))
})

cat("#### Test for plotPvalues.asreml41\n")
test_that("plotPvalues.asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(dae)
  LeafSucculence.diff <- readRDS("./data/LeafSucculence.diff")
  LeafSucculence.diff <- LeafSucculence.diff[[1]]
  
  pdata <- plotPvalues(LeafSucculence.diff, gridspacing = 3, show.sig = TRUE, 
                       axis.labels = TRUE)
  testthat::expect_equal(nrow(pdata$pvalue), 144)
  testthat::expect_equal(ncol(pdata$pvalues), 3)
  testthat::expect_true(all(c("Rows","Columns","p") %in% names(pdata$pvalues)))
  
  pdata <- plotPvalues(LeafSucculence.diff, factors.per.grid = 2, show.sig = TRUE, 
                       axis.labels = TRUE)
  testthat::expect_equal(nrow(pdata$pvalues), 144)
  testthat::expect_equal(ncol(pdata$pvalues), 3)
  testthat::expect_true(all(c("Rows","Columns","p") %in% names(pdata$pvalues)))
  
  pdata <- plotPvalues(LeafSucculence.diff, sections = c("Depths","Slope"), 
                       show.sig = TRUE, axis.labels = TRUE)
  testthat::expect_equal(nrow(pdata$pvalues), 144)
  testthat::expect_equal(ncol(pdata$pvalues), 5)
  testthat::expect_true(all(c("Rows","Columns","p","sections1","sections2") %in% names(pdata$pvalues)))
  
})

cat("#### Test for factor combinations asreml41\n")
test_that("factor.combinations.asreml41", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(dae)
  LeafSucculence.diff <- readRDS("./data/LeafSucculence.diff")
  LeafSucculence.diff <- LeafSucculence.diff[[1]]
  
  LeafSucculence.diff <- recalcLSD(LeafSucculence.diff, LSDtype = "factor.combinations", 
                                   LSDby = "Species")
  testthat::expect_warning(LeafSucculence.diff <- redoErrorIntervals(LeafSucculence.diff, 
                                                                    error.intervals = "half"))
  testthat::expect_equal(nrow(LeafSucculence.diff$LSD), 3)
  testthat::expect_equal(ncol(LeafSucculence.diff$LSD), 8)
  testthat::expect_true(all(c("P1","P2","P3") %in% rownames(LeafSucculence.diff$LSD)))
  testthat::expect_false("lower.halfLeastSignificant.limit" %in% names(LeafSucculence.diff$predictions))
  testthat::expect_true(names(LeafSucculence.diff$predictions)[length(names(
    LeafSucculence.diff$predictions))] == "est.status")
  
})

cat("#### Test for recalcLSD.alldiffs4\n")
test_that("recalcLSD.alldiffs4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)
  asreml.options(keep.order = TRUE) #required for asreml4 only
  testthat::expect_output(current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                                                random = ~ Benches:MainPlots,
                                                data= WaterRunoff.dat))
  current.asrt <- as.asrtests(current.asr, NULL, NULL)
  diffs <- predictPlus.asreml(classify = "Sources:Type", 
                              asreml.obj = current.asr, tables = "none", 
                              wald.tab = current.asrt$wald.tab, 
                              present = c("Type","Species","Sources"))
  testthat::expect_is(diffs, "alldiffs")
  
  diffs <- recalcLSD.alldiffs(diffs, LSDtype = "factor.combinations", LSDby = "Sources")
  testthat::expect_equal(nrow(diffs$LSD), 6)
  testthat::expect_equal(ncol(diffs$LSD), 8)
  testthat::expect_warning(diffs <- redoErrorIntervals(diffs, 
                                                       error.intervals = "halfLeastSignificant"))
  testthat::expect_false("upper.halfLeastSignificant.limit" %in% names(diffs$predictions))
  
})


cat("#### Test for LSDby4\n")
test_that("LSDby4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml, lib.loc = asr41.lib)
  library(asremlPlus)
  library(dae)
  #example 9-1 from Montgomery 5 edn

  #Set up data.frame
  Pressure.lev <- c(10,15,20)
  Speed.lev <- c(100,120,140)
  Nozzle.lev <- c("A", "B", "C")
  Fac3Syrup.dat <- fac.gen(generate=list(Nozzle = Nozzle.lev, 
                                         Pressure = Pressure.lev, Speed = Speed.lev),
                           each=2)
  Fac3Syrup.dat <- within(Fac3Syrup.dat, 
                          {
                            SpeedPress <- fac.combine(list(Speed,Pressure), 
                                                       combine.levels = TRUE)
                            WSpeedPress <- fac.nested(SpeedPress)
                          })
  Fac3Syrup.dat <- data.frame(Test = factor(1:54), Fac3Syrup.dat)
  Fac3Syrup.dat$Loss <- c(-35,-25,-45,-60,-40,15, 110,75,-10,30,80,54,
                          4,5,-40,-30,31,36, 17,24,-65,-58,20,4, 
                          55,120,-55,-44,110,44, -23,-5,-64,-62,-20,-31,
                          -39,-35,-55,-67,15,-30, 90,113,-28,-26,110,135,
                          -30,-55,-61,-52,54,4)+70
  Fac3Syrup.dat <- with(Fac3Syrup.dat, Fac3Syrup.dat[order(SpeedPress, WSpeedPress),])
  #Analysis
  interaction.ABC.plot(Loss, Pressure, Speed, Nozzle, data=Fac3Syrup.dat)
  Fac3Syrup.aov <- aov(Loss ~ Nozzle * Pressure * Speed + Error(Test), Fac3Syrup.dat)
  summary(Fac3Syrup.aov)
  
  
  m1 <- do.call("asreml",
                args = list(Loss ~ Nozzle * Pressure * Speed, 
                            residual = ~idh(SpeedPress):WSpeedPress,
                            data = Fac3Syrup.dat))
  testthat::expect_true(abs(summary(m1)$varcomp$component[2] - 27.5) < 1e-05)
  wald.tab <- wald.asreml(m1, denDF = "numeric")$Wald
  testthat::expect_equal(nrow(wald.tab), 8)
  diffs <- predictPlus(m1, classify = "Nozzle:Pressure:Speed", 
                       #linear.transformation = ~(Nozzle + Pressure):Speed,
                       wald.tab = wald.tab,
                       tables = "none")
  testthat::expect_true("upper.Confidence.limit" %in% names(diffs$predictions))
  testthat::expect_true(all(c( "LSDtype", "LSDstatistic") %in% names(attributes(diffs))))
  testthat::expect_true(is.null(attr(diffs, which = "LSDby")))
  testthat::expect_true((attr(diffs, which = "LSDtype") == "overall"))
  
    
  #Calculate LSD, but leave as CIs
  diffs.LSD <- recalcLSD(diffs, LSDtype = "factor",
                         LSDby = c("Speed","Pressure"))
  testthat::expect_equal(nrow(diffs.LSD$LSD), 9)
  testthat::expect_true(abs(diffs.LSD$LSD$minLSD[1]- 11.92550) < 1e-05)
  testthat::expect_true(all(abs(diffs.LSD$LSD$minLSD- diffs.LSD$LSD$maxLSD) < 1e-05))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD))))
  testthat::expect_true((attr(diffs.LSD, which = "LSDtype") == "factor.combinations"))
  testthat::expect_true("upper.Confidence.limit" %in% names(diffs$predictions))
  
  #Convert from CI to LSI
  diffs.LSI <- redoErrorIntervals(diffs.LSD, error.intervals = "half")
  testthat::expect_true("upper.halfLeastSignificant.limit" %in% names(diffs.LSI$predictions))
  testthat::expect_equal(nrow(diffs.LSI$LSD), 9)
  diffs <- redoErrorIntervals(diffs, error.intervals = "half", LSDtype = "factor",
                       LSDby = c("Speed","Pressure"), wald.tab = wald.tab,
                       tables = "none")
  testthat::expect_true("upper.halfLeastSignificant.limit" %in% names(diffs$predictions))
  testthat::expect_equal(nrow(diffs$LSD), 9)
  testthat::expect_true(abs(diffs$LSD$minLSD[1]- 11.92550) < 1e-05)
  testthat::expect_true(all(abs(diffs$LSD$minLSD- diffs$LSD$maxLSD) < 1e-05))
  
  #Test changing the LSDby
  testthat::expect_warning(diff.Press <- 
                             redoErrorIntervals(diffs, error.intervals = "half", 
                                                LSDtype = "factor",
                                                LSDby = "Pressure", wald.tab = wald.tab,
                                                tables = "none"))
  diff.Press$LSD
  testthat::expect_equal(nrow(diff.Press$LSD), 3)
  testthat::expect_true(abs(diff.Press$LSD$minLSD[1]- 11.92550) < 1e-05)
  testthat::expect_true(abs(diff.Press$LSD$meanLSD[1]- 41.13342) < 1e-05)
  testthat::expect_true(abs(diff.Press$LSD$maxLSD[1]- 67.62672) < 1e-05)
  
  #No LSDtype
  testthat::expect_warning(diff.Press <- 
                             redoErrorIntervals(diffs, error.intervals = "half", 
                                                LSDby = "Pressure", wald.tab = wald.tab,
                                                tables = "none"))
  testthat::expect_equal(nrow(diff.Press$LSD), 3)
  testthat::expect_true(abs(diff.Press$LSD$minLSD[1]- 11.92550) < 1e-05)
  testthat::expect_true(abs(diff.Press$LSD$meanLSD[1]- 41.13342) < 1e-05)
  testthat::expect_true(abs(diff.Press$LSD$maxLSD[1]- 67.62672) < 1e-05)
  
  testthat::expect_warning(diff.all <- 
                             redoErrorIntervals(diffs, error.intervals = "half", 
                                                LSDtype = "overall",
                                                LSDby = NULL, wald.tab = wald.tab,
                                                tables = "none"))
  testthat::expect_equal(nrow(diff.all$LSD), 1)
  testthat::expect_true(rownames(diff.all$LSD) == "overall")
  testthat::expect_true(abs(diff.all$LSD$minLSD[1]- 11.92550) < 1e-05)

  #LSDtype = overall only
  testthat::expect_warning(diff.all <- 
                             redoErrorIntervals(diffs, error.intervals = "half", 
                                                LSDtype = "overall",
                                                wald.tab = wald.tab,
                                                tables = "none"))
  testthat::expect_equal(nrow(diff.all$LSD), 1)
  testthat::expect_true(abs(diff.all$LSD$minLSD[1]- 11.92550) < 1e-05)

  #Test predictPlus with LSD options
  #With linear transformation and LSDtype = "factor combinations"
  diffs.LSD <- predictPlus(m1, classify = "Nozzle:Pressure:Speed", 
                       linear.transformation = ~(Nozzle + Pressure):Speed,
                       error.intervals = "half", LSDtype = "factor", LSDby = c("Speed", "Pressure"),
                       wald.tab = wald.tab,
                       tables = "none")
  testthat::expect_true("upper.halfLeastSignificant.limit" %in% names(diffs.LSD$predictions))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic") %in% names(attributes(diffs.LSD))))
  testthat::expect_true((attr(diffs.LSD, which = "LSDtype") == "factor.combinations"))
  testthat::expect_true(all(c( "LSDtype", "LSDby", "LSDstatistic", "LSDvalues") %in% 
                              names(attributes(diffs.LSD$predictions))))
  testthat::expect_true(attr(diffs.LSD$predictions, which = "LSDtype") == "factor.combinations")
  testthat::expect_true(attr(diffs.LSD$predictions, which = "LSDstatistic") == "mean")

})

