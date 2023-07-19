#devtools::test("asremlPlus")
context("lme4_alldiffs")

cat("#### Test for predictions.frame on Oats with lme4\n")
test_that("PredictionsFrame_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  ## Use lmerTest and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
    Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
    Var.preds <- summary(Var.emm)
    Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
  
  ## Check the class and validity of the predictions.frame
  testthat::expect_true(is.predictions.frame(Var.preds))
  testthat::expect_true(validPredictionsFrame(Var.preds))
  }
  
})

cat("#### Test for LSD.frame on Oats with lme4\n")
test_that("PredictionsFrame_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  ## Use lme4 and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
    ## Set up a wald.tab
    int <- as.data.frame(rbind(rep(NA,4)))
    rownames(int) <- "(Intercept)"
    wald.tab <- anova(m1.lmer, ddf = "Kenward", type = 1)[,3:6]
    names(wald.tab) <- names(int) <- c("Df", "denDF", "F.inc", "Pr")
    wald.tab <- rbind(int, wald.tab)
    #Get predictions
    Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
    Var.preds <- summary(Var.emm)
    ## Modify Var.preds to be compatible with a predictions.frame
    Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
    Var.vcov <- vcov(Var.emm)
    Var.sed <- NULL
    den.df <- wald.tab[match("Nitrogen:Variety", rownames(wald.tab)), "denDF"]
    testthat::expect_equal(den.df, 45)
    
    #Test for LSDs in allDifference
    Var.LSD.diffs <- allDifferences(predictions = Var.preds, classify = "Variety:Nitrogen", 
                                     sed = Var.sed, vcov = Var.vcov, tdf = den.df,
                                     sortFactor = "Variety", decreasing = TRUE)
    testthat::expect_true(all(abs(Var.LSD.diffs$LSD -  
                                    c(66,15.47425, 18.54065, 19.56706, 18.54065, 0.1653881,2,4)) < 1e-05))
    testthat::expect_true(setequal(names(Var.LSD.diffs$LSD), c("c", "minLSD", "meanLSD", "maxLSD", 
                                                              "assignedLSD", "accuracyLSD", 
                                                              "falsePos", "falseNeg")))
    testthat::expect_equal(rownames(Var.LSD.diffs$LSD), "overall")
    
    #Test for LSDby in allDifference
    Var.LSD.diffs <- allDifferences(predictions = Var.preds, classify = "Variety:Nitrogen", 
                                    sed = Var.sed, vcov = Var.vcov, tdf = den.df, 
                                    LSDtype = "factor.combinations", LSDby = "Variety",
                                    sortFactor = "Variety", decreasing = TRUE)
    testthat::expect_true(all(abs(Var.LSD.diffs$LSD[c("minLSD", "meanLSD", "maxLSD", 
                                                      "assignedLSD")] -  15.47425) < 1e-05))
    testthat::expect_true(all(Var.LSD.diffs$LSD["accuracyLSD"] < 100 * .Machine$double.eps))
    testthat::expect_true(setequal(names(Var.LSD.diffs$LSD), c("c", "minLSD", "meanLSD", "maxLSD", 
                                                               "assignedLSD", "accuracyLSD", 
                                                               "falsePos", "falseNeg")))
    testthat::expect_true(setequal(rownames(Var.LSD.diffs$LSD), c("Marvellous", "Golden Rain", "Victory")))
    
    #test for LSDby recalcLSD
    Var.LSD.diffs <- recalcLSD(Var.LSD.diffs, LSDtype = "factor.combinations", LSDby = "Variety")
    testthat::expect_true(is.alldiffs(Var.LSD.diffs))
    testthat::expect_true(validAlldiffs(Var.LSD.diffs))
    testthat::expect_true(all(abs(Var.LSD.diffs$LSD[c("minLSD", "meanLSD", "maxLSD", 
                                                      "assignedLSD")] -  15.47425) < 1e-05))
    testthat::expect_true(all(Var.LSD.diffs$LSD["accuracyLSD"] < 100 * .Machine$double.eps))
    testthat::expect_true(setequal(names(Var.LSD.diffs$LSD), c("c", "minLSD", "meanLSD", "maxLSD", 
                                                               "assignedLSD", "accuracyLSD", 
                                                               "falsePos", "falseNeg")))
    testthat::expect_true(setequal(rownames(Var.LSD.diffs$LSD), c("Marvellous", "Golden Rain", "Victory")))
  }
})

cat("#### Test for allDifferences.data.frame on Oats with lme4\n")
test_that("alldiffs_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  ## Use lme4 and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
    ## Set up a wald.tab
    int <- as.data.frame(rbind(rep(NA,4)))
    rownames(int) <- "(Intercept)"
    wald.tab <- anova(m1.lmer, ddf = "Kenward", type = 1)[,3:6]
    names(wald.tab) <- names(int) <- c("Df", "denDF", "F.inc", "Pr")
    wald.tab <- rbind(int, wald.tab)
    #Get predictions
    Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
    Var.preds <- summary(Var.emm)
    ## Modify Var.preds to be compatible with a predictions.frame
    Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
    Var.vcov <- vcov(Var.emm)
    Var.sed <- NULL
    den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  
  #Test for as.alldiffs
  Var.diffs <- as.alldiffs(predictions = Var.preds, 
                           sed = Var.sed, vcov = Var.vcov, 
                           classify = "Nitrogen:Variety", response = "Yield", tdf = den.df)
  testthat::expect_true(is.alldiffs(Var.diffs))
  testthat::expect_true(validAlldiffs(Var.diffs))
  testthat::expect_equal(nrow(Var.diffs$predictions),12)
  testthat::expect_equal(ncol(Var.diffs$predictions),8)
  testthat::expect_true(is.null(attr(Var.diffs, which = "sortOrder")))
  testthat::expect_true(all(names(Var.diffs$predictions)[1:2] ==  c("Nitrogen","Variety")))
  testthat::expect_true(as.character(Var.diffs$predictions$Variety[1]) == "Victory" & 
                          as.character(Var.diffs$predictions$Variety[2]) == "Victory")
  
  #Test for allDifferences
  Var.sort.diffs <- allDifferences(predictions = Var.preds, classify = "Variety:Nitrogen", 
                                   sed = Var.sed, vcov = Var.vcov, tdf = den.df,
                                   sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(is.alldiffs(Var.sort.diffs))
  testthat::expect_true(validAlldiffs(Var.sort.diffs))
  testthat::expect_true(as.character(Var.sort.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.sort.diffs$predictions$Variety[2]) == "Marvellous")
  testthat::expect_equal(length(attr(Var.sort.diffs, which = "sortOrder")),3)
  testthat::expect_true(all(names(Var.sort.diffs$predictions)[1:2] ==  c("Variety","Nitrogen")))
  
  ## Change the order of the factors in the alldiffs object and reorder components
  Var.reord.diffs <- allDifferences(predictions = Var.preds,
                                    classify = "Variety:Nitrogen", 
                                    sed = Var.sed, vcov = Var.vcov, tdf = den.df)
  testthat::expect_true(is.alldiffs(Var.reord.diffs))
  testthat::expect_true(validAlldiffs(Var.reord.diffs))
  testthat::expect_true(as.character(Var.reord.diffs$predictions$Variety[1]) == "Victory" &
                          as.character(Var.reord.diffs$predictions$Variety[2]) == "Victory")
  
  #Test for re-order factors with renewClassify
  Var.reord.diffs <- renewClassify(Var.diffs, newclassify = "Variety:Nitrogen")
  testthat::expect_true(as.character(Var.reord.diffs$predictions$Variety[1]) == "Victory" &
                          as.character(Var.reord.diffs$predictions$Variety[2]) == "Victory")
  
  #Test for re-order factors and sort
  Var.both.diffs <- allDifferences(predictions = Var.preds,
                                   classify = "Variety:Nitrogen", 
                                   sed = Var.sed, vcov = Var.vcov, tdf = den.df, 
                                   sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(as.character(Var.both.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.both.diffs$predictions$Variety[2]) == "Marvellous")
  Var.both.diffs <- renewClassify(Var.diffs, newclassify = "Variety:Nitrogen", 
                                    sortFactor = "Variety", decreasing = TRUE)
  testthat::expect_true(as.character(Var.both.diffs$predictions$Variety[1]) == "Marvellous" & 
                          as.character(Var.both.diffs$predictions$Variety[2]) == "Marvellous")
  }
})  

cat("#### Test for alldiffs validity on Oats with lme4\n")
test_that("validity_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data(Oats.dat)
  
  ## Use lme4 and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
    ## Set up a wald.tab
    int <- as.data.frame(rbind(rep(NA,4)))
    rownames(int) <- "(Intercept)"
    wald.tab <- anova(m1.lmer, ddf = "Kenward", type = 1)[,3:6]
    names(wald.tab) <- names(int) <- c("Df", "denDF", "F.inc", "Pr")
    wald.tab <- rbind(int, wald.tab)
    #Get predictions
    Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
    Var.preds <- summary(Var.emm)
    ## Modify Var.preds to be compatible with a predictions.frame
    Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
    Var.vcov <- vcov(Var.emm)
    Var.sed <- NULL
    den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  
  #Test for as.alldiffs
  Var.diffs <- as.alldiffs(predictions = Var.preds, 
                           sed = Var.sed, vcov = Var.vcov, 
                           response = "Yield", tdf = den.df)
  Var.diffs$vcov <- NULL
  Var.diffs$predictions <- Var.diffs$predictions[, -match("est.status", 
                                                          names(Var.diffs$predictions))]
  testthat::expect_true(is.alldiffs(Var.diffs))
  testthat::expect_warning(msg <- validAlldiffs(Var.diffs))
  testthat::expect_equal(length(msg),5)
  testthat::expect_true(all(unlist(lapply(msg, nchar)) ==  c(41, 63, 118, 62, 119)))
  testthat::expect_true(grepl("Predictions.frame", msg[2], fixed = TRUE))
  testthat::expect_true(grepl("classify", msg[5], fixed = TRUE))
  
  #Test to show that variables in the classify must be in same order as in the predictions.frame
  Var.diffs <- as.alldiffs(predictions = Var.preds, 
                           sed = Var.sed, vcov = Var.vcov, 
                           classify = "Variety:Nitrogen", response = "Yield", tdf = den.df)
  testthat::expect_equal(length(attributes(Var.diffs)),6)
  testthat::expect_true(is.null(attr(Var.diffs, which = "sortOrder")))
  testthat::expect_true(is.character(msg <- validAlldiffs(Var.diffs))) #invalid classify
  testthat::expect_true(grepl("classify var", msg[2], fixed = TRUE))
  testthat::expect_false(all(names(Var.diffs$predictions)[1:2] ==  c("Variety","Nitrogen")))
  }
})  

cat("#### Test for linTransform on WaterRunoff with lme4\n")
test_that("linTransform_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data("WaterRunoff.dat")
  
  ## Use lmeTest and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(pH ~ Benches + (Sources * Species) + 
                                (1|Benches:MainPlots),
                              data=na.omit(WaterRunoff.dat))
    SS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Species)
    SS.preds <- summary(SS.emm)
    den.df <- min(SS.preds$df, na.rm = TRUE)
    ## Modify SS.preds to be compatible with a predictions.frame
    SS.preds <- as.predictions.frame(SS.preds, predictions = "emmean", 
                                     se = "SE", interval.type = "CI", 
                                     interval.names = c("lower.CL", "upper.CL"))
    
    ## Form an all.diffs object and check its validity
    SS.vcov <- vcov(SS.emm)
    SS.diffs <- allDifferences(predictions = SS.preds, classify = "Sources:Species", 
                               vcov = SS.vcov, tdf = den.df)
    testthat::expect_true(validAlldiffs(SS.diffs))

    #Get additive predictions
    diffs.sub <- linTransform(SS.diffs, classify = "Sources:Species", 
                              linear.transformation = ~ Sources + Species,
                              Vmatrix = TRUE, tables = "none")
    testthat::expect_true(validAlldiffs(diffs.sub))
    testthat::expect_equal(diffs.sub$predictions$predicted.value[1] - 
                             diffs.sub$predictions$predicted.value[6],
                           diffs.sub$predictions$predicted.value[7] - 
                             diffs.sub$predictions$predicted.value[12])
    
    #Contrast matrix for differences between each species and non-planted for the last source
    L <- cbind(matrix(rep(0,7*32), nrow = 7, ncol = 32),
               diag(1, nrow = 7), 
               matrix(rep(-1, 7), ncol = 1))
    rownames(L) <- as.character(diffs.sub$predictions$Species[33:39])
    diffs.L <- linTransform(diffs.sub, 
                            classify = "Sources:Species",
                            linear.transformation = L,
                            tables = "predictions")
    testthat::expect_true(abs(diffs.L$predictions$predicted.value[1] + 0.04097328) < 1e-04)
    testthat::expect_true(diffs.L$predictions$Combination[1] == "S. iqscjbogxah")
  }
})


cat("#### Test for alldiffs functions on WaterRunoff with lme4\n")
test_that("alldiffs_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data("WaterRunoff.dat")
  
  ## Use lmeTest and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(pH ~ Benches + (Sources * (Type + Species)) + 
                                (1|Benches:MainPlots),
                              data=na.omit(WaterRunoff.dat))
    TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Type)
    TS.preds <- summary(TS.emm)
    den.df <- min(TS.preds$df, na.rm = TRUE)
    ## Modify TS.preds to be compatible with a predictions.frame
    TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                     se = "SE", interval.type = "CI", 
                                     interval.names = c("lower.CL", "upper.CL"))
    
    ## Form an all.diffs object and check its validity
    els <- as.numeric(rownames(TS.preds))
    TS.vcov <- vcov(TS.emm)[els,els]
    TS.diffs <- allDifferences(predictions = TS.preds, classify = "Sources:Type", 
                               vcov = TS.vcov, tdf = den.df)
    testthat::expect_true(validAlldiffs(TS.diffs))
  }  
  
  
  ##Plot p-values for predictions obtained using asreml or lme4  
  if (exists("TS.diffs"))
  {
    ##Plot p-values based on diffs obtained from predictions using asreml or lme4  
    pdata <- plotPvalues(TS.diffs, gridspacing = rep(c(3,4), c(4,2)), show.sig = TRUE)
    testthat::expect_equal(nrow(pdata$pvalues), 400)
    testthat::expect_equal(ncol(pdata$pvalues), 3)
    testthat::expect_true(all(c("Rows","Columns","p") %in% names(pdata$pvalues)))
    testthat::expect_equal(sum(!is.na(pdata$pvalues$p)), 380)
    testthat::expect_equal(length(pdata), 2)
    testthat::expect_equal(length(pdata$plots), 1)
    
    pdata <- plotPvalues(TS.diffs, sections = "Sources", show.sig = TRUE, axis.labels = TRUE)
    testthat::expect_equal(nrow(pdata$pvalues), 400)
    testthat::expect_equal(ncol(pdata$pvalues), 5)
    testthat::expect_true(all(c("Rows","Columns","p","sections1","sections2") %in% names(pdata$pvalues)))
    testthat::expect_equal(sum(!is.na(pdata$pvalues$p)), 380)
    
    p <- TS.diffs$p.differences
    testthat::expect_silent(p <- within(reshape2::melt(p), 
                                         { 
                                           Var1 <- factor(Var1, levels=dimnames(TS.diffs$p.differences)[[1]])
                                           Var2 <- factor(Var2, levels=levels(Var1))
                                         }))
  names(p) <- c("Rows","Columns","p")
  testthat::expect_silent(plotPvalues(p, x = "Rows", y = "Columns", 
                                      gridspacing = rep(c(3,4), c(4,2)), 
                                      show.sig = TRUE))
    
    
    ##Recalculate the LSD values for predictions obtained using asreml or lme4  
    TS.diffs <- recalcLSD(TS.diffs, LSDtype = "factor.combinations", 
                          LSDby = "Sources")
    testthat::expect_equal(nrow(TS.diffs$LSD), 6)
    testthat::expect_equal(ncol(TS.diffs$LSD), 8)
    testthat::expect_warning(TS.diffs <- redoErrorIntervals(TS.diffs, 
                                                            error.intervals = "halfLeast"))
    testthat::expect_false("upper.halfLeastSignificant.limit" %in% names(TS.diffs$predictions))
    
    ##Use subset.alldiffs to select a subset of the alldiffs object
    TS.diffs.subs <- subset(TS.diffs, 
                            subset = grepl("R", Sources, fixed = TRUE) & 
                              Type %in% c("Control","Medicinal"))
    testthat::expect_is(TS.diffs.subs, "alldiffs")
    testthat::expect_true(validAlldiffs(TS.diffs.subs))
    testthat::expect_equal(nrow(TS.diffs.subs$predictions),10)
    testthat::expect_equal(ncol(TS.diffs.subs$predictions),8)
    testthat::expect_false(any(TS.diffs.subs$predictions$Sources == "Tap water"))
    testthat::expect_false(any(TS.diffs.subs$predictions$B %in% c("Landscape","Culinary")))
    testthat::expect_equal(length(attributes(TS.diffs.subs)),9)
  }
})

cat("#### Test for exploreLSD on WaterRunoff with lme4\n")
test_that("explore_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data("WaterRunoff.dat")
  
  ## Use lmeTest and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(pH ~ Benches + (Sources * (Type + Species)) + 
                                (1|Benches:MainPlots),
                              data=na.omit(WaterRunoff.dat))
    TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Type)
    TS.preds <- summary(TS.emm)
    den.df <- min(TS.preds$df, na.rm = TRUE)
    ## Modify TS.preds to be compatible with a predictions.frame
    TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                     se = "SE", interval.type = "CI", 
                                     interval.names = c("lower.CL", "upper.CL"))
    
    ## Form an all.diffs object and check its validity
    els <- as.numeric(rownames(TS.preds))
    TS.vcov <- vcov(TS.emm)[els,els]
    TS.diffs <- allDifferences(predictions = TS.preds, classify = "Sources:Type", 
                               vcov = TS.vcov, tdf = den.df)
    testthat::expect_true(validAlldiffs(TS.diffs))
  }  
  
  
  ##Explore the LSDs #the values are not the same as those from asreml
  if (exists("TS.diffs"))
  {
    ##Explore the LSD values for predictions obtained using asreml or lmerTest  
    LSDstat <- exploreLSDs(TS.diffs, LSDtype = "factor.combinations", 
                           LSDby = "Sources")
    testthat::expect_equal(names(LSDstat), c("frequencies", "distinct.vals", "statistics", "accuracy", 
                                             "false.pos", "false.neg", "per.pred.accuracy", "LSD"))
    testthat::expect_true(all(lapply(c("frequencies", "statistics", "accuracy"), 
                                     function(k, LSDstat) nrow(LSDstat[[k]]), LSDstat = LSDstat) == 6))
    testthat::expect_true(all(unlist(lapply(c("frequencies", "statistics", "accuracy"), function(k, LSDstat, dat) 
      all(rownames(LSDstat[[k]])== levels(WaterRunoff.dat$Sources)), LSDstat = LSDstat, dat = WaterRunoff.dat))))
    testthat::expect_equal(names(LSDstat$frequencies), as.character(seq(0.19, 0.41, 0.02)))
    testthat::expect_equal(LSDstat$distinct.vals$`Rain+Basalt`, c(0.204,0.306,0.319))
    testthat::expect_true(all(abs(LSDstat$statistics[1,] - 
                                    c(3, 0.2049027,0.2252329,0.2557283,0.2814212,
                                      0.3065539,0.3126769,0.3163506,0.3187998)) < 1e-05))
    testthat::expect_true(all(abs(LSDstat$accuracy[1,] - 
                                    c(3,0.5558598,0.4154228,0.2466349,0.2719005,
                                      0.3315934,0.3446824,0.3522925,0.3572686)) < 1e-05))
    testthat::expect_true(all(lapply(c("per.pred.accuracy", "LSD"), 
                                     function(k, LSDstat) nrow(LSDstat[[k]]), LSDstat = LSDstat) == 20))
    testthat::expect_equal(rownames(LSDstat$per.pred.accuracy), 
                           as.character(fac.combine(as.list(TS.diffs$predictions[c("Sources","Type")]), 
                                                    combine.levels = TRUE)))
    testthat::expect_true(all(abs(LSDstat$per.pred.accuracy[1,] - 
                                    c(0.4960954,0.361053,0.1987486,0.2719005,
                                      0.3315934,0.3446824,0.3522925,0.3572686)) < 1e-05))
  }
})

cat("#### Test for sort.alldiffs on Smarthouse with lme4\n")
test_that("sort.alldiffs_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data(Smarthouse.dat)
  
 
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    #Analyse the first response
    m1.lmer <- lmerTest::lmer(y1 ~ Genotype*A*B + (1|Replicate/Mainplot),
                              data=na.omit(Smarthouse.dat))
    GAB.emm <- emmeans::emmeans(m1.lmer, specs = ~ Genotype:A:B)
    GAB.preds <- summary(GAB.emm)
    den.df <- min(GAB.preds$df, na.rm = TRUE)
    ## Modify GAB.preds to be compatible with a predictions.frame
    GAB.preds <- as.predictions.frame(GAB.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
    
    ## Form an all.diffs object and check its validity
    GAB.vcov <- vcov(GAB.emm)
    GAB.diffs <- allDifferences(predictions = GAB.preds, classify = "Genotype:A:B", 
                                vcov = GAB.vcov, tdf = den.df)
    testthat::expect_true(is.alldiffs(GAB.diffs))
    testthat::expect_true(validAlldiffs(GAB.diffs))
    testthat::expect_equal(nrow(GAB.diffs$predictions),120)
    testthat::expect_equal(ncol(GAB.diffs$predictions),9)
    testthat::expect_equal(as.character(GAB.diffs$predictions$Genotype[1]),"Axe")
    testthat::expect_true(is.null(attr(GAB.diffs, which = "sortOrder")))
    
    #Test reodering of the classify
    GAB.diffs.reord <- renewClassify(GAB.diffs, newclassify = "A:B:Genotype")
    testthat::expect_equal(as.character(GAB.diffs.reord$predictions$Genotype[1]),"Axe")
    testthat::expect_equal(as.character(GAB.diffs.reord$predictions$Genotype[2]),"Espada")
    testthat::expect_true(abs(GAB.diffs.reord$predictions$predicted.value[2] - -0.2265723017) < 1e-06)
    
    
    #Use sort.alldiffs and save order for use with other response variables
    GAB.diffs.sort <- sort(GAB.diffs, sortFactor = "Genotype")
    sort.order <- attr(GAB.diffs.sort, which = "sortOrder")
    testthat::expect_is(GAB.diffs.sort, "alldiffs")
    testthat::expect_true(validAlldiffs(GAB.diffs.sort))
    testthat::expect_equal(nrow(GAB.diffs.sort$predictions),120)
    testthat::expect_equal(ncol(GAB.diffs.sort$predictions),9)
    testthat::expect_equal(as.character(GAB.diffs.sort$predictions$Genotype[1]),"Gladius")
    testthat::expect_equal(length(attributes(GAB.diffs.sort)),9)
    testthat::expect_equal(length(attr(GAB.diffs.sort, which = "sortOrder")),10)
    
    
    #Analyse the second response
    m2.lmer <- lmerTest::lmer(y2 ~ Genotype*A*B + (1|Replicate/Mainplot),
                              data=na.omit(Smarthouse.dat))
    GAB.emm <- emmeans::emmeans(m2.lmer, specs = ~ Genotype:A:B)
    GAB.preds <- summary(GAB.emm)
    den.df <- min(GAB.preds$df, na.rm = TRUE)
    ## Modify GAB.preds to be compatible with a predictions.frame
    GAB.preds <- as.predictions.frame(GAB.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
    
    ## Form an all.diffs object, sorting it using the y1 sort.order and check its validity
    GAB.vcov <- vcov(GAB.emm)
    GAB.diffs2.sort <- allDifferences(predictions = GAB.preds, classify = "Genotype:A:B", 
                                      vcov = GAB.vcov, tdf = den.df, 
                                      sortFactor = "Genotype", 
                                      sortOrder = sort.order)
    testthat::expect_true(validAlldiffs(GAB.diffs2.sort))
    testthat::expect_equal(as.character(GAB.diffs.sort$predictions$Genotype[1]),
                           as.character(GAB.diffs2.sort$predictions$Genotype[1]))
    testthat::expect_equal(attr(GAB.diffs.sort, which = "sortOrder"),
                           attr(GAB.diffs2.sort, which = "sortOrder"))
  } 
})

cat("#### Test for sort.alldiffs on WateRunoff with lme4\n")
test_that("sort.alldiffsWater_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)

  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    #Analyse first response  
    m1.lmer <- lmerTest::lmer(pH ~ Benches + (Sources * (Type + Species)) + 
                                (1|Benches:MainPlots),
                              data=na.omit(WaterRunoff.dat))
    TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Type)
    TS.preds <- summary(TS.emm)
    den.df <- min(TS.preds$df, na.rm = TRUE)
    ## Modify TS.preds to be compatible with a predictions.frame
    TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                     se = "SE", interval.type = "CI", 
                                     interval.names = c("lower.CL", "upper.CL"))
    
    ## Form an all.diffs object and check its validity
    els <- as.numeric(rownames(TS.preds))
    TS.vcov <- vcov(TS.emm)[els,els]
    TS.diffs <- allDifferences(predictions = TS.preds, 
                               classify = "Sources:Type", 
                               vcov = TS.vcov, tdf = den.df)
    
    testthat::expect_true(is.alldiffs(TS.diffs))
    testthat::expect_true(validAlldiffs(TS.diffs))
    testthat::expect_equal(nrow(TS.diffs$predictions),20)
    testthat::expect_equal(ncol(TS.diffs$predictions),8)
    testthat::expect_equal(as.character(TS.diffs$predictions$Type[1]),"Landscape")
    testthat::expect_true(is.null(attr(TS.diffs, which = "sortOrder")))
    
    #Test reodering of the classify
    TS.diffs.reord <- renewClassify(TS.diffs, newclassify = "Type:Sources")
    testthat::expect_equal(as.character(TS.diffs.reord$predictions$Sources[1]),"Rainwater")
    testthat::expect_equal(as.character(TS.diffs.reord$predictions$Sources[2]),"Recycled water")
    testthat::expect_true(abs(TS.diffs.reord$predictions$predicted.value[2] - 7.646389) < 1e-06)
    
    #Test sort.alldiffs and save order for use with other response variables
    TS.diffs.sort <- sort(TS.diffs, sortFactor = "Sources", 
                          sortParallelToCombo = list(Type = "Control"))
    sort.order <- attr(TS.diffs.sort, which = "sortOrder")
    testthat::expect_is(TS.diffs.sort, "alldiffs")
    testthat::expect_true(validAlldiffs(TS.diffs.sort))
    testthat::expect_equal(nrow(TS.diffs.sort$predictions),20)
    testthat::expect_equal(ncol(TS.diffs.sort$predictions),8)
    testthat::expect_equal(as.character(TS.diffs.sort$predictions$Sources[1]),"Recycled water")
    testthat::expect_equal(length(attributes(TS.diffs.sort)),9)
    testthat::expect_equal(length(attr(TS.diffs.sort, which = "sortOrder")),6)
    
    #Test sort.alldiffs with supplied sortOrder
    m2.lmer <- lmerTest::lmer(Turbidity ~ Benches + (Sources * (Type + Species)) + 
                                (1|Benches:MainPlots),
                              data=na.omit(WaterRunoff.dat))
    TS.emm <- emmeans::emmeans(m2.lmer, specs = ~ Sources:Type)
    TS.preds <- summary(TS.emm)
    den.df <- min(TS.preds$df, na.rm = TRUE)
    ## Modify TS.preds to be compatible with a predictions.frame
    TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                     se = "SE", interval.type = "CI", 
                                     interval.names = c("lower.CL", "upper.CL"))
    
    ## Form an all.diffs object, sorting it using the pH sort.order and check its validity
    els <- as.numeric(rownames(TS.preds))
    TS.vcov <- vcov(TS.emm)[els,els]
    TS.diffs2.sort <- allDifferences(predictions = TS.preds, 
                                     classify = "Sources:Type", 
                                     vcov = TS.vcov, tdf = den.df,
                                     sortFactor = "Sources", 
                                     sortOrder = sort.order)
    validAlldiffs(TS.diffs2.sort)
    testthat::expect_equal(as.character(TS.diffs.sort$predictions$Sources[1]),
                           as.character(TS.diffs2.sort$predictions$Sources[1]))
    testthat::expect_equal(attr(TS.diffs.sort, which = "sortOrder"),
                           attr(TS.diffs2.sort, which = "sortOrder"))
  }
})

cat("#### Test for facCombine.alldiffs on Ladybird with lme4\n")
test_that("facCombine.alldiffs_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data("Ladybird.dat")
  
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
  
  if (exists("HCL.preds"))
  {
    ## Form an all.diffs object with predictions obtained with either asreml or lmerTest
    HCL.diffs <- as.alldiffs(predictions = HCL.preds, classify = "Host:Cadavers:Ladybird", 
                             sed = HCL.sed, vcov = HCL.vcov, tdf = den.df)
    
    ## check the class and validity of the alldiffs object
    testthat::expect_true(is.alldiffs(HCL.diffs))
    testthat::expect_true(validAlldiffs(HCL.diffs))
    testthat::expect_equal(nrow(HCL.diffs$predictions),12)
    testthat::expect_equal(ncol(HCL.diffs$predictions),9)
    
    ## Combine Cadavers and Ladybird
    HCL.Comb.diffs <- facCombine(HCL.diffs, factors = c("Cadavers","Ladybird"))
    
    ## check the validity of HCL.Comb.diffs
    validAlldiffs(HCL.Comb.diffs)
    testthat::expect_equal(nrow(HCL.Comb.diffs$predictions),12)
    testthat::expect_equal(ncol(HCL.Comb.diffs$predictions),8)
    testthat::expect_true(all(c("Host", "Cadavers_Ladybird", "predicted.value") %in% 
                                names(HCL.Comb.diffs$predictions)))
  }
  
})

cat("#### Test for as.predictions.frame on Ladybird with lme4\n")
test_that("as.predictions.frame_lme4", {
  #  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data("Ladybird.dat")
  
  ## Use lmerTest and emmmeans to get the predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(logitP ~ Host*Cadavers*Ladybird + (1|Run),
                              data=Ladybird.dat)
    HCL.emm <- emmeans::emmeans(m1.lmer, specs = ~ Host:Cadavers:Ladybird)
    HCL.preds <- summary(HCL.emm)
    den.df <- min(HCL.preds$df)
    ## Modify HCL.preds to be compatible with a predictions.frame
    HCL.preds <- as.predictions.frame(HCL.preds, predictions = "emmean", 
                                     se = "SE", interval.type = "CI", 
                                     interval.names = c("lower.CL", "upper.CL"))
    testthat::expect_true(validPredictionsFrame(HCL.preds))
    testthat::expect_silent(plotPredictions(HCL.preds, y = "predicted.value", "Host:Cadavers:Ladybird"))
  } 
})

cat("#### Test for addBacktransforms on WaterRunoff with lme4\n")
test_that("addBacktransforms_WaterRunoff_lme4", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asremlPlus)
  library(dae)
  data(WaterRunoff.dat)
  

  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(log.Turbidity ~ Benches + (Sources * (Type + Species)) + 
                                (1|Benches:MainPlots),
                              data=na.omit(WaterRunoff.dat))
    TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Species)
    TS.preds <- summary(TS.emm)
    den.df <- min(TS.preds$df, na.rm = TRUE)
    ## Modify TS.preds to be compatible with a predictions.frame
    TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                     se = "SE", interval.type = "CI", 
                                     interval.names = c("lower.CL", "upper.CL"))
    
    ## Form an all.diffs object and check its validity
    els <- as.numeric(rownames(TS.preds))
    TS.vcov <- vcov(TS.emm)[els,els]
    TS.diffs <- allDifferences(predictions = TS.preds, classify = "Sources:Species", 
                               vcov = TS.vcov, tdf = den.df)
    validAlldiffs(TS.diffs)
  }  
  
  ## Plot p-values for predictions obtained using asreml or lmerTest
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
