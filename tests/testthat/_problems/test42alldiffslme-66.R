# Extracted from test42alldiffslme.r:66

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("lme4_alldiffs")
cat("#### Test for predictions.frame on Oats with lme4\n")
cat("#### Test for LSD.frame on Oats with lme4\n")

# test -------------------------------------------------------------------------
library(asremlPlus)
library(dae)
data(Oats.dat)
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
