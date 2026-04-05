# Extracted from test42alldiffsasr.r:2212

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "asremlPlus", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("prediction_alldiffs")
cat("#### Test for allDifferences.data.frame sort.alldiffs on Oats with asreml42\n")
cat("#### Test for LSDs and halfLSIs on system data with asreml42\n")
cat("#### Test for LSD on Oats with asreml42\n")
cat("#### Test for sort.alldiffs on Smarthouse with asreml42\n")
cat("#### Test for LSD with sort.alldiffs on Smarthouse with asreml42\n")
cat("#### Test for LSDsupplied on Oats with asreml42\n")
cat("#### Test for single-prediction LSDs in 821 Barley with asreml42\n")
cat("#### Test for LSD on WaterRunoff with asreml42\n")
cat("#### Test for exploreLSDs on WaterRunoff with asreml42\n")
cat("#### Test for exploreLSDs on Oats with asreml42\n")
cat("#### Test for sort.alldiffs on WaterRunoff with asreml42\n")
cat("#### Test for sort.alldiffs on Oats with asreml42\n")
cat("#### Test for subset.alldiffs on Smarthouse with asreml42\n")
cat("#### Test for facCombine.alldiffs on Ladybird with asreml42\n")
cat("#### Test for facRecast.alldiffs on Ladybird with asreml42\n")
cat("#### Test for linear.transformation on Oats with asreml42\n")
cat("#### Test for linear.transformation on dat699 with asreml42\n")
cat("#### Test for linear.transformation on WaterRunoff with asreml42\n")

# test -------------------------------------------------------------------------
skip_if_not_installed("asreml")
skip_on_cran()
library(asreml)
library(asremlPlus)
library(dae)
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE)
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        data= WaterRunoff.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
diffs <- predictPlus(classify = "Sources:Species", Vmatrix = TRUE, 
                       asreml.obj = current.asr, tables = "none", 
                       wald.tab = current.asrt$wald.tab, 
                       present = c("Type","Species","Sources"))
diffs.sub <- linTransform(diffs, classify = "Sources:Species", Vmatrix = TRUE,
                            linear.transformation = ~ Sources + Species,
                            tables = "none")
testthat::expect_equal(diffs.sub$predictions$predicted.value[1] - 
                           diffs.sub$predictions$predicted.value[6],
                         diffs.sub$predictions$predicted.value[7] - 
                           diffs.sub$predictions$predicted.value[12])
L <- kronecker(diag(1, nrow = 4), 
                 cbind(diag(1, nrow = 5), matrix(rep(-1, 5), ncol = 1)))
L <- mat.dirsum(list(L, 
                       kronecker(diag(1, nrow = 2), 
                                 cbind(diag(1, nrow = 7), 
                                       matrix(rep(-1, 7), ncol = 1)))))
testthat::expect_silent(diffs.L.EGLS <- linTransform(diffs.sub, 
                                                   classify = "Sources:Species",
                                                   linear.transformation = L,
                                                   tables = "none"))
testthat::expect_silent(diffs.L <- linTransform(diffs.sub, 
                                                  classify = "Sources:Species",
                                                  linear.transformation = L, 
                                                  EGLS.linTransform = FALSE,
                                                  tables = "none"))
testthat::expect_false(isCompoundSymmetric(diffs.sub$vcov))
ksed <- diffs.L$sed
ksed <- na.omit(ksed[upper.tri(ksed)])
ksed <- ksed[!(ksed < 1e-08)]
testthat::expect_true(length(ksed) == diffs.L$LSD["c"] && diffs.L$LSD["c"] == 484)
testthat::expect_true(all(abs(diffs.L$predictions$predicted.value[c(1,6,11,16,21,28)] - 
                                  (diffs.sub$predictions$predicted.value[1] - 
                                     diffs.sub$predictions$predicted.value[6])) < 1e-06))
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE)
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                        random = ~ Benches:MainPlots,
                        data= WaterRunoff.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
diffs.sub <- predictPlus.asreml(classify = "Sources:Species", Vmatrix = TRUE, 
                                  linear.transformation = ~ Sources + Species,
                                  asreml.obj = current.asr, tables = "none", 
                                  wald.tab = current.asrt$wald.tab, 
                                  present = c("Type","Species","Sources"))
L <- cbind(matrix(rep(0,7*32), nrow = 7, ncol = 32),
             diag(1, nrow = 7), 
             matrix(rep(-1, 7), ncol = 1))
rownames(L) <- as.character(diffs.sub$predictions$Species[33:39])
diffs.L <- linTransform(diffs.sub, 
                          classify = "Sources:Species",
                          linear.transformation = L,
                          tables = "predictions")
testthat::expect_true(abs(diffs.L$predictions$predicted.value[1] + 0.04270763) < 1e-04)
testthat::expect_true(diffs.L$predictions$Combination[1] == "S. iqscjbogxah")
L1 <- matrix(L[1,1:40], nrow = 1)
diffs.L1 <- linTransform(diffs.sub, 
                           classify = "Sources:Species",
                           linear.transformation = L1,
                           tables = "predictions")
twoway <- data.frame(A = factor(rep(1:2, c(4,5))),
                       B = factor(rep(c(1:3,1:3), c(2,1,1,2,1,2))),
                       y = c(6,4,3,3,3,5,5,5,7))
twoway <- twoway[-4,]
mod <- do.call("asreml", 
                 list(y ~ A*B, data = twoway))
pred.full <- predict(mod, classify = "A:B")
testthat::expect_true(any(is.na(pred.full$pvals$predicted.value)))
diffs.sub <- predictPlus.asreml(asreml.obj = mod, classify = "A:B", 
                                  linear.transformation = ~ A + B,
                                  tables = "predictions", 
                                  error.intervals = "Stand")
testthat::expect_equal(nrow(diffs.sub$predictions), 5)
mod.add <- asreml(y ~ A+B, data = twoway)
pred.add <- predict(mod.add, classify = "A:B")
testthat::expect_true(!any(is.na(pred.add$pvals$predicted.value)))
testthat::expect_true(abs(diffs.sub$predictions$predicted.value[1] - 
                              pred.add$pvals$predicted.value[1]) < 1e-04)
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
asreml.options(keep.order = TRUE, ai.sing = TRUE, extra = 5)
for (k in 1:nresp)
  {
    fix <- paste(responses.lRGR[k], 
                 " ~ Smarthouse/(xMainPosn + xLane) + Genotype.ID*Treatment.1",
                 sep= "")
    HEB25.asr <- do.call("asreml",
                         args = list(fixed = as.formula(fix), 
                                     random = ~ Smarthouse:Zones:Mainplots, 
                                    # residual = ~ idh(Treat.Smarthouse):Zones:Mainplots, 
                                     data = cart.dat, workspace="1gb", 
                                     na.action=na.method(y="include", x="include"),
                                     maxit=50))
    summary(HEB25.asr)$varcomp
    current.asrt <- as.asrtests(HEB25.asr)
    current.asrt <- rmboundary(current.asrt)
    current.asr <- current.asrt$asreml.obj
    wald.tab <- recalcWaldTab(asrtests.obj=current.asrt, dDF.fault="residual")
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
                                              error.intervals = "Conf", 
                                              tables = "none")
  }
testthat::expect_equal(sum(unlist(lapply(save$lRGR_sm_32_42, is.null))),1)
