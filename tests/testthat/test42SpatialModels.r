#devtools::test("asremlPlus")
context("spatial_modelling")

cat("#### Test for makeTPPSplineMats with both wheat datasets with asreml42\n")
test_that("Wheat_makeTPPSplineMats_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data("wheat94.dat")
  wheat94.dat$R <- as.factor(wheat94.dat$row)
  wheat94.dat$C <- as.factor(wheat94.dat$col)
  
  #'## Use makeTPPSplineMats to produce design matrices
  tps.mbf <- makeTPPSplineMats(wheat94.dat, row.covar = "row", col.covar = "col", 
                               asreml.option = "mbf")
  tps.grp <- makeTPPSplineMats(wheat94.dat, row.covar = "row", col.covar = "col", 
                               asreml.option = "grp")
  
  #'## Generate the full length columns from tps.mbf 
  df.mbf <- tps.mbf[[1]]$data.plus
  BrZxx.df$TP.row <- factor(BrZxx.df$TP.row)
  tmp <- dplyr::left_join(df.mbf, BrZxx.df)
  TP.C.1 <- tmp[paste0("V",1:22)]*tmp$TP.C.1
  names(TP.C.1) <- paste0("TP.C.1_frow_",1:22)
  TP.C.2 <- tmp[paste0("V",1:22)]*tmp$TP.C.2
  names(TP.C.2) <- paste0("TP.C.2_frow_",1:22)
  BcZxx.df$TP.col <- factor(BcZxx.df$TP.col)
  tmp <- dplyr::left_join(df.mbf, BcZxx.df)
  TP.R.1 <- tmp[paste0("V",1:15)]*tmp$TP.R.1
  names(TP.R.1) <- paste0("TP.R.1_fcol_",1:15)
  TP.R.2 <- tmp[paste0("V",1:15)]*tmp$TP.R.2
  names(TP.R.2) <- paste0("TP.R.2_fcol_",1:15)
  df.mbf <- cbind(df.mbf, TP.C.1, TP.C.2, TP.R.1, TP.R.2)
  
  #'## Show that both sets of fixed terms are equal
  fix.terms <- c("TP.C.1","TP.C.2","TP.R.1","TP.R.2")
  testthat::expect_true(all.equal(df.mbf[fix.terms], tps.grp[[1]]$data.plus[fix.terms]))
  
  #'## Test whether columns generated from the mbf data.frames are the same as the grp data.frame columns 
  cols <- c(paste0("TP.C.1_frow_",1:10), paste0("TP.C.2_frow_",1:10), 
            paste0("TP.R.1_fcol_",1:15), paste0("TP.R.2_fcol_",1:15))
  testthat::expect_true(all.equal(df.mbf[cols], tps.grp[[1]]$data.plus[cols]))
  
  #'## Test whether random interaction columns are the smae
  cols <- paste0("TP_fcol_frow_", 1:330)
  names(BcrZxx.df)[1:330] <- cols
  testthat::expect_true(all.equal(BcrZxx.df[cols], tps.grp[[1]]$data.plus[cols]))
  
  #'## Repeat for wheat data from asremlPlus
  data(Wheat.dat)
  
  #Add row and column covariates
  tmp.dat <- within(Wheat.dat, 
                    {
                      cColumn <- dae::as.numfac(Column)
                      cColumn <- cColumn  - mean(unique(cColumn))
                      cRow <- dae::as.numfac(Row)
                      cRow <- cRow - mean(unique(cRow))
                    })
  
  #'## Use makeTPPSplineMats to produce design matrices
  rm(list = fix.terms)
  rm(BcZxx.df, BrZxx.df, BcrZxx.df)
  tps.mbf <- makeTPPSplineMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn", 
                                 asreml.option = "mbf")
  tps.grp <- makeTPPSplineMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn", 
                               asreml.option = "grp")
  
  #Check makeTPPSplineMats with grp
  testthat::expect_true(all(names(tps.grp[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                     "BcrZ.df","dim","trace","grp","data.plus")))
  testthat::expect_true(all(names(tps.grp[[1]]$data.plus[,1:19]) == 
                              c("Rep","Row","Column", "WithinColPairs","Variety","yield",
                                "cRow","cColumn","TP.col","TP.row",
                                "TP.CxR","TP.C.1","TP.C.2","TP.R.1","TP.R.2", 
                                "TP.CR.1","TP.CR.2","TP.CR.3","TP.CR.4")))
  testthat::expect_true(all(grepl("TP\\.",names(tps.grp[[1]]$data.plus[,20:50]))))
  testthat::expect_true(all(grepl("TP\\_",names(tps.grp[[1]]$data.plus)[81:ncol(tps.grp[[1]]$data.plus)])))

  #Check makeTPPSplineMats with mbf
  testthat::expect_true(exists("tps.mbf"))
  testthat::expect_true(exists("BcZxx.df"))
  testthat::expect_true(all(names(tps.mbf[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                     "BcrZ.df","dim","trace","data.plus")))
  testthat::expect_true(all(names(tps.mbf[[1]]$data.plus) == 
                              c("cRow","cColumn","Rep","Row","Column", 
                                "WithinColPairs","Variety","yield","TP.col","TP.row",
                                "TP.CxR","TP.C.1","TP.C.2","TP.R.1","TP.R.2", 
                                "TP.CR.1","TP.CR.2","TP.CR.3","TP.CR.4")))

  ## Compare tps.mbf and tps.grp
  # Generate the full length columns from tps.mbf 
  df.mbf <- tps.mbf[[1]]$data.plus
  BrZxx.df$TP.row <- factor(BrZxx.df$TP.row)
  tmp <- dplyr::left_join(df.mbf, BrZxx.df)
  TP.C.1 <- tmp[paste0("V",1:10)]*tmp$TP.C.1
  names(TP.C.1) <- paste0("TP.C.1_frow_",1:10)
  TP.C.2 <- tmp[paste0("V",1:10)]*tmp$TP.C.2
  names(TP.C.2) <- paste0("TP.C.2_frow_",1:10)
  BcZxx.df$TP.col <- factor(BcZxx.df$TP.col)
  tmp <- dplyr::left_join(df.mbf, BcZxx.df)
  TP.R.1 <- tmp[paste0("V",1:15)]*tmp$TP.R.1
  names(TP.R.1) <- paste0("TP.R.1_fcol_",1:15)
  TP.R.2 <- tmp[paste0("V",1:15)]*tmp$TP.R.2
  names(TP.R.2) <- paste0("TP.R.2_fcol_",1:15)
  df.mbf <- cbind(df.mbf, TP.C.1, TP.C.2, TP.R.1, TP.R.2)
  
  #'## Show that both sets of fixed terms are equal
  fix.terms <- c("TP.C.1","TP.C.2","TP.R.1","TP.R.2")
  testthat::expect_true(all.equal(df.mbf[fix.terms], tps.grp[[1]]$data.plus[fix.terms]))
  
  #'## Test whether columns generated from the mbf data.frames are the same as the grp data.frame columns 
  cols <- c(paste0("TP.C.1_frow_",1:10), paste0("TP.C.2_frow_",1:10), 
            paste0("TP.R.1_fcol_",1:15), paste0("TP.R.2_fcol_",1:15))
  testthat::expect_true(all.equal(df.mbf[cols], tps.grp[[1]]$data.plus[cols]))
  
  #'## Test whether random interaction columns are the same
  cols <- paste0("TP_fcol_frow_", 1:150)
  names(BcrZxx.df)[1:150] <- cols
  testthat::expect_true(all.equal(BcrZxx.df[cols], tps.grp[[1]]$data.plus[cols]))
  
  #'## Use makeTPPSplineMats to produce design matrices for theta = c(0,60)
  rm(list = fix.terms)
  rm(BcZxx.df, BrZxx.df, BcrZxx.df)
  testthat::expect_warning(
    tps.mbf_0_60 <- makeTPPSplineMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn", 
                                       asreml.option = "mbf", theta = c(0,60)), 
    regexp = "The following objects are being overwritten: BcZxx.df, BrZxx.df, BcrZxx.df")
  tps.grp_0_60 <- makeTPPSplineMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn", 
                                     asreml.option = "grp", theta = c(0,60))
  
  
  #'## Generate the full length columns from tps.mbf_0_60 (used to cause problems in asreml)
  df.mbf <- tps.mbf_0_60[[1]]$data.plus
  BrZxx.df$TP.row <- factor(BrZxx.df$TP.row)
  tmp <- dplyr::left_join(df.mbf, BrZxx.df)
  TP.C.1 <- tmp[paste0("V",1:10)]*tmp$TP.C.1
  names(TP.C.1) <- paste0("TP.C.1_frow_",1:10)
  TP.C.2 <- tmp[paste0("V",1:10)]*tmp$TP.C.2
  names(TP.C.2) <- paste0("TP.C.2_frow_",1:10)
  BcZxx.df$TP.col <- factor(BcZxx.df$TP.col)
  tmp <- dplyr::left_join(df.mbf, BcZxx.df)
  TP.R.1 <- tmp[paste0("V",1:15)]*tmp$TP.R.1
  names(TP.R.1) <- paste0("TP.R.1_fcol_",1:15)
  TP.R.2 <- tmp[paste0("V",1:15)]*tmp$TP.R.2
  names(TP.R.2) <- paste0("TP.R.2_fcol_",1:15)
  df.mbf <- cbind(df.mbf, TP.C.1, TP.C.2, TP.R.1, TP.R.2)
  
  #'## Show that both sets of fixed terms are equal
  fix.terms <- c("TP.C.1","TP.C.2","TP.R.1","TP.R.2")
  testthat::expect_true(all.equal(df.mbf[fix.terms], tps.grp_0_60[[1]]$data.plus[fix.terms]))
  
  #'## Test whether columns generated from the mbf data.frames are the same as the grp data.frame columns 
  cols <- c(paste0("TP.C.1_frow_",1:10), paste0("TP.C.2_frow_",1:10), 
            paste0("TP.R.1_fcol_",1:15), paste0("TP.R.2_fcol_",1:15))
  testthat::expect_true(all.equal(df.mbf[cols], tps.grp_0_60[[1]]$data.plus[cols]))
  
  #'## Test whether random interaction columns are the same
  cols <- paste0("TP_fcol_frow_", 1:150)
  names(BcrZxx.df)[1:150] <- cols
  testthat::expect_true(all.equal(BcrZxx.df[cols], tps.grp_0_60[[1]]$data.plus[cols]))
  
  #Test mbf.env = NULL
  if (exists("BcZxx.df")) rm("BcZxx.df")
  testthat::expect_silent(
    tps.mbf <- makeTPPSplineMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn", mbf.env = NULL))
  testthat::expect_true(exists("tps.mbf"))
  testthat::expect_true(exists("BcZxx.df"))
  
  #Test trapping of illegal nsect argument
  testthat::expect_error(
    tps <- makeTPPSplineMats(tmp.dat, row.covar = "cRow", col.covar = "cColumn", nsect = 2, 
                             asreml.option = "grp"),
    regexp = "the argument\\(s\\) nsect are not legal arguments for 'tpsmmb'")
})

cat("#### Test spatial modelling for chick pea example with asreml42\n")
test_that("chickpea_spatial_mod_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  

  data(chkpeadat)
  tmp.dat <- within(chkpeadat, 
                    {
                      vMPosn <- as.numfac(fac.recast(Mainplot, newlevels = rep(1:11, times = 4)))
                      vMPosn <- vMPosn - mean(unique(vMPosn))
                    })
  #Remove some Lanes in SE Smarthouse so have different grids in the Smarthouses
  tmp.dat <- tmp.dat[-(1013:1056), ]
  (table(tmp.dat$Smarthouse))
  
  asreml.options(design = TRUE)
  current.asr <- do.call(asreml, 
                         list(fixed = Biomass.plant ~ Smarthouse + Lines * TRT, 
                              random = ~Smarthouse:Zone/Mainplot, 
                              data = tmp.dat))
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Lane and Position effects")
  init.asrt <- rmboundary(init.asrt)
  
  #Test makeTPPSplineMats with sections and grp
  tps.grp <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                               row.covar = "vLanes", col.covar = "vMPosn",
                               asreml.option = "grp")
  testthat::expect_true(all(names(tps.grp) == c("SW","SE")))
  testthat::expect_true(all(names(tps.grp[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                 "BcrZ.df","dim","trace","grp","data.plus")))
  testthat::expect_true(all(names(tps.grp[[1]]$data.plus[,1:19]) == 
                              c("Smarthouse","Lane","Position","Zone","vLanes","vPos",
                                "Mainplot","Subplot","Lines","TRT","Rep",
                                "X100.SW","Biomass.plant","Pods.plant","Filled.pods.plant", 
                                "Empty.pods.plant","Seed.No.plant","Seed.weight.plant","vMPosn")))
  testthat::expect_true(all(grepl("TP\\.",names(tps.grp[[1]]$data.plus[,20:100]))))
  testthat::expect_true(all(grepl("TP\\_",names(tps.grp[[1]]$data.plus)[101:ncol(tps.grp[[1]]$data.plus)])))
  testthat::expect_equal(tps.grp[[1]]$grp$TP.C.1_frow[1], tps.grp[[1]]$grp$All[1])
  testthat::expect_equal(length(tps.grp[[1]]$grp$All), 334)
  
  tps.mbf <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                               row.covar = "vLanes", col.covar = "vMPosn",
                               asreml.option = "mbf")
  testthat::expect_true(all(names(tps.mbf) == c("SW","SE")))
  testthat::expect_true(all(names(tps.mbf[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                     "BcrZ.df","dim","trace","data.plus")))
  testthat::expect_equal(ncol(tps.mbf[[1]]$data.plus), 30)
  testthat::expect_true(all(grepl("TP\\.",names(tps.mbf[[1]]$data.plus[,20:20]))))
  

  #'# Fit SE first using mbf
  dat <- tps.mbf[["SE"]]$data.plus
  mbf.lis <- tps.mbf[["SE"]]$mbflist
  asr.mbf.SE <- do.call(asreml, 
                        list(fixed = Biomass.plant ~ Smarthouse + Lines * TRT +  
                               at(Smarthouse,  'SE'):TP.CR.2 + 
                               at(Smarthouse,  'SE'):TP.CR.3 + at(Smarthouse,  'SE'):TP.CR.4,
                             random = ~ Smarthouse:Zone/Mainplot +  
                               at(Smarthouse,  'SE'):TP.C.1:mbf(TP.row) + 
                               at(Smarthouse,  'SE'):TP.C.2:mbf(TP.row) + 
                               at(Smarthouse,  'SE'):TP.R.1:mbf(TP.col) + 
                               at(Smarthouse,  'SE'):TP.R.2:mbf(TP.col) + 
                               at(Smarthouse,  'SE'):mbf(TP.CxR),
                             mbf = mbf.lis, data = dat))
  
  des.mbf.SE <- asr.mbf.SE$design
  des.mbf.SE <- as.data.frame(as.matrix(des.mbf.SE))
  cols.mbf.SE <- c(paste0("at(Smarthouse, 'SE'):TP.C.1:mbf(TP.row)_V",1:22), 
                   paste0("at(Smarthouse, 'SE'):TP.C.2:mbf(TP.row)_V",1:22), 
                   paste0("at(Smarthouse, 'SE'):TP.R.1:mbf(TP.col)_V",1:11), 
                   paste0("at(Smarthouse, 'SE'):TP.R.2:mbf(TP.col)_V",1:11))
  
  #'## Add SW using mbf
  dat <- tps.mbf[["SW"]]$data.plus
  mbf.lis <- tps.mbf[["SW"]]$mbflist
  asr.mbf.SW <- update(asr.mbf.SE, 
                       fixed. = . ~ . + at(Smarthouse,  'SW'):TP.CR.2 + 
                         at(Smarthouse,  'SW'):TP.CR.3 + at(Smarthouse,  'SW'):TP.CR.4,
                       random. = ~ . + at(Smarthouse,  'SW'):TP.C.1:mbf(TP.row) + 
                         at(Smarthouse,  'SW'):TP.C.2:mbf(TP.row) + 
                         at(Smarthouse,  'SW'):TP.R.1:mbf(TP.col) + 
                         at(Smarthouse,  'SW'):TP.R.2:mbf(TP.col) + 
                         at(Smarthouse,  'SW'):mbf(TP.CxR),
                       mbf = mbf.lis, data = dat)
  des.mbf.SW <- asr.mbf.SW$design
  des.mbf.SW <- as.data.frame(as.matrix(des.mbf.SW))
  cols.mbf.SW <- c(paste0("TP.C.1:mbf(TP.row)_V",1:24,":at(Smarthouse, 'SW')"), 
                   paste0("mbf(TP.row)_V",1:24,":TP.C.2:at(Smarthouse, 'SW')"), 
                   paste0("TP.R.1:mbf(TP.col)_V",1:11,":at(Smarthouse, 'SW')"), 
                   paste0("mbf(TP.col)_V",1:11,":TP.R.2:at(Smarthouse, 'SW')"))
  
  #'# Fit SE first using grp
  dat <- tps.grp[["SE"]]$data.plus
  grp.lis <- tps.grp[["SE"]]$grp
  asr.grp.SE <- do.call(asreml, 
                        list(fixed = Biomass.plant ~ Smarthouse + Lines * TRT +  
                               at(Smarthouse,  'SE'):TP.CR.2 + 
                               at(Smarthouse,  'SE'):TP.CR.3 + at(Smarthouse,  'SE'):TP.CR.4,
                             random = ~ Smarthouse:Zone/Mainplot +  
                               at(Smarthouse,  'SE'):grp(TP.C.1_frow) + 
                               at(Smarthouse,  'SE'):grp(TP.C.2_frow) + 
                               at(Smarthouse,  'SE'):grp(TP.R.1_fcol) + 
                               at(Smarthouse,  'SE'):grp(TP.R.2_fcol) + 
                               at(Smarthouse,  'SE'):grp(TP_fcol_frow),
                             group = grp.lis, data = dat))
  
  des.grp.SE <- asr.grp.SE$design
  des.grp.SE <- as.data.frame(as.matrix(des.grp.SE))
  cols.grp.SE <- c(paste0("at(Smarthouse, 'SE'):grp(TP.C.1_frow)_TP.C.1_frow_", 1:22), 
                   paste0("at(Smarthouse, 'SE'):grp(TP.C.2_frow)_TP.C.2_frow_", 1:22), 
                   paste0("at(Smarthouse, 'SE'):grp(TP.R.1_fcol)_TP.R.1_fcol_", 1:11), 
                   paste0("at(Smarthouse, 'SE'):grp(TP.R.2_fcol)_TP.R.2_fcol_", 1:11))
  
  #'## Add SW using mbf
  dat <- tps.grp[["SW"]]$data.plus
  grp.lis <- tps.grp[["SW"]]$grp
  asr.grp.SW <- update(asr.mbf.SE, 
                       fixed. = . ~ .+ at(Smarthouse,  'SW'):TP.CR.2 + 
                         at(Smarthouse,  'SW'):TP.CR.3 + at(Smarthouse,  'SW'):TP.CR.4,
                       random. = ~ Smarthouse:Zone/Mainplot +  
                         at(Smarthouse,  'SW'):grp(TP.C.1_frow) + 
                         at(Smarthouse,  'SW'):grp(TP.C.2_frow) + 
                         at(Smarthouse,  'SW'):grp(TP.R.1_fcol) + 
                         at(Smarthouse,  'SW'):grp(TP.R.2_fcol) + 
                         at(Smarthouse,  'SW'):grp(TP_fcol_frow),
                       group = grp.lis, data = dat)
  des.grp.SW <- asr.grp.SW$design
  des.grp.SW <- as.data.frame(as.matrix(des.grp.SW))
  cols.grp.SW <- c(paste0("at(Smarthouse, 'SW'):grp(TP.C.1_frow)_TP.C.1_frow_", 1:24), 
                   paste0("at(Smarthouse, 'SW'):grp(TP.C.2_frow)_TP.C.2_frow_", 1:24), 
                   paste0("at(Smarthouse, 'SW'):grp(TP.R.1_fcol)_TP.R.1_fcol_", 1:11), 
                   paste0("at(Smarthouse, 'SW'):grp(TP.R.2_fcol)_TP.R.2_fcol_", 1:11))

  #### The following shows (i) that mbf  and grp give the same design matrices, 
  ##   (ii) that the design matrices in asreml are the same for mbf and grp, and 
  ##   (iii) that the design matrices in asreml are not equivalent to those generated in tps
  
  #Compare design matrices for SE in asreml.obj - TRUE
  des.mbf.SE <- des.mbf.SE[529:1012,cols.mbf.SE]
  des.grp.SE <- des.grp.SE[529:1012,cols.grp.SE]
  names(des.grp.SE) <- cols.mbf.SE
  testthat::expect_true(all.equal(des.mbf.SE, des.grp.SE))
  
  #Compare design matrices for SW in asreml.obj - TRUE
  des.mbf.SW <- des.mbf.SW[1:528,cols.mbf.SW]
  des.grp.SW <- des.grp.SW[1:528,cols.grp.SW]
  names(des.grp.SW) <- cols.mbf.SW
  testthat::expect_true(all.equal(des.mbf.SE, des.grp.SE))
  
  #Get the design matrices from the tps objects and test if equal to those from asreml.obj
  #Neither is TRUE
  dat.grp.SE <- tps.grp[["SE"]]$data.plus[529:1012,c(31:52, 55:76, 79:100)]
  colnames(dat.grp.SE) <- colnames(des.grp.SE) <- cols.grp.SE
  testthat::show_failure(all.equal(dat.grp.SE, des.grp.SE))
  dat.grp.SW <- tps.grp[["SW"]]$data.plus[1:528,c(31:100)]
  colnames(dat.grp.SW) <- colnames(des.grp.SW) <- cols.grp.SW
  testthat::show_failure(all.equal(dat.grp.SW, des.grp.SW))
  
  #'## Generate the full length columns from tps.mbf 
  df.mbf <- tps.mbf[[1]]$data.plus
  BrZSE.df$TP.row <- factor(BrZSE.df$TP.row)
  tmp <- dplyr::left_join(df.mbf, BrZSE.df)
  TP.C.1 <- tmp[paste0("V",1:22)]*tmp$TP.C.1
  names(TP.C.1) <- paste0("at(Smarthouse, 'SE'):TP.C.1:mbf(TP.row)_V",1:22)
  TP.C.2 <- tmp[paste0("V",1:22)]*tmp$TP.C.2
  names(TP.C.2) <- paste0("at(Smarthouse, 'SE'):TP.C.2:mbf(TP.row)_V",1:22)
  BcZSE.df$TP.col <- factor(BcZSE.df$TP.col)
  tmp <- dplyr::left_join(df.mbf, BcZSE.df)
  TP.R.1 <- tmp[paste0("V",1:11)]*tmp$TP.R.1
  names(TP.R.1) <- paste0("at(Smarthouse, 'SE'):TP.R.1:mbf(TP.col)_V",1:11)
  TP.R.2 <- tmp[paste0("V",1:11)]*tmp$TP.R.2
  names(TP.R.2) <- paste0("at(Smarthouse, 'SE'):TP.R.2:mbf(TP.col)_V",1:11)
  df.mbf <- cbind(df.mbf, TP.C.1, TP.C.2, TP.R.1, TP.R.2)
  cols.df.SE <- c(names(TP.C.1), names(TP.C.2), names(TP.R.1), names(TP.R.2))
 
   #check that cols generated from mbf df are the same as the tps.grp cols - TRUE
  colnames(dat.grp.SE) <- colnames(des.grp.SE) <- cols.df.SE
  all.equal(df.mbf[529:1012,cols.df.SE], dat.grp.SE) #Only the names differ
  #col generated from mbf are not the same as those in the asreml design matrix - FALSE
  colnames(des.mbf.SE) <- cols.df.SE
  testthat::show_failure(all.equal(df.mbf[529:1012,cols.df.SE], des.mbf.SE)) #Only the names differ

  #'## Test whether random interaction columns are the same in desing for mbf and grp = TRUE
  des.grp.SE.CxR <- asr.grp.SE$design
  des.grp.SE.CxR <- des.grp.SE.CxR[, grepl("TP_fcol_frow", colnames(des.grp.SE.CxR))]
  des.mbf.SE.CxR <- asr.mbf.SE$design
  des.mbf.SE.CxR <- des.mbf.SE.CxR[, grepl("TP.CxR", colnames(des.mbf.SE.CxR))]
  colnames(des.grp.SE.CxR) <- colnames(des.mbf.SE.CxR)
  testthat::expect_true(all.equal(des.grp.SE.CxR, des.mbf.SE.CxR))
})

cat("#### Test for wheat76 spatial models with asreml42\n")
test_that("Wheat_spatial_models_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(Wheat.dat)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  #Add row and column covariates
  tmp.dat <- within(Wheat.dat, 
                    {
                      cColumn <- dae::as.numfac(Column)
                      cColumn <- cColumn  - mean(unique(cColumn))
                      cRow <- dae::as.numfac(Row)
                      cRow <- cRow - mean(unique(cRow))
                    })
  
  #Fit initial model - Row and column random
  asreml.options(design = TRUE)
  current.asr <- do.call(asreml, 
                         list(yield ~ Rep + WithinColPairs + Variety, 
                              random = ~ Row + Column,
                              data=tmp.dat, maxit = 50))
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1720.891), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Row and Column effects")
  
  
  # Check for and remove any boundary terms
  init.asrt <- rmboundary(init.asrt, IClikelihood = "full")
  testthat::expect_lt(abs(init.asrt$test.summary$AIC - 1720.891), 0.50)

  # Try call with illegal argument
  testthat::expect_error(
    current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                        row.covar = "cRow", col.covar = "cColumn",
                                        dropRowterm = "Row", dropColterm = "Column",
                                        nsect = 2,
                                        asreml.option = "grp"), 
    regexp = "the argument\\(s\\) nsect are not legal arguments for 'changeModelOnIC', 'asreml'")
  
  # Try TPPS model with grp
  grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(grp.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1643.467), 0.10)
  testthat::expect_equal(rownames(summary(grp.asrt$asreml.obj)$varcomp), 
                         c("grp(TP.C.2_frow)", "dev(cRow)", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)!cor", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)_1", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)_2", 
                           "grp(TP_fcol_frow)", "units!R"))
  testthat::expect_equal(rownames(grp.asrt$wald.tab), c("(Intercept)", "Rep", "WithinColPairs", 
                                                        "Variety", 
                                                        "TP.CR.2", "TP.CR.3", "TP.CR.4"))
  
  
  #Repeat to make sure no carry-over effects 
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1643.467), 0.10)
  
  # Try TPPS model with mbf
  mbf.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  dropRowterm = "Row", dropColterm = "Column",
                                  asreml.option = "mbf")
  info <- infoCriteria(list(grp.asrt$asreml.obj, mbf.asrt$asreml.obj), IClikelihood = "full")
  testthat::expect_true(all.equal(info[1,], info[2,], tolerance = 1e-06, check.attributes = FALSE)) #mbf & grp are same
  
  
  #Rotate the penalty matrix with grp, testing if full rotated model is better than full, unrotated model
  grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  dropRowterm = "Row", dropColterm = "Column", 
                                  rotateX = TRUE, ngridangles = c(9,9), 
                                  asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1643.467), 0.10) #this unrotated AIC
  
  
  #Rotate the penalty matrix with grp
  grp.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  dropRowterm = "Row", dropColterm = "Column", 
                                  rotateX = TRUE, ngridangles = c(9,9), 
                                  asreml.option = "grp")
  info <- infoCriteria(grp.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1650.192), 0.10)
  testthat::expect_equal(rownames(summary(grp.asrt$asreml.obj)$varcomp), 
                         c("grp(TP.C.2_frow)", "dev(cRow)", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)!cor", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)_1", 
                           "grp(TP.R.1_fcol)+grp(TP.R.2_fcol)!corh(2)_2", 
                           "grp(TP_fcol_frow)", "units!R"))
  testthat::expect_equal(rownames(grp.asrt$wald.tab), c("(Intercept)", "Rep", "WithinColPairs", 
                                                        "Variety", 
                                                        "TP.CR.2", "TP.CR.3", "TP.CR.4"))
  testthat::expect_true(all(attr(grp.asrt, which = "theta.opt")[[1]] == c(20,60)))
  
  #Rotate the penalty matrix with mbf
  mbf.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  dropRowterm = "Row", dropColterm = "Column", 
                                  rotateX = TRUE, ngridangles = c(9,9), 
                                  asreml.option = "mbf")
  info <- infoCriteria(list(grp.asrt$asreml.obj, mbf.asrt$asreml.obj), IClikelihood = "full")
  testthat::expect_true(all.equal(info[1,], info[2,], check.attributes = FALSE))
  info <- infoCriteria(mbf.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1650.192), 0.10)
  testthat::expect_equal(rownames(summary(mbf.asrt$asreml.obj)$varcomp), 
                         c("dev(cRow)", "mbf(TP.row):TP.C.2", 
                           "TP.R.1:mbf(TP.col)+mbf(TP.col):TP.R.2!corh(2)!cor", 
                           "TP.R.1:mbf(TP.col)+mbf(TP.col):TP.R.2!corh(2)_1", 
                           "TP.R.1:mbf(TP.col)+mbf(TP.col):TP.R.2!corh(2)_2", 
                           "mbf(TP.CxR)", "units!R"))
  testthat::expect_equal(rownames(mbf.asrt$wald.tab), c("(Intercept)", "Rep", "WithinColPairs", 
                                                        "Variety", 
                                                        "TP.CR.2", "TP.CR.3", "TP.CR.4"))
  testthat::expect_true(all(attr(mbf.asrt, which = "theta.opt")[[1]] == c(20,60)))
  
  #Compare the data.frames on which the rotated fits are based
  grp.dat <- grp.asrt$asreml.obj$call$data
  mbf.dat <- mbf.asrt$asreml.obj$call$data
  testthat::expect_true(all.equal(grp.dat[c("TP.col","TP.row","TP.CxR",
                                            "TP.C.1","TP.C.2","TP.R.1","TP.R.2", 
                                            "TP.CR.1","TP.CR.2","TP.CR.3","TP.CR.4")], 
                                  mbf.dat[c("TP.col","TP.row","TP.CxR",
                                            "TP.C.1","TP.C.2","TP.R.1","TP.R.2", 
                                            "TP.CR.1","TP.CR.2","TP.CR.3","TP.CR.4")]))
  
  
  # Try TPNCSS model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPNCSS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1639.792), 0.10)
  
  # Try corr model - at the moment this fit fails because the addition of a unit terms 
  # clashes with the addition of a variance term for cRow:exp(cColumn)
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column", 
                                      IClikelihood = "full")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1653.096), 0.10)
  
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn", 
                                      row.factor = "Row", col.factor = "Column", 
                                      row.corrFitfirst = FALSE,
                                      IClikelihood = "full")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1653.096), 0.10)

  #Return all of the models
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 4)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS", "TPP1LS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == 
                              c("nonspatial", "corr", "TPNCSS", "TPPCS", "TPP1LS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1720.891, 1653.096, 1639.792, 1643.467, 1710.282)) < 0.10))
  testthat::expect_equal(spatial.asrts$best.spatial.mod, "TPNCSS")
  
  #Fit two models and return both
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 2)
  testthat::expect_equal(names(spatial.asrts$asrts), c("TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - c(1720.891, 1639.792, 1643.467)) < 0.10))
  
  #Fit all models with Row and Column random and return all
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr", "TPN", "TPPC", "TPP1"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 4)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS", "TPP1LS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == 
                              c("nonspatial", "corr", "TPNCSS", "TPPCS", "TPP1LS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1720.891, 1653.096, 1639.792, 1643.467, 1710.282)) < 0.10))
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                                 row.covar = "cRow", col.covar = "cColumn",
                                                 row.factor = "Row", col.factor = "Column")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   dropRowterm = "Row", dropColterm = "Column")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "cRow", col.covar = "cColumn",
                                                  dropRowterm = "Row", dropColterm = "Column",
                                                  asreml.option = "grp")
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   dropRowterm = "Row", dropColterm = "Column",
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
  infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
  testthat::expect_true(all.equal(spatial.asrts$spatial.IC[-1,], infoEach[,-3], 
                                  tolerance = 0.5))

  #Test rotateX with grp and no parallel processing
  spatial.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  dropRowterm = "Row", dropColterm = "Column",
                                  rotateX = TRUE, ngridangles = c(9,9),
                                  asreml.option = "grp")
  testthat::expect_true(abs(infoCriteria(spatial.asrt$asreml.obj, IClikelihood = "full")$AIC - 1650.192) < 1e-03)
  
  #Test rotateX with grp and parallel processing
  spatial.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  dropRowterm = "Row", dropColterm = "Column",
                                  rotateX = TRUE, ngridangles = c(9,9), 
                                  which.rotacriterion = "AIC", 
                                  nrotacores = parallel::detectCores(), 
                                  asreml.option = "grp")
  testthat::expect_true(abs(infoCriteria(spatial.asrt$asreml.obj, IClikelihood = "full")$AIC - 1650.192) < 1e-03)
  
  #Test rotateX with mbf and no parallel processing
  spatial.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  dropRowterm = "Row", dropColterm = "Column",
                                  rotateX = TRUE, ngridangles = c(9,9),
                                  asreml.option = "mbf")
  testthat::expect_true(abs(infoCriteria(spatial.asrt$asreml.obj, IClikelihood = "full")$AIC - 1650.192) < 1e-03)

  #Test rotateX with mbf and parallel processing - gives error
  testthat::expect_error(
    spatial.asrt <- addSpatialModel(init.asrt, spatial.model = "TPPS",
                                    row.covar = "cRow", col.covar = "cColumn",
                                    dropRowterm = "Row", dropColterm = "Column",
                                    rotateX = TRUE, ngridangles = c(9,9),
                                    which.rotacriterion = "AIC",
                                    nrotacores = parallel::detectCores(),
                                    asreml.option = "mbf"), 
    regexp = paste("Parallel processing has not been implemented for asreml.option set to mbf;",
                   "nrotacores must be one"))
  # testthat::expect_true(abs(infoCriteria(spatial.asrt$asreml.obj, IClikelihood = "full")$AIC - 1650.192) < 1e-03)
  
  #Fit initial model - Row and column fixed
  current.asr <- do.call(asreml, 
                         list(yield ~ Rep + WithinColPairs + Row + Column + Variety, 
                              data=tmp.dat, maxit = 50))
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 1)
  testthat::expect_lt(abs(info$AIC - 1690.964), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Row and Column effects")
  init.asrt <- rmboundary(init.asrt)
  
  # Try a TPNCSS model with fixed Row and Column
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPNCSS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 6)
  testthat::expect_lt(abs(info$AIC - 1639.792), 0.10)
  facs <- c("Row", "Column")
  #Check Row and COlumn terms not in model
  testthat::expect_false(any(facs %in% rownames(current.asrt$wald.tab)) &&
                           any(facs %in% names(current.asrt$asreml.obj$vparameters)))
  
  # Try TPPS model with fixed Row and Column
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "grp")
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1643.467), 0.10)
  #Check Row and COlumn terms not in model
  facs <- c("Row", "Column")
  testthat::expect_false(any(facs %in% rownames(current.asrt$wald.tab)) &&
                           any(facs %in% names(current.asrt$asreml.obj$vparameters)))
  
  #Fit all models with Row and Column fixed and return all
  #NB chooseSpatialModel returns that spatialICs for the best fitting correlation model
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = "all", 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column",
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 4)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS", "TPP1LS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == 
                              c("nonspatial", "corr", "TPNCSS", "TPPCS", "TPP1LS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1690.964, 1653.978, 1639.792, 1643.467, 1690.964)) < 0.10))
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                                     row.covar = "cRow", col.covar = "cColumn",
                                                     row.factor = "Row", col.factor = "Column")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   dropRowterm = "Row", dropColterm = "Column")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "cRow", col.covar = "cColumn",
                                                  dropRowterm = "Row", dropColterm = "Column",
                                                  asreml.option = "grp")
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS",
                                                   row.covar = "cRow", col.covar = "cColumn",
                                                   dropRowterm = "Row", dropColterm = "Column",
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
  infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
  testthat::expect_true(all.equal(spatial.asrts$spatial.IC[2:4,], infoEach[1:3 ,-3], 
                                  tolerance = 0.5))
})

cat("#### Test for wheat76 spatial models using mbf with asreml42\n")
test_that("Wheat_spatial_models_mbf_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(Wheat.dat)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  #Add row and column covariates
  tmp.dat <- within(Wheat.dat, 
                    {
                      cColumn <- dae::as.numfac(Column)
                      cColumn <- cColumn  - mean(unique(cColumn))
                      cRow <- dae::as.numfac(Row)
                      cRow <- cRow - mean(unique(cRow))
                    })
  
  #Fit initial model - Row and column random
  current.asr <- do.call(asreml, 
                         list(yield ~ Rep + WithinColPairs + Variety, 
                              random = ~ Row + Column,
                              data=tmp.dat, maxit = 50))
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1720.891), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Row and Column effects")
  
  
  # Check for and remove any boundary terms
  init.asrt <- rmboundary(init.asrt, IClikelihood = "full")
  testthat::expect_lt(abs(init.asrt$test.summary$AIC - 1720.891), 0.50)
  
  # Try TPPS model using mbf - does not fit the same model as grp
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      dropRowterm = "Row", dropColterm = "Column",
                                      asreml.option = "mbf", 
                                      update = FALSE)
  
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 7)
  testthat::expect_lt(abs(info$AIC - 1643.467), 0.10)
})

cat("#### Test for wheat703 corr spatial models with asreml42\n")
test_that("Wheat703_corr_models_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(indv703.dat)
  data(summ703)
  
# summ <- list()
  for (kresp in responses.test)
  { 
    mod.ch <- paste(kresp, "~ Block + Line")
    mod <- as.formula(mod.ch)
    cat("\n\n#### ",mod.ch,"\n\n")
    asreml.options(keep.order = TRUE)
    current.asr <- do.call(asreml, 
                           args=list(fixed = mod,
                                     random = ~ SubBlock/Block, 
                                     data = indv703.dat, maxiter=50))
    current.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                                label = "Initial model")
    current.asrt <- rmboundary(current.asrt)
    corr.asrt <- chooseSpatialModelOnIC(current.asrt, 
                                        row.covar = "cLane", col.covar = "cPosn", 
                                        row.factor = "Lane", col.factor = "Position", 
                                        trySpatial = "corr")
    ksumm <- summary(corr.asrt$asrts[[1]]$asreml.obj)$varcomp
#    print(all.equal(ksumm, summ[[kresp]], tolerance = 1e-05))
    testthat::expect_true(all.equal(ksumm, summ[[kresp]], tolerance = 1e-05))
#    summ <- c(summ, list(summary(corr.asrt$asrts[[1]]$asreml.obj)$varcomp))
  }
#  names(summ) <- responses.test
#  save(summ, file = "./data/summ703.rda")
})

cat("#### Test for wheat76 corr spatial models with asreml42\n")
test_that("Wheat76_corr_models_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data(Wheat.dat)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  #Add row and column covariates
  tmp.dat <- within(Wheat.dat, 
                    {
                      cColumn <- dae::as.numfac(Column)
                      cColumn <- cColumn  - mean(unique(cColumn))
                      cRow <- dae::as.numfac(Row)
                      cRow <- cRow - mean(unique(cRow))
                    })
  
  #Fit initial model - Row and column random
  current.asr <- do.call(asreml, 
                         list(yield ~ Rep + WithinColPairs + Variety, 
                              random = ~ Row + Column,
                              data=tmp.dat, maxit = 50))
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1720.891), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Row and Column effects")
  init.asrt <- rmboundary(init.asrt)
  
  # Try Row:ar1(Column) model
  current.asrt <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  row.factor = "Row", col.factor = "Column",
                                  corr.funcs = c("", "ar1"))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 4)
  testthat::expect_lt(abs(info$AIC - 1669.928), 0.10)
  
  # Try exp(cRow):Column model
  current.asrt <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  row.factor = "Row", col.factor = "Column",
                                  corr.funcs = c("exp", ""))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1714.379), 0.10)
  
  # Try exp(cRow):ar1(Column) model
  current.asrt <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                  row.covar = "cRow", col.covar = "cColumn",
                                  row.factor = "Row", col.factor = "Column",
                                  corr.funcs = c("exp", "ar1"))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1714.379), 0.10)
  testthat::expect_equal(names(current.asrt$asreml.obj$vparameters), 
                         c("Row", "Column", "Column:cRow", "Column:cRow!cRow!pow", "units!R"))
  
  #Compare lvr and TPP1LS models
  spatial.asrts <- list()
  spatial.asrts[["lvr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                            row.covar = "cRow", col.covar = "cColumn",
                                            row.factor = "Row", col.factor = "Column",
                                            corr.funcs = c("lvr", "lvr"))
  spatial.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS",
                                               row.covar = "cRow", col.covar = "cColumn",
                                               dropRowterm = "Row", dropColterm = "Column",
                                               degree = c(1,1), difforder = c(1,1),
                                               asreml.option = "grp")
  infoEach <- do.call(rbind, 
                      lapply(spatial.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
  
  testthat::expect_true(all.equal(infoEach$AIC, c(1714.861, 1710.282), tolerance = 1e-05))

  #Check trap for all id 
  testthat::expect_error(
    current.asrt <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                    row.covar = "cRow", col.covar = "cColumn",
                                    row.factor = "Row", col.factor = "Column",
                                    corr.funcs = c("idv", "")), 
    regexp = "Both correlation functions are id or equivalent")
  
  # Try id(Row):ar1(Column) model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      corr.funcs = c("id", "ar1"))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - 1667.54), 0.10)
  testthat::expect_equal(names(current.asrt$asreml.obj$vparameters), 
                         c("Column", "Row:Column", "Row:Column!Column!cor", "units!R"))
  tests <- current.asrt$test.summary
  testthat::expect_equal(nrow(tests), 4)
  testthat::expect_false(all(grepl("Try row", tests$terms)))
  
  # Try ar1(Row):id(Column) model
  current.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                      row.covar = "cRow", col.covar = "cColumn",
                                      row.factor = "Row", col.factor = "Column",
                                      corr.funcs = c("ar1", "id"))
  info <- infoCriteria(current.asrt$asreml.obj, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 5)
  testthat::expect_lt(abs(info$AIC - 1709.001), 0.10)
  testthat::expect_equal(names(current.asrt$asreml.obj$vparameters), 
                         c("Row", "Column", "Column:Row", "Column:Row!Row!cor", "units!R"))
  
  #Fit a correlation model with id to check spatial.IC
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr",  "TPPC"), 
                                          row.covar = "cRow", col.covar = "cColumn",
                                          row.factor = "Row", col.factor = "Column", 
                                          corr.funcs = c("id", "ar1"),
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 2)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPPCS"))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(1720.891, 1667.540, 1643.467)) < 0.10))
})  

cat("#### Test for barley03 spatial models with asreml42\n")
test_that("barely_spatial_models_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(asreml)
  library(asremlPlus)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  #This data is for the Durban 2003 barley data cited in Piepho, Boer and Williams (2022)
  data("barley.dat")
  
  #Fit initial model - Row and column random
  current.asr <- do.call(asreml, 
                         list(yield ~ rep + gen, 
                              random = ~ row + col,
                              data=barley.dat, maxit = 50))
  info <- infoCriteria(current.asr, IClikelihood = "full")
  testthat::expect_equal(info$varDF, 3)
  testthat::expect_lt(abs(info$AIC - -484.1135), 0.10)
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random row and col effects")
  init.asrt <- rmboundary(init.asrt)
  
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModel(init.asrt, spatial.model = "corr", 
                                                 row.covar = "crow", col.covar = "ccol",
                                                 row.factor = "row", col.factor = "col")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "crow", col.covar = "ccol",
                                                   dropRowterm = "row", dropColterm = "col")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "crow", col.covar = "ccol",
                                                  dropRowterm = "row", dropColterm = "col",
                                                  asreml.option = "grp")
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                   row.covar = "crow", col.covar = "ccol",
                                                   dropRowterm = "row", dropColterm = "col",
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
  infoEach <- lapply(spatialEach.asrts, function(asrt) infoCriteria(asrt$asreml.obj, 
                                                                    IClikelihood = "full"))
  (infoEach <- do.call(rbind, infoEach))
  testthat::expect_true(all.equal(infoEach$AIC, c(-641.2598, -611.8811, -616.8260, -646.7571), 
                                  tolerance = 1e-02))
  infoEach <- lapply(spatialEach.asrts, function(asrt) infoCriteria(asrt$asreml.obj, 
                                                                    IClikelihood = "REML"))
  (infoEach <- do.call(rbind, infoEach))
  testthat::expect_true(all.equal(infoEach$AIC, c(-230.4462, -191.8063, -226.5424, -230.1942), 
                                  tolerance = 1e-02))
})

cat("#### Test for nonfitting spatial models with asreml42\n")
test_that("nonfit_spatial_models_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  data("gw.dat")  
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  gw.dat <- within(gw.dat, 
                   {
                     cRow <- as.numfac(Row)
                     cRow <- cRow - mean(unique(cRow))
                     cCol <- as.numfac(Column)
                     cCol <- cCol - mean(unique(cCol))
                   })
  
  #Fit initial model
  current.asr <- do.call(asreml, 
                         args = list(y ~ Species:Substrate:Irrigation + cRow +cCol, 
                                     data = gw.dat, maxit = 50))
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Row and Column trends")
  init.asrt <- rmboundary(init.asrt)
  
  #Test for trySpatial = "none"
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = "none")
  testthat::expect_true(all(names(spatial.asrts) == 
                              c("asrts","spatial.IC","best.spatial.mod","best.spatial.IC")))
  testthat::expect_equal(names(spatial.asrts$asrts), "nonspatial")
  testthat::expect_equal(spatial.asrts$best.spatial.mod, "nonspatial")
  testthat::expect_true(abs(spatial.asrts$best.spatial.IC - 892.861) < 1e-04)
  testthat::expect_true(abs(spatial.asrts$spatial.IC$AIC - 892.861) < 1e-04)
  
  #Fit two models and return both - neither fits
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cCol",
                                          row.factor = "Row", col.factor = "Column", 
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 2)
  testthat::expect_equal(names(spatial.asrts$asrts), c("TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "TPNCSS", "TPPCS")))
  testthat::expect_true(all(abs(spatial.asrts$spatial.IC$AIC - 
                                  c(892.861, 892.861, 892.861)) < 0.10))
  
  #Fit all models and return all - corr fits
  spatial.asrts <- chooseSpatialModelOnIC(init.asrt, trySpatial = c("corr", "TPN", "TPPC"), 
                                          row.covar = "cRow", col.covar = "cCol",
                                          row.factor = "Row", col.factor = "Column", 
                                          dropRowterm = "Row", dropColterm = "Column",
                                          asreml.option = "grp", return.asrts = "all")
  testthat::expect_equal(length(spatial.asrts$asrts), 3)
  testthat::expect_equal(names(spatial.asrts$asrts), c("corr", "TPNCSS", "TPPCS"))
  testthat::expect_true(all(rownames(spatial.asrts$spatial.IC) == c("nonspatial", "corr", "TPNCSS", "TPPCS")))
  #The spline models AIC are greater than the nonspatial model and so chooseModelOnIC returns the nonspatial AIC for them
  testthat::expect_true(all(abs(na.omit(spatial.asrts$spatial.IC$AIC) - 
                                  c(892.861, 888.5976, 892.861, 892.861)) < 0.10))
  testthat::expect_equal(spatial.asrts$best.spatial.mod, "corr")
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModelOnIC(init.asrt, spatial.model = "corr", 
                                                     row.covar = "cRow", col.covar = "cCol",
                                                     row.factor = "Row", col.factor = "Column")
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(init.asrt, spatial.model = "TPN", 
                                                   row.covar = "cRow", col.covar = "cCol")
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS", 
                                                  row.covar = "cRow", col.covar = "cCol",
                                                  dropRowterm = "Row", dropColterm = "Column",
                                                  asreml.option = "grp")
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(init.asrt, spatial.model = "TPPS",
                                                   row.covar = "cRow", col.covar = "cCol",
                                                   dropRowterm = "Row", dropColterm = "Column",
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
  infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, , IClikelihood = "full")))
  testthat::expect_true(all(abs(infoEach$AIC - c(888.5976,897.4360,899.2390,895.8357)) < 0.001))
  #The spline models AIC are greater than the nonspatial model and so chooseModelOnIC returns the nonspatial AIC for them
  testthat::expect_false(all(abs(infoEach$AIC[1:3] - spatial.asrts$spatial.IC$AIC[2:4]) < 0.001))
  #testthat::expect_true(all.equal(spatial.asrts$spatial.IC[2:4,], infoEach[-4,-3], tolerance = 1e-01))
})

cat("#### Test spatial modelling for chick pea example with asreml42\n")
test_that("chickpea_spatial_mod_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  data(chkpeadat)
  tmp.dat <- within(chkpeadat, 
                    {
                      vMPosn <- as.numfac(fac.recast(Mainplot, newlevels = rep(1:11, times = 4)))
                      vMPosn <- vMPosn - mean(unique(vMPosn))
                    })
  asreml.options(design = TRUE)
  current.asr <- do.call(asreml, 
                         list(fixed = Biomass.plant ~ Smarthouse + Lines * TRT, 
                              random = ~Smarthouse:Zone/Mainplot, 
                              data = tmp.dat, maxit = 50))
  
  #Create an asrtests object, removing boundary terms
  init.asrt <- as.asrtests(current.asr, NULL, NULL, IClikelihood = "full", 
                           label = "Random Lane and Position effects")
  init.asrt <- rmboundary(init.asrt)
  
  # Try TPPS model with Mainplots and two Smarthouses
  TPPS.Main.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                            sections = "Smarthouse", 
                                            row.covar = "vLanes", col.covar = "vMPosn",
                                            dropRowterm = "Lane", dropColterm = NULL,
                                            asreml.option = "grp")
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.Main.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 4001.819)) < 0.10))
  
  # Try TPPS model with Lanes x Positions and two Smarthouses
  TPPS.LP.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                          sections = "Smarthouse", 
                                          row.covar = "vLanes", col.covar = "vPos",
                                          dropRowterm = NULL, dropColterm = NULL,
                                          asreml.option = "grp")
  
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.LP.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 3999.176)) < 0.10))
  
  # Try TPPS model with rotation for Mainplots and two Smarthouses
  TPPSRot.Main.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                               sections = "Smarthouse", 
                                               row.covar = "vLanes", col.covar = "vMPosn",
                                               dropRowterm = "Lane", dropColterm = NULL,
                                               rotateX = TRUE, ngridangles = c(3,3),
                                               asreml.option = "grp")
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPSRot.Main.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,8)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 3981.618)) < 0.10))
  theta.opt <- attr(TPPSRot.Main.grp.asrt$asreml.obj, which = "theta.opt")
  testthat::expect_true(all(theta.opt$SW == 90))
  testthat::expect_true(all(theta.opt$SE == c(30,0)))
  
  # Try TPPS model with rotation for Lanes x Positions and two Smarthouses
  TPPS.LP.grp.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                          sections = "Smarthouse", 
                                          row.covar = "vLanes", col.covar = "vPos",
                                          dropRowterm = NULL, dropColterm = NULL,
                                          rotateX = TRUE, ngridangles = c(3,3),
                                          asreml.option = "grp")
  
  info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.LP.grp.asrt$asreml.obj), 
                       IClikelihood = "full")
  testthat::expect_true(all(info$varDF == c(3,8)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 3979.302)) < 0.10))
  
  #Test makeTPPSplineMats with sections with Lanes x Positions  and mbf
  tps.mbf <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                               row.covar = "vLanes", col.covar = "vPos",
                               asreml.option = "mbf")
  testthat::expect_true(all(names(tps.mbf) == c("SW","SE")))
  testthat::expect_true(all(names(tps.mbf[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                     "BcrZ.df","dim","trace","data.plus")))
  testthat::expect_true(all(names(tps.mbf[[1]]$data.plus[,1:19]) == 
                              c("Smarthouse", "vLanes","vPos","Lane","Position","Zone",
                                "Mainplot","Subplot","Lines","TRT","Rep",
                                "X100.SW","Biomass.plant","Pods.plant","Filled.pods.plant", 
                                "Empty.pods.plant","Seed.No.plant","Seed.weight.plant","vMPosn")))
  testthat::expect_true(all(grepl("TP\\.",names(tps.mbf[[1]]$data.plus[,20:30]))))
  testthat::expect_equal(nrow(tps.mbf[[1]]$data.plus), 1056)
  testthat::expect_equal(ncol(tps.mbf[[1]]$data.plus), 30)
  testthat::expect_true(all(sapply(tps.mbf[[1]]$mbflist, function(x) grepl("SW.df", x$cov))))
  testthat::expect_equal(sapply(list(BcZSE.df,BrZSE.df,BcrZSE.df), nrow), c(22, 24, 528))
  testthat::expect_equal(sapply(list(BcZSE.df,BrZSE.df,BcrZSE.df), ncol), c(23, 25, 529))
  
  # Try TPPS model with Lanes x Positions and two Smarthouses and mbf
  TPPS.LP.mbf.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                          sections = "Smarthouse", 
                                          row.covar = "vLanes", col.covar = "vPos",
                                          dropRowterm = NULL, dropColterm = NULL,
                                          asreml.option = "mbf")
  
  print(info <- infoCriteria(list(split = init.asrt$asreml.obj, TPPS = TPPS.LP.mbf.asrt$asreml.obj), 
                       IClikelihood = "full"))
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 3999.176 )) < 0.10))
  
 
  #Test makeTPPSplineMats with sections with Lanes x MainPosns  and mbf
  tps.mbf <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                               row.covar = "vLanes", col.covar = "vMPosn",
                               asreml.option = "mbf")
  testthat::expect_true(all(names(tps.mbf) == c("SW","SE")))
  testthat::expect_true(all(names(tps.mbf[[1]]) == c("data","mbflist","BcZ.df","BrZ.df",
                                                     "BcrZ.df","dim","trace","data.plus")))
  testthat::expect_true(all(names(tps.mbf[[1]]$data.plus[,1:19]) == 
                              c("Smarthouse", "vLanes","vMPosn","Lane","Position","Zone",
                                "vPos", "Mainplot","Subplot","Lines","TRT","Rep",
                                "X100.SW","Biomass.plant","Pods.plant","Filled.pods.plant", 
                                "Empty.pods.plant","Seed.No.plant","Seed.weight.plant")))
  testthat::expect_true(all(grepl("TP\\.",names(tps.mbf[[1]]$data.plus[,20:30]))))
  testthat::expect_equal(nrow(tps.mbf[[1]]$data.plus), 1056)
  testthat::expect_equal(ncol(tps.mbf[[1]]$data.plus), 30)
  testthat::expect_equal(sapply(list(BcZSE.df,BrZSE.df,BcrZSE.df), nrow), c(11, 24, 264))
  testthat::expect_equal(sapply(list(BcZSE.df,BrZSE.df,BcrZSE.df), ncol), c(12, 25, 265))
  
  # Try TPPS model with Lanes x Positions and two Smarthouses and supplying tpps object
  TPPS.Main.mbf.asrt <- addSpatialModelOnIC(init.asrt, spatial.model = "TPPS", 
                                            sections = "Smarthouse", 
                                            row.covar = "vLanes", col.covar = "vMPosn",
                                            dropRowterm = NULL, dropColterm = NULL,
                                            asreml.option = "mbf", tpps4mbf.obj = tps.mbf)
  
  print(info <- infoCriteria(list(split = init.asrt$asreml.obj, 
                                  TPPS = TPPS.Main.mbf.asrt$asreml.obj), 
                             IClikelihood = "full"))
  testthat::expect_true(all(info$varDF == c(3,11)))
  testthat::expect_true(all(abs(info$AIC - c(4289.513, 4001.819)) < 0.10))
})

cat("#### Test hetero variances for HEB25 with asreml42\n")
test_that("HEB25_heterovar_asreml42", {
  skip_if_not_installed("asreml")
  skip_on_cran()
  library(dae)
  library(asreml)
  library(asremlPlus)
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  
  #Re-arrange and re-order the data.frame and add factors 
  data(cart.dat)
  tmp.dat <- within(cart.dat, 
                    { 
                      Smarthouse.Treat <- fac.combine(list(Smarthouse, Treatment.1))
                      Lanes <- factor(Lanes)
                      xPosition <- dae::as.numfac(Positions)
                      xPosition <- xPosition - mean(unique(xPosition))
                    })
  tmp.dat <- tmp.dat[c("Snapshot.ID.Tag", "Smarthouses", "Lanes", "Positions", 
                       "Genotype.ID", "Lines.nos", "Check", "Treatment.1", "Conditions", 
                       "Smarthouse", "Treat.Smarthouse", "Smarthouse.Treat", 
                       "Zones", "Rows", "Mainplots", "Subplots", 
                       "xLane", "xPosition", "xMainPosn", "MainCol", 
                       "Fresh.Weight", "Dry.Weight", "Number.Tillers.correct" , 
                       "Plant.Length", "ratio", "Water_Amount", "SSA", "ASA", 
                       "Caliper.Length", "Convex.Hull.Area", "Height", "SCR", 
                       "WUE.2", "agrOST", "rgrOST", "linOST", "logOST", "linm4OST", 
                       "logm4OST", "linm5OST", "logm5OST", 
                       "agrm4OST", "rgrm4OST", "agrm5OST", "rgrm5OST", 
                       "WUE100", "CHA10000", "ASA10000", "SSA10000", 
                       "ShootArea_sm", "AGR_sm_32_42", "RGR_sm_32_42", 
                       "AGR_sm_42_50", "RGR_sm_42_50", "AGR_sm_50_59", "RGR_sm_50_59")]
  names(tmp.dat)[match(c("xLane", "xPosition", "xMainPosn", "MainCol"), names(tmp.dat))] <- 
    c("cLane", "cPosition", "cMainPosn", "MainPosn")
  tmp.dat <- with(tmp.dat, tmp.dat[order(Treat.Smarthouse, Zones, Mainplots), ])
  
  
  #Fit an initial model that includes the random term us(Treatment):Genotype
  asreml.options(keep.order = TRUE) #required for asreml4 only
  HEB25.asr <- do.call(asreml, 
                       list(fixed = Dry.Weight ~ Smarthouse + Check + Treatment.1 + 
                              Check:Treatment.1, 
                            random = ~ us(Treatment.1):Genotype.ID + 
                              (at(Smarthouse, 'NW') + at(Smarthouse, 'NE')):Zones:Mainplots, 
                            residual = ~idh(Treat.Smarthouse):Zones:Mainplots, 
                            data = tmp.dat, na.action=na.method(y="include", x="include"), 
                            maxit = 100, trace = FALSE))
  summ <- summary(HEB25.asr)$varcomp
  testthat::expect_equal(nrow(summ), 10)
  testthat::expect_equal(summ$bound, c("P","P","P","P","P","F","P","P","P","P"))
  
  HEB25.idh.asrt <- as.asrtests(HEB25.asr, NULL, NULL, label = "Nonspatial model", 
                                IClikelihood = "full")
  suppressWarnings(
    testthat::expect_true(all(abs(infoCriteria(HEB25.idh.asrt$asreml.obj)[c("AIC","BIC")] - 
                                    c(539.034, 583.5741)) < 1e-03)))
  #print(HEB25.idh.asrt)
  
  #Test spatial models on Lanes x MainPosn
  #Check makeTPPSplineMats - must be ordered for Smarthouse then Treatment.1
  tpsLM.mat <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                                 row.covar = "cLane", col.covar = "cMainPosn",
                                 asreml.option = "grp")
  testthat::expect_equal(names(tpsLM.mat), c("NW", "NE"))
  testthat::expect_equal(c(nrow(tpsLM.mat$NW$data), nrow(tpsLM.mat$NE$data)), c(264, 264))
  testthat::expect_equal(c(ncol(tpsLM.mat$NW$data), ncol(tpsLM.mat$NE$data)), c(348, 348))
  testthat::expect_equal(c(nrow(tpsLM.mat$NW$data.plus), nrow(tpsLM.mat$NE$data.plus)), c(1056, 1056))
  testthat::expect_equal(c(ncol(tpsLM.mat$NW$data.plus), ncol(tpsLM.mat$NE$data.plus)), c(401, 401))
  #Check that order of dat.plus is the same as in the original tmp.dat
  testthat::expect_equal(tpsLM.mat$NW$data.plus$Snapshot.ID.Tag, tmp.dat$Snapshot.ID.Tag)
  
  #Choose best model for L x M spatial variation
  HEB25.spatialLM.asrts <- chooseSpatialModelOnIC(HEB25.idh.asrt, 
                                                  sections = "Smarthouse", 
                                                  row.covar = "cLane", col.covar = "cMainPosn",
                                                  row.factor = "Lanes", col.factor = "MainPosn",
                                                  asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25.spatialLM.asrts$spatial.IC$AIC - 
                                  c(525.5955, 488.4491, 474.0910, 473.2412, 479.7447) < 1e-03)))
  testthat::expect_equal(names(HEB25.spatialLM.asrts$asrts), c("corr",  "TPNCSS", "TPPCS",  "TPP1LS"))
  summ <- summary(HEB25.spatialLM.asrts$asrts$TPPCS$asreml.obj)$varcomp
  summ$bound[summ$bound == " "] <- "P" #hack to overcome asreml returning spaces
  testthat::expect_equal(nrow(summ), 18)
  testthat::expect_true(all((summ$bound[-14] == "P")))
  testthat::expect_true(all((summ$bound[14] == "F")))
  summ <- summary(HEB25.spatialLM.asrts$asrts$corr$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ), 13)
  testthat::expect_equal(summ$bound, c("P","U","U","P","U","P","P","P",
                                       "F","P","P","P","P"))
  
  #Check that calculated spatial.IC is the same as those for models fitted using addSpatialModel
  spatialEach.asrts <- list()
  spatialEach.asrts[["corr"]] <- addSpatialModelOnIC(HEB25.idh.asrt, spatial.model = "corr", 
                                                     sections = "Smarthouse", 
                                                     row.covar = "cLane", col.covar = "cMainPosn",
                                                     row.factor = "Lanes", col.factor = "MainPosn")
  testthat::expect_true(any(grepl("NW", spatialEach.asrts[["corr"]]$test.summary$terms)) & 
                          any(grepl("NE", spatialEach.asrts[["corr"]]$test.summary$terms)))
  spatialEach.asrts[["TPNCSS"]] <- addSpatialModel(HEB25.idh.asrt, spatial.model = "TPN", 
                                                   sections = "Smarthouse", 
                                                   row.covar = "cLane", col.covar = "cMainPosn")
  testthat::expect_true(any(grepl("NW", spatialEach.asrts[["TPNCSS"]]$test.summary$terms)) & 
                          any(grepl("NE", spatialEach.asrts[["TPNCSS"]]$test.summary$terms)))
  spatialEach.asrts[["TPPCS"]] <- addSpatialModel(HEB25.idh.asrt, spatial.model = "TPPS", 
                                                  sections = "Smarthouse", 
                                                  row.covar = "cLane", col.covar = "cMainPosn",
                                                  asreml.option = "grp")
  testthat::expect_true(any(grepl("NW", spatialEach.asrts[["TPPCS"]]$test.summary$terms)) & 
                          any(grepl("NE", spatialEach.asrts[["TPPCS"]]$test.summary$terms)))
  spatialEach.asrts[["TPP1LS"]] <- addSpatialModel(HEB25.idh.asrt, spatial.model = "TPPS",
                                                   sections = "Smarthouse", 
                                                   row.covar = "cLane", col.covar = "cMainPosn",
                                                   degree = c(1,1), difforder = c(1,1),
                                                   asreml.option = "grp")
  testthat::expect_true(any(grepl("NW", spatialEach.asrts[["TPP1LS"]]$test.summary$terms)) & 
                          any(grepl("NE", spatialEach.asrts[["TPP1LS"]]$test.summary$terms)))
  infoEach <- do.call(rbind, 
                      lapply(spatialEach.asrts, 
                             function(asrt) infoCriteria(asrt$asreml.obj, IClikelihood = "full")))
  testthat::expect_true(all.equal(HEB25.spatialLM.asrts$spatial.IC[c(2,4:5),], infoEach[c(1,3:4), -3], 
                                  tolerance = 1e-05))
  
  #Test spatial models on Lanes x Positions
  #Check makeTPPSplineMats - must be ordered for Smarthouse then Treatment.1
  tpsLP.mat <- makeTPPSplineMats(tmp.dat, sections = "Smarthouse", 
                                 row.covar = "cLane", col.covar = "cPosition",
                                 asreml.option = "grp")
  testthat::expect_equal(names(tpsLP.mat), c("NW", "NE"))
  testthat::expect_equal(c(nrow(tpsLP.mat$NW$data), nrow(tpsLP.mat$NE$data)), c(528, 528))
  testthat::expect_equal(c(ncol(tpsLP.mat$NW$data), ncol(tpsLP.mat$NE$data)), c(634, 634))
  testthat::expect_equal(c(nrow(tpsLP.mat$NW$data.plus), nrow(tpsLP.mat$NE$data.plus)), c(1056, 1056))
  testthat::expect_equal(c(ncol(tpsLP.mat$NW$data.plus), ncol(tpsLP.mat$NE$data.plus)), c(687, 687))
  
  #Choose best model for L x P spatial variation
  HEB25.spatialLP.asrts <- chooseSpatialModelOnIC(HEB25.idh.asrt,  
                                                  sections = "Smarthouse", 
                                                  row.covar = "cLane", col.covar = "cPosition",
                                                  row.factor = "Lanes", col.factor = "Positions",
                                                  asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25.spatialLP.asrts$spatial.IC$AIC - 
                                  c(525.5955, 486.4039, 471.5088, 472.8215, 476.6325) < 0.1)))
  testthat::expect_equal(names(HEB25.spatialLP.asrts$asrts), c("corr",  "TPNCSS", "TPPCS",  "TPP1LS"))
  summ <- summary(HEB25.spatialLP.asrts$asrts$TPPCS$asreml.obj)$varcomp
  summ$bound[summ$bound == " "] <- "P" #hack to overcome asreml returning spaces
  testthat::expect_equal(nrow(summ), 19)
  testthat::expect_true(all((summ$bound[-15] == "P")))
  testthat::expect_true(all((summ$bound[15] == "F")))
  summ <- summary(HEB25.spatialLP.asrts$asrts$corr$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ), 14)
  testthat::expect_equal(summ$bound, c("P","P","U","U","P","U",
                                       "P","P","P","F","P","B","P","P"))
  
  #Return two P-spline models with rotation  for L x P spatial variation
  HEB25Rot.spatialLP.asrts <- chooseSpatialModelOnIC(HEB25.idh.asrt,  trySpatial = c("TPPCS", "TPP1LS"),
                                                     sections = "Smarthouse", 
                                                     row.covar = "cLane", col.covar = "cPosition",
                                                     row.factor = "Lanes", col.factor = "Positions",
                                                     rotateX = TRUE, ngridangles = c(2,2),
                                                     asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25Rot.spatialLP.asrts$spatial.IC$AIC - 
                                  c(525.5955, 467.2578, 476.6325) < 0.1)))
  testthat::expect_equal(names(HEB25Rot.spatialLP.asrts$asrts), c("TPPCS",  "TPP1LS"))
  summ <- summary(HEB25Rot.spatialLP.asrts$asrts$TPPCS$asreml.obj)$varcomp
  summ$bound[summ$bound == " "] <- "P" #hack to overcome asreml returning spaces
  testthat::expect_equal(nrow(summ), 17)
  testthat::expect_true(all((summ$bound[-c(13)] == "P")))
  testthat::expect_true(all((summ$bound[13] == "F")))
  summ <- summary(HEB25Rot.spatialLP.asrts$asrts$TPP1LS$asreml.obj)$varcomp
  summ$bound[summ$bound == " "] <- "P" #hack to overcome asreml returning spaces
  testthat::expect_equal(nrow(summ), 17)
  testthat::expect_true(all((summ$bound[-13] == "P")))
  testthat::expect_true(all((summ$bound[13] == "F")))
  theta.opt <- attr(HEB25Rot.spatialLP.asrts$asrts$TPPCS$asreml.obj, which = "theta.opt")
  testthat::expect_true(all(theta.opt$NW == c(45,90)))
  testthat::expect_true(all(theta.opt$NE == c(90,90)))
  
  #Test dsum 
  HEB25.asr <- do.call(asreml,
                       list(fixed = Dry.Weight ~ Smarthouse + Check + Treatment.1 + 
                              Check:Treatment.1, 
                            random = ~ us(Treatment.1):Genotype.ID + 
                              (at(Smarthouse, 'NW') + at(Smarthouse, 'NE')):Zones:Mainplots, 
                            residual = ~ dsum(~ Zones:Mainplots | Treat.Smarthouse), 
                            data = tmp.dat, na.action=na.method(y="include", x="include"), 
                            maxit = 100, trace = FALSE))
  
  HEB25.ds.asrt <- as.asrtests(HEB25.asr, NULL, NULL, label = "Nonspatial model", 
                               IClikelihood = "full")
  suppressWarnings(
    testthat::expect_true(all(abs(infoCriteria(HEB25.ds.asrt$asreml.obj)[c("AIC","BIC")] - 
                                    c(539.034, 583.5741)) < 1e-03)))
  summ.idh <- summary(HEB25.idh.asrt$asreml.obj)$varcomp
  summ.ds <- summary(HEB25.ds.asrt$asreml.obj)$varcomp
  #Check that varcomp is the same for idh and dsum
  testthat::expect_true(all.equal(summ.idh[-6,], summ.ds, tolerance = 1e-03, 
                                  check.attributes = FALSE))
  #print(HEB25.ds.asrt)
  
  #Choose best model for L x M spatial variation when variance specified using dsum
  HEB25.spatialLM.ds.asrts <- chooseSpatialModelOnIC(HEB25.ds.asrt, 
                                                     sections = "Smarthouse", 
                                                     row.covar = "cLane", col.covar = "cMainPosn",
                                                     row.factor = "Lanes", col.factor = "MainPosn",
                                                     asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25.spatialLM.ds.asrts$spatial.IC$AIC - 
                                  c(525.5954, 488.4492, 474.0911, 473.2411, 479.7447) < 1e-03)))
  testthat::expect_equal(names(HEB25.spatialLM.ds.asrts$asrts), c("corr",  "TPNCSS", "TPPCS",  "TPP1LS"))
  #Check TPPCS
  summ.idh <- summary(HEB25.spatialLM.asrts$asrts$TPPCS$asreml.obj)$varcomp
  summ.idh$bound[summ.idh$bound == " "] <- "P" #hack to overcome asreml returning spaces
  summ.ds <- summary(HEB25.spatialLM.ds.asrts$asrts$TPPCS$asreml.obj)$varcomp
  summ.ds$bound[summ.ds$bound == " "] <- "P" #hack to overcome asreml returning spaces
  testthat::expect_equal(nrow(summ.ds), 17)
  testthat::expect_true(all((summ.ds$bound[-15] == "P")))
  testthat::expect_true(all((summ.ds$bound[15] == "P")))
  testthat::expect_equal(rownames(summ.idh)[1:13], rownames(summ.ds)[1:13])
  testthat::expect_true(all.equal(summ.idh[-14,-5], summ.ds[,-5], tolerance = 1e-02, 
                                  check.attributes = FALSE))
  #Check corr
  summ.idh <- summary(HEB25.spatialLM.asrts$asrts$corr$asreml.obj)$varcomp
  summ.ds <- summary(HEB25.spatialLM.ds.asrts$asrts$corr$asreml.obj)$varcomp
  testthat::expect_equal(nrow(summ.ds), 12)
  testthat::expect_equal(summ.ds$bound, c("P","U","U","P","U","P",
                                          "P","P","P","P","P","P"))
  #Two components are missing from dsum
  testthat::expect_equal(rownames(summ.idh)[c(1:8)], rownames(summ.ds)[1:8])
  # testthat::expect_true(all.equal(summ.idh[-c(6,7,11), 1:2], summ.ds[, 1:2], tolerance = 0.1, 
  #                                 check.attributes = FALSE))
  
  #Choose best model for L x M spatial variation when variance specified using dsum
  HEB25.spatialLP.ds.asrts <- chooseSpatialModelOnIC(HEB25.ds.asrt,  
                                                     sections = "Smarthouse", 
                                                     row.covar = "cLane", col.covar = "cPosition",
                                                     row.factor = "Lanes", col.factor = "Positions",
                                                     asreml.option = "grp", return.asrts = "all")
  testthat::expect_true(all(abs(HEB25.spatialLP.ds.asrts$spatial.IC$AIC - 
                                  c(525.5954, 480.4911, 471.5088, 472.8214, 476.6324) < 0.1)))
  testthat::expect_equal(names(HEB25.spatialLP.ds.asrts$asrts), c("corr",  "TPNCSS", "TPPCS",  "TPP1LS"))
  #Check TPPCS
  summ.idh <- summary(HEB25.spatialLP.asrts$asrts$TPPCS$asreml.obj)$varcomp
  summ.idh$bound[summ$bound == " "] <- "P" #hack to overcome asreml returning spaces
  summ.ds <- summary(HEB25.spatialLP.ds.asrts$asrts$TPPCS$asreml.obj)$varcomp
  summ.ds$bound[summ$bound == " "] <- "P" #hack to overcome asreml returning spaces
  testthat::expect_equal(rownames(summ.idh)[1:13], rownames(summ.ds)[1:13])
  testthat::expect_true(all.equal(summ.idh[-15,], summ.ds, tolerance = 1e-03, 
                                  check.attributes = FALSE))
  #Check corr
  summ.idh <- summary(HEB25.spatialLP.asrts$asrts$corr$asreml.obj)$varcomp
  summ.ds <- summary(HEB25.spatialLP.ds.asrts$asrts$corr$asreml.obj)$varcomp
  testthat::expect_equal(rownames(summ.idh)[c(1:4,7:9)], rownames(summ.ds)[c(1:4,8:10)])
  #idh and ds do not give equivalent answers
  #testthat::expect_true(all.equal(summ.idh[-11,"component"], summ.ds[,"component"], tolerance = 1e-03, 
  #                                 check.attributes = FALSE))
  infoAIC <- infoCriteria(list(idh = HEB25.spatialLP.asrts$asrts$corr$asreml.obj, 
                            ds = HEB25.spatialLP.ds.asrts$asrts$corr$asreml.obj))["AIC"]
  testthat::expect_true((infoAIC$AIC[1] - infoAIC$AIC[2]) >5)
})

