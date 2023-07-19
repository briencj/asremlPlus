#The following functions have been developed for single and multiple section experiments

checkTrySpatial <- function(trySpatial)
{
  trySpat.opts <- c("none", "corr", "TPNCSS", "TPPCS", "TPP1LS", "all")
  trySpatial <- trySpat.opts[unlist(lapply(trySpatial, check.arg.values, options=trySpat.opts))]
  
  if (length(intersect(trySpatial, trySpat.opts)) == 0)
    stop("trySpatial must be one of ", paste0(trySpat.opts, collapse = ", "))
  if ("all" %in% trySpatial)
    trySpatial <- c("corr", "TPNCSS", "TPPCS", "TPP1LS")
  if ("none" %in% trySpatial && length(trySpatial) > 1)
    trySpatial <= "none"
  return(trySpatial)
}

calc.nsect <- function(dat, sections)
{
  if (!is.null(sections))
  {    
    if (!is.factor(dat[[sections]]))
      stop("sections, if not NULL, must be a factor")
    else
    {
      nsect <- nlevels(dat[[sections]]) 
    }
  } else #Sections is NULL
    nsect <- 1
  
  return(nsect)
}

addSpatialModel.asrtests <- function(asrtests.obj, spatial.model = "TPPS", 
                                     sections = NULL, 
                                     row.covar = "cRow", col.covar = "cCol", 
                                     row.factor = "Row", col.factor = "Col", 
                                     corr.funcs = c("ar1", "ar1"), 
                                     row.corrFitfirst = TRUE, 
                                     dropRowterm = NULL, dropColterm = NULL, 
                                     nsegs = NULL, nestorder = c(1, 1), 
                                     degree = c(3,3), difforder = c(2,2), 
                                     asreml.option = "mbf", tpps4mbf.obj = NULL,  
                                     allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                                     checkboundaryonly = FALSE, update = FALSE, 
                                     IClikelihood = "full", ...)
{    
  #Deal with arguments for tpsmmb and changeModelOnIC
  inargs <- list(...)
  checkEllipsisArgs(c("changeModelOnIC", "asreml"), inargs)
  checkEllipsisArgs("tpsmmb", inargs, pkg = "asremlPlus")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check IClikelihood options
  options <- c("REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  
  #Check nsegs
  if (length(nsegs) > 2)
    stop("nsegs must specify no more than 2 values")
  
  #Check spatial.model options
  options <- c("corr", "TPNCSS", "TPPS")
  spatial.mod <- options[check.arg.values(spatial.model, options)]
  
  #Check asreml.option
  options <- c("mbf", "grp")
  asreml.opt <- options[check.arg.values(asreml.option, options)]
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  spatial.asrts <- list()
  
  #Fit a local spatial model involving correlated effects
  if ("corr" %in% spatial.mod)
    spatial.asrt <- fitCorrMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               row.factor = row.factor, col.factor = col.factor, 
                               corr.funcs = corr.funcs, 
                               row.corrFitfirst = row.corrFitfirst, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, chooseOnIC = FALSE, 
                               IClikelihood = ic.lik, ...)
  #Fit a local spatial model involving TPNCSS
  if ("TPNCSS" %in% spatial.mod)
    spatial.asrt <- fitTPNCSSMod(asrtests.obj, sections = sections, 
                                 row.covar = row.covar, col.covar = col.covar, 
                                 dropRowterm = dropRowterm, dropColterm = dropColterm, 
                                 allow.unconverged = allow.unconverged, 
                                 allow.fixedcorrelation = allow.fixedcorrelation,
                                 checkboundaryonly = checkboundaryonly, 
                                 update = update, chooseOnIC = FALSE, 
                                 IClikelihood = "full", ...)
  
  #Fit a residual spatial model involving TPPS
  if ("TPPS" %in% spatial.mod)
    spatial.asrt <- fitTPPSMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               dropRowterm = dropRowterm, dropColterm = dropColterm, 
                               nsegs = nsegs, nestorder = nestorder, 
                               degree = degree, difforder = difforder,
                               asreml.opt = asreml.opt, 
                               tpps4mbf.obj = tpps4mbf.obj,
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, chooseOnIC = FALSE, 
                               IClikelihood = ic.lik, ...)
  
  return(spatial.asrt)  
}

addSpatialModelOnIC.asrtests <- function(asrtests.obj, spatial.model = "TPPS", 
                                         sections = NULL, 
                                         row.covar = "cRow", col.covar = "cCol", 
                                         row.factor = "Row", col.factor = "Col", 
                                         corr.funcs = c("ar1", "ar1"), 
                                         row.corrFitfirst = TRUE, 
                                         dropRowterm = NULL, dropColterm = NULL, 
                                         nsegs = NULL, nestorder = c(1, 1), 
                                         degree = c(3,3), difforder = c(2,2), 
                                         asreml.option = "mbf", tpps4mbf.obj = NULL,  
                                         allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                                         checkboundaryonly = FALSE, update = FALSE, 
                                         IClikelihood = "full", which.IC = "AIC", 
                                         ...)
{    
  #Deal with arguments for tpsmmb and changeModelOnIC
  inargs <- list(...)
  checkEllipsisArgs(c("changeModelOnIC", "asreml"), inargs)
  checkEllipsisArgs("tpsmmb", inargs, pkg = "asremlPlus")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check nsegs
  if (length(nsegs) > 2)
    stop("nsegs must specify no more than 2 values")
  
  #Check IClikelihood options
  options <- c("REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  
  options <- c("AIC", "BIC") #, "both")
  ic.type <- options[check.arg.values(which.IC, options)]
  
  #Check spatial.model options
  options <- c("corr", "TPNCSS", "TPPS")
  spatial.mod <- options[check.arg.values(spatial.model, options)]
  
  #Check asreml.option
  options <- c("mbf", "grp")
  asreml.opt <- options[check.arg.values(asreml.option, options)]
  
  asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
  spatial.asrts <- list()
  
  #Fit a local spatial model involving correlated effects
  if ("corr" %in% spatial.mod)
    spatial.asrt <- fitCorrMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               row.factor = row.factor, col.factor = col.factor, 
                               corr.funcs = corr.funcs, 
                               row.corrFitfirst = row.corrFitfirst, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, chooseOnIC = TRUE, 
                               IClikelihood = ic.lik, which.IC = ic.type, 
                               ...)
  #Fit a local spatial model involving TPNCSS
  if ("TPNCSS" %in% spatial.mod)
    spatial.asrt <- fitTPNCSSMod(asrtests.obj, sections = sections, 
                                 row.covar = row.covar, col.covar = col.covar, 
                                 dropRowterm = dropRowterm, dropColterm = dropColterm, 
                                 allow.unconverged = allow.unconverged, 
                                 allow.fixedcorrelation = allow.fixedcorrelation,
                                 checkboundaryonly = checkboundaryonly, 
                                 update = update, chooseOnIC = TRUE, 
                                 IClikelihood = ic.lik, which.IC = ic.type, 
                                 ...)
  
  #Fit a residual spatial model involving TPPS
  if ("TPPS" %in% spatial.mod)
    spatial.asrt <- fitTPPSMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               dropRowterm = dropRowterm, dropColterm = dropColterm, 
                               nsegs = nsegs, nestorder = nestorder, 
                               degree = degree, difforder = difforder,
                               asreml.opt = asreml.opt, 
                               tpps4mbf.obj = tpps4mbf.obj,
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, chooseOnIC = TRUE, 
                               IClikelihood = ic.lik, which.IC = ic.type, 
                               ...)
  
  return(spatial.asrt)  
}

#This function calculates the IC statistics for a fitted spatial model, 
#irrespective of whether the finally fitted model was the better than the nospatial model
calcSpatialICs <- function(spatial.asrt, spatial.mod, corr.funcs = c("ar1", "ar1"), 
                           spatial.IC)
{
  tests.cols <- c("DF","denDF","AIC","BIC")
  tests <- spatial.asrt$test.summary
  if (spatial.mod == "corr")
  { 
    corr.funcs <- corr.funcs[corr.funcs != ""]
    tests.sp <- do.call(rbind,
                        lapply(c(corr.funcs), 
                               function(func, tests) 
                                 tests[grepl(paste0(" ", func), tests$terms), ], 
                               tests = tests))
    tests.sp <- unique(tests.sp)
    if (nrow(tests.sp > 0))
      tests.sp <- tests.sp[tests.sp$action == "Swapped", ]
    #Take into account a test for a nugget term so that have IC statistics for the best fitting corr model
    if (spatial.mod == "corr" && any(grepl("Try fixed residual variance", tests$terms) & 
                                     grepl("Swapped", tests$action)))
    {  
      tests.unit <- as.data.frame(tests[grepl("Try fixed residual variance", tests$terms) & 
                                          grepl("Swapped", tests$action), ])
      tests.sp <- rbind(tests.sp, tests.unit)
    }
  }
  else #get tensor spline tests 
    tests.sp <- tests[grepl("Try tensor", tests$terms) &  grepl("wapped", tests$action), ]
  
  #Form spatial.IC, calculating ICs and loglik
  if (nrow(tests.sp) < 1 || !any(grepl("wapped", tests.sp$action))) #spatial not fitted
  {
    tests.sp <- as.data.frame(matrix(rep(NA, 5), nrow = 1))
    names(tests.sp) <- c(tests.cols, "loglik")
  } else
  {
    tests.sp <- as.data.frame(tests.sp[ ,tests.cols])
    tests.sp$loglik <- 0
  }
  names(tests.sp)[1:2] <- c("fixedDF","varDF")
  tests.sp <- as.data.frame(t(as.matrix(colSums(rbind(spatial.IC[1,],tests.sp)), nrow = 1)))
  rownames(tests.sp) <- spatial.mod
  tests.sp$loglik <- with(tests.sp, -0.5 * (AIC - 2 * (fixedDF + varDF)))
  spatial.IC <- rbind(spatial.IC,tests.sp)
  return(spatial.IC)
}

chooseSpatialModelOnIC.asrtests <- function(asrtests.obj, trySpatial = "all", 
                                            sections = NULL, 
                                            row.covar = "cRow", col.covar = "cCol", 
                                            row.factor = "Row", col.factor = "Col", 
                                            corr.funcs = c("ar1", "ar1"), 
                                            row.corrFitfirst = TRUE, 
                                            dropRowterm = NULL, dropColterm = NULL, 
                                            nsegs = NULL, nestorder = c(1, 1), 
                                            asreml.option = "mbf", tpps4mbf.obj = NULL, 
                                            allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                                            checkboundaryonly = FALSE, update = FALSE, 
                                            IClikelihood = "full", which.IC = "AIC", 
                                            return.asrts = "best", ...)
{    
  #Deal with arguments for tpsmmb and changeModelOnIC
  inargs <- list(...)
  checkEllipsisArgs(c("changeModelOnIC", "asreml"), inargs)
  checkEllipsisArgs("tpsmmb", inargs, pkg = "asremlPlus")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  trySpatial <- checkTrySpatial(trySpatial)
  
  #Check nsegs
  if (length(nsegs) > 2)
    stop("nsegs must specify no more than 2 values")
  
  #Check IClikelihood options
  options <- c("REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  
  options <- c("AIC", "BIC") #, "both")
  ic.type <- options[check.arg.values(which.IC, options)]
  
  #Check spatial.model options
  options <- c("mbf", "grp")
  asreml.opt <- options[check.arg.values(asreml.option, options)]
  
  #Check return.asrts options
  options <- c("best", "all")
  return.opt <- options[check.arg.values(return.asrts, options)]
  
  if ("none" %in% trySpatial)
  {
    spatial.asrts <- list(nonspatial = asrtests.obj)
    spatial.IC <- infoCriteria(asrtests.obj$asreml.obj, IClikelihood = ic.lik)
    spatial.IC <- spatial.IC[-match("NBound", names(spatial.IC))]
    rownames(spatial.IC) <- "nonspatial"
    #Find min AIC and, if multiple mins, select in specified order
    spatial.comp <- round(spatial.IC[[which.IC]], digits = 3)
    names(spatial.comp) <- rownames(spatial.IC)
    min.asrt <- which.min(spatial.comp)
  } else #fit a spatial model
  {
    asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
    spatial.asrts <- list()
    spatial.IC <- infoCriteria(asrtests.obj$asreml.obj, IClikelihood = ic.lik)
    spatial.IC <- spatial.IC[-match("NBound", names(spatial.IC))]
    rownames(spatial.IC) <- "nonspatial"
    #Fit a local spatial model involving correlated effects
    if ("corr" %in% trySpatial)
    { 
      spatial.asrts[["corr"]] <- fitCorrMod(asrtests.obj, sections = sections, 
                                            row.covar = row.covar, col.covar = col.covar, 
                                            row.factor = row.factor, col.factor = col.factor, 
                                            corr.funcs = corr.funcs, 
                                            row.corrFitfirst = row.corrFitfirst, 
                                            allow.unconverged = allow.unconverged, 
                                            allow.fixedcorrelation = allow.fixedcorrelation,
                                            checkboundaryonly = checkboundaryonly, 
                                            update = update, chooseOnIC = TRUE, 
                                            IClikelihood = ic.lik, which.IC = ic.type, 
                                            ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["corr"]], spatial.mod = "corr", 
                                   corr.funcs = corr.funcs, spatial.IC = spatial.IC)
    }
    #Fit a local spatial model involving TPNCSS
    if ("TPNCSS" %in% trySpatial)
    { 
      spatial.asrts[["TPNCSS"]] <- fitTPNCSSMod(asrtests.obj, sections = sections, 
                                                row.covar = row.covar, col.covar = col.covar, 
                                                dropRowterm = dropRowterm, dropColterm = dropColterm, 
                                                allow.unconverged = allow.unconverged, 
                                                allow.fixedcorrelation = allow.fixedcorrelation,
                                                checkboundaryonly = checkboundaryonly, 
                                                update = update, chooseOnIC = TRUE, 
                                                IClikelihood = ic.lik, which.IC = ic.type, 
                                                ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["TPNCSS"]] , spatial.mod = "TPNCSS", 
                                   spatial.IC = spatial.IC)
    }
    
    #Fit a residual spatial model involving TPPCS
    if ("TPPCS" %in% trySpatial)
    { 
      spatial.asrts[["TPPCS"]] <- fitTPPSMod(asrtests.obj, sections = sections, 
                                             row.covar = row.covar, col.covar = col.covar, 
                                             dropRowterm = dropRowterm, dropColterm = dropColterm, 
                                             nsegs = nsegs, nestorder = nestorder, 
                                             degree = c(3,3), difforder = c(2,2), 
                                             asreml.opt = asreml.opt, 
                                             tpps4mbf.obj = tpps4mbf.obj,
                                             allow.unconverged = allow.unconverged, 
                                             allow.fixedcorrelation = allow.fixedcorrelation,
                                             checkboundaryonly = checkboundaryonly, 
                                             update = update, chooseOnIC = TRUE, 
                                             IClikelihood = ic.lik, which.IC = ic.type, 
                                             ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["TPPCS"]] , spatial.mod = "TPPCS", 
                                   spatial.IC = spatial.IC)
    }
    
    #Fit a residual spatial model involving TPP1LS
    if ("TPP1LS" %in% trySpatial)
    { 
      spatial.asrts[["TPP1LS"]] <- fitTPPSMod(asrtests.obj, sections = sections, 
                                              row.covar = row.covar, col.covar = col.covar, 
                                              dropRowterm = dropRowterm, dropColterm = dropColterm, 
                                              nsegs = nsegs, nestorder = nestorder, 
                                              degree = c(1,1), difforder = c(1,1), 
                                              asreml.opt = asreml.opt, 
                                              tpps4mbf.obj = tpps4mbf.obj,
                                              allow.unconverged = allow.unconverged, 
                                              allow.fixedcorrelation = allow.fixedcorrelation,
                                              checkboundaryonly = checkboundaryonly, 
                                              update = update, chooseOnIC = TRUE, 
                                              IClikelihood = ic.lik, which.IC = ic.type, 
                                              ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["TPP1LS"]] , spatial.mod = "TPP1LS", 
                                   spatial.IC = spatial.IC)
    }
    
    #Find min AIC and, if multiple mins, select in specified order
    spatial.comp <- round(spatial.IC[[which.IC]], digits = 3)
    names(spatial.comp) <- rownames(spatial.IC)
    min.asrt <- which.min(spatial.comp)
    if (length(min.asrt) > 1)
    {
      #pick one in the order given below
      if ("nonspatial" %in% names(min.asrt)) min.asrt <- min.asrt["nonspatial"]
      if ("TPPCS" %in% names(min.asrt)) min.asrt <- min.asrt["TPPCS"]
      if ("TPNCSS" %in% names(min.asrt)) min.asrt <- min.asrt["TPNCSS"]
      if ("TPP1LS" %in% names(min.asrt)) min.asrt <- min.asrt["TPP1LS"]
      if ("corr" %in% names(min.asrt)) min.asrt <- min.asrt["corr"]
    }
    #If return only best, get the best asrtests.obj
    if (return.opt == "best")
      spatial.asrts <- c(list(nonspatial = asrtests.obj), spatial.asrts)[names(min.asrt)]
    
  } 
  return(list(asrts = spatial.asrts, spatial.IC = spatial.IC, 
              best.spatial.mod = names(min.asrt), 
              best.spatial.IC = spatial.comp[min.asrt]))
}

#This function assumes that row.factor and col.factor are in the data
makeCorrSpec1D <- function(corr.funcs, dimension, 
                           row.covar, col.covar, row.factor, col.factor, 
                           met.funcs, unimpl.funcs)
{  
  if (corr.funcs[dimension] %in% met.funcs)
    corr1D <- paste0(corr.funcs[dimension],"(",
                     ifelse(dimension == 1, row.covar, col.covar),")")
  else
  { 
    if (corr.funcs[dimension] != "")
      corr1D <- paste0(corr.funcs[dimension],"(",
                       ifelse(dimension == 1, row.factor, col.factor),")")
    else
      corr1D <- ifelse(dimension == 1, row.factor, col.factor)
  }
  return(corr1D)
}  

fitCorrMod <- function(asrtests.obj, sections = NULL,
                       row.covar = "cRow", col.covar = "cCol", 
                       row.factor = "Row", col.factor = "Col", 
                       corr.funcs = c("ar1", "ar1"), 
                       row.corrFitfirst = TRUE, 
                       allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                       checkboundaryonly = FALSE, update = TRUE, 
                       chooseOnIC = TRUE, 
                       IClikelihood = "full", which.IC = "AIC", 
                       ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  if (!asr4)
    stop(paste("Fitting spatial models using correlation/variance models", 
               "has not been implemented for asreml version less than 4.0"))
  
  #Check that named columns are in the data
  dat.in <- asrtests.obj$asreml.obj$call$data
  if (is.symbol(dat.in))
    dat.in <- eval(dat.in)
  
  #Check the correlation functions and set them up
  met.funcs <- c("exp", "gau", "lvr")
  met.funcs <- c(met.funcs, sapply(met.funcs, function(f) paste0(f, c("v","h"))))
  unimpl.funcs <- c("iexp", "igau", "ieuc", "sph", "cir", "aexp", "agau", "mtrn")
  unimpl.funcs <- c(unimpl.funcs, sapply(unimpl.funcs, function(f) paste0(f, c("v","h"))))
  if (any(unimpl.funcs %in% corr.funcs))
    stop("Some of the following corr.funcs ar not implemented for spatial modelling: ", 
         paste(unimpl.funcs[unimpl.funcs %in% corr.funcs], collapse = ","))
  if (all(corr.funcs %in% c("", "id", "idv")))
    stop("Both correlation functions are id or equivalent")

  #Check the grid covars and factors
  #row.factor and col.factor must be in dat.in so that each dimension can be fitted independently
  #(so can have a term without a corr function when fitting a dimension)
  grid.cols <- c(sections, row.factor, col.factor)
  if (any(corr.funcs %in% met.funcs))
    grid.cols <- c(grid.cols, row.covar, col.covar)
  checkNamesInData(c(sections, row.factor, col.factor), dat.in)
  if (!all(sapply(dat.in[c(row.factor, col.factor)], is.factor)))
    stop("Both row.factor and col.factor must be factors in the data stored in the asreml.obj")
  
  nsect <- calc.nsect(dat.in, sections)
  
  #Get row and col corr models
  row.corr <- makeCorrSpec1D(corr.funcs = corr.funcs, dimension = 1, 
                             row.covar = row.covar, col.covar = col.covar, 
                             row.factor = row.factor, col.factor = col.factor, 
                             met.funcs = met.funcs, unimpl.funcs = unimpl.funcs)
  col.corr <- makeCorrSpec1D(corr.funcs = corr.funcs, dimension = 2, 
                             row.covar = row.covar, col.covar = col.covar, 
                             row.factor = row.factor, col.factor = col.factor, 
                             met.funcs = met.funcs, unimpl.funcs = unimpl.funcs)
  
  #Remove units if in model
  if (grepl("units", as.character(getFormulae(asrtests.obj$asreml.obj)$random)[2]))
    asrtests.obj <- changeTerms(asrtests.obj, 
                                dropRandom = "units", label = "Remove units term", 
                                allow.unconverged = allow.unconverged, 
                                allow.fixedcorrelation = allow.fixedcorrelation,
                                checkboundaryonly = TRUE, 
                                update = update, 
                                IClikelihood = IClikelihood, ...)
  
  #Check if correlations already included in a term 
  t <- mapply(function(func, corr, asr)
  { 
    if (!(func %in% c("", "id", "idv")) && 
        any(sapply(as.character(getFormulae(asr)), 
                   function(mod, corr) grepl(corr, mod, fixed = TRUE), 
                   corr = corr)))
      warning("The correlation function ", corr, " is already in the model")
    invisible()
  }, func = corr.funcs, corr = c(row.corr, col.corr), 
  MoreArgs = (list(asr = asrtests.obj$asreml.obj)))
  
  #Prepare to fit
  facs <- c(row.factor, col.factor)
  rfuncs <- corr.funcs
  rterms <- c(row.corr, col.corr)
  
  if (!row.corrFitfirst)
  {
    facs <- facs[c(2,1)]
    rfuncs <- rfuncs[c(2,1)]
    rterms <- rterms[c(2,1)]
  }
  
  #Loop over the sections
  corr.asrt <- asrtests.obj
  for (i in 1:nsect)
  {
    corr.term <- FALSE
    #Check have a corr func
    if (any(rfuncs[1] == c("", "id", "idv")))
      result1 <- "Unswapped"
    else
    { 
      #Try first correl in current section
      ran.term1 <- paste0(rterms[1], ":", facs[2])
      lab1 <- paste0("Try ", rterms[1])
      if (nsect > 1)
      {  
        ran.term1 <- paste0("at(", sections, ", ",i, "):", ran.term1)
        lab1 <- paste0(lab1, " for ", sections, " ",i)
      }
      if (chooseOnIC)
        corr.asrt <- changeModelOnIC(corr.asrt, 
                                     addRandom = ran.term1, label = lab1, 
                                     allow.unconverged = allow.unconverged, 
                                     allow.fixedcorrelation = allow.fixedcorrelation,
                                     checkboundaryonly = TRUE, 
                                     update = update, 
                                     IClikelihood = IClikelihood, 
                                     which.IC = which.IC, 
                                     ...)
      else
        corr.asrt <- changeTerms(corr.asrt, 
                                 addRandom = ran.term1, label = lab1, 
                                 allow.unconverged = allow.unconverged, 
                                 allow.fixedcorrelation = allow.fixedcorrelation,
                                 checkboundaryonly = TRUE, 
                                 update = update, 
                                 IClikelihood = IClikelihood, ...)
      
      if (largeVparChange(corr.asrt$asreml.obj, 0.75))
        corr.asrt <- iterate(corr.asrt)
      result1 <- getTestEntry(corr.asrt, label = lab1)$action
    }
    
    #Try 2nd correl in current section
    if (!any(rfuncs[2] == c("", "id", "idv")))
    {  
      lab <- paste0("Try ", rterms[2])
      if (nsect > 1)
        lab <- paste0(lab, " for ", sections, " ",i)
      if (!grepl("Unswapped", result1) && !grepl("Unchanged", result1)) #first fac ar1 fitted
      { 
        corr.term <- TRUE
        last.term <- ran.term1
        #Check for ran.term1 in random formula and if absent check for different order
        last.term <- chk4TermInFormula(corr.asrt$asreml.obj$call$random, term = last.term, 
                                       asreml.obj = corr.asrt$asreml.obj)
        ran.term <- paste0(rterms[1], ":", rterms[2])
        if (nsect > 1)
          ran.term <- paste0("at(", sections, ", ",i, "):", ran.term)
        
        if (chooseOnIC)
          corr.asrt <- changeModelOnIC(corr.asrt, 
                                       addRandom = ran.term, 
                                       dropRandom = last.term, label = lab, 
                                       allow.unconverged = allow.unconverged, 
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = TRUE, 
                                       update = update, 
                                       IClikelihood = IClikelihood, 
                                       which.IC = which.IC, 
                                       ...)
        else
          corr.asrt <- changeTerms(corr.asrt, 
                                   addRandom = ran.term, 
                                   dropRandom = last.term, label = lab, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = TRUE, 
                                   update = update, 
                                   IClikelihood = IClikelihood, ...)
        
        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
        if (!(grepl("Unswapped", getTestEntry(corr.asrt, label = lab)$action)) && 
            !(grepl("Unchanged", getTestEntry(corr.asrt, label = lab)$action)))
          last.term <- ran.term
      } else #no first fac ar1
      { 
        ran.term <- paste0(facs[1], ":", rterms[2])
        if (nsect > 1)
          ran.term <- paste0("at(", sections, ", ",i, "):", ran.term)
        if (chooseOnIC)
          corr.asrt <- changeModelOnIC(corr.asrt, 
                                       addRandom = ran.term, label = lab, 
                                       allow.unconverged = allow.unconverged, 
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = TRUE, 
                                       update = update, 
                                       IClikelihood = IClikelihood, 
                                       which.IC = which.IC, 
                                       ...)
        else
          corr.asrt <- changeTerms(corr.asrt, 
                                   addRandom = ran.term, label = lab, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = TRUE, 
                                   update = update, 
                                   IClikelihood = IClikelihood, ...)

        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
        result <- getTestEntry(corr.asrt, label = lab)$action
        if (!grepl("Unswapped", result) && !grepl("Unchanged", result))
        { 
          corr.term <- TRUE
          last.term <- ran.term
          #Check for ran.term1 in random formula and if absent check for different order
          last.term <- chk4TermInFormula(corr.asrt$asreml.obj$call$random, term = last.term, 
                                         asreml.obj = corr.asrt$asreml.obj)
          
          #Check for unchanged because unconverged or fixed correlation for first correlation fitted 
          #- if found, retry first correlation
          if (any(sapply(c("Unchanged"), #could add other actions for when the first fit needs retrying
                         function(act, result1) grepl(act, result1),
                         result1 = result1)))
          {
            ran.term1 <- paste0(rterms[1], ":", rterms[2])
            if (nsect > 1)
              ran.term1 <- paste0("at(", sections, ", ",i, "):", ran.term1)
            if (chooseOnIC)
              corr.asrt <- changeModelOnIC(corr.asrt, 
                                           dropRandom = last.term,
                                           addRandom = ran.term1, label = lab1, 
                                           allow.unconverged = allow.unconverged, 
                                           allow.fixedcorrelation = allow.fixedcorrelation,
                                           checkboundaryonly = TRUE, 
                                           update = update, 
                                           IClikelihood = IClikelihood, 
                                           which.IC = which.IC, 
                                           ...)
            else
              corr.asrt <- changeTerms(corr.asrt, 
                                       dropRandom = last.term,
                                       addRandom = ran.term, label = lab1, 
                                       allow.unconverged = allow.unconverged, 
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = TRUE, 
                                       update = update, 
                                       IClikelihood = IClikelihood, ...)
            
            if (largeVparChange(corr.asrt$asreml.obj, 0.75))
              corr.asrt <- iterate(corr.asrt)
            result1 <- getTestEntry(corr.asrt, label = lab1)$action
          }
        }
      }
    }

    #Test the nugget residual variance
    if (chooseOnIC && corr.term)
    {
      if (largeVparChange(corr.asrt$asreml.obj, 0.75))
        corr.asrt <- iterate(corr.asrt)
      if (asr4)
      {
        if (asr4.2)
        { 
          vpc <- corr.asrt$asreml.obj$vparameters.con
          vpt <- corr.asrt$asreml.obj$vparameters.type
        } else
        {
          vpc <- vpc.char(corr.asrt$asreml.obj)
          vpt <- vpt.char(corr.asrt$asreml.obj)
        }
        vpc <- vpc[grepl("!R$", names(vpc))]
        vpt <- vpt[names(vpc)]
      }
      if (length(vpc) == 1 && !(vpc %in% c("F", "B") && vpt == "V"))
      {
        lab <- "Try fixed residual variance"
        resvar <- names(vpc)
        resmod <- as.character(getFormulae(corr.asrt$asreml.obj)$residual)
        if (length(resmod) == 0)
          resmod <- NULL
        else
          resmod <- resmod[2]

        if (chooseOnIC)
          corr.asrt <- changeModelOnIC(corr.asrt, 
                                       newResidual = resmod, label = lab, 
                                       set.terms = resvar, initial.values = 1, 
                                       bounds = "F", ignore.suffices = FALSE, 
                                       allow.unconverged = allow.unconverged, 
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = TRUE, 
                                       update = update, 
                                       IClikelihood = IClikelihood, 
                                       which.IC = which.IC, 
                                       ...)
        else
          corr.asrt <- changeTerms(corr.asrt, 
                                   newResidual = resmod, label = lab, 
                                   set.terms = resvar, initial.values = 1, 
                                   bounds = "F", ignore.suffices = FALSE, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = TRUE, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC, 
                                   ...)
        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
      }
    }
  }

  #Having made all model changes with checkboundaryonly = TRUE, update for checkboundaryonly set to FALSE
  if (!checkboundaryonly)
    corr.asrt <- rmboundary(corr.asrt, checkboundaryonly = checkboundaryonly, update = TRUE)
  
  return(corr.asrt)
}

fitTPNCSSMod <- function(asrtests.obj, sections = NULL, 
                         row.covar = "cRow", col.covar = "cCol", 
                         dropRowterm = NULL, dropColterm = NULL, 
                         nsegs = NULL, 
                         allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                         checkboundaryonly = FALSE, update = TRUE, 
                         chooseOnIC = TRUE, 
                         IClikelihood = "full", which.IC = "AIC", 
                         ...)
{ 
  #Check that named columns are in the data
  dat.in <- asrtests.obj$asreml.obj$call$data
  if (is.symbol(dat.in))
    dat.in <- eval(dat.in)
  checkNamesInData(c(sections, row.covar, col.covar, dropRowterm, dropColterm), dat.in)
  
  #Check conformability of covars and factors
  if (!is.null(dropRowterm) && nlevels(dat.in[[dropRowterm]]) != length(unique(dat.in[[row.covar]])))
    stop(dropRowterm, " does not have the same number of levels as there are values of ", row.covar)
  if (!is.null(dropColterm) && nlevels(dat.in[[dropColterm]]) != length(unique(dat.in[[col.covar]])))
    stop(dropColterm, " does not have the same number of levels as there are values of ", col.covar)
  
  #Are dropRowterm and dropColterm already in the model?
  facs <- c(dropRowterm, dropColterm)
  drop.fix <- NULL
  if (any(facs %in% rownames(asrtests.obj$wald.tab)))
    drop.fix <- paste(facs[facs %in% rownames(asrtests.obj$wald.tab)], collapse = " + ")
  drop.ran <- NULL
  if (any(facs %in% names(asrtests.obj$asreml.obj$vparameters)))
    drop.ran <- paste(facs[facs %in% names(asrtests.obj$asreml.obj$vparameters)], 
                      collapse = " + ")
  
  nsect <- calc.nsect(dat.in, sections)
  
  #spatial using tensor NCS splines
  if (nsect == 1)
  { 
    sect.fac <- NULL
    fix.terms <- paste0(row.covar, " + ", col.covar, " + ", row.covar, ":", col.covar)
  }
  else
  { 
    
    sect.fac <- paste0("at(", sections, "):")
    fix.terms <- paste(paste0(sect.fac, c(row.covar, col.covar)), collapse = " + ")
  }
  
  #Construct random terms without sections
  spl.row <- paste0("spl(", row.covar, ")")
  spl.col <- paste0("spl(", col.covar, ")")
  spl.terms <- c(spl.row, spl.col, 
                 paste0("dev(", row.covar, ")"), paste0("dev(", col.covar, ")"), 
                 paste0(spl.row, ":", col.covar), paste0(row.covar, ":", spl.col), 
                 paste0(spl.row, ":", spl.col))
  
  for (i in 1:nsect)
  {
    if (nsect > 1) 
    { 
      sect.fac <- paste0("at(", sections, ", ", i, "):")
      lab <- paste0("Try tensor NCS splines for ", sections, " ",i)
    } else
      lab <- paste0("Try tensor NCS splines")
    #Fit TPNCSS to a section 
    if (chooseOnIC)
      tspl.asrt <- changeModelOnIC(asrtests.obj, 
                                   addFixed = fix.terms, 
                                   dropFixed = drop.fix, 
                                   dropRandom = drop.ran, 
                                   addRandom = paste(sect.fac, spl.terms, 
                                                     collapse = " + "), 
                                   label = lab, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = checkboundaryonly, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC, 
                                   ...)
    else
      tspl.asrt <- changeTerms(asrtests.obj, 
                               addFixed = fix.terms, 
                               dropFixed = drop.fix, 
                               dropRandom = drop.ran, 
                               addRandom = paste(sect.fac, spl.terms, 
                                                 collapse = " + "), 
                               label = lab, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, 
                               IClikelihood = IClikelihood, ...)
    
    tspl.asrt <- rmboundary(tspl.asrt)
  }
  
  return(tspl.asrt)
}

#Creates a list with the tpps bits for asreml.options = "grp"
addPSdesign.mat <- function(dat, sections = NULL, nsect = 1, 
                            row.coords, col.coords, 
                            nsegs = NULL, nestorder = c(1, 1), 
                            degree = c(3,3), difforder = c(2,2), 
                            asreml.opt = "grp", stub = "xx", 
                            ...)
{
  if (nsect != 1)
  {  
    tmp <- split(dat, f = dat[[sections]])
    sect.levs <- levels(dat[[sections]])
    tps.XZmat  <- lapply(tmp, 
                         function(data, columncoordinates, rowcoordinates, 
                                  sects, nsegments, asreml.opt)
                         { 
                           if (all(sapply(nsegments, is.null)))
                             nsegments <- c(length(unique(data[[columncoordinates]]))-1,
                                            length(unique(data[[rowcoordinates]]))-1)
                           #stub <- data[[sects]][1]
                           XZ.mat <- tpsmmb(columncoordinates = columncoordinates, 
                                            rowcoordinates = rowcoordinates, 
                                            data = data, 
                                            stub = stub, nsegments = nsegments, 
                                            nestorder = nestorder, 
                                            degree = degree, difforder = difforder,
                                            asreml = asreml.opt, 
                                            ...)
                           return(XZ.mat)
                         }, columncoordinates = col.coords, rowcoordinates = row.coords, 
                         sects = sect.levs, nsegments = nsegs, 
                         asreml.opt = asreml.opt)
  } else
  {
    
    if (all(sapply(nsegs, is.null)))
      nsegs <- c(length(unique(dat[[col.coords]]))-1,
                 length(unique(dat[[row.coords]]))-1)
    tps.XZmat <- list(tpsmmb(columncoordinates = col.coords, 
                             rowcoordinates = row.coords, 
                             data = dat, 
                             stub = stub, nsegments = nsegs, 
                             nestorder = nestorder, 
                             degree = degree, difforder = difforder,
                             asreml = asreml.opt, 
                             ...))
  }
  
  return(tps.XZmat)
}

#makes sure that the data components of tps.XZmat are conformable between sections for the grp asreml option 
conformTPSSections <- function(tps.XZmat)
{
  #Check if data is conformable
  dat.cols <- sapply(tps.XZmat, function(mat) ncol(mat$data))
  if (!all(dat.cols == dat.cols[1]))
  {
    #Check same group names
    grp.names <- lapply(tps.XZmat, function(mat) names(mat$grp))
    if (!all(sapply(grp.names, function(nam, grp1) all(nam == grp1), grp1 = grp.names[[1]])))
      stop("There are not the same groups for each section")
    if (!all(sapply(tps.XZmat, function(mat) all(diff(mat$grp[["all"]]) == 1))))
      stop("Not all sections occupy sequential columns")
    if (!all(diff(sapply(tps.XZmat, function(mat) mat$grp[[1]][1])) == 0))
      stop("Not all sections start in the groups in the same column")
    grp.names <- grp.names[[1]]
    
    #Find start and max lengths
    start <- tps.XZmat$grp[[1]][1]
    new.grps.lens <- sapply(grp.names, 
                            function(agrp, tps.XZmat)
                              max.len <- max(sapply(tps.XZmat, function(mat) length(mat$grp[[agrp]]))),
                            tps.XZmat = tps.XZmat)
    ngrps <- length(new.grps.lens)
    if (sum(new.grps.lens[1:(ngrps-1)]) != new.grps.lens[ngrps])
      stop("THe number of columns in the all group does not match the sum of the other groups")
    starts <- c(0,cumsum(new.grps.lens[1:(ngrps-2)])) + start
    ends <- cumsum(new.grps.lens)[1:(ngrps-1)] + start - 1
    names(starts) <- names(starts) <- grp.names[1:(ngrps-1)]
    new.col.names <- unlist(mapply(function(grp, len) {paste(grp, 1:len, sep = "_")}, 
                                   grp = grp.names[-ngrps], len = new.grps.lens[-ngrps]))
    new.mat <- lapply(tps.XZmat, function(mat, grp.names, starts, ends, new.col.names)
    { 
      ngrps = length(starts)
      new.grps.lens <- ends - starts + 1
      olddat <- mat$data
      oldgrps <- mat$grp
      newdat <- olddat[,1:(starts[1]-1)]
      tpsbdat <- as.data.frame(matrix(0, nrow = nrow(newdat), ncol = (ends[ngrps] - starts[1] + 1)))
      names(tpsbdat) <- new.col.names
      newdat <- cbind(newdat, tpsbdat)
      mat$grp <- mapply(function(grp, start) #change grps
      {
        grp.len <- length(grp)
        grp <- start + 0:(grp.len-1)
        return(grp)
        
      }, grp = mat$grp[-length(mat$grp)], start = starts, SIMPLIFY = FALSE)
      #Copy across columns from old to new dat
      for (nam in grp.names)
        newdat[,mat$grp[[nam]]] <- olddat[,oldgrps[[nam]]]
      mat$data <- newdat
      mat$grp$All <- unlist(mat$grp)
      names(mat$grp$All) <- NULL
      return (mat)
    }, grp.names = grp.names[-length(grp.names)], starts = starts, ends = ends, new.col.names = new.col.names)
  } else #data already conformable
    new.mat <- tps.XZmat
  return(new.mat)
}

#'### Function to create spline basis matrices and data for TPS splines
makeTPPSplineMats.data.frame <- function(data, sections = NULL, 
                                           row.covar, col.covar, 
                                           nsegs = NULL, nestorder = c(1, 1), 
                                           degree = c(3,3), difforder = c(2,2),
                                           asreml.option = "mbf", 
                                           ...)
{
  #Check that named columns are in the data
  checkNamesInData(c(sections, row.covar, col.covar), data)
  
  #Make sure that row covar is not centred so that tpsmmb does not create a faulty index
  row.min <- min(data[[row.covar]], na.rm = TRUE)
  if (row.min < 0)
    data[row.covar] <- data[[row.covar]] - row.min + 1
  
  col.min <- min(data[[col.covar]], na.rm = TRUE)
  if (col.min < 0)
    data[col.covar] <- data[[col.covar]] - col.min + 1
  
  #Deal with arguments for tpsmmb
  inargs <- list(...)
  checkEllipsisArgs("tpsmmb", inargs, pkg = "asremlPlus")
  
  #Check nsegs, nestorder, degree and difforder
  if (length(nsegs) > 2)
    stop("nsegs must specify no more than 2 values, one for each of the row and column dimensions")
  if (length(nestorder) != 2)
    stop("nestorder must specify 2 values, one for each of the row and column dimensions")
  if (length(degree) != 2)
    stop("degree must specify 2 values, one for each of the row and column dimensions")
  if (length(difforder) != 2)
    stop("difforder must specify 2 values, one for each of the row and column dimensions")
  
  #Check asreml.option
  options <- c("mbf","grp")
  asreml.opt <- options[check.arg.values(asreml.option, options)]
  
  #Remove any previous Tensor Spline basis columns from the data
  if (any(grepl("TP\\.", names(data))) || any(grepl("TP\\_", names(data))))
    data <- data[,-c(grep("TP\\.", names(data)), grep("TP\\_", names(data)))]
  
  #Spatial local spatial model using tensor P-splines with the grp option 
  rc.cols <- c(sections, row.covar, col.covar)
  dat.rc <- data[rc.cols]
  dat.rc <- unique(dat.rc)
  nsect <- calc.nsect(data, sections)

  #Get data for mbf 
  if (asreml.opt == "mbf")
  {
    mbf.env <- parent.frame(1)
    #stop('Sorry, but the mbf setting of asreml.opt is not functioning yet - use asreml.opt = "grp".')
    if (nsect != 1)
    {  
      tmp <- split(dat.rc, f = dat.rc[[sections]])
      sect.levs <- levels(dat.rc[[sections]])
      tps.XZmat  <- lapply(tmp, 
                           function(dat, columncoordinates, rowcoordinates, 
                                    sects, nsegments, nestorder, 
                                    degree, difforder, asreml.opt)
                           { 
                             stub <- as.character(dat[[sects]][1])
                             dat <- unique(dat)
                             rownames(dat) <- NULL
                             XZ.mat <- tpsmmb(columncoordinates = columncoordinates, 
                                              rowcoordinates = rowcoordinates, 
                                              data = dat, 
                                              stub = stub, 
                                              nsegments = nsegments, 
                                              nestorder = nestorder, 
                                              degree = degree, difforder = difforder,
                                              asreml = asreml.opt, 
                                              ...)
                             return(XZ.mat)
                           }, columncoordinates = col.covar, rowcoordinates = row.covar, 
                           sects = sections, nsegments = nsegs, nestorder = nestorder, 
                           degree = degree, difforder = difforder, asreml.opt = asreml.opt)
      #Put the mbf data.frames into the parent frame
      sect.levs <- levels(dat.rc[[sections]])
      for (sect.lev in sect.levs)
      {
        Zmat.names <- paste0(paste0(c("BcZ", "BrZ", "BcrZ"),sect.lev), ".df")
        assign(Zmat.names[1], tps.XZmat[[sect.lev]]$BcZ.df, envir = mbf.env)
        assign(Zmat.names[2], tps.XZmat[[sect.lev]]$BrZ.df, envir = mbf.env)
        assign(Zmat.names[3], tps.XZmat[[sect.lev]]$BcrZ.df, envir = mbf.env)
      }        
    } else
    {
      stub = "xx"
      tps.XZmat <- list(tpsmmb(columncoordinates = col.covar, 
                               rowcoordinates = row.covar, 
                               data = dat.rc, 
                               stub = stub, nsegments = nsegs, 
                               nestorder = nestorder, 
                               degree = degree, difforder = difforder,
                               asreml = asreml.opt, 
                               ...))
      Zmat.names <- paste0(paste0(c("BcZ", "BrZ", "BcrZ"),stub), ".df")
      if (any(sapply(Zmat.names, exists, envir = parent.frame(1))))
        warning("THe following objects are being overwritten: ", 
                paste(Zmat.names[sapply(Zmat.names, exists, envir = parent.frame(2))], 
                      collapse = ", "))
      tps.XZmat[[1]]$BcZ.df$TP.col <- factor(tps.XZmat[[1]]$BcZ.df$TP.col,
                                             levels = tps.XZmat[[1]]$BcZ.df$TP.col)
      tps.XZmat[[1]]$BrZ.df$TP.row <- factor(tps.XZmat[[1]]$BrZ.df$TP.row,
                                             levels = tps.XZmat[[1]]$BrZ.df$TP.row)
      tps.XZmat[[1]]$BcrZ.df$TP.CxR <- factor(tps.XZmat[[1]]$BcrZ.df$TP.CxR,
                                              levels = tps.XZmat[[1]]$BcrZ.df$TP.CxR)
      assign(Zmat.names[1], tps.XZmat[[1]]$BcZ.df, envir = mbf.env)
      assign(Zmat.names[2], tps.XZmat[[1]]$BrZ.df, envir = mbf.env)
      assign(Zmat.names[3], tps.XZmat[[1]]$BcrZ.df, envir = mbf.env)
    }
    
    #Build the data.frame to be used in the analysis
    kextra <- ncol(data) - length(rc.cols)
    if (nsect == 1)
      dat <- tps.XZmat[[1]]$data
    else
    { 
      dat <- lapply(tps.XZmat, function(mat) mat$data)
      dat <- do.call(rbind, dat)
    }
    dat$TP.col <- factor(dat$TP.col)
    dat$TP.row <- factor(dat$TP.row)
    dat$TP.CxR <- factor(dat$TP.CxR)
    dat <- merge(data, dat, all = TRUE, sort = FALSE)
    dat <- dat[c(rc.cols, setdiff(names(dat), rc.cols))]
    attr(dat, which = "mbf.env") <- parent.frame(1)
  } else #doing grp
  {  
    tps.XZmat <- addPSdesign.mat(dat.rc, sections = sections, nsect = nsect, 
                                 row.coords = row.covar, col.coords = col.covar, 
                                 nsegs = nsegs, nestorder = nestorder, 
                                 degree = degree, difforder = difforder,
                                 asreml.opt = "grp", 
                                 stub = "xx", ...)
    
    #Build the data.frame to be used in the analysis
    kextra <- ncol(data) - length(rc.cols)
    if (nsect == 1)
      dat <- tps.XZmat[[1]]$data
    else
    { 
      #Check, and if necessary, make data from different sections conformable
      tps.XZmat <- conformTPSSections(tps.XZmat)
      dat <- lapply(tps.XZmat, function(mat) mat$data)
      dat <- do.call(rbind, dat)
    }
    if (nrow(dat.rc) < nrow(data))
    { 
      tmp <- within(data, .row.ord <- 1:nrow(data))
      dat <- merge(tmp, dat, all.x = TRUE, by = rc.cols, sort = FALSE)
      dat <- dat[order(dat$.row.ord), ]
      dat <- dat[, -match(".row.ord", names(dat))]
    } else    
      dat <- merge(data, dat, all.x = TRUE, by = rc.cols, sort = FALSE)
    
    #Adjust the grp columns for merging
    if (nsect == 1)
    { 
      tps.XZmat[[1]]$grp <- lapply(tps.XZmat[[1]]$grp, 
                                   function(grp, kextra) grp <- grp + kextra, 
                                   kextra = kextra)
    }
    else
      tps.XZmat <- lapply(tps.XZmat, 
                          function(mat, kextra) 
                          {
                            mat$grp <- lapply(mat$grp, 
                                              function(grp, kextra) grp <- grp + kextra, 
                                              kextra = kextra)
                            return(mat)
                          }, kextra = kextra)
  }
  
  #Add to returned object
  tps.XZmat <- lapply(tps.XZmat, function(mat, data = dat) c(mat, list(data.plus = dat)))
  
  return(tps.XZmat)
}

fitTPSModSect <- function(tspl.asrt, mat, sect.fac, dropRowterm, dropColterm, 
                          row.covar, col.covar, lab, 
                          degree = c(3,3), difforder = c(2,2), 
                          asreml.opt = "mbf", stub = stub, 
                          allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                          chooseOnIC = TRUE, 
                          checkboundaryonly = FALSE, update = TRUE, 
                          IClikelihood = "full", which.IC = "AIC", ...)
{
  inargs <- list(...)
  
  #Are dropRowterm and dropColterm already in the model?
  facs <- c(dropRowterm, dropColterm)
  drop.fix <- NULL
  drop.ran <- NULL
  if (!is.null(facs))
  {
    if (any(facs %in% rownames(tspl.asrt$wald.tab)))
      drop.fix <- paste(facs[facs %in% rownames(tspl.asrt$wald.tab)], collapse = " + ")
    if (any(facs %in% names(tspl.asrt$asreml.obj$vparameters)))
      drop.ran <- paste(facs[facs %in% names(tspl.asrt$asreml.obj$vparameters)], 
                        collapse = " + ")
  }

  nfixterms <- difforder[1] * difforder[2] 
  if (nfixterms > 1)
    fix.ch <- paste(paste0(sect.fac, paste0("TP.CR.", 2:nfixterms)), collapse = " + ")
  else
    fix.ch <- NULL

  if (chooseOnIC)
    fitfunc <- changeModelOnIC
  else
    fitfunc <- changeTerms

  if (asreml.opt == "mbf")
  {
    mbf.lis <- mat$mbflist
    ran.ch <- paste(paste0(sect.fac,  
                           c(paste0("TP.C.",1:difforder[1],":mbf(TP.row)"),  
                             paste0("TP.R.",1:difforder[2],":mbf(TP.col)"), 
                             "mbf(TP.CxR)", 
                             paste0("dev(",row.covar,")"), 
                             paste0("dev(",col.covar,")"))), 
                    collapse = " + ")
    tspl.asrt <- do.call(fitfunc, 
                         args = c(list(tspl.asrt,
                                       addFixed = fix.ch,
                                       dropFixed = drop.fix, 
                                       addRandom = ran.ch,
                                       dropRandom = drop.ran, 
                                       mbf = mbf.lis,
                                       label = lab,
                                       allow.unconverged = allow.unconverged, 
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = checkboundaryonly, 
                                       update = update, 
                                       IClikelihood = IClikelihood, 
                                       which.IC = which.IC), 
                                  inargs))
  } else #grp
  {
    grp <- mat$grp
    ran.ch <- paste(paste0(sect.fac,  
                           c(paste0("grp(TP.C.",1:difforder[1],"_frow)"),  
                             paste0("grp(TP.R.",1:difforder[2],"_fcol)"), 
                             "grp(TP_fcol_frow)", 
                             paste0("dev(",row.covar,")"), 
                             paste0("dev(",col.covar,")"))), 
                    collapse = " + ")
    tspl.asrt <- do.call(fitfunc, 
                         args = c(list(tspl.asrt, 
                                       addFixed = fix.ch,
                                       dropFixed = drop.fix, 
                                       addRandom = ran.ch,
                                       dropRandom = drop.ran, 
                                       group = grp,
                                       label = lab, 
                                       allow.unconverged = allow.unconverged, 
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = checkboundaryonly, 
                                       update = update, 
                                       IClikelihood = IClikelihood, 
                                       which.IC = which.IC), 
                                  inargs))
  }
  return(tspl.asrt)
}



#Fit a tensor-spline spatial model
fitTPPSMod <- function(asrtests.obj, sections = NULL, 
                       row.covar = "cRow", col.covar = "cCol", 
                       dropRowterm = NULL, dropColterm = NULL, 
                       nsegs = NULL, nestorder = c(1, 1), 
                       degree = c(3,3), difforder = c(2,2), 
                       asreml.opt = "mbf", 
                       tpps4mbf.obj = NULL, 
                       allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                       checkboundaryonly = FALSE, update = TRUE, 
                       chooseOnIC = TRUE, 
                       IClikelihood = "full", which.IC = "AIC",
                       ...)
{ 

  inargs <- list(...)
  checkEllipsisArgs("makeTPPSplineMats.data.frame", inargs)
  checkEllipsisArgs("tpsmmb", inargs, pkg = "asremlPlus")
  
  #Get data for mbf 
  if (asreml.opt == "mbf")
  {
    #stop('Sorry, but the mbf setting of asreml.opt is not functioning yet - use asreml.opt = "grp".')
    if (is.null(tpps4mbf.obj))
      stop("For asreml.option set to mbf, the argument tpps4mbf must be set to a list",
           " produced using makeTPPSSplineMats.data.frame also with asreml.option set to mbf")
    assign("tps.XZmat", tpps4mbf.obj)
    
    #Get the data
    dat <- tpps4mbf.obj[[1]]$data.plus
    checkNamesInData(c(sections, row.covar, col.covar, dropRowterm, dropColterm), dat)
    #Update the asreml.obj for the new data.frame
    asreml.obj  <- asrtests.obj$asreml.obj
    asreml.obj <- asreml::update.asreml(asreml.obj, data = dat)
    #Set the mbf.env in the model.frame attribute 
    if (is.null(asreml.obj$mf)) #mf is not in mf component and so must be an RDS file
    {
      mf.file <- asreml.obj$call$model.frame
      mf <- readRDS(file = mf.file)
      attr(mf, which = "mbf.env") <- attr(tps.XZmat[[1]]$data.plus, which = "mbf.env")
      saveRDS(mf, file = mf.file)
    } else
    {
      if (inherits(asreml.obj$mf, what = "asr.model.frame") || 
          inherits(asreml.obj$mf, what = "asreml.model.frame"))
        attr(asreml.obj$mf, which = "mbf.env") <- attr(tps.XZmat[[1]]$data.plus, which = "mbf.env")
      else
        stop("For asreml.option set to mbf, cannot find the asreml model frame to set the mbf environment")
    }
  } else #doing grp
  {  
    #Check that named columns are in the data
    dat.in <- asrtests.obj$asreml.obj$call$data
    if (is.symbol(dat.in))
      dat.in <- eval(dat.in)
    checkNamesInData(c(sections, row.covar, col.covar, dropRowterm, dropColterm), dat.in)
    #Check conformability of covars and factors
    if (!is.null(dropRowterm) && nlevels(dat.in[[dropRowterm]]) != length(unique(dat.in[[row.covar]])))
      stop(dropRowterm, " does not have the same number of levels as there are values of ", row.covar)
    if (!is.null(dropColterm) && nlevels(dat.in[[dropColterm]]) != length(unique(dat.in[[col.covar]])))
      stop(dropColterm, " does not have the same number of levels as there are values of ", col.covar)
    
    tps.XZmat <- makeTPPSplineMats(dat.in, sections = sections, 
                                    row.covar = row.covar, col.covar = col.covar,
                                    nsegs = nsegs, nestorder = nestorder,
                                    degree = degree, difforder = difforder,
                                    asreml.opt = "grp", ...)
    dat <- tps.XZmat[[1]]$data.plus
}

  #Update the asreml.obj for the new data.frame
  asreml.obj  <- asrtests.obj$asreml.obj
  asreml.obj <- asreml::update.asreml(asreml.obj, data = dat)
  tspl.asrt <- as.asrtests(asreml.obj = asreml.obj, NULL, NULL, 
                           IClikelihood = "full", label = "Change to new data.frame with TPS bits")

  #Fit spatial TPPS to sections
  nsect <- calc.nsect(dat, sections)
  for (i in 1:nsect)
  {
    if (nsect == 1)
    { 
      stub = "xx"
      sect.fac <- NULL
      lab <- paste0("Try tensor P-splines")
    }
    else
    { 
      stub <- levels(dat[[sections]])[i]
      sect.fac <- paste0("at(", sections, ",  ", i, "):")
      lab <- paste0("Try tensor P-splines for ", sections, " ",i)
    }
    tspl.asrt <- fitTPSModSect(tspl.asrt, mat = tps.XZmat[[i]], sect.fac = sect.fac, 
                               dropRowterm = dropRowterm, dropColterm = dropColterm, 
                               row.covar = row.covar, col.covar = col.covar, 
                               degree = degree, difforder = difforder,
                               lab = lab, asreml.opt = asreml.opt, stub = stub, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, 
                               chooseOnIC = chooseOnIC, 
                               IClikelihood = IClikelihood, 
                               which.IC = which.IC, 
                               ...)
  }
  
  return(tspl.asrt)
}

