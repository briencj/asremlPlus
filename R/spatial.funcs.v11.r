#The following functions have been developed for single and multiple section experiments

addSpatialModel.asrtests <- function(asrtests.obj, spatial.model = "TPPS", 
                                     sections = NULL, 
                                     row.covar = "cRow", col.covar = "cCol", 
                                     row.factor = "Row", col.factor = "Col", 
                                     corr.funcs = c("ar1", "ar1"), corr.orders = c(0, 0), 
                                     row.corrFitfirst = TRUE, 
                                     allow.corrsJointFit = TRUE, nugget.variance = TRUE, 
                                     dropFixed = NULL, dropRandom = NULL, 
                                     nsegs = NULL, nestorder = c(1, 1), 
                                     degree = c(3,3), difforder = c(2,2), 
                                     usRandLinCoeffs = TRUE, 
                                     rotateX = FALSE, ngridangles = NULL, 
                                     which.rotacriterion = "AIC", nrotacores = 1, 
                                     asreml.option = "grp", tpps4mbf.obj = NULL,  
                                     allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                                     checkboundaryonly = FALSE, update = TRUE, trace = FALSE, 
                                     maxit = 30, IClikelihood = "full", which.IC = "AIC", 
                                     ...)
{    
  #Deal with arguments for tpsmmb and changeModelOnIC
  inargs <- list(...)
  if (any(c("dropRowterm", "dropColterm") %in% names(inargs)))
    stop(paste("The arguments dropRowterm and dropColterm have been deprecated;",
               "use dropFixed and dropRandom instead"))
  
  checkEllipsisArgs_tpsmmb(c("changeTerms.asrtests", "asreml"), inargs)
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
 
  #Check if have separate section random and residual terms
  checkSections4RanResTerms(asrtests.obj, sections = sections, asr4 = asr4)
   
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

  #Fit a local spatial model involving correlated effects
  if ("corr" %in% spatial.mod)
    spatial.asrt <- fitCorrMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               row.factor = row.factor, col.factor = col.factor, 
                               corr.funcs = corr.funcs, corr.orders = corr.orders, 
                               row.corrFitfirst = row.corrFitfirst, 
                               allow.corrsJointFit = allow.corrsJointFit, 
                               nugget.variance = nugget.variance, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, trace = trace, chooseOnIC = FALSE, 
                               maxit = maxit, IClikelihood = ic.lik, 
                               which.IC = ic.type, ...)
  #Fit a local spatial model involving TPNCSS
  if ("TPNCSS" %in% spatial.mod)
    spatial.asrt <- fitTPNCSSMod(asrtests.obj, sections = sections, 
                                 row.covar = row.covar, col.covar = col.covar, 
                                 dropFixed = dropFixed, dropRandom = dropRandom, 
                                 allow.unconverged = allow.unconverged, 
                                 allow.fixedcorrelation = allow.fixedcorrelation,
                                 checkboundaryonly = checkboundaryonly, 
                                 update = update, trace = trace, chooseOnIC = FALSE, 
                                 maxit = maxit, IClikelihood = "full", 
                                 which.IC = ic.type, ...)
  
  #Fit a residual spatial model involving TPPS
  if ("TPPS" %in% spatial.mod)
    spatial.asrt <- fitTPPSMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               dropFixed = dropFixed, dropRandom = dropRandom, 
                               nsegs = nsegs, nestorder = nestorder, 
                               degree = degree, difforder = difforder, 
                               usRandLinCoeffs = usRandLinCoeffs, 
                               rotateX = rotateX, ngridangles = ngridangles, 
                               which.rotacriterion = which.rotacriterion, 
                               nrotacores = nrotacores, 
                               asreml.opt = asreml.opt, 
                               tpps4mbf.obj = tpps4mbf.obj,
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, trace = trace, chooseOnIC = FALSE, 
                               maxit = maxit, IClikelihood = ic.lik, 
                               which.IC = ic.type, ...)
  
  return(spatial.asrt)  
}

addSpatialModelOnIC.asrtests <- function(asrtests.obj, spatial.model = "TPPS", 
                                         sections = NULL, 
                                         row.covar = "cRow", col.covar = "cCol", 
                                         row.factor = "Row", col.factor = "Col", 
                                         corr.funcs = c("ar1", "ar1"), corr.orders = c(0, 0), 
                                         row.corrFitfirst = TRUE, 
                                         allow.corrsJointFit = TRUE, nugget.variance = TRUE, 
                                         dropFixed = NULL, dropRandom = NULL, 
                                         nsegs = NULL, nestorder = c(1, 1), 
                                         degree = c(3,3), difforder = c(2,2), 
                                         usRandLinCoeffs = TRUE, 
                                         rotateX = FALSE, ngridangles = NULL, 
                                         which.rotacriterion = "AIC", 
                                         nrotacores = 1, 
                                         asreml.option = "grp", tpps4mbf.obj = NULL,  
                                         allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                                         checkboundaryonly = FALSE, update = TRUE, trace = FALSE, 
                                         maxit = 30, IClikelihood = "full", which.IC = "AIC", 
                                         ...)
{    
  #Deal with arguments for tpsmmb and changeModelOnIC
  inargs <- list(...)
  if (any(c("dropRowterm", "dropColterm") %in% names(inargs)))
    stop(paste("The arguments dropRowterm and dropColterm have been deprecated;",
               "use dropFixed and dropRandom instead"))
  checkEllipsisArgs_tpsmmb(c("changeModelOnIC.asrtests", "asreml"), inargs)

  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check if have separate section random and residual terms
  checkSections4RanResTerms(asrtests.obj, sections = sections, asr4 = asr4)
  
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

  #Fit a local spatial model involving correlated effects
  if ("corr" %in% spatial.mod)
    spatial.asrt <- fitCorrMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               row.factor = row.factor, col.factor = col.factor, 
                               corr.funcs = corr.funcs, corr.orders = corr.orders, 
                               row.corrFitfirst = row.corrFitfirst, 
                               allow.corrsJointFit = allow.corrsJointFit, 
                               nugget.variance = nugget.variance, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, trace = trace, chooseOnIC = TRUE, 
                               maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                               ...)
  #Fit a local spatial model involving TPNCSS
  if ("TPNCSS" %in% spatial.mod)
    spatial.asrt <- fitTPNCSSMod(asrtests.obj, sections = sections, 
                                 row.covar = row.covar, col.covar = col.covar, 
                                 dropFixed = dropFixed, dropRandom = dropRandom, 
                                 allow.unconverged = allow.unconverged, 
                                 allow.fixedcorrelation = allow.fixedcorrelation,
                                 checkboundaryonly = checkboundaryonly, 
                                 update = update, trace = trace, chooseOnIC = TRUE, 
                                 maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                                 ...)

  #Fit a residual spatial model involving TPPS
  if ("TPPS" %in% spatial.mod)
    spatial.asrt <- fitTPPSMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               dropFixed = dropFixed, dropRandom = dropRandom, 
                               nsegs = nsegs, nestorder = nestorder, 
                               degree = degree, difforder = difforder, 
                               usRandLinCoeffs = usRandLinCoeffs, 
                               rotateX = rotateX, ngridangles = ngridangles, 
                               asreml.opt = asreml.opt, 
                               tpps4mbf.obj = tpps4mbf.obj,
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation, 
                               which.rotacriterion = which.rotacriterion, 
                               nrotacores = nrotacores, 
                               checkboundaryonly = checkboundaryonly, 
                               update = update, trace = trace, chooseOnIC = TRUE, 
                               maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                               ...)
  return(spatial.asrt)  
}

#This function calculates the IC statistics for a fitted spatial model, 
#irrespective of whether the finally fitted model was better than the nonspatial model
calcSpatialICs <- function(spatial.asrt, spatial.mod, IClikelihood = "full", 
                           spatial.IC)
{
  tests.sp <- infoCriteria(spatial.asrt$asreml.obj, IClikelihood = IClikelihood)
  tests.sp <- tests.sp[-match("NBound", names(tests.sp))]
  rownames(tests.sp) <- spatial.mod
  spatial.IC <- rbind(spatial.IC,tests.sp)
  return(spatial.IC)
}

chooseSpatialModelOnIC.asrtests <- function(asrtests.obj, trySpatial = "all", 
                                            sections = NULL, 
                                            row.covar = "cRow", col.covar = "cCol", 
                                            row.factor = "Row", col.factor = "Col", 
                                            corr.funcs = c("ar1", "ar1"), corr.orders = c(0, 0), 
                                            row.corrFitfirst = TRUE, 
                                            allow.corrsJointFit = TRUE, nugget.variance = TRUE, 
                                            dropFixed = NULL, dropRandom = NULL, 
                                            nsegs = NULL, nestorder = c(1, 1), 
                                            usRandLinCoeffs = TRUE, 
                                            rotateX = FALSE, ngridangles = NULL, 
                                            which.rotacriterion = "AIC", nrotacores = 1, 
                                            asreml.option = "grp", tpps4mbf.obj = NULL, 
                                            allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                                            checkboundaryonly = FALSE, update = TRUE, trace = FALSE, 
                                            maxit = 30, IClikelihood = "full", which.IC = "AIC", 
                                            return.asrts = "best", ...)
{    
  #Deal with arguments for tpsmmb and changeModelOnIC
  inargs <- list(...)
  if (any(c("dropRowterm", "dropColterm") %in% names(inargs)))
    stop(paste("The arguments dropRowterm and dropColterm have been deprecated;",
               "use dropFixed and dropRandom instead"))
  checkEllipsisArgs_tpsmmb(c("changeModelOnIC.asrtests", "asreml"), inargs)

  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check if have separate section random and residual terms
  checkSections4RanResTerms(asrtests.obj, sections = sections, asr4 = asr4)
  
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
                                            corr.funcs = corr.funcs, corr.orders = corr.orders, 
                                            row.corrFitfirst = row.corrFitfirst, 
                                            allow.corrsJointFit = allow.corrsJointFit, 
                                            nugget.variance = nugget.variance, 
                                            allow.unconverged = allow.unconverged, 
                                            allow.fixedcorrelation = allow.fixedcorrelation,
                                            checkboundaryonly = checkboundaryonly, 
                                            update = update, trace = trace, chooseOnIC = TRUE, 
                                            maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                                            ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["corr"]], spatial.mod = "corr", 
                                   IClikelihood = ic.lik, spatial.IC = spatial.IC)
    }
    
    #Fit a local spatial model involving TPNCSS
    if ("TPNCSS" %in% trySpatial)
    { 
      spatial.asrts[["TPNCSS"]] <- fitTPNCSSMod(asrtests.obj, sections = sections, 
                                                row.covar = row.covar, col.covar = col.covar, 
                                                dropFixed = dropFixed, dropRandom = dropRandom, 
                                                allow.unconverged = allow.unconverged, 
                                                allow.fixedcorrelation = allow.fixedcorrelation,
                                                checkboundaryonly = checkboundaryonly, 
                                                update = update, trace = trace, chooseOnIC = TRUE, 
                                                maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                                                ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["TPNCSS"]] , spatial.mod = "TPNCSS", 
                                   IClikelihood = ic.lik, spatial.IC = spatial.IC)
    }
    
    #Fit a residual spatial model involving TPPSC2
    if ("TPPSC2" %in% trySpatial)
    { 
      spatial.asrts[["TPPSC2"]] <- fitTPPSMod(asrtests.obj, sections = sections, 
                                             row.covar = row.covar, col.covar = col.covar, 
                                             dropFixed = dropFixed, dropRandom = dropRandom, 
                                             nsegs = nsegs, nestorder = nestorder, 
                                             degree = c(3,3), difforder = c(2,2), 
                                             rotateX = rotateX, ngridangles = ngridangles, 
                                             usRandLinCoeffs = usRandLinCoeffs, 
                                             which.rotacriterion = which.rotacriterion, 
                                             nrotacores = nrotacores, 
                                             asreml.opt = asreml.opt, 
                                             tpps4mbf.obj = tpps4mbf.obj,
                                             allow.unconverged = allow.unconverged, 
                                             allow.fixedcorrelation = allow.fixedcorrelation,
                                             checkboundaryonly = checkboundaryonly, 
                                             update = update, trace = trace, chooseOnIC = TRUE, 
                                             maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                                             ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["TPPSC2"]] , spatial.mod = "TPPSC2", 
                                   IClikelihood = ic.lik, spatial.IC = spatial.IC)
    }
    
    #Fit a residual spatial model involving TPPSL1
    if ("TPPSL1" %in% trySpatial)
    { 
      spatial.asrts[["TPPSL1"]] <- fitTPPSMod(asrtests.obj, sections = sections, 
                                              row.covar = row.covar, col.covar = col.covar, 
                                              dropFixed = dropFixed, dropRandom = dropRandom, 
                                              nsegs = nsegs, nestorder = nestorder, 
                                              degree = c(1,1), difforder = c(1,1), 
                                              usRandLinCoeffs = FALSE, rotateX = FALSE, 
                                              asreml.opt = asreml.opt, 
                                              tpps4mbf.obj = tpps4mbf.obj,
                                              allow.unconverged = allow.unconverged, 
                                              allow.fixedcorrelation = allow.fixedcorrelation,
                                              checkboundaryonly = checkboundaryonly, 
                                              update = update, trace = trace, chooseOnIC = TRUE, 
                                              maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                                              ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["TPPSL1"]] , spatial.mod = "TPPSL1", 
                                   IClikelihood = ic.lik, spatial.IC = spatial.IC)
    }
    
    #Find min AIC and, if multiple mins, select in specified order
    spatial.comp <- round(spatial.IC[[which.IC]], digits = 3)
    names(spatial.comp) <- rownames(spatial.IC)
    min.asrt <- which.min(spatial.comp)
    if (length(min.asrt) > 1)
    {
      #pick one in the order given below
      if ("nonspatial" %in% names(min.asrt)) min.asrt <- min.asrt["nonspatial"]
      if ("TPPSC2" %in% names(min.asrt)) min.asrt <- min.asrt["TPPSC2"]
      if ("TPNCSS" %in% names(min.asrt)) min.asrt <- min.asrt["TPNCSS"]
      if ("TPPSL1" %in% names(min.asrt)) min.asrt <- min.asrt["TPPSL1"]
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

fitCorrMod <- function(asrtests.obj, sections = NULL,
                       row.covar = "cRow", col.covar = "cCol", 
                       row.factor = "Row", col.factor = "Col", 
                       corr.funcs = c("ar1", "ar1"), corr.orders = c(0, 0), 
                       row.corrFitfirst = TRUE, 
                       allow.corrsJointFit = TRUE, nugget.variance = TRUE,
                       allow.unconverged = allow.unconverged, 
                       allow.fixedcorrelation = allow.fixedcorrelation,
                       checkboundaryonly = checkboundaryonly, update = update, 
                       trace = trace, chooseOnIC = TRUE, maxit = 30, 
                       IClikelihood = "full", which.IC = "AIC", 
                       ...)
{
  inargs <- list(...)

  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  if (!asr4)
    stop(paste("Fitting spatial models using correlation/variance models", 
               "has not been implemented for asreml version less than 4.0"))
  
  #Save ai.sing setting so can make sure that it is resored on exit (reset in corb at moment)
  ksing <-   get("asr_options", envir = getFromNamespace(".asremlEnv", "asreml"))$ai.sing

  
  sing.excl <- c("S", "?")
  bounds.excl <- c("B", sing.excl)
  all.bounds.excl <- c(bounds.excl, "F")
  corr.types <- c("R", "P")
  
  #Check that named columns are in the data
  dat.in <- asrtests.obj$asreml.obj$call$data
  if (is.symbol(dat.in))
    dat.in <- eval(dat.in)
  
  #Check the correlation functions
  if (any(get.specials("unimpl.specials") %in% corr.funcs))
    stop("Some of the following corr.funcs ar not implemented for spatial modelling: ", 
         paste(get.specials("unimpl.specials")[unimpl.funcs %in% corr.funcs], collapse = ","))
  if (all(corr.funcs %in% get.specials("id.specials")))
    stop("Both correlation functions are id or equivalent")

  #Check the grid covars and factors
  #row.factor and col.factor must be in dat.in so that each dimension can be fitted independently
  #(so can have a term without a corr function when fitting a dimension)
  grid.cols <- c(sections, row.factor, col.factor)
  if (any(corr.funcs %in% get.specials("met.specials")))
    grid.cols <- c(grid.cols, row.covar, col.covar)
  checkNamesInData(c(sections, row.factor, col.factor), dat.in)
  if (!all(sapply(dat.in[c(row.factor, col.factor)], is.factor)))
    stop("Both row.factor and col.factor must be factors in the data stored in the asreml.obj")
  
  nsect <- calc.nsect(dat.in, sections)
  
  #Get row and col corr models
  row.corr <- makeCorrSpec1D(corr.funcs = corr.funcs, corr.orders = corr.orders, dimension = 1, 
                             row.covar = row.covar, col.covar = col.covar, 
                             row.factor = row.factor, col.factor = col.factor)
  col.corr <- makeCorrSpec1D(corr.funcs = corr.funcs, corr.orders = corr.orders, dimension = 2, 
                             row.covar = row.covar, col.covar = col.covar, 
                             row.factor = row.factor, col.factor = col.factor)
  
  #Remove units if in random model
  if (grepl("units", as.character(getFormulae(asrtests.obj$asreml.obj)$random)[2]))
    asrtests.obj <- changeTerms(asrtests.obj, 
                                dropRandom = "units", label = "Remove random units term", 
                                maxit = maxit, 
                                allow.unconverged = allow.unconverged, 
                                allow.fixedcorrelation = allow.fixedcorrelation,
                                checkboundaryonly = TRUE, 
                                update = update, 
                                IClikelihood = IClikelihood, ...)
  
  #Check if correlations already included in a term 
  t <- mapply(function(func, corr, asr)
  { 
    if (!(func %in% get.specials("id.specials")) && 
        any(sapply(as.character(getFormulae(asr)), 
                   function(mod, corr) grepl(corr, mod, fixed = TRUE), 
                   corr = func)))#corr)))
      warning("The correlation function ", corr, " is already in the model")
    invisible()
  }, func = corr.funcs, corr = c(row.corr, col.corr), 
  MoreArgs = (list(asr = asrtests.obj$asreml.obj)))
  
  #Prepare to fit
  facs <- c(row.factor, col.factor)
  rfuncs <- corr.funcs
  rorders <- corr.orders
  rterms <- c(row.corr, col.corr)
  sterms <- c(ifelse(corr.funcs[1] %in% get.specials("met.specials"), 
                     row.covar, row.factor),
              ifelse(corr.funcs[2] %in% get.specials("met.specials"), 
                     col.covar, col.factor))
  
  if (!row.corrFitfirst)
  {
    facs <- facs[c(2,1)]
    rfuncs <- rfuncs[c(2,1)]
    rorders <- rorders[c(2,1)]
    rterms <- rterms[c(2,1)]
    sterms <- sterms[c(2,1)]
  }
  
  #Loop over the sections
  nuggsOK <- nugget.variance
  corr.asrt <- asrtests.obj
  for (i in 1:nsect)
  {
    init.asrt <- corr.asrt #This remains unchanged to retain the supplied model
    if (trace)  {cat("\n#### Initial fit for section", i, "\n\n"); print(init.asrt)}
    if (chooseOnIC)
    { fitfunc <- "changeModelOnIC"
    } else 
    {  fitfunc <- "changeTerms"}
    
    spat.var <- paste0(sterms, collapse = ":")
    if (nsect > 1)
    { 
      stub <- levels(dat.in[[sections]])[i]
      spat.var <- paste0("at(", sections, ", '",stub, "'):", spat.var)
    } else
      stub <- NULL

    startFixRes <- FALSE
    #### For corb, add random variance with fixed Residual
    if (any(grepl("corb", corr.funcs)) && nuggsOK)
    {  
      if (trace) 
      {cat("\n#### Try to fit nugget variance before fitting correlations\n\n"); print(corr.asrt)}

      #Do not allow singularities with corb functions because crashes R
      if (any(grepl("corb", corr.funcs)))
          asreml::asreml.options(ai.sing = FALSE)
      
      lab0 <- paste("Add random",  paste0(sterms, collapse = ":"), "and fix residual")
      if (nsect > 1)
        lab0 <- paste0(lab0, " for ", sections, " ",stub)
      #Get random and residual terms for correlation model for the current section
      vpc.corr <- getSectionVpars(corr.asrt$asreml.obj, 
                                  sections = sections, stub = stub, 
                                  sterm.facs = sterms, 
                                  asr4.2 = asr4.2)
      #Check that residual does not have heterogeneous terms unrelated to sections 
      if (is.null(vpc.corr$res)) #implies multiple residual terms, none involving sections
      {  nuggsOK <- FALSE
      } else
      {
        #Try fixing either the  single residual variance term or that for the current section 
        old.inargs <- inargs
        #add set.terms arguments to inargs so promulgated to further fitting functions
        inargs <- addSetterms2inargs(setterms = list(set.terms = names(vpc.corr$res), 
                                                     ignore.suffices = FALSE, 
                                                     bounds = "F", initial.values = 1),
                                     inargs)
        tmp.asrt <- tryCatchLog(
          do.call(changeTerms,
                  c(list(corr.asrt, 
                         addRandom = spat.var, 
                         label = lab0, 
                         allow.fixedcorrelation = allow.fixedcorrelation, 
                         allow.unconverged = allow.unconverged, 
                         IClikelihood = IClikelihood, 
                         checkboundaryonly = TRUE), #need so don't remove spat.var = B when res is F
                    inargs)),
          error = function(e) 
          {print(paste("Failed attempting to fit correlations to both dimensions;",
                       "continued analysis without them")); NULL}, 
          include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
        
        if (!is.allnull(tmp.asrt))
        {       
          vpc.corr <- getSectionVpars(tmp.asrt$asreml.obj, 
                                      sections = sections, stub = stub, 
                                      sterm.facs = sterms, 
                                      asr4.2 = asr4.2)
          lasttest <- tail(tmp.asrt$test.summary, 1)
          if (grepl("and fix residual", lasttest$terms) && 
              !(grepl("Unswapped", lasttest$action) || grepl("Unchanged", lasttest$action)) &&
              !is.allnull(vpc.corr$ran) && !(vpc.corr$ran[spat.var] %in% bounds.excl))
          {  
            startFixRes <- TRUE
            corr.asrt <- tmp.asrt
          } else
          {  
            if (lasttest$action == "Changed random" && 
                (vpc.corr$ran[spat.var] %in% bounds.excl))
            test.summary <- addtoTestSummary(corr.asrt$test.summary, terms = lab0, 
                                             DF = NA, denDF = NA, p = NA, 
                                             AIC = NA, BIC = NA, 
                                             action = "Unchanged - bound")
            corr.asrt$test.summary <- test.summary
            inargs <- old.inargs
          }
        } else
        {  
          test.summary <- addtoTestSummary(corr.asrt$test.summary, terms = lab0, 
                                           DF = NA, denDF = NA, p = NA, 
                                           AIC = NA, BIC = NA, 
                                           action = "Unchanged - Singular")
          corr.asrt$test.summary <- test.summary
          inargs <- old.inargs
        }
      }
    } #End fix nugget section 

    #Start fitting the first correlation
    corr.term <- FALSE
    #Check have a corr func
    if (any(rfuncs[1] == get.specials("id.specials")))
    { result1 <- "Unswapped"
    } else
    { 
      vpc <- getSectionVpars(corr.asrt$asreml.obj, 
                             sections = sections, stub = stub, 
                             sterm.facs = sterms, 
                             asr4.2 = asr4.2)
      #Check if residual terms for each section
      if (!length(vpc$res))
        warning("Could not find a residual term for ", sections, " ", stub)

      #### Try first correl in current section
      if (trace) {cat("\n#### Fit first correlation\n\n"); print(corr.asrt)}
      ran.term1 <- paste0(rterms[1], ":", facs[2])
      lab1 <- paste0("Try ", rterms[1])
      if (nsect > 1)
      {  
        ran.term1 <- paste0("at(", sections, ", '",stub, "'):", ran.term1)
        lab1 <- paste0(lab1, " for ", sections, " ",stub)
      }
      #Determine if spat.var is a fitted random term
      vpc.ran <-vpc$ran
      spat.term <- findterm(spat.var, names(vpc.ran)) #allows for changed order
      if (length(spat.term) == 1 && spat.term != 0) #have got a single spat.term
      { drop.spatvar <- names(spat.term)
      } else
      {  drop.spatvar <- NULL}
      tmp.asrt <- tryCatchLog(
        do.call(fitfunc, 
                c(list(corr.asrt, label = lab1, 
                       addRandom = ran.term1, dropRandom = drop.spatvar, 
                       allow.unconverged = allow.unconverged, 
                       allow.fixedcorrelation = allow.fixedcorrelation,
                       maxit = maxit, 
                       checkboundaryonly = TRUE, 
                       update = update, 
                       IClikelihood = IClikelihood, 
                       which.IC = which.IC), 
                  inargs)),
        error = function(e) {print("Analysis continued"); NULL}, 
        include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
      if (largeVparChange(corr.asrt$asreml.obj, 0.75))
        corr.asrt <- iterate(corr.asrt)
      
      #Check for singular (S) spatial terms - only change model if var == S can be made F and 
      #none of the correlations are S
      corr.asrt <- chk4SingularSpatTerms(tmp.asrt, corr.asrt,  label = lab1, 
                                         sections = sections, stub = stub, 
                                         sterm.facs = sterms, 
                                         asr4 = asr4, asr4.2 = asr4.2, 
                                         maxit = maxit, 
                                         allow.unconverged = allow.unconverged, 
                                         allow.fixedcorrelation = allow.fixedcorrelation,
                                         checkboundaryonly = checkboundaryonly, 
                                         update = update, 
                                         IClikelihood = IClikelihood, 
                                         which.IC = which.IC,
                                         bounds.excl =  bounds.excl, 
                                         sing.excl = sing.excl)
      if (is.allnull(tmp.asrt))
      {
        test.summary <- addtoTestSummary(corr.asrt$test.summary, terms = lab1, 
                                         DF = NA, denDF = NA, p = NA, 
                                         AIC = NA, BIC = NA, 
                                         action = "Unchanged - Singular")
        corr.asrt$test.summary <- test.summary
      }
      result1 <- getTestEntry(corr.asrt, label = lab1)$action
      
      #If corb and rorder == 0, try to fit corb up to order 10
      corr.lis <- do.call(fitCorbPlus1, 
                          c(list(corr.asrt, ran.term = ran.term1, rorder = rorders[1], 
                                 lab = lab1, result = result1, dimension = 0, 
                                 IClikelihood = IClikelihood, trace = trace), 
                            inargs))
      corr.asrt <- corr.lis$asrt
      if (ran.term1 != corr.lis$last.term)
      {
        lab1 <- corr.lis$last.lab
        result1 <- corr.lis$result
        ran.term1 <- corr.lis$last.term
      }
      #If all correlation model terms are bound, reinstate initial model
      lab.drop <- "Dropped correlations"
      corr.asrt <- allBoundSectionVpars(corr.asrt, init.asrt, 
                                        lab = lab.drop, #lab1, 
                                        sections = sections, stub = stub, 
                                        sterm.facs = sterms, 
                                        all.bounds.excl  = all.bounds.excl, 
                                        asr4.2 = asr4.2)
      result.drop <- getTestEntry(corr.asrt, label = lab.drop, error.absent = FALSE)$action
      if (!is.null(result.drop)) 
         result1 <- result.drop

      #Determine if any correlations have been fitted
      corr.term <- setCorrTerm(corr.asrt$asreml.obj, spat.var = spat.var, 
                               asr4.2 = asr4.2)
      
    } #End of first correlation section
    
    #### Try 2nd correl in current section
    if (!any(rfuncs[2] == get.specials("id.specials")))
    {  
      first.asrt <- corr.asrt
      lab <- paste0("Try ", rterms[2])
      if (nsect > 1)
        lab <- paste0(lab, " for ", sections, " ",stub)
      # Has first fac corr.func been fitted
      if (corr.term) #yes
      { 
        if (trace) 
          {cat("\n#### Add second correlation to first correlation\n\n"); print(corr.asrt)}
        last.term <- ran.term1
        
        #Check for ran.term1 in random formula, and if absent, check for different order
        last.term <- chk4TermInFormula(corr.asrt$asreml.obj$call$random, term = last.term, 
                                       asreml.obj = corr.asrt$asreml.obj)
        ran.term <- paste0(rterms[1], ":", rterms[2])
        if (nsect > 1)
          ran.term <- paste0("at(", sections, ", '",stub, "'):", ran.term)
        
        tmp.asrt <- tryCatchLog(
          do.call(fitfunc, 
                  c(list(corr.asrt, 
                         addRandom = ran.term, 
                         dropRandom = last.term, label = lab, 
                         maxit = maxit, 
                         allow.unconverged = allow.unconverged, 
                         allow.fixedcorrelation = allow.fixedcorrelation,
                         checkboundaryonly = checkboundaryonly, 
                         update = update, 
                         IClikelihood = IClikelihood, 
                         which.IC = which.IC),
                    inargs)),
          error = function(e) {print("Analysis continued"); NULL}, 
          include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
        #Check for singular (S) spatial terms - only change model if var == S can be made F and 
        #none of the correlations are S
        corr.asrt <- chk4SingularSpatTerms(tmp.asrt, corr.asrt,  label = lab, 
                                           sections = sections, stub = stub, 
                                           sterm.facs = sterms, 
                                           asr4 = asr4, asr4.2 = asr4.2, 
                                           maxit = maxit, 
                                           allow.unconverged = allow.unconverged, 
                                           allow.fixedcorrelation = allow.fixedcorrelation,
                                           checkboundaryonly = checkboundaryonly, 
                                           update = update, 
                                           IClikelihood = IClikelihood, 
                                           which.IC = which.IC,
                                           bounds.excl =  bounds.excl, 
                                           sing.excl = sing.excl)
        if (is.allnull(tmp.asrt))
        {
          test.summary <- addtoTestSummary(corr.asrt$test.summary, terms = lab, 
                                           DF = NA, denDF = NA, p = NA, 
                                           AIC = NA, BIC = NA, 
                                           action = "Unchanged - Singular")
          corr.asrt$test.summary <- test.summary
        }
        result <- getTestEntry(corr.asrt, label = lab)$action
        if (!(grepl("Unswapped", result)) && !(grepl("Unchanged", result)))
          last.term <- ran.term
       
        #If corb and rorder == 0, try to fit corb up to order 10
        corr.lis <- do.call(fitCorbPlus1, 
                            c(list(corr.asrt, ran.term = ran.term, rorder = rorders[2], 
                                   lab = lab, result = result, dimension = 2, 
                                   IClikelihood = IClikelihood, trace = trace), 
                              inargs))
        corr.asrt <- corr.lis$asrt
        if (ran.term != corr.lis$last.term)
        {
          lab <- corr.lis$last.lab
          result <- corr.lis$result
          ran.term <- corr.lis$last.term
        }
        #If all correlation model terms are bound, reinstate first correlation model
        lab.drop <- "Dropped seconf correlation"
        corr.asrt <- allBoundSectionVpars(corr.asrt, first.asrt, 
                                          lab = lab.drop, #lab1, 
                                          sections = sections, stub = stub, 
                                          sterm.facs = sterms, 
                                          all.bounds.excl  = all.bounds.excl, 
                                          asr4.2 = asr4.2)
        result.drop <- getTestEntry(corr.asrt, label = lab.drop, error.absent = FALSE)$action
        if (!is.null(result.drop)) 
          result <- result.drop
      } else #### no first fac corr
      { 
        if (trace) 
          {cat("\n### Fit second correlation - no first correlation\n\n"); print(corr.asrt)}
        ran.term <- paste0(facs[1], ":", rterms[2])
        if (nsect > 1)
          ran.term <- paste0("at(", sections, ", '",stub, "'):", ran.term)
        #Determine if spat.var is a fitted random term
        vpc.ran <- getSectionVpars(corr.asrt$asreml.obj, 
                                   sections = sections, stub = stub, 
                                   sterm.facs = sterms, 
                                   asr4.2 = asr4.2)$ran
        spat.term <- findterm(spat.var, names(vpc.ran)) #allows for changed order
        if (length(spat.term) == 1 && spat.term != 0) #have got a single spat.term in fit
        {   drop.spatvar <- names(spat.term)
        } else
        {  drop.spatvar <- NULL}

        tmp.asrt <- tryCatchLog(
          do.call(fitfunc, 
                  c(list(corr.asrt, label = lab, 
                         addRandom = ran.term, dropRandom = drop.spatvar, 
                         maxit = maxit, 
                         allow.unconverged = allow.unconverged, 
                         allow.fixedcorrelation = allow.fixedcorrelation,
                         checkboundaryonly = checkboundaryonly, 
                         update = update, 
                         IClikelihood = IClikelihood, 
                         which.IC = which.IC), 
                    inargs)),
          error = function(e) {print("Analysis continued"); NULL}, 
          include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
        #Check for singular (S) spatial terms - only change model if var == S can be made F and 
        #none of the correlations are S
        corr.asrt <- chk4SingularSpatTerms(tmp.asrt, corr.asrt, label = lab, 
                                           sections = sections, stub = stub, 
                                           sterm.facs = sterms, 
                                           asr4 = asr4, asr4.2 = asr4.2, 
                                           maxit = maxit, 
                                           allow.unconverged = allow.unconverged, 
                                           allow.fixedcorrelation = allow.fixedcorrelation,
                                           checkboundaryonly = checkboundaryonly, 
                                           update = update, 
                                           IClikelihood = IClikelihood, 
                                           which.IC = which.IC,
                                           bounds.excl =  bounds.excl, 
                                           sing.excl = sing.excl)
        
        if (is.allnull(tmp.asrt))
        {
          test.summary <- addtoTestSummary(corr.asrt$test.summary, terms = lab, 
                                           DF = NA, denDF = NA, p = NA, 
                                           AIC = NA, BIC = NA, 
                                           action = "Unchanged - Singular")
          corr.asrt$test.summary <- test.summary
        }
        #Determine if any correlations have been fitted
        corr.term <- setCorrTerm(corr.asrt$asreml.obj, spat.var = spat.var, 
                                 asr4.2 = asr4.2)
        
        
        result <- getTestEntry(corr.asrt, label = lab)$action
        if (corr.term)
        { 
          #If corb and rorder == 0, try to fit corb up to order 10
          corr.lis <- do.call(fitCorbPlus1, 
                              c(list(corr.asrt, ran.term = ran.term, rorder = rorders[2], 
                                     lab = lab, result = result, dimension = 2, 
                                     IClikelihood = IClikelihood), 
                                inargs))
          corr.asrt <- corr.lis$asrt
          if (ran.term != corr.lis$last.term)
          {
            lab <- corr.lis$last.lab
            result <- corr.lis$result
            ran.term <- corr.lis$last.term
          }
         if (!grepl("Unswapped", result) && !grepl("Unchanged", result))
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
            if (trace) 
            {cat(paste("\n#### Try to fit first correlation again,",
                       "after fitting second correlation\n\n")); print(corr.asrt)}
            ran.term1 <- paste0(rterms[1], ":", rterms[2])
            if (nsect > 1)
              ran.term1 <- paste0("at(", sections, ", '",stub, "'):", ran.term1)
            tmp.asrt <- tryCatchLog(
              do.call(fitfunc, 
                      c(list(corr.asrt, 
                             dropRandom = last.term,
                             addRandom = ran.term1, label = lab1, 
                             maxit = maxit, 
                             allow.unconverged = FALSE, 
                             allow.fixedcorrelation = FALSE,
                             checkboundaryonly = FALSE, 
                             update = update, 
                             IClikelihood = IClikelihood, 
                             which.IC = which.IC), 
                        inargs)), 
              error = function(e) {print("Analysis continued"); NULL}, 
              include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
            if (largeVparChange(corr.asrt$asreml.obj, 0.75))
              corr.asrt <- iterate(corr.asrt)
            #Check for singular (S) spatial terms - only change model if var == S can be made F and 
            #none of the correlations are S
            corr.asrt <- chk4SingularSpatTerms(tmp.asrt, corr.asrt,  label = lab1, 
                                               sections = sections, stub = stub, 
                                               sterm.facs = sterms, 
                                               asr4 = asr4, asr4.2 = asr4.2, 
                                               maxit = maxit, 
                                               allow.unconverged = allow.unconverged, 
                                               allow.fixedcorrelation = allow.fixedcorrelation,
                                               checkboundaryonly = checkboundaryonly, 
                                               update = update, 
                                               IClikelihood = IClikelihood, 
                                               which.IC = which.IC,
                                               bounds.excl =  bounds.excl, 
                                               sing.excl = sing.excl)
            
            if (is.allnull(tmp.asrt))
            {
              test.summary <- addtoTestSummary(corr.asrt$test.summary, terms = lab1, 
                                               DF = NA, denDF = NA, p = NA, 
                                               AIC = NA, BIC = NA, 
                                               action = "Unchanged - Singular")
              corr.asrt$test.summary <- test.summary
            }
            result1 <- getTestEntry(corr.asrt, label = lab1)$action
            if (largeVparChange(corr.asrt$asreml.obj, 0.75))
              corr.asrt <- iterate(corr.asrt)
            #If corb and rorder == 0, try to fit corb up to order 10
            corr.lis <- do.call(fitCorbPlus1, 
                                c(list(corr.asrt, ran.term = ran.term1, rorder = rorders[1], 
                                       lab = lab1, result = result1, dimension = 1, 
                                       IClikelihood = IClikelihood), 
                                  inargs))
            corr.asrt <- corr.lis$asrt
            if (ran.term1 != corr.lis$last.term)
            {
              lab1 <- corr.lis$last.lab
              result1 <- corr.lis$result
              ran.term1 <- corr.lis$last.term
            }
          }
          #If all correlation model terms are bound, reinstate initial model
          lab.drop <- "Dropped correlations"
          corr.asrt <- allBoundSectionVpars(corr.asrt, init.asrt, 
                                            lab = lab.drop, #lab1, 
                                            sections = sections, stub = stub, 
                                            sterm.facs = sterms, 
                                            all.bounds.excl  = all.bounds.excl, 
                                            asr4.2 = asr4.2)
          result.drop <- getTestEntry(corr.asrt, label = lab.drop, error.absent = FALSE)$action
          if (!is.null(result.drop)) 
            result1 <- result.drop
        }
      }
      
      #Determine if any correlations have been fitted
      corr.term <- setCorrTerm(corr.asrt$asreml.obj, spat.var = spat.var, 
                               asr4.2 = asr4.2)

      #### If no  correlation fitted and both rows and cols have corr funcs, try fitting them together
      if (!corr.term && (!any(rfuncs %in% get.specials("id.specials"))) 
          && allow.corrsJointFit)
      {
        if (trace) {cat("\n### Try joint correlation\n\n"); print(corr.asrt)}
        #Try first correl in current section
        ran.term2 <- paste0(rterms[1], ":", rterms[2])
        lab2 <- paste0("Try ", rterms[1], " and ", rterms[2])
        if (nsect > 1)
        {  
          ran.term2 <- paste0("at(", sections, ", '",stub, "'):", ran.term2)
          lab2 <- paste0(lab2, " for ", sections, " ",stub)
        }
        #Determine if spat.var is a fitted random term
        vpc.ran <- getSectionVpars(corr.asrt$asreml.obj, 
                                   sections = sections, stub = stub, 
                                   sterm.facs = sterms, 
                                   asr4.2 = asr4.2)$ran
        spat.term <- findterm(spat.var, names(vpc.ran)) #allows for changed order
        if (length(spat.term) == 1 && spat.term != 0) #have got a single spat.term in fit
          drop.spatvar <- names(spat.term)
        else
          drop.spatvar <- NULL
        tmp.asrt <- tryCatchLog(
          do.call(fitfunc, 
                  c(list(corr.asrt, label = lab2,
                         addRandom = ran.term2, dropRandom = drop.spatvar, 
                         maxit = maxit, 
                         allow.unconverged = FALSE, 
                         allow.fixedcorrelation = FALSE,
                         checkboundaryonly = TRUE, 
                         update = update, 
                         IClikelihood = IClikelihood, 
                         which.IC = which.IC), 
                    inargs)),
          error = function(e) {print(paste("Failed attempting to fit correlations to both dimensions;",
                                           "continued analysis without them")); NULL}, 
          include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
        
        #Check for singular (S) spatial terms - only change model if var == S can be made F and 
        #none of the correlations are S
        corr.asrt <- chk4SingularSpatTerms(tmp.asrt, corr.asrt, label = lab2, 
                                           sections = sections, stub = stub, 
                                           sterm.facs = sterms, 
                                           asr4 = asr4, asr4.2 = asr4.2, 
                                           maxit = maxit, 
                                           allow.unconverged = allow.unconverged, 
                                           allow.fixedcorrelation = allow.fixedcorrelation,
                                           checkboundaryonly = checkboundaryonly, 
                                           update = update, 
                                           IClikelihood = IClikelihood, 
                                           which.IC = which.IC,
                                           bounds.excl =  bounds.excl, 
                                           sing.excl = sing.excl)
        if (is.allnull(tmp.asrt))
        {
          test.summary <- addtoTestSummary(corr.asrt$test.summary, terms = lab2, 
                                           DF = NA, denDF = NA, p = NA, 
                                           AIC = NA, BIC = NA, 
                                           action = "Unchanged - Singular")
          corr.asrt$test.summary <- test.summary
        }
        
        result2 <- getTestEntry(corr.asrt, label = lab2)$action

        #Determine if any correlations have been fitted
        corr.term <- setCorrTerm(corr.asrt$asreml.obj, spat.var = spat.var, 
                                 asr4.2 = asr4.2)
        if (corr.term) #two-factor corr fitted
        { 
          #If corb and rorder == 0, try to fit corb up to order 10
          corr.lis <- do.call(fitCorbPlus1, 
                              c(list(corr.asrt, ran.term = ran.term2, rorder = rorders[1], 
                                     lab = lab2, result = result2, dimension = 1, 
                                     IClikelihood = IClikelihood), 
                                inargs))
          corr.asrt <- corr.lis$asrt
          if (ran.term2 != corr.lis$last.term)
          {
            lab2 <- corr.lis$last.lab
            result2 <- corr.lis$result
            ran.term2 <- corr.lis$last.term
          }

          #If corb and rorder == 0, try to fit corb up to order 10
          corr.lis <- do.call(fitCorbPlus1, 
                              c(list(corr.asrt, ran.term = ran.term2, rorder = rorders[2], 
                                     lab = lab2, result = result2, dimension = 2, 
                                     IClikelihood = IClikelihood), 
                                inargs))
          corr.asrt <- corr.lis$asrt
          if (ran.term2 != corr.lis$last.term)
          {
            lab2 <- corr.lis$last.lab
            result2 <- corr.lis$result
            ran.term2 <- corr.lis$last.term
          }
          last.term <- ran.term2
          #If all correlation model terms are bound, reinstate initial model
          lab.drop <- "Dropped correlations"
          corr.asrt <- allBoundSectionVpars(corr.asrt, init.asrt, 
                                            lab = lab.drop, #lab1, 
                                            sections = sections, stub = stub, 
                                            sterm.facs = sterms, 
                                            all.bounds.excl  = all.bounds.excl, 
                                            asr4.2 = asr4.2)
          result.drop <- getTestEntry(corr.asrt, label = lab.drop, error.absent = FALSE)$action
          if (!is.null(result.drop)) 
            result2 <- result.drop
        }
      }
    } #end of 2nd correl in current section
    
    #Determine if any correlations have been fitted
    corr.term <- setCorrTerm(corr.asrt$asreml.obj, spat.var = spat.var, 
                             asr4.2 = asr4.2)
    ##### Test for nugget variance, 
    # - only if the residual model is a variance model related to sections
    # - chooseOnIC is TRUE 
    # - if chooseOnIC and startFixRes are FALSE, then starting model has nugget variance
    # - if chooseOnIC is FALSE and startFixRes is TRUE, then need to try 
    # -    P for Residual bound provided nuggsOK is TRUE
    if ((chooseOnIC || (!chooseOnIC && startFixRes)) && corr.term && nuggsOK)
    {
      if (trace) {cat("\n#### Testing nugget variance\n\n"); print(corr.asrt)}
      #Get random and residual terms for correlation model for the current section
      vpc.corr <- getSectionVpars(corr.asrt$asreml.obj, 
                                  sections = sections, stub = stub, 
                                  sterm.facs = sterms, 
                                  asr4.2 = asr4.2)
      
      #Determine the correlation terms, if any
      vpt.corr <- getVpars(corr.asrt$asreml.obj, asr4.2)$vpt
      vpt.ran <- vpt.corr[names(vpc.corr$ran)]
      vpt.r <- vpt.ran[vpt.ran %in% c("R", "P", "C")]
      vpc.r <- vpc.corr$ran[names(vpt.r)]
      vpt.ran <- vpc.corr$ran[names(vpt.ran[vpt.ran == "G"])]
      #Is there (i) vpc.corr$res == NULL, implying multiple residual terms, none with 
      #  sections, or (ii) no ran terms, or (iii) all ran terms are bound
      if (is.null(vpc.corr$res) || length(vpc.r) == 0 || all(vpc.r %in% all.bounds.excl))
        nuggsOK <- FALSE
      #if have Fixed Residual, only if startFixRes then try positive 
      if (nuggsOK && length(vpc.corr$res) && 
          (startFixRes || (vpc.corr$res !=  "F" && 
                           (length(vpt.ran) == 1 && !(vpt.ran %in% all.bounds.excl)))))
      {
        #Try fixing either the  single residual variance term or that for the current section 
        tmp.asrt <- do.call(chgResTermBound, 
                            c(list(corr.asrt, sections = sections, stub = stub, 
                                   asr4 = asr4, asr4.2 = asr4.2, 
                                   fitfunc = fitfunc, 
                                   sterm.facs = sterms, vpc.res = vpc.corr$res, 
                                   maxit = maxit, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = TRUE, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC, bounds.excl =  bounds.excl), 
                              inargs))
        lasttest <- tail(tmp.asrt$test.summary, 1)
        if (grepl("nugget", lasttest$terms) && 
            (grepl("Unswapped", lasttest$action) || grepl("Unchanged", lasttest$action)))
          corr.asrt <- tmp.asrt #Change so have test.summary entry
        else
        { 
          new.vpc.corr <- getSectionVpars(tmp.asrt$asreml.obj, 
                                          sections = sections, stub = stub, 
                                          sterm.facs = sterms, 
                                          asr4.2 = asr4.2)
          bound.res <- ifelse(startFixRes, "P", "F")
          #Change residual to fixed if either 
          #     (i) all random correlation model terms are unbound 
          #   or (ii) none have changed
          n.new <- length(new.vpc.corr$ran)
          n.old <- length(vpc.corr$ran)
          if (new.vpc.corr$res == bound.res && 
              (n.new > 0 && (!any(new.vpc.corr$ran %in% bounds.excl) || 
                             (n.new == n.old && all(new.vpc.corr$ran == vpc.corr$ran)))))
            corr.asrt <- tmp.asrt
          else
          { 
            test.summary <- tmp.asrt$test.summary
            test.summary$action[nrow(test.summary)] <- "Unchanged residual"
            corr.asrt$test.summary <- test.summary
          }
        }
        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
      }
    } #end of nugget variance test

    
    #Check for either singular spatial variance ot residual variance
    corr.asrt <- chk4SingularSpatResVarTerms(corr.asrt, init.asrt, 
                                             corr.term = corr.term, 
                                             spat.var = spat.var, nuggsOK = nuggsOK, 
                                             sections = sections, stub = stub, 
                                             sterms = sterms, 
                                             allow.unconverged = allow.unconverged, 
                                             sing.excl = sing.excl, 
                                             bounds.excl = bounds.excl,
                                             all.bounds.excl = all.bounds.excl, 
                                             trace = trace, asr4.2 = asr4.2)

    #### Having made all model changes with checkboundaryonly = TRUE, 
    #### update for checkboundaryonly set to FALSE
    if (!checkboundaryonly)
    { 
      if (trace) {cat("\n#### Remove bound components\n\n"); print(corr.asrt)}
      corr.asrt <- do.call(rmboundary,
                           c(list(corr.asrt, checkboundaryonly = checkboundaryonly, 
                                  update = update, IClikelihood = IClikelihood), 
                             inargs))
      if (trace) 
      {cat("\n#### Have updated with checkboundary set to FALSE\n\n"); print(corr.asrt)}
    }
    
    #If all correlation model terms are bound, reinstate initial model
    lab.drop <- "Dropped correlations"
    corr.asrt <- allBoundSectionVpars(corr.asrt, init.asrt, 
                                      lab = lab.drop, #lab1, 
                                      sections = sections, stub = stub, 
                                      sterm.facs = sterms, 
                                      all.bounds.excl  = all.bounds.excl, 
                                      asr4.2 = asr4.2)
    result.drop <- getTestEntry(corr.asrt, label = lab.drop, error.absent = FALSE)$action
    if (!is.null(result.drop)) 
      result1 <- result.drop
    
    #### Further attempts to deal with bound random and residual terms when 
    # (i) checkboundary only is FALSE and (ii) there are correlation terms
    corr.term <- setCorrTerm(corr.asrt$asreml.obj, spat.var = spat.var, asr4.2 = asr4.2)
    if (corr.term && !checkboundaryonly)
    {
      if (trace) {cat("\n#### Further investigation of bound terms\n\n"); print(corr.asrt)}
      for (j in i:1)
      {
        if (nsect > 1)
        {
          stub <- levels(dat.in[[sections]])[j]
        } else
          stub <- NULL
        #Get random and residual terms for the current section
        vpc.ran <- getSectionVpars(corr.asrt$asreml.obj, which = "ran",
                                   sections = sections, stub = stub, 
                                   sterm.facs = sterms, 
                                   asr4.2 = asr4.2)$ran
        vpc.res <- getVpars(corr.asrt$asreml.obj, asr4.2)$vpc
        vpc.res <- vpc.res[grepl("!R$", names(vpc.res))]
        vpc.corr <- c(vpc.ran, vpc.res)
 
        if (any(unlist(vpc.corr) %in% all.bounds.excl))  #changed from bounds.excl
        {
          #Get vpc for this section
          vpc.corr <- getSectionVpars(corr.asrt$asreml.obj, 
                                      sections = sections, stub = stub, 
                                      sterm.facs = sterms, 
                                      asr4.2 = asr4.2)
          #### Only process if have bound residual and/or random corr model terms
          if (any(unlist(vpc.corr) %in% all.bounds.excl))
          {
            ### get bound random terms
            vpc.bran <- vpc.corr$ran[vpc.corr$ran %in% all.bounds.excl]
            if (length(vpc.bran) > 0)
            { 
              vpt.bran <- getVpars(corr.asrt$asreml.obj, asr4.2)$vpt[names(vpc.bran)]
              #If any random correlations bound, remove the corresponding correlation term
              if (any(vpt.bran %in% corr.types))
              {
                #Get bound corr vpars and remove
                vpc.bC <- vpc.bran[vpt.bran %in% corr.types] 
                corr.asrt <- rmboundCorrVpar(corr.asrt, vpcbound = names(vpc.bC), 
                                             maxit = maxit, 
                                             allow.unconverged = allow.unconverged,
                                             allow.fixedcorrelation = allow.fixedcorrelation,
                                             checkboundaryonly = FALSE,
                                             update = update,
                                             IClikelihood = IClikelihood,
                                             inargs = inargs)

                if (largeVparChange(corr.asrt$asreml.obj, 0.75))
                  corr.asrt <- iterate(corr.asrt)
                
                #update constraints for corr.bound under the new model
                vpc.corr <- getSectionVpars(corr.asrt$asreml.obj, 
                                            sections = sections, stub = stub, 
                                            sterm.facs = sterms, 
                                            asr4.2 = asr4.2)
              }
            }
          } #end vpars section
        } #end dealing with bound vpars
      } #end bounds within a sections
    } #end checking bounds
    
    if (largeVparChange(corr.asrt$asreml.obj, 0.75))
      corr.asrt <- iterate(corr.asrt)
    
    #Final check for either singular spatial variance o5 residual variance
    corr.asrt <- chk4SingularSpatResVarTerms(corr.asrt, init.asrt, 
                                             corr.term = corr.term, 
                                             spat.var = spat.var, nuggsOK = nuggsOK, 
                                             sections = sections, stub = stub, 
                                             sterms = sterms, 
                                             allow.unconverged = allow.unconverged, 
                                             sing.excl = sing.excl, 
                                             bounds.excl = bounds.excl,
                                             all.bounds.excl = all.bounds.excl, 
                                             trace = trace, asr4.2 = asr4.2)
  } #end of sections loop
  
  #Check for uncoverged analysis when allow.unconverged is FALSE
  if (!allow.unconverged && !corr.asrt$asreml.obj$converge)
    corr.asrt <- revert2previousFit(corr.asrt, init.asrt, 
                                    terms = "Unconverged spatial model", 
                                    action = "Revert to initial fit")
    
  if (trace) {cat("\n#### Exiting corr model fitting\n\n"); print(corr.asrt)}
  
  #Ensure setting of ai.sing is reinstated to the value on entry (for corb)
  asreml::asreml.options(ai.sing = ksing)

  return(corr.asrt)
}

fitTPNCSSMod <- function(asrtests.obj, sections = NULL, 
                         row.covar = "cRow", col.covar = "cCol", 
                         dropFixed = dropFixed, dropRandom = dropRandom, 
                         nsegs = NULL, 
                         allow.unconverged = allow.unconverged, 
                         allow.fixedcorrelation = allow.fixedcorrelation,
                         checkboundaryonly = checkboundaryonly, update = update, 
                         trace = trace, chooseOnIC = TRUE, maxit = 30, 
                         IClikelihood = "full", which.IC = "AIC", 
                         ...)
{ 
  #Check that named columns are in the data
  dat.in <- asrtests.obj$asreml.obj$call$data
  if (is.symbol(dat.in))
    dat.in <- eval(dat.in)
  checkNamesInData(c(sections, row.covar, col.covar), dat.in)
  
  
  #Check conformability of covars and factors
  # if (!is.null(dropRowterm) && nlevels(dat.in[[dropRowterm]]) != length(unique(dat.in[[row.covar]])))
  #   stop(dropRowterm, " does not have the same number of levels as there are values of ", row.covar)
  # if (!is.null(dropColterm) && nlevels(dat.in[[dropColterm]]) != length(unique(dat.in[[col.covar]])))
  #   stop(dropColterm, " does not have the same number of levels as there are values of ", col.covar)
  
  #Are dropRowterm and dropColterm already in the model?
  # facs <- c(dropRowterm, dropColterm)
  # drop.fix <- NULL
  # if (any(facs %in% rownames(asrtests.obj$wald.tab)))
  #   drop.fix <- paste(facs[facs %in% rownames(asrtests.obj$wald.tab)], collapse = " + ")
  # drop.ran <- NULL
  # if (any(facs %in% names(asrtests.obj$asreml.obj$vparameters)))
  #   drop.ran <- paste(facs[facs %in% names(asrtests.obj$asreml.obj$vparameters)], 
  #                     collapse = " + ")
  
  nsect <- calc.nsect(dat.in, sections)
  #Check dropFixed and dropRandom for length
  ndropF <- length(dropFixed)
  ndropR <- length(dropRandom)
  if (!is.null(dropFixed) && !(ndropF %in% c(1,nsect)))
    stop("The length of dropFixed must 1 or the number of levels in ", sections, " (",  nsect, ")")
  if (!is.null(dropRandom) && !(ndropR %in% c(1,nsect)))
    stop("The length of dropRandom must 1 or the number of levels in ", sections, " (",  nsect, ")")
  if (ndropF == 1 && nsect != 1) 
    dropFixed <- c(dropFixed, rep(NA, nsect - ndropF))
  if (ndropR == 1 && nsect != 1) 
    dropRandom <- c(dropRandom, rep(NA, nsect - ndropR))
  
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
  
  tspl.asrt <- asrtests.obj
  for (i in 1:nsect)
  {
    if (nsect > 1) 
    { 
      stub <- levels(dat.in[[sections]])[i]
      sect.fac <- paste0("at(", sections, ", '", stub, "'):")
      lab <- paste0("Try tensor NCS splines for ", sections, " ",stub)
    } else
      lab <- paste0("Try tensor NCS splines")
    #Fit TPNCSS to a section 
    drop.fix <- dropFixed[i]; if (!is.null(drop.fix) && is.na(drop.fix)) drop.fix <- NULL
    drop.ran <- dropRandom[i]; if (!is.null(drop.ran) && is.na(drop.ran)) drop.ran <- NULL
    if (chooseOnIC)
      tspl.asrt <- changeModelOnIC(tspl.asrt, 
                                   addFixed = fix.terms, 
                                   dropFixed = drop.fix[i], 
                                   dropRandom = drop.ran[i], 
                                   addRandom = paste(sect.fac, spl.terms, 
                                                     collapse = " + "), 
                                   allow.absentDropTerms = TRUE, 
                                   label = lab, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = checkboundaryonly, 
                                   update = update, 
                                   maxit = maxit, IClikelihood = IClikelihood, 
                                   which.IC = which.IC, 
                                   ...)
    else
      tspl.asrt <- changeTerms(tspl.asrt, 
                               addFixed = fix.terms, 
                               dropFixed = drop.fix[i], 
                               dropRandom = drop.ran[i], 
                               addRandom = paste(sect.fac, spl.terms, 
                                                 collapse = " + "), 
                               label = lab, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, 
                               maxit = maxit, IClikelihood = IClikelihood, ...)
    
  }
  
  #Final check for boundary terms
  if (!checkboundaryonly)
    tspl.asrt <- rmboundary(tspl.asrt, checkboundaryonly = checkboundaryonly, 
                            update = update, IClikelihood = IClikelihood)
  
  #Check for uncoverged analysis when allow.unconverged is FALSE
  if (!allow.unconverged && !tspl.asrt$asreml.obj$converge)
    corr.asrt <- revert2previousFit(tspl.asrt, asrtests.obj, 
                                    terms = "Unconverged spatial model", 
                                    action = "Revert to initial fit")
  
  return(tspl.asrt)
}

#Creates a list with the tpps bits for asreml.option = "grp"
addPSdesign.mat <- function(dat, sections = NULL, nsect = 1, 
                            row.coords, col.coords, 
                            nsegs = NULL, nestorder = c(1, 1), 
                            degree = c(3,3), difforder = c(2,2), 
                            rotateX = FALSE, theta = c(0, 0), 
                            asreml.opt = "grp", stub = "xx", 
                            ...)
{
  tpsmmb.args <- getTpsmmb.args(list(...))
  
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
                           stub <- data[[sects]][1]
                           XZ.mat <- tpsmmb(columncoordinates = columncoordinates, 
                                            rowcoordinates = rowcoordinates, 
                                            data = data, 
                                            stub = stub, nsegments = nsegments, 
                                            nestorder = nestorder, 
                                            degree = degree, difforder = difforder,
                                            rotateX = rotateX, theta = theta, 
                                            asreml = asreml.opt, ...)
                           return(XZ.mat)
                         }, columncoordinates = col.coords, rowcoordinates = row.coords, 
                         sects = sections, nsegments = nsegs, 
                         asreml.opt = asreml.opt)
  } else
  {
    
    if (all(sapply(nsegs, is.null)))
      nsegs <- c(length(unique(dat[[col.coords]]))-1,
                 length(unique(dat[[row.coords]]))-1)
    tps.XZmat <- list(
      do.call(tpsmmb, 
              c(list(columncoordinates = col.coords, 
                     rowcoordinates = row.coords, 
                     data = dat, 
                     stub = stub, nsegments = nsegs, 
                     nestorder = nestorder, 
                     degree = degree, difforder = difforder,
                     rotateX = rotateX, theta = theta, 
                     asreml = asreml.opt), 
                tpsmmb.args)))
  }
  attr(tps.XZmat, which = "nsegs") <- nsegs
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
    start <- tps.XZmat[[1]]$grp[[1]][1]
    new.grps.lens <- sapply(grp.names, 
                            function(agrp, tps.XZmat)
                              max.len <- max(sapply(tps.XZmat, function(mat) length(mat$grp[[agrp]]))),
                            tps.XZmat = tps.XZmat)
    ngrps <- length(new.grps.lens)
    if (sum(new.grps.lens[1:(ngrps-1)]) != new.grps.lens[ngrps])
      stop("The number of columns in the all group does not match the sum of the other groups")
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
#'
#' * The arguments degree, difforder and theta are in the order column followed by row dimension. 
makeTPPSplineMats.data.frame <- function(data, sections = NULL, 
                                         row.covar, col.covar, 
                                         nsegs = NULL, nestorder = c(1, 1), 
                                         degree = c(3,3), difforder = c(2,2), 
                                         rotateX = FALSE, theta = c(0, 0), 
                                         asreml.option = "grp", mbf.env = sys.frame(), 
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
  
  #Check nsegs, nestorder, degree, difforder and ngridangles
  if (length(nsegs) != 2 && !is.null(nsegs))
    stop("nsegs must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(nestorder) != 2)
    stop("nestorder must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(degree) != 2)
    stop("degree must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(difforder) != 2)
    stop("difforder must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(theta) != 2)
    stop("theta must specify exactly 2 values, one for each of the column and row dimensions")

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

  rotateX <- rotateX
  theta <- theta
  tps.XZmat <- addPSdesign.mat(dat.rc, sections = sections, nsect = nsect, 
                               row.coords = row.covar, col.coords = col.covar, 
                               nsegs = nsegs, nestorder = nestorder, 
                               degree = degree, difforder = difforder,
                               rotateX = rotateX, theta = theta, 
                               asreml.opt = asreml.opt, stub = "xx", ...)
  #Get data for mbf 
  if (asreml.opt == "mbf")
  {
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
    dat <- suppressMessages(dplyr::left_join(data, dat))
    dat <- dat[c(rc.cols, setdiff(names(dat), rc.cols))]
    if (!is.null(mbf.env)) #assign the spline-basis data.frames
    { 
      attr(dat, which = "mbf.env") <- mbf.env
      if (nsect != 1)
      {  
        #Put the mbf data.frames into the mbf.env
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
        Zmat.names <- paste0(paste0(c("BcZ", "BrZ", "BcrZ"),stub), ".df")
        if (any(sapply(Zmat.names, exists, envir = parent.frame(1))))
          warning("The following objects are being overwritten: ", 
                  paste(Zmat.names[sapply(Zmat.names, exists, envir = parent.frame(2))], 
                        collapse = ", "))
        assign(Zmat.names[1], tps.XZmat[[1]]$BcZ.df, envir = mbf.env)
        assign(Zmat.names[2], tps.XZmat[[1]]$BrZ.df, envir = mbf.env)
        assign(Zmat.names[3], tps.XZmat[[1]]$BcrZ.df, envir = mbf.env)
      }
    }
  } else #doing grp
  {  
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
      dat <- suppressMessages(dplyr::left_join(tmp, dat, by = rc.cols))
      dat <- dat[order(dat$.row.ord), ]
      dat <- dat[, -match(".row.ord", names(dat))]
    } else    
      dat <- suppressMessages(dplyr::left_join(data, dat, by = rc.cols))
    dat$TP.col <- factor(dat$TP.col)
    dat$TP.row <- factor(dat$TP.row)
    dat$TP.CxR <- factor(dat$TP.CxR)
    
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
  attr(tps.XZmat, which = "nsect") <- nsect
  attr(tps.XZmat, which = "nestorder") <- nestorder
  attr(tps.XZmat, which = "degree") <- degree
  attr(tps.XZmat, which = "difforder") <- difforder
  attr(tps.XZmat, which = "theta_opt") <- theta

  return(tps.XZmat)
}

fitTPSModSect <- function(tspl.asrt, data, mat, ksect, sect.fac, 
                          drop.fix, drop.ran, 
                          sections = NULL, 
                          row.covar, col.covar, lab, 
                          nsegs = NULL, nestorder = c(1, 1), 
                          degree = c(3,3), difforder = c(2,2), 
                          rotateX = FALSE, theta = c(0.0), 
                          usRandLinCoeffs = TRUE, 
                          asreml.opt = "grp", stub = "xx", 
                          allow.unconverged = allow.unconverged, 
                          allow.fixedcorrelation = allow.fixedcorrelation,
                          checkboundaryonly = checkboundaryonly, update = update, 
                          trace = trace, chooseOnIC = TRUE, maxit = 30, 
                          IClikelihood = "full", which.IC = "AIC", ...)
{
  inargs <- list(...)
  
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  sing.excl <- c("B","S","?")
  dorotate <- rotateX && any(difforder == 2)
  
  #Determine terms specified by dropFixed and dropRandom to remove from the model?
  if (!is.null(drop.fix) && is.na(drop.fix)) drop.fix <- NULL
  if (!is.null(drop.ran) && is.na(drop.ran)) drop.ran <- NULL
  
  nfixterms <- difforder[1] * difforder[2] 
  if (nfixterms > 1)
    fix.ch <- paste(paste0(sect.fac, paste0("TP.CR.", 2:nfixterms)), collapse = " + ")
  else
    fix.ch <- NULL

  if (chooseOnIC)
  { 
    fitfunc <- changeModelOnIC
    allow.absentDropTerms <- TRUE
  }
  else
  { 
    fitfunc <- changeTerms
    allow.absentDropTerms <- NULL
  }
  
  init.asrt <- tspl.asrt
  theta.opt <- theta
  if (dorotate)
  { 
    labunrot <- gsub("tensor", "unrotated tensor", lab)
    labrot <- gsub("tensor", 
                   paste0("rotated ", paste(round(theta.opt,0.5), collapse = ","), 
                          " tensor"), lab)
  }
  else
  { 
    labunrot <- lab
  }
  
  if (asreml.opt == "mbf")
  {
    #Set the mbf.env in asreml.obj to the current environment
    mbf.env <- sys.frame()
    asreml.obj <- tspl.asrt$asreml.obj
    asreml.obj <- setmbfenv(asreml.obj, dat = asreml.obj$call$data, mbf.env = mbf.env)
    
    #Assign basis data.frames to the current environment
    Zmat.names <- paste0(paste0(c("BcZ", "BrZ", "BcrZ"), stub), ".df")
    if (any(sapply(Zmat.names, exists, envir = mbf.env)))
      warning("THe following objects are being overwritten: ", 
              paste(Zmat.names[sapply(Zmat.names, exists, envir = parent.frame(2))], 
                    collapse = ", "))
    assign(Zmat.names[1], mat$BcZ.df, envir = mbf.env)
    assign(Zmat.names[2], mat$BrZ.df, envir = mbf.env)
    assign(Zmat.names[3], mat$BcrZ.df, envir = mbf.env)
    
    mbf.lis <- mat$mbflist
    
    ran.ch <- paste(paste0(sect.fac,  
                           c(paste0("TP.C.",1:difforder[1],":mbf(TP.row)"), 
                             paste0("TP.R.",1:difforder[2],":mbf(TP.col)"), 
                             "mbf(TP.CxR)", 
                             paste0("dev(",row.covar,")"), 
                             paste0("dev(",col.covar,")"))), 
                    collapse = " + ")

    #Fit the full P-spline model, without rotation
    #tspl.asrt.unrot 
    tspl.asrt <- do.call(fitfunc, 
                         args = c(list(tspl.asrt,
                                       addFixed = fix.ch,
                                       dropFixed = drop.fix, 
                                       addRandom = ran.ch,
                                       dropRandom = drop.ran, 
                                       allow.absentDropTerms = allow.absentDropTerms, 
                                       mbf = mbf.lis,
                                       label = labunrot,
                                       allow.unconverged = allow.unconverged, 
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = checkboundaryonly, 
                                       update = update, 
                                       maxit = maxit, 
                                       IClikelihood = IClikelihood, 
                                       which.IC = which.IC), 
                                  inargs))

    #Fit the model using the optimal theta for rotating the penalty eigenvectors
    if (dorotate)
    {
      #Fit the P-splines for the optimal theta
      rm(list = Zmat.names, envir = mbf.env)
      mat <- makeTPPSplineMats(data, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar,
                               nsegs = nsegs, nestorder = nestorder,
                               degree = degree, difforder = difforder,
                               rotateX = rotateX, theta = theta.opt, 
                               asreml.opt = "grp", mbf.env = NULL)[[ksect]]
      if (any(sapply(Zmat.names, exists, envir = mbf.env)))
        warning("THe following objects are being overwritten: ", 
                paste(Zmat.names[sapply(Zmat.names, exists, envir = parent.frame(2))], 
                      collapse = ", "))
      assign(Zmat.names[1], mat$BcZ.df, envir = mbf.env)
      assign(Zmat.names[2], mat$BrZ.df, envir = mbf.env)
      assign(Zmat.names[3], mat$BcrZ.df, envir = mbf.env)
      mbf.lis <- mat$mbflist
      dat <- mat$data.plus
      if (chooseOnIC)
      { 
        tspl.asrt <- updateOnIC.asrtests(tspl.asrt, data = dat, 
                                         mbf = mbf.lis, maxit = maxit, 
                                         label = labrot, 
                                         maxit = maxit, 
                                         IClikelihood = IClikelihood, 
                                         which.IC = which.IC)
        if (grepl("Unswapped", getTestEntry(tspl.asrt, label = labrot)$action))
          theta.opt <- c(0,0)
      } else
        tspl.asrt <- update.asrtests(tspl.asrt, data = mat$data.plus, 
                                     mbf = mbf.lis, maxit = maxit, 
                                     label = labrot, IClikelihood = IClikelihood)
      #Check criteria
      # print(infoCriteria(list(old = tspl.asrt$asreml.obj, new = new.asr), IClikelihood = "full"))
      # new.dev <- deviance.asr(new.asr) 
      # print(c(dev, new.dev))
    }
  } else #grp
  {    
    grp <- mat$grp
    #Following needed to ensure that group information is present in the asreml.obj
    tspl.asrt$asreml.obj <-  newfit(tspl.asrt$asreml.obj, group = grp)

    ran.ch <- paste(paste0(sect.fac,  
                           c(paste0("grp(TP.C.",1:difforder[1],"_frow)"), 
                             paste0("grp(TP.R.",1:difforder[2],"_fcol)"), 
                             "grp(TP_fcol_frow)", 
                             paste0("dev(",row.covar,")"), 
                             paste0("dev(",col.covar,")"))), 
                    collapse = " + ")
    
    #Fit the full P-spline model, without rotation
    #tspl.asrt.unrot 
    tspl.asrt <- do.call(fitfunc, 
                         args = c(list(tspl.asrt, 
                                       addFixed = fix.ch,
                                       dropFixed = drop.fix, 
                                       addRandom = ran.ch,
                                       dropRandom = drop.ran, 
                                       allow.absentDropTerms = allow.absentDropTerms, 
                                       group = grp,
                                       label = labunrot, 
                                       allow.unconverged = allow.unconverged, 
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = checkboundaryonly, 
                                       update = update, 
                                       maxit = maxit, 
                                       IClikelihood = IClikelihood, 
                                       which.IC = which.IC), 
                                  inargs))

    #Fit the model using the optimal theta for rotating the penalty eigenvectors
    if (dorotate)
    {
      #Fit the P-splines for the optimal theta
      mat <- makeTPPSplineMats(data, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar,
                               nsegs = nsegs, nestorder = nestorder,
                               degree = degree, difforder = difforder,
                               rotateX = rotateX, theta = theta.opt, 
                               asreml.opt = "grp")[[ksect]]
      grp <- mat$grp
      #Following needed to ensure that group information is present in the asreml.obj
      tspl.asrt$asreml.obj <-  newfit(tspl.asrt$asreml.obj, group = grp)
      if (chooseOnIC)
      { 
        tspl.asrt <- updateOnIC.asrtests(tspl.asrt, data = mat$data.plus, 
                                         group = grp, maxit = maxit, 
                                         label = labrot, IClikelihood = IClikelihood, 
                                         which.IC = which.IC)
        if (grepl("Unswapped", getTestEntry(tspl.asrt, label = labrot)$action))
          theta <- c(0,0)
      } else
        tspl.asrt <- update.asrtests(tspl.asrt, data = mat$data.plus, 
                                     group = grp, maxit = maxit, 
                                     label = labrot, IClikelihood = IClikelihood)
    }
  }
  
  #Prepare for fitting unstructured model to row and col marginal terms
  if (usRandLinCoeffs)
  { 
    vpars.all <- names(tspl.asrt$asreml.obj$vparameters)
    repln <- 1
    
    #Do not allow singularities in this section of the code
    # ksing <-   get("asr_options", envir = getFromNamespace(".asremlEnv", "asreml"))$ai.sing
    # asreml::asreml.options(ai.sing = FALSE)
    
    if (!is.null(sections))
    {   
      klev <- levels(data[[sections]])[ksect]
      vpars.all <- vpars.all[grepl(klev, vpars.all)]
      if (!asr4.2)
        vpars.all <- gsub(paste0("'",klev,"'"), ksect, vpars.all)
      repln <- length(levels(data[[sections]]))
    }
    
    #If more than one col variable in the marginal random row term in this section, try unstructured model
    rowmarg.vpar <- vpars.all[grepl("TP\\.C\\.", vpars.all)]
    nr <- length(rowmarg.vpar)
    if (nr > 1)
    {
      us.func <- ifelse(nr > 2, "corgh", "corh")
      drop.ran <-paste(rowmarg.vpar, collapse = " + ") 
      add.ran <- paste0("str( ~ ", drop.ran, ", ~ ", 
                        us.func, "(", length(rowmarg.vpar), "):id(", mat$dim['nbr']*repln,"))")
      lab.r <- "Try column-parameters covariance for random row terms"
      if (asreml.opt == "mbf")
        tmp.asrt <- do.call(fitfunc, 
                            args = c(list(tspl.asrt, 
                                          addRandom = add.ran,
                                          dropRandom = drop.ran, 
                                          mbf = mbf.lis,
                                          label = lab.r, 
                                          allow.unconverged = allow.unconverged, 
                                          allow.fixedcorrelation = allow.fixedcorrelation,
                                          checkboundaryonly = TRUE, #remove using bespoke code
                                          update = FALSE, #to ensure clean refit
                                          maxit = maxit, 
                                          IClikelihood = IClikelihood, 
                                          which.IC = which.IC), 
                                     inargs))
      else
        tmp.asrt <- do.call(fitfunc, 
                            args = c(list(tspl.asrt, 
                                          addRandom = add.ran,
                                          dropRandom = drop.ran, 
                                          group = grp,
                                          label = lab.r, 
                                          allow.unconverged = allow.unconverged, 
                                          allow.fixedcorrelation = allow.fixedcorrelation,
                                          checkboundaryonly = TRUE, #remove using bespoke code
                                          update = FALSE, #to ensure clean refit
                                          maxit = maxit, 
                                          IClikelihood = IClikelihood, 
                                          which.IC = which.IC), 
                                     inargs))
      #Check that no marginal random col parameters are bound and, if they are, remove all corh (corgh) parameters
      result <- getTestEntry(tmp.asrt, label = lab.r)
      if (!grepl("Unswapped", result$action) && !grepl("Unchanged", result$action))
      {
        vpars <- getVpars(tmp.asrt$asreml.obj, asr4.2 = asr4.2)
        vpc <- vpars$vpc
        names(vpc) <- gsub('\"', "\'", names(vpc))
        vpc.col <- names(vpc)[grepl(gsub(" \\+ ", "+", drop.ran), names(vpc), fixed = TRUE)]
        if (length(vpc.col) > 0)
        {
          vpc.bound <- vpc[vpc.col]
          if (any(vpc.bound %in% sing.excl))
          {
            test.summary <- addtoTestSummary(tmp.asrt$test.summary, terms = drop.ran, 
                                             DF = result$DF, denDF = NA, p = NA, 
                                             AIC = result$AIC, BIC = result$BIC, 
                                             action = "Unchanged - Boundary")
            tspl.asrt$test.summary <- test.summary
          } else
            tspl.asrt <- tmp.asrt #no bound terms
        } else
          tspl.asrt <- tmp.asrt #no terms found
      } else
        tspl.asrt <- tmp.asrt #model remained unchanged
    }
    
    #If more than one row variable in the marginal random col term in this section, try unstructured model
    vpars.all <- vpars.all[vpars.all %in% names(tspl.asrt$asreml.obj$vparameters)]
    colmarg.vpar <- vpars.all[grepl("TP\\.R\\.", vpars.all)]
    nc <- length(colmarg.vpar)
    lab.c <- "Try row-parameters covariance for random column terms"
    if (length(colmarg.vpar) > 1)
    {
      us.func <- ifelse(nc > 2, "corgh", "corh")
      drop.ran <-paste(colmarg.vpar, collapse = " + ") 
      add.ran <- paste0("str( ~ ", drop.ran, 
                        ", ~ ", us.func, "(", nc, "):id(", mat$dim['nbc']*repln,"))")
      if (asreml.opt == "mbf")
        tmp.asrt <- do.call(fitfunc, 
                            args = c(list(tspl.asrt, 
                                          addRandom = add.ran,
                                          dropRandom = drop.ran, 
                                          mbf = mbf.lis,
                                          label = lab.c, 
                                          allow.unconverged = allow.unconverged, 
                                          allow.fixedcorrelation = allow.fixedcorrelation,
                                          checkboundaryonly = TRUE, #remove using bespoke code
                                          update = FALSE, #to ensure clean refit
                                          maxit = maxit, 
                                          IClikelihood = IClikelihood, 
                                          which.IC = which.IC), 
                                     inargs))
      else
        tmp.asrt <- do.call(fitfunc, 
                            args = c(list(tspl.asrt, 
                                          addRandom = add.ran,
                                          dropRandom = drop.ran, 
                                          group = grp,
                                          label = lab.c, 
                                          allow.unconverged = allow.unconverged, 
                                          allow.fixedcorrelation = allow.fixedcorrelation,
                                          checkboundaryonly = TRUE, #remove using bespoke code
                                          update = FALSE, #to ensure clean refit
                                          maxit = maxit, 
                                          IClikelihood = IClikelihood, 
                                          which.IC = which.IC), 
                                     inargs))
      #Check that no marginal random col parameters are bound and, if they are, remove all corh (corgh) parameters
      result <- getTestEntry(tmp.asrt, label = lab.c)
      if (!grepl("Unswapped", result$action) && !grepl("Unchanged", result$action))
      {
        vpars <- getVpars(tmp.asrt$asreml.obj, asr4.2 = asr4.2)
        vpc <- vpars$vpc
        names(vpc) <- gsub('\"', "\'", names(vpc))
        vpc.col <- names(vpc)[grepl(gsub(" \\+ ", "+", drop.ran), names(vpc), fixed = TRUE)]
        if (length(vpc.col) > 0)
        {
          vpc.bound <- vpc[vpc.col]
          if (any(vpc.bound %in% sing.excl))
          {
            test.summary <- addtoTestSummary(tmp.asrt$test.summary, terms = drop.ran, 
                                             DF = result$DF, denDF = NA, p = NA, 
                                             AIC = result$AIC, BIC = result$BIC, 
                                             action = "Unchanged - Boundary")
            tspl.asrt$test.summary <- test.summary
          } else
            tspl.asrt <- tmp.asrt #no bound terms
        } else
          tspl.asrt <- tmp.asrt #no terms found
      } else
        tspl.asrt <- tmp.asrt #model remained unchanged
      
      #Ensure setting of ai.sing is reinstated to the value on entry
      # asreml::asreml.options(ai.sing = ksing)
    }
  }
  
  #Check for singular fixed P-spline terms and, if any, remove
  
  asreml.obj <- tspl.asrt$asreml.obj
  coefF <- summary(asreml.obj, coef=TRUE)$coef.fixed
  which.cF <- rownames(coefF)[is.na(coefF[, "z.ratio"])]
  if (any(grepl("TP.CR.", which.cF)))
  {
    which.cF <- which.cF[grepl("TP.CR.", which.cF)]
    asreml.obj <- newfit(asreml.obj, 
                         fixed. = as.formula(paste("~ . - ", 
                                                   paste(which.cF, collapse = " - "))))
    tspl.asrt$asreml.obj <- asreml.obj
    test.summary <- addtoTestSummary(tspl.asrt$test.summary, terms = "Aliased fixed spline term(s)", 
                                     DF = NA, denDF = NA, p = NA, 
                                     AIC = NA, BIC = NA, 
                                     action = "Dropped")
    tspl.asrt$test.summary <- test.summary
  }
  
  tspl.asrt <- rmboundary(tspl.asrt, checkboundaryonly = checkboundaryonly, 
                          update = update, IClikelihood = IClikelihood)
  attr(tspl.asrt$asreml.obj, which = "theta.opt") <- theta.opt

  return(tspl.asrt)
}

#Fit a tensor-spline spatial model
fitTPPSMod <- function(asrtests.obj, sections = NULL, 
                       row.covar = "cRow", col.covar = "cCol", 
                       dropFixed = NULL, dropRandom = NULL, 
                       nsegs = NULL, nestorder = c(1, 1), 
                       degree = c(3,3), difforder = c(2,2), 
                       usRandLinCoeffs = TRUE, 
                       rotateX = FALSE, ngridangles = NULL,
                       which.rotacriterion = "AIC", nrotacores = 1, 
                       asreml.opt = "grp", 
                       tpps4mbf.obj = NULL, 
                       allow.unconverged = allow.unconverged, 
                       allow.fixedcorrelation = allow.fixedcorrelation,
                       checkboundaryonly = checkboundaryonly, update = update, 
                       trace = trace, chooseOnIC = TRUE, maxit = 30, 
                       IClikelihood = "full", which.IC = "AIC",
                       ...)
{ 
  inargs <- list(...)
  tpsmmb.args <- checkEllipsisArgs_tpsmmb("makeTPPSplineMats.data.frame", inargs)

  #Check which.criterion options
  options <- c("deviance", "likelihood", "AIC", "BIC")
  which.rotacriterion <- options[check.arg.values(which.rotacriterion, options)]
  
  #Stop parallel processing for mbf
  if (asreml.opt == "mbf" && nrotacores > 1)
    stop(paste("Parallel processing has not been implemented for asreml.option set to mbf;",
               "nrotacores must be one"))
  
  if (rotateX && !is.null(tpps4mbf.obj))
    warning("The supplied tpps4mbf.obj will not be appropriate for model fitted using a rotated penalty matrix")

  #Check nsegs, nestorder, degree, difforder and ngridangles
  if (length(nsegs) != 2 && !is.null(nsegs))
    stop("nsegs must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(nestorder) != 2)
    stop("nestorder must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(degree) != 2)
    stop("degree must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(difforder) != 2)
    stop("difforder must specify exactly 2 values, one for each of the column and row dimensions")
  if (rotateX && any(difforder != 2))
    stop("Rotation of the penalty matrix is only implements for difforder equal to 2")
  if (!any(is.null(ngridangles)) && length(ngridangles) != 2)
    stop(paste("ngridangles must ether be NULL or specify exactly 2 values,",
               "one for each of the column and row dimensions"))

  #Get the data from the original call and check that named columns are in the data
  dat.in <- asrtests.obj$asreml.obj$call$data
  if (is.symbol(dat.in))
    dat.in <- eval(dat.in)
  checkNamesInData(c(sections, row.covar, col.covar), dat.in)
  #Check conformability of covars and factors
  # if (!is.null(dropRowterm) && nlevels(dat.in[[dropRowterm]]) != length(unique(dat.in[[row.covar]])))
  #   stop(dropRowterm, " does not have the same number of levels as there are values of ", row.covar)
  # if (!is.null(dropColterm) && nlevels(dat.in[[dropColterm]]) != length(unique(dat.in[[col.covar]])))
  #   stop(dropColterm, " does not have the same number of levels as there are values of ", col.covar)
  
  if (is.null(tpps4mbf.obj))
  {   
    #Create spline basis functions
    #do not set to NULL so that the mbf df will be assigned in the current environment
    # - needed for the rmboundary call at the end
    tps.XZmat <- do.call(makeTPPSplineMats, 
                         c(list(dat.in, sections = sections, 
                                row.covar = row.covar, col.covar = col.covar,
                                nsegs = nsegs, nestorder = nestorder,
                                degree = degree, difforder = difforder, 
                                rotateX = rotateX,
                                asreml.opt = asreml.opt), 
                           tpsmmb.args))
  }
  else #user supplied
    tps.XZmat <- tpps4mbf.obj
  dat <- tps.XZmat[[1]]$data.plus

  #Update the asreml.obj for the new data.frame
  asreml.obj  <- asrtests.obj$asreml.obj
  asreml.obj <- newfit(asreml.obj, data = dat)
  tspl.asrt <- as.asrtests(asreml.obj = asreml.obj, NULL, 
                           test.summary = asrtests.obj$test.summary, 
                           IClikelihood = "full", 
                           label = "Change to new data.frame with TPS bits")

  #Prepare for model fitting
  nsect <- calc.nsect(dat.in, sections)
  #Check dropFixed and dropRandom for length
  ndropF <- length(dropFixed)
  ndropR <- length(dropRandom)
  if (!is.null(dropFixed) && !(ndropF %in% c(1,nsect)))
    stop("The length of dropFixed must 1 or the number of levels in ", sections, " (",  nsect, ")")
  if (!is.null(dropRandom) && !(ndropR %in% c(1,nsect)))
    stop("The length of dropRandom must 1 or the number of levels in ", sections, " (",  nsect, ")")
  if (ndropF == 1 && nsect != 1) 
    dropFixed <- c(dropFixed, rep(NA, nsect - ndropF))
  if (ndropR == 1 && nsect != 1) 
    dropRandom <- c(dropRandom, rep(NA, nsect - ndropR))
  
  rotated <- rotateX & any(difforder > 1)
  if (rotated) 
    theta.opt <- getRotationThetas(tspl.asrt, data = dat.in, mat = tps.XZmat, sections = sections, 
                                   dropFixed = dropFixed, dropRandom = dropRandom, 
                                   row.covar = row.covar, col.covar = col.covar, 
                                   nsegs = nsegs, nestorder = nestorder, 
                                   degree = degree, difforder = difforder, 
                                   rotateX = rotateX, ngridangles =ngridangles, 
                                   which.rotacriterion = which.rotacriterion, 
                                   nrotacores = nrotacores, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = checkboundaryonly, 
                                   asreml.opt = asreml.opt, 
                                   update = update, 
                                   maxit = maxit, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC)
  else
  {
    theta.opt <- lapply(1:nsect, function(ksect) c(0,0))
  }

  #Fit spatial TPPS to sections
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
      sect.fac <- paste0("at(", sections, ",  '", stub, "'):")
      lab <- paste0("Try tensor P-splines for ", sections, " ",stub)
    }
    if (length(inargs))
      other.args <- inargs[setdiff(names(inargs), names(getTpsmmb.args(inargs)))]
    else 
      other.args <- NULL
    tspl.asrt <- do.call(fitTPSModSect, 
                         c(list(tspl.asrt, data = dat.in, mat = tps.XZmat[[i]], 
                                ksect = i, sect.fac = sect.fac, 
                                drop.fix = dropFixed[i], drop.ran = dropRandom[i], 
                                sections = sections, 
                                row.covar = row.covar, col.covar = col.covar, 
                                nsegs = nsegs, nestorder = nestorder, 
                                degree = degree, difforder = difforder, 
                                usRandLinCoeffs = usRandLinCoeffs,
                                rotateX = rotateX, theta = theta.opt[[i]], 
                                lab = lab, asreml.opt = asreml.opt, stub = stub, 
                                allow.unconverged = allow.unconverged, 
                                allow.fixedcorrelation = allow.fixedcorrelation,
                                checkboundaryonly = checkboundaryonly, 
                                update = update, 
                                maxit = maxit, 
                                chooseOnIC = chooseOnIC, 
                                IClikelihood = IClikelihood, 
                                which.IC = which.IC), 
                           other.args))
  }
  
  #Set the mbf.env of the asreml.obj in tspl.asrt to the current environment
  asreml.obj <- tspl.asrt$asreml.obj
  asreml.obj <- setmbfenv(asreml.obj, dat = asreml.obj$call$data)
  tspl.asrt$asreml.obj <- asreml.obj
  
  #Final check for boundary terms
  if (!checkboundaryonly)
    tspl.asrt <- rmboundary(tspl.asrt, checkboundaryonly = checkboundaryonly, 
                            update = update, IClikelihood = IClikelihood)
  
  #Add theta.opt attribute
  if (rotated)
  {
    if (nsect > 1)
      names(theta.opt) <- levels(dat.in[[sections]])
    attr(tspl.asrt$asreml.obj, which = "theta.opt") <- theta.opt
  }
  
  #Check for uncoverged analysis when allow.unconverged is FALSE
  if (!allow.unconverged && !tspl.asrt$asreml.obj$converge)
    corr.asrt <- revert2previousFit(tspl.asrt, asrtests.obj, 
                                    terms = "Unconverged spatial model", 
                                    action = "Revert to initial fit")
  
  return(tspl.asrt)
}

