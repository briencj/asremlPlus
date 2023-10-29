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

checkSections4RanResTerms <- function(asrtests.obj, sections = NULL, asr4)
{
  #Check that separate residual terms have been fitted for multiple sections 
  if (!is.null(sections))
  {
    ran <- getFormulae(asrtests.obj$asreml.obj)$random
    if (!is.null(ran))
    { 
      ran <- getTerms.formula(ran)
      ran <- ran[ran != sections] #Remove sections term, if present
      if (!any((grepl("idh\\(", ran) | grepl("at\\(", ran)) & grepl(sections, ran)))
        warning(paste0("There are no 'idh' or 'at' terms for ", sections, 
                       " in the supplied random model; the fitting may be improved ",
                       "if separate random block terms are specified for each section"))
    }
    if (asr4)
    {  
      res <- getFormulae(asrtests.obj$asreml.obj)$residual
      if (!is.null(res))
      { 
        res <- getTerms.formula(res)
        if (!any((grepl("idh\\(", res) | grepl("dsum\\(", res)) & grepl(sections, res)))
          warning(paste0("There are no 'idh' or 'dsum' terms involving ", sections, 
                         " in the supplied residual model; the fitting may be improved ",
                         "if separate units terms are specified for each section"))
      }
    }  
    else
    {
      res <- getFormulae(asrtests.obj$asreml.obj)$rcov
      if (!is.null(res))
      { 
        res <- getTerms.formula(res)
        if (!any((grepl("idh\\(", res) | grepl("at\\(", res)) & grepl(sections, res)))
          warning(paste0("There are no 'idh' or 'at' terms for ", sections, 
                         " in the supplied rcov model; the fitting may be improved ",
                         "if separate units terms are specified for each section"))
      }
    }
  }

  invisible()
}

addSpatialModel.asrtests <- function(asrtests.obj, spatial.model = "TPPS", 
                                     sections = NULL, 
                                     row.covar = "cRow", col.covar = "cCol", 
                                     row.factor = "Row", col.factor = "Col", 
                                     corr.funcs = c("ar1", "ar1"), 
                                     row.corrFitfirst = TRUE, allow.corrsJointFit = TRUE, 
                                     dropRowterm = NULL, dropColterm = NULL, 
                                     nsegs = NULL, nestorder = c(1, 1), 
                                     degree = c(3,3), difforder = c(2,2), 
                                     rotateX = FALSE, ngridangles = c(18, 18), 
                                     which.rotacriterion = "AIC", nrotacores = 1, 
                                     asreml.option = "mbf", tpps4mbf.obj = NULL,  
                                     allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                                     checkboundaryonly = FALSE, update = FALSE, 
                                     maxit = 30, IClikelihood = "full", ...)
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
 
  #Check if have separate section random and residual terms
  checkSections4RanResTerms(asrtests.obj, sections = sections, asr4 = asr4)
   
  #Check IClikelihood options
  options <- c("REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  
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
                               allow.corrsJointFit = allow.corrsJointFit, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, chooseOnIC = FALSE, 
                               maxit = maxit, IClikelihood = ic.lik, ...)
  #Fit a local spatial model involving TPNCSS
  if ("TPNCSS" %in% spatial.mod)
    spatial.asrt <- fitTPNCSSMod(asrtests.obj, sections = sections, 
                                 row.covar = row.covar, col.covar = col.covar, 
                                 dropRowterm = dropRowterm, dropColterm = dropColterm, 
                                 allow.unconverged = allow.unconverged, 
                                 allow.fixedcorrelation = allow.fixedcorrelation,
                                 checkboundaryonly = checkboundaryonly, 
                                 update = update, chooseOnIC = FALSE, 
                                 maxit = maxit, IClikelihood = "full", ...)
  
  #Fit a residual spatial model involving TPPS
  if ("TPPS" %in% spatial.mod)
    spatial.asrt <- fitTPPSMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               dropRowterm = dropRowterm, dropColterm = dropColterm, 
                               nsegs = nsegs, nestorder = nestorder, 
                               degree = degree, difforder = difforder, 
                               rotateX = rotateX, ngridangles = ngridangles, 
                               which.rotacriterion = which.rotacriterion, 
                               nrotacores = nrotacores, 
                               asreml.opt = asreml.opt, 
                               tpps4mbf.obj = tpps4mbf.obj,
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, chooseOnIC = FALSE, 
                               maxit = maxit, IClikelihood = ic.lik, ...)
  
  return(spatial.asrt)  
}

addSpatialModelOnIC.asrtests <- function(asrtests.obj, spatial.model = "TPPS", 
                                         sections = NULL, 
                                         row.covar = "cRow", col.covar = "cCol", 
                                         row.factor = "Row", col.factor = "Col", 
                                         corr.funcs = c("ar1", "ar1"), 
                                         row.corrFitfirst = TRUE, allow.corrsJointFit = TRUE, 
                                         dropRowterm = NULL, dropColterm = NULL, 
                                         nsegs = NULL, nestorder = c(1, 1), 
                                         degree = c(3,3), difforder = c(2,2), 
                                         rotateX = FALSE, ngridangles = c(18, 18), 
                                         which.rotacriterion = "AIC", 
                                         nrotacores = 1, 
                                         asreml.option = "mbf", tpps4mbf.obj = NULL,  
                                         allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                                         checkboundaryonly = FALSE, update = FALSE, 
                                         maxit = 30, IClikelihood = "full", which.IC = "AIC", 
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
  spatial.asrts <- list()
  
  #Fit a local spatial model involving correlated effects
  if ("corr" %in% spatial.mod)
    spatial.asrt <- fitCorrMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               row.factor = row.factor, col.factor = col.factor, 
                               corr.funcs = corr.funcs, 
                               row.corrFitfirst = row.corrFitfirst, 
                               allow.corrsJointFit = allow.corrsJointFit, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, chooseOnIC = TRUE, 
                               maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
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
                                 maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                                 ...)

  #Fit a residual spatial model involving TPPS
  if ("TPPS" %in% spatial.mod)
    spatial.asrt <- fitTPPSMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               dropRowterm = dropRowterm, dropColterm = dropColterm, 
                               nsegs = nsegs, nestorder = nestorder, 
                               degree = degree, difforder = difforder, 
                               rotateX = rotateX, ngridangles = ngridangles, 
                               asreml.opt = asreml.opt, 
                               tpps4mbf.obj = tpps4mbf.obj,
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation, 
                               which.rotacriterion = which.rotacriterion, 
                               nrotacores = nrotacores, 
                               checkboundaryonly = checkboundaryonly, 
                               update = update, chooseOnIC = TRUE, 
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
                                            corr.funcs = c("ar1", "ar1"), 
                                            row.corrFitfirst = TRUE, allow.corrsJointFit = TRUE, 
                                            dropRowterm = NULL, dropColterm = NULL, 
                                            nsegs = NULL, nestorder = c(1, 1), 
                                            rotateX = FALSE, ngridangles = c(18, 18), 
                                            which.rotacriterion = "AIC", nrotacores = 1, 
                                            asreml.option = "mbf", tpps4mbf.obj = NULL, 
                                            allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                                            checkboundaryonly = FALSE, update = FALSE, 
                                            maxit = 30, IClikelihood = "full", which.IC = "AIC", 
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
                                            corr.funcs = corr.funcs, 
                                            row.corrFitfirst = row.corrFitfirst, 
                                            allow.corrsJointFit = allow.corrsJointFit, 
                                            allow.unconverged = allow.unconverged, 
                                            allow.fixedcorrelation = allow.fixedcorrelation,
                                            checkboundaryonly = checkboundaryonly, 
                                            update = update, chooseOnIC = TRUE, 
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
                                                dropRowterm = dropRowterm, dropColterm = dropColterm, 
                                                allow.unconverged = allow.unconverged, 
                                                allow.fixedcorrelation = allow.fixedcorrelation,
                                                checkboundaryonly = checkboundaryonly, 
                                                update = update, chooseOnIC = TRUE, 
                                                maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                                                ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["TPNCSS"]] , spatial.mod = "TPNCSS", 
                                   IClikelihood = ic.lik, spatial.IC = spatial.IC)
    }
    
    #Fit a residual spatial model involving TPPCS
    if ("TPPCS" %in% trySpatial)
    { 
      spatial.asrts[["TPPCS"]] <- fitTPPSMod(asrtests.obj, sections = sections, 
                                             row.covar = row.covar, col.covar = col.covar, 
                                             dropRowterm = dropRowterm, dropColterm = dropColterm, 
                                             nsegs = nsegs, nestorder = nestorder, 
                                             degree = c(3,3), difforder = c(2,2), 
                                             rotateX = rotateX, ngridangles = ngridangles, 
                                             which.rotacriterion = which.rotacriterion, 
                                             nrotacores = nrotacores, 
                                             asreml.opt = asreml.opt, 
                                             tpps4mbf.obj = tpps4mbf.obj,
                                             allow.unconverged = allow.unconverged, 
                                             allow.fixedcorrelation = allow.fixedcorrelation,
                                             checkboundaryonly = checkboundaryonly, 
                                             update = update, chooseOnIC = TRUE, 
                                             maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                                             ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["TPPCS"]] , spatial.mod = "TPPCS", 
                                   IClikelihood = ic.lik, spatial.IC = spatial.IC)
    }
    
    #Fit a residual spatial model involving TPP1LS
    if ("TPP1LS" %in% trySpatial)
    { 
      spatial.asrts[["TPP1LS"]] <- fitTPPSMod(asrtests.obj, sections = sections, 
                                              row.covar = row.covar, col.covar = col.covar, 
                                              dropRowterm = dropRowterm, dropColterm = dropColterm, 
                                              nsegs = nsegs, nestorder = nestorder, 
                                              degree = c(1,1), difforder = c(1,1), 
                                              rotateX = FALSE, ngridangles = c(0, 0), 
                                              which.rotacriterion = which.rotacriterion, 
                                              nrotacores = nrotacores, 
                                              asreml.opt = asreml.opt, 
                                              tpps4mbf.obj = tpps4mbf.obj,
                                              allow.unconverged = allow.unconverged, 
                                              allow.fixedcorrelation = allow.fixedcorrelation,
                                              checkboundaryonly = checkboundaryonly, 
                                              update = update, chooseOnIC = TRUE, 
                                              maxit = maxit, IClikelihood = ic.lik, which.IC = ic.type, 
                                              ...)
      spatial.IC <- calcSpatialICs(spatial.asrt = spatial.asrts[["TPP1LS"]] , spatial.mod = "TPP1LS", 
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

getVpars <- function(asreml.obj, asr4, asr4.2)
{
  if (asr4.2)
  { 
    vpc <- asreml.obj$vparameters.con
    vpt <- asreml.obj$vparameters.type
  } else
  {
    vpc <- vpc.char(asreml.obj)
    vpt <- vpt.char(asreml.obj)
  }
  return(list(vpc = vpc, vpt = vpt))
}

#Function to identify residual and correlation model terms that are currently fitted
getSectionVpars <- function( asreml.obj, sections, stub, corr.facs, which = c("res", "ran"), 
                             asr4, asr4.2)
{
  vpar <- getVpars(asreml.obj, asr4 = asr4, asr4.2 = asr4.2)
  vpc <- vpar$vpc
  vpt <- vpar$vpt

  #Get the resdiual term
  if ("res" %in% which)
  {
    #Get residual variance term using vpc and vpt
    vpc.res <- vpc[grepl("!R$", names(vpc))]
    vpt.res <- vpt[names(vpc.res)]
    vpt.res <- vpt.res[vpt.res == "V"]
    vpc.res <- vpc.res[names(vpt.res)]
    if (length(vpc.res) > 1) 
    { 
      if (!is.null(sections) && 
          #Check all residual variances include the sections name 
          all(sapply(names(vpc.res), 
                     function(vp) 
                     {
                       vp <- strsplit(vp, "\\_")[[1]][1]
                       vp <- vp == sections
                     })))
      { 
        kres.term <- grepl(sections, names(vpc.res)) & grepl(stub, names(vpc.res))
        vpc.res <- vpc.res[kres.term]
      } else
      { 
        warning("There are multiple residual terms, but at least some of these do not involve ",sections,
                "; consequeently the model cannot accommodate nugget variances.")
        vpc.res <- NULL
      }
    }
    
    if (length(vpc.res) != 1)
      warning("Could not find a residual term for ", sections, " ", stub) 
  } else
    vpc.res <- NULL
  
  #Get the random correlation terms (if sections, from all sections)
  if ("ran" %in%  which)
  { 
    vpc.ran <- vpc[!grepl("!R$", names(vpc))]
    if (length(vpc.ran) > 0)
      vpc.ran <- vpc.ran[sapply(names(vpc.ran),
                                function(term, corr.facs)
                                {
                                  term <- fac.getinTerm(rmTermDescription(term), asr4.2 = asr4.2)
                                  ran.corr <- 
                                    {
                                      if (length(term) == length(corr.facs))
                                        all(sapply(corr.facs,
                                                   {
                                                     function(fac, term)
                                                       any(grepl(fac, term))
                                                   }, term = term))
                                      else
                                        FALSE
                                    }
                                  return(ran.corr)
                                }, corr.facs = c(sections, corr.facs))]
    if (length(vpc.ran) == 0)
      vpc.ran <- NULL
  } else
    vpc.ran <- NULL
  return(list(ran = vpc.ran, res = vpc.res))
}

#Function to fix a single residual term or the residual variance for current level of section 
fixResTerm <- function(corr.asrt, sections, stub, asr4, asr4.2, 
                       fitfunc = "changeTerms", 
                       vpc.res, initial.values = 1,
                       maxit = 30, 
                       allow.unconverged = allow.unconverged, 
                       allow.fixedcorrelation = allow.fixedcorrelation,
                       checkboundaryonly = TRUE, 
                       update = update, 
                       IClikelihood = IClikelihood, 
                       which.IC = which.IC, ...)
{
  inargs <- list(...)
  
  bounds.excl <- c("S","B")
  
  #If already fixed, don't process
  if (vpc.res != "F")
  {
    #A single residual variance
    if (is.null(sections))
    { 
      resmod <- as.character(getFormulae(corr.asrt$asreml.obj)$residual)
      if (length(resmod) == 0)
        resmod <- NULL
      else
        resmod <- resmod[2]
      
      if (vpc.res %in% bounds.excl)
        fitfunc <- "changeTerms"
      lab <- "Try fixed residual variance"
      if (fitfunc == "changeTerms")
        lab <- "Force fixed residual variance"
      corr.asrt <- do.call(fitfunc, 
                           c(list(corr.asrt, 
                                  newResidual = resmod, label = lab, 
                                  set.terms = names(vpc.res), 
                                  initial.values = initial.values, 
                                  bounds = "F", ignore.suffices = FALSE, 
                                  maxit = maxit, 
                                  allow.unconverged = allow.unconverged, 
                                  allow.fixedcorrelation = allow.fixedcorrelation,
                                  checkboundaryonly = TRUE, 
                                  update = update, 
                                  IClikelihood = IClikelihood, 
                                  which.IC = which.IC), 
                             inargs))
      
      if (largeVparChange(corr.asrt$asreml.obj, 0.75))
        corr.asrt <- iterate(corr.asrt)
    } else #sections with multiple residuals - fix the one for this section
    {
      kres.term <- names(vpc.res)
      lab <- paste("Try fixed", names(kres.term))
      fitfunc <- "changeModelOnIC"
      if (vpc.res %in% bounds.excl)
      { 
        fitfunc <- "changeTerms"
        lab <- paste("Force fixed", kres.term)
      }
        corr.asrt <- do.call(fitfunc, 
                             c(list(corr.asrt, 
                                    label = lab, 
                                    set.terms = kres.term, 
                                    initial.values = initial.values, 
                                    bounds = "F", ignore.suffices = FALSE, 
                                    maxit = maxit, 
                                    allow.unconverged = allow.unconverged, 
                                    allow.fixedcorrelation = allow.fixedcorrelation,
                                    checkboundaryonly = TRUE, 
                                    update = update, 
                                    IClikelihood = IClikelihood, 
                                    which.IC = which.IC), 
                               inargs))
    }
  }
  return(corr.asrt)
}

rmRanTerm <- function(corr.asrt, vpbound, 
                      maxit = 30, 
                      allow.unconverged, allow.fixedcorrelation,
                      checkboundaryonly, update,
                      IClikelihood, which.IC, 
                      inargs)
{ 
  for (bound in vpbound)
  { 
    #Find term in random model formula to delete
    ran.terms <- getFormulae(corr.asrt$asreml.obj)$random
    ran.terms <- getTerms.formula(ran.terms)
    facs.bound <- fac.getinTerm(rmTermDescription(bound))
    facs.bound <- gsub('\\\"', "'", facs.bound)
    facs.bound <- gsub('\\(', "\\\\(", facs.bound)
    facs.bound <- gsub('\\)', "\\\\)", facs.bound)
    terms.bound <-  sapply(ran.terms,
                           function(term, facs.bound)
                           {
                             all(sapply(facs.bound,
                                        {
                                          function(fac, term)
                                            grepl(fac, term)
                                        }, term = term))
                           }, facs.bound = facs.bound)
    terms.bound <- names(terms.bound)[terms.bound]
    terms.bound <- paste(terms.bound, collapse = " + ")

    #Remove a bound term
    if (!is.null(terms.bound) && terms.bound != "")
      corr.asrt <- do.call(changeTerms,
                           c(list(corr.asrt,
                                  dropRandom = terms.bound,
                                  label = paste("Drop bound",terms.bound),
                                  maxit = maxit, 
                                  allow.unconverged = allow.unconverged,
                                  allow.fixedcorrelation = allow.fixedcorrelation,
                                  checkboundaryonly = FALSE,
                                  update = update,
                                  IClikelihood = IClikelihood,
                                  which.IC = which.IC),
                             inargs))
  }
  return(corr.asrt)
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

chk4SingularCorrTerms <- function(asrtests.obj, corr.asrt, label, 
                                  sections, stub, corr.facs, asr4, asr4.2)
{
  vpc.corr <- getSectionVpars(asrtests.obj$asreml.obj, 
                              sections = sections, stub = stub, 
                              corr.facs = corr.facs, 
                              asr4 = asr4, asr4.2 = asr4.2)
  
  
  #Determine the correlation terms, if any
  vpt.corr <- getVpars(asrtests.obj$asreml.obj, asr4, asr4.2)$vpt
  vpt.ran <- vpt.corr[names(vpc.corr$ran)]
  vpt.r <- vpt.ran[vpt.ran %in% c("R", "P", "C")]
  vpc.r <- vpc.corr$ran[names(vpt.r)]
  #Are there singular r terms
  if (length(vpc.r) > 0 && any(unlist(vpc.r) %in% "S"))
  {
    entry <- getTestEntry(asrtests.obj, label = label)
    entry$action <- "Unchanged - singular term(s)"
    corr.asrt$test.summary <- rbind(corr.asrt$test.summary, entry) 
  } else #no S terms
    corr.asrt <- asrtests.obj
  return(corr.asrt)
}

fitCorrMod <- function(asrtests.obj, sections = NULL,
                       row.covar = "cRow", col.covar = "cCol", 
                       row.factor = "Row", col.factor = "Col", 
                       corr.funcs = c("ar1", "ar1"), 
                       row.corrFitfirst = TRUE, allow.corrsJointFit = TRUE, 
                       allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                       checkboundaryonly = FALSE, update = TRUE, 
                       chooseOnIC = TRUE, 
                       maxit = 30, IClikelihood = "full", which.IC = "AIC", 
                       ...)
{
  inargs <- list(...)

  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  if (!asr4)
    stop(paste("Fitting spatial models using correlation/variance models", 
               "has not been implemented for asreml version less than 4.0"))
  
  bounds.excl <- c("B", "S")
  all.bounds.excl <- c(bounds.excl, "F")
  
  #Check that named columns are in the data
  dat.in <- asrtests.obj$asreml.obj$call$data
  if (is.symbol(dat.in))
    dat.in <- eval(dat.in)
  
  #Check the correlation functions and set them up
  id.funcs <- c("", "id", "idv")
  cor.funcs <- c("ar1", "ar2", "ar3", "sar","sar2",
                 "ma1", "ma2", "arma", "cor", "corb", "corg")
  cor.funcs <- c(cor.funcs, sapply(cor.funcs, function(f) paste0(f, c("v","h"))))
  met.funcs <- c("exp", "gau", "lvr")
  met.funcs <- c(met.funcs, sapply(met.funcs, function(f) paste0(f, c("v","h"))))
  unimpl.funcs <- c("iexp", "igau", "ieuc", "sph", "cir", "aexp", "agau", "mtrn")
  unimpl.funcs <- c(unimpl.funcs, sapply(unimpl.funcs, function(f) paste0(f, c("v","h"))))
  if (any(unimpl.funcs %in% corr.funcs))
    stop("Some of the following corr.funcs ar not implemented for spatial modelling: ", 
         paste(unimpl.funcs[unimpl.funcs %in% corr.funcs], collapse = ","))
  if (all(corr.funcs %in% id.funcs))
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
    if (!(func %in% id.funcs) && 
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
  nuggsOK <- TRUE
  corr.asrt <- asrtests.obj
  for (i in 1:nsect)
  {
    if (chooseOnIC)
      fitfunc <- "changeModelOnIC"
    else
      fitfunc <- "changeTerms"
    
    if (nsect > 1)
      stub <- levels(dat.in[[sections]])[i]
    else
      stub <- NULL
    corr.term <- FALSE
    #Check have a corr func
    if (any(rfuncs[1] == id.funcs))
      result1 <- "Unswapped"
    else
    { 
      #Try first correl in current section
      ran.term1 <- paste0(rterms[1], ":", facs[2])
      lab1 <- paste0("Try ", rterms[1])
      if (nsect > 1)
      {  
        ran.term1 <- paste0("at(", sections, ", '",stub, "'):", ran.term1)
        lab1 <- paste0(lab1, " for ", sections, " ",stub)
      }
      tmp.asrt <- do.call(fitfunc, 
                          c(list(corr.asrt, 
                                 addRandom = ran.term1, label = lab1, 
                                 allow.unconverged = allow.unconverged, 
                                 allow.fixedcorrelation = allow.fixedcorrelation,
                                 maxit = maxit, 
                                 checkboundaryonly = TRUE, 
                                 update = update, 
                                 IClikelihood = IClikelihood, 
                                 which.IC = which.IC), 
                            inargs))
      #Check for singular (S) correlation model terms and only change model if none
      corr.asrt <- chk4SingularCorrTerms(tmp.asrt, corr.asrt,  label = lab1, 
                                         sections = sections, stub = stub, 
                                         corr.facs = facs, 
                                         asr4 = asr4, asr4.2 = asr4.2)
      if (largeVparChange(corr.asrt$asreml.obj, 0.75))
        corr.asrt <- iterate(corr.asrt)
      result1 <- getTestEntry(corr.asrt, label = lab1)$action
    }
    
    #Try 2nd correl in current section
    if (!any(rfuncs[2] == id.funcs))
    {  
      lab <- paste0("Try ", rterms[2])
      if (nsect > 1)
        lab <- paste0(lab, " for ", sections, " ",stub)
      if (!grepl("Unswapped", result1) && !grepl("Unchanged", result1)) #first fac ar1 fitted
      { 
        corr.term <- TRUE
        last.term <- ran.term1
        #Check for ran.term1 in random formula and if absent check for different order
        last.term <- chk4TermInFormula(corr.asrt$asreml.obj$call$random, term = last.term, 
                                       asreml.obj = corr.asrt$asreml.obj)
        ran.term <- paste0(rterms[1], ":", rterms[2])
        if (nsect > 1)
          ran.term <- paste0("at(", sections, ", '",stub, "'):", ran.term)
        tmp.asrt <- do.call(fitfunc, 
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
                              inargs))
        #Check for singular (S) correlation model terms and only change model if none
        corr.asrt <- chk4SingularCorrTerms(tmp.asrt, corr.asrt,  label = lab, 
                                           sections = sections, stub = stub, 
                                           corr.facs = facs, 
                                           asr4 = asr4, asr4.2 = asr4.2)
        
        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
        if (!(grepl("Unswapped", getTestEntry(corr.asrt, label = lab)$action)) && 
            !(grepl("Unchanged", getTestEntry(corr.asrt, label = lab)$action)))
          last.term <- ran.term
      } else #no first fac corr
      { 
        ran.term <- paste0(facs[1], ":", rterms[2])
        if (nsect > 1)
          ran.term <- paste0("at(", sections, ", '",stub, "'):", ran.term)
        tmp.asrt <- do.call(fitfunc, 
                            c(list(corr.asrt, 
                                   addRandom = ran.term, label = lab, 
                                   maxit = maxit, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = checkboundaryonly, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC), 
                              inargs))
        #Check for singular (S) correlation model terms and only change model if none
        corr.asrt <- chk4SingularCorrTerms(tmp.asrt, corr.asrt, label = lab, 
                                           sections = sections, stub = stub, 
                                           corr.facs = facs, 
                                           asr4 = asr4, asr4.2 = asr4.2)
        
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
              ran.term1 <- paste0("at(", sections, ", '",stub, "'):", ran.term1)
            tmp.asrt <- do.call(fitfunc, 
                                c(list(corr.asrt, 
                                       dropRandom = last.term,
                                       addRandom = ran.term1, label = lab1, 
                                       maxit = maxit, 
                                       allow.unconverged = allow.unconverged, 
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = FALSE, 
                                       update = update, 
                                       IClikelihood = IClikelihood, 
                                       which.IC = which.IC), 
                                  inargs))
            #Check for singular (S) correlation model terms and only change model if none
            corr.asrt <- chk4SingularCorrTerms(tmp.asrt, corr.asrt, label = lab1, 
                                               sections = sections, stub = stub, 
                                               corr.facs = facs, 
                                               asr4 = asr4, asr4.2 = asr4.2)
            
            result1 <- getTestEntry(corr.asrt, label = lab1)$action
          }
        }
      }
      if (largeVparChange(corr.asrt$asreml.obj, 0.75))
        corr.asrt <- iterate(corr.asrt)

      #If no  correlation fitted and both rows and cols have corr funcs, try fitting them together
      if (!corr.term && (!any(rfuncs %in% id.funcs)) && allow.corrsJointFit)
      {
        #Try first correl in current section
        ran.term2 <- paste0(rterms[1], ":", rterms[2])
        lab2 <- paste0("Try ", rterms[1], " and ", rterms[2])
        if (nsect > 1)
        {  
          ran.term2 <- paste0("at(", sections, ", '",stub, "'):", ran.term2)
          lab2 <- paste0(lab2, " for ", sections, " ",stub)
        }
        tmp.asrt <- do.call(fitfunc, 
                            c(list(corr.asrt, 
                                   addRandom = ran.term2, label = lab2, 
                                   maxit = maxit, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = TRUE, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC), 
                              inargs))
        #Check for singular (S) correlation model terms and only change model if none
        corr.asrt <- chk4SingularCorrTerms(tmp.asrt, corr.asrt, label = lab2, 
                                           sections = sections, stub = stub, 
                                           corr.facs = facs, 
                                           asr4 = asr4, asr4.2 = asr4.2)
        
        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
        result2 <- getTestEntry(corr.asrt, label = lab2)$action
        if (!grepl("Unswapped", result2) && !grepl("Unchanged", result2)) #two-factor corr fitted
        { 
          corr.term <- TRUE
          last.term <- ran.term2
        }
      }
    } #end of 2nd correl in current section
    
    if (largeVparChange(corr.asrt$asreml.obj, 0.75))
      corr.asrt <- iterate(corr.asrt)

    #Test for nugget variance, only if the residual model is a variance model related to sections
    if (chooseOnIC && corr.term && nuggsOK)
    {
      #Get random and residual terms for correlation model for the current section
      vpc.corr <- getSectionVpars(corr.asrt$asreml.obj, 
                                  sections = sections, stub = stub, 
                                  corr.facs = facs, 
                                  asr4 = asr4, asr4.2 = asr4.2)

      #Determine the correlation terms, if any
      vpt.corr <- getVpars(corr.asrt$asreml.obj, asr4, asr4.2)$vpt
      vpt.ran <- vpt.corr[names(vpc.corr$ran)]
      vpt.r <- vpt.ran[vpt.ran %in% c("R", "P", "C")]
      vpc.r <- vpc.corr$ran[names(vpt.r)]
      #Are there no r terms or all r terms are bound or there is no residual term
      if (length(vpc.r) == 0 || all(vpc.r %in% all.bounds.excl) || is.null(vpc.corr$res))
        nuggsOK <- FALSE
      if (nuggsOK)
      {
        #Try fixing either the  single residual variance term or that for the current section 
        tmp.asrt <- do.call(fixResTerm, 
                            c(list(corr.asrt, sections = sections, stub = stub, 
                                   asr4 = asr4, asr4.2 = asr4.2, 
                                   fitfunc = "changeModelOnIC", vpc.res = vpc.corr$res, 
                                   maxit = maxit, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = TRUE, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC), 
                              inargs))
        lasttest <- tail(tmp.asrt$test.summary, 1)
        if (grepl("fixed residual variance", lasttest$terms) && 
            (grepl("Unswapped", lasttest$action) || grepl("Unchanged", lasttest$action)))
          corr.asrt <- tmp.asrt
        else
        { 
          new.vpc.corr <- getSectionVpars(tmp.asrt$asreml.obj, 
                                          sections = sections, stub = stub, 
                                          corr.facs = facs, 
                                          asr4 = asr4, asr4.2 = asr4.2)
          #Change if residual is fixed and either all random correlation model terms are unbound 
          #   or none have changed
          n.new <- length(new.vpc.corr$ran)
          n.old <- length(vpc.corr$ran)
          if (new.vpc.corr$res == "F" && 
              (n.new > 0 && (!any(new.vpc.corr$ran %in% bounds.excl) || 
                             (n.new == n.old && all(new.vpc.corr$ran == vpc.corr$ran)))))
            corr.asrt <- tmp.asrt
        }
        if (largeVparChange(corr.asrt$asreml.obj, 0.75))
          corr.asrt <- iterate(corr.asrt)
      }
    } #end of nugget variance test

    #Having made all model changes with checkboundaryonly = TRUE, update for checkboundaryonly set to FALSE
    if (!checkboundaryonly)
      corr.asrt <- rmboundary(corr.asrt, checkboundaryonly = checkboundaryonly, 
                              update = update, IClikelihood = IClikelihood)
    
    #Determine if there is a correlation term
    vpt.corr <- getVpars(corr.asrt$asreml.obj, asr4, asr4.2)$vpt
    vpt.corr <- vpt.corr[vpt.corr %in% c("R", "P", "C")]
    corr.term <- length(vpt.corr) > 0 

    #Further atttempts to deal with bound random and residual terms when 
    # (i) checkboundary only is FALSE and (ii) there are correlation terms
    if (corr.term && !checkboundaryonly)
    {
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
                                   corr.facs = facs, 
                                   asr4 = asr4, asr4.2 = asr4.2)$ran
        vpc.res <- getVpars(corr.asrt$asreml.obj, asr4, asr4.2)$vpc
        vpc.res <- vpc.res[grepl("!R$", names(vpc.res))]
        vpc.corr <- c(vpc.ran, vpc.res)
        
        if (any(unlist(vpc.corr) %in% all.bounds.excl))  #changed from bounds.excl
        {
          #Get vpc for this section
          vpc.corr <- getSectionVpars(corr.asrt$asreml.obj, 
                                      sections = sections, stub = stub, 
                                      corr.facs = facs, 
                                      asr4 = asr4, asr4.2 = asr4.2)
          #Only process if have bound residual and/or random corr model terms
          if (any(unlist(vpc.corr) %in% all.bounds.excl))
          {
            #get bound random terms
            vpc.bran <- vpc.corr$ran[vpc.corr$ran %in% all.bounds.excl]
            if (length(vpc.bran) > 0)
            { 
              vpt.bran <- getVpars(corr.asrt$asreml.obj, asr4, asr4.2)$vpt[names(vpc.bran)]
              #If any random correlations bound, remove corresponding term
              if (any(vpt.bran %in% c("R", "P", "C")))
              {
                #Get bound corr vpars and remove
                vpc.bC <- vpc.bran[vpt.bran %in% c("R", "P", "C")] 
                corr.asrt <- rmRanTerm(corr.asrt, vpbound = names(vpc.bC), 
                                       maxit = maxit, 
                                       allow.unconverged = allow.unconverged,
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = FALSE,
                                       update = update,
                                       IClikelihood = IClikelihood,
                                       which.IC = which.IC, 
                                       inargs = inargs)
                if (largeVparChange(corr.asrt$asreml.obj, 0.75))
                  corr.asrt <- iterate(corr.asrt)
                
                #update constraints for corr.bound under the new model
                vpc.corr <- getSectionVpars(corr.asrt$asreml.obj, 
                                            sections = sections, stub = stub, 
                                            corr.facs = facs, 
                                            asr4 = asr4, asr4.2 = asr4.2)
              }
            }
            
            #If Res is F or B and a ran variance is bound, remove bound ran term
            if (!is.null(vpc.corr$res) && vpc.corr$res %in% all.bounds.excl)
            {
              vpc.ran <- vpc.corr$ran
              vpt.ran <- getVpars(corr.asrt$asreml.obj, asr4, asr4.2)$vpt[names(vpc.ran)]
              vpc.Vran <- vpc.ran[vpt.ran %in% c("V", "G")]
              if (length(vpc.Vran) > 0 && any(vpc.Vran %in% bounds.excl))
              {
                #Get bound vars vpars
                vpc.bV <- vpc.Vran[vpc.Vran %in% bounds.excl] 
                if (length(vpc.bV) > 0)
                { 
                  #Remove the terms
                  corr.asrt <- rmRanTerm(corr.asrt, vpbound = names(vpc.bV),
                                         maxit = maxit, 
                                         allow.unconverged = allow.unconverged,
                                         allow.fixedcorrelation = allow.fixedcorrelation,
                                         checkboundaryonly = FALSE,
                                         update = update,
                                         IClikelihood = IClikelihood,
                                         which.IC = which.IC, 
                                         inargs = inargs)
                  
                  if (largeVparChange(corr.asrt$asreml.obj, 0.75))
                    corr.asrt <- iterate(corr.asrt)
                  
                  #update constraints for corr.bound under the new model
                  vpc.corr <- getSectionVpars(corr.asrt$asreml.obj, 
                                              sections = sections, stub = stub, 
                                              corr.facs = facs, 
                                              asr4 = asr4, asr4.2 = asr4.2)
                }              
              }
            }
            
            
            if (length(vpc.corr$ran) > 0)
            { 
              vpt.ran <- getVpars(corr.asrt$asreml.obj, asr4, asr4.2)$vpt[names(vpc.corr$ran)]
              vpc.Vran <- vpc.corr$ran[vpt.ran %in% c("V","G")]
              vpc.bVran <- vpc.Vran[vpc.Vran %in% all.bounds.excl]
            } else
              vpc.bVran <- vpc.Vran <- vpt.ran <- NULL
            if (xor(length(vpc.bVran) > 0, 
                    length(vpc.corr$res) > 0 && (vpc.corr$res %in% c("B","S"))))
            { 
              if (vpc.corr$res %in% c("B","S"))
              { 
                #loop until the section residual variance is not in c("B","S") 
                kloop <- 0
                while (vpc.corr$res %in% c("B","S") && kloop < 3)
                { 
                  tmp.asrt <- do.call(fixResTerm, 
                                      c(list(corr.asrt, sections = sections, stub = stub, 
                                             asr4 = asr4, asr4.2 = asr4.2, 
                                             fitfunc = "changeTerms", vpc.res = vpc.corr$res, 
                                             maxit = maxit, 
                                             allow.unconverged = allow.unconverged, 
                                             allow.fixedcorrelation = TRUE,
                                             checkboundaryonly = TRUE, 
                                             update = update, 
                                             IClikelihood = IClikelihood, 
                                             which.IC = which.IC), 
                                        inargs))
                  lab <- "Force fixed residual variance"
                  if (!is.null(sections))
                    lab <- paste("Force fixed", names(vpc.corr$res))
                  result <- getTestEntry(tmp.asrt, lab)$action
                  #update the constraints under the new model
                  vpc.corr <- getSectionVpars(tmp.asrt$asreml.obj, 
                                              sections = sections, stub = stub, 
                                              corr.facs = facs, 
                                              asr4 = asr4, asr4.2 = asr4.2)
                  #Check if a corr has gone bound and, it is has, remove it
                  vpc.bran <- vpc.corr$ran[vpc.corr$ran %in% all.bounds.excl]
                  if (length(vpc.bran) > 0)
                  { 
                    vpt.bran <- getVpars(tmp.asrt$asreml.obj, asr4, asr4.2)$vpt[names(vpc.bran)]
                    
                    #If any random correlations bound, remove corresponding term
                    if (any(vpt.bran %in% c("R", "P", "C")))
                    {
                      #Get bound corr vpars and remove
                      vpc.bC <- vpc.bran[vpt.bran %in% c("R", "P", "C")] 
                      corr.asrt <- rmRanTerm(corr.asrt, vpbound = names(vpc.bC),
                                             maxit = maxit, 
                                             allow.unconverged = allow.unconverged,
                                             allow.fixedcorrelation = allow.fixedcorrelation,
                                             checkboundaryonly = FALSE,
                                             update = update,
                                             IClikelihood = IClikelihood,
                                             which.IC = which.IC, 
                                             inargs = inargs)
                    } 
                  } else 
                  { 
                    if(grepl("Unchanged", result))
                    {
                      vpc.ran <- vpc.corr$ran[1]
                      corr.asrt <- rmRanTerm(corr.asrt, vpbound = names(vpc.ran),
                                             maxit = maxit, 
                                             allow.unconverged = allow.unconverged,
                                             allow.fixedcorrelation = allow.fixedcorrelation,
                                             checkboundaryonly = FALSE,
                                             update = update,
                                             IClikelihood = IClikelihood,
                                             which.IC = which.IC, 
                                             inargs = inargs)
                    } else
                      corr.asrt <- tmp.asrt
                  }
                  
                  #update the constraints under the new model
                  vpc.corr <- getSectionVpars(corr.asrt$asreml.obj, 
                                              sections = sections, stub = stub, 
                                              corr.facs = facs, 
                                              asr4 = asr4, asr4.2 = asr4.2)
                  kloop <- kloop + 1
                } #end while loop
              } else
              {
                for (bound in names(vpc.bVran))
                {
                  lab <- paste("Force fixed", bound)
                  tmp.asrt <- do.call(changeTerms, 
                                      c(list(corr.asrt, 
                                             label = lab, 
                                             set.terms = bound, initial.values = 1, 
                                             bounds = "F", ignore.suffices = FALSE, 
                                             maxit = maxit, 
                                             allow.unconverged = allow.unconverged, 
                                             allow.fixedcorrelation = allow.fixedcorrelation,
                                             checkboundaryonly = TRUE, 
                                             update = update, 
                                             IClikelihood = IClikelihood, 
                                             which.IC = which.IC), 
                                        inargs))
                  #Check for singular (S) correlation model terms and only change model if none
                  corr.asrt <- chk4SingularCorrTerms(tmp.asrt, corr.asrt, label = lab,  
                                                     sections = sections, stub = stub, 
                                                     corr.facs = facs, 
                                                     asr4 = asr4, asr4.2 = asr4.2)
                  result <- getTestEntry(corr.asrt, label = lab)$action
                  if (grepl("Unchanged", result))
                  {
                    corr.asrt <- rmRanTerm(corr.asrt, vpbound = bound,
                                           maxit = maxit, 
                                           allow.unconverged = allow.unconverged,
                                           allow.fixedcorrelation = allow.fixedcorrelation,
                                           checkboundaryonly = FALSE,
                                           update = update,
                                           IClikelihood = IClikelihood,
                                           which.IC = which.IC, 
                                           inargs = inargs)
                  }
                }   
              }
            } else #end xor(res, V ran bound)
            {
              #If both Res and Vran are in all.bounds.excl, remove the V ran terms
              vpc.Vran <- vpc.corr$ran[vpt.ran %in% c("V","G")]
              vpc.FBVran <- vpc.Vran[vpc.Vran %in% all.bounds.excl]
              if (length(vpc.FBVran) > 0 && vpc.corr$res %in% c("B","S"))
                corr.asrt <- rmRanTerm(corr.asrt, vpbound = names(vpc.FBVran),
                                       maxit = maxit, 
                                       allow.unconverged = allow.unconverged,
                                       allow.fixedcorrelation = allow.fixedcorrelation,
                                       checkboundaryonly = FALSE,
                                       update = update,
                                       IClikelihood = IClikelihood,
                                       which.IC = which.IC, 
                                       inargs = inargs)
            } #end one variance bound section
            #} #end one variance bound loop - don't iterate as can unfix a term
          } #end vpars section
        } #end dealing with bound vpars
      } #end bounds within a sections
    } #end checking bounds
  } #end of sections loop
  
  return(corr.asrt)
}

fitTPNCSSMod <- function(asrtests.obj, sections = NULL, 
                         row.covar = "cRow", col.covar = "cCol", 
                         dropRowterm = NULL, dropColterm = NULL, 
                         nsegs = NULL, 
                         allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                         checkboundaryonly = FALSE, update = TRUE, 
                         chooseOnIC = TRUE, 
                         maxit = 30, IClikelihood = "full", which.IC = "AIC", 
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
    if (chooseOnIC)
      tspl.asrt <- changeModelOnIC(tspl.asrt, 
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
                                   maxit = maxit, IClikelihood = IClikelihood, 
                                   which.IC = which.IC, 
                                   ...)
    else
      tspl.asrt <- changeTerms(tspl.asrt, 
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
                               maxit = maxit, IClikelihood = IClikelihood, ...)
    
  }
  
  #Final check for boundary terms
  if (!checkboundaryonly)
    tspl.asrt <- rmboundary(tspl.asrt, checkboundaryonly = checkboundaryonly, 
                            update = update, IClikelihood = IClikelihood)
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
                                            asreml = asreml.opt, 
                                            ...)
                           return(XZ.mat)
                         }, columncoordinates = col.coords, rowcoordinates = row.coords, 
                         sects = sections, nsegments = nsegs, 
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
                             rotateX = rotateX, theta = theta, 
                             asreml = asreml.opt, 
                             ...))
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
                                         asreml.option = "mbf", mbf.env = sys.frame(), 
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
                          dropRowterm, dropColterm, 
                          sections = NULL, 
                          row.covar, col.covar, lab, 
                          nsegs = NULL, nestorder = c(1, 1), 
                          degree = c(3,3), difforder = c(2,2), 
                          rotateX = FALSE, ngridangles = c(18,18), 
                          which.rotacriterion = "AIC", nrotacores = 1, 
                          asreml.opt = "mbf", stub = "xx", 
                          allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                          chooseOnIC = TRUE, 
                          checkboundaryonly = FALSE, update = TRUE, 
                          maxit = 30, IClikelihood = "full", which.IC = "AIC", ...)
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
  
  init.asrt <- tspl.asrt
  theta.opt <- NULL
  
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
    if (rotateX)
      labunrot <- gsub("tensor", "unrotated tensor", lab)
    else
      labunrot <- lab
    #tspl.asrt.unrot 
    tspl.asrt <- do.call(fitfunc, 
                         args = c(list(tspl.asrt,
                                       addFixed = fix.ch,
                                       dropFixed = drop.fix, 
                                       addRandom = ran.ch,
                                       dropRandom = drop.ran, 
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

    #Find the optimal theta for rotating the penalty eigenvectors, fit model with the rotation
    if (rotateX && any(difforder == 2))
    {
      ran.rot.ch <- paste(paste0(sect.fac,  
                                 c(paste0("TP.C.",1:difforder[1],":mbf(TP.row)"), 
                                   paste0("TP.R.",1:difforder[2],":mbf(TP.col)")),
                                 collapse = " + ")) 
      #Fit the reduced random model
      rot.asrt <- do.call(changeTerms, 
                          args = list(init.asrt, 
                                      addFixed = fix.ch,
                                      dropFixed = drop.fix, 
                                      addRandom = ran.rot.ch,
                                      dropRandom = drop.ran, 
                                      mbf = mbf.lis,
                                      label = "Fit model for rotation gridding", 
                                      allow.unconverged = TRUE, 
                                      allow.fixedcorrelation = TRUE,
                                      checkboundaryonly = TRUE, 
                                      update = update, 
                                      maxit = maxit, 
                                      IClikelihood = IClikelihood, 
                                      which.IC = which.IC))
      rot.asr <- rot.asrt$asreml.obj
      dev <- deviance.asr(rot.asr)
      
      #FInd the optimal thetas
      theta_opt <- rotate.penalty.U(rot.asr, data, sections = sections, ksect = ksect, 
                                    row.covar = row.covar, col.covar = col.covar,
                                    nsegs = nsegs, nestorder = nestorder,
                                    degree = degree, difforder = difforder,
                                    rotateX = rotateX, ngridangles = ngridangles, 
                                    which.rotacriterion = which.rotacriterion, 
                                    nrotacores = nrotacores, maxit = maxit, 
                                    asreml.opt = "mbf", mbf.env = sys.frame(), 
                                    stub = stub)
      theta.opt <- theta_opt$theta.opt
      cat("\n\n#### Optimal thetas:", paste(theta.opt, collapse = ", "), "\n\n")
      
      #Fit the P-splines for the optimal theta
      rm(list = Zmat.names, envir = mbf.env)
      mat <- makeTPPSplineMats(data, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar,
                               nsegs = nsegs, nestorder = nestorder,
                               degree = degree, difforder = difforder,
                               rotateX = rotateX, theta = theta.opt, 
                               asreml.opt = "mbf", mbf.env = NULL)[[ksect]]
      if (any(sapply(Zmat.names, exists, envir = mbf.env)))
        warning("THe following objects are being overwritten: ", 
                paste(Zmat.names[sapply(Zmat.names, exists, envir = parent.frame(2))], 
                      collapse = ", "))
      assign(Zmat.names[1], mat$BcZ.df, envir = mbf.env)
      assign(Zmat.names[2], mat$BrZ.df, envir = mbf.env)
      assign(Zmat.names[3], mat$BcrZ.df, envir = mbf.env)
      mbf.lis <- mat$mbflist
      dat <- mat$data.plus
      labrot <- gsub("tensor", 
                     paste0("rotated ", paste(round(theta.opt,0.5), collapse = ","), 
                            " tensor"), lab)
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
    
    ran.ch <- paste(paste0(sect.fac,  
                           c(paste0("grp(TP.C.",1:difforder[1],"_frow)"), 
                             paste0("grp(TP.R.",1:difforder[2],"_fcol)"), 
                             "grp(TP_fcol_frow)", 
                             paste0("dev(",row.covar,")"), 
                             paste0("dev(",col.covar,")"))), 
                    collapse = " + ")
    
    #Fit the full P-spline model, without rotation
    if (rotateX)
      labunrot <- gsub("tensor", "unrotated tensor", lab)
    else
      labunrot <- lab
    #tspl.asrt.unrot 
    tspl.asrt <- do.call(fitfunc, 
                         args = c(list(tspl.asrt, 
                                       addFixed = fix.ch,
                                       dropFixed = drop.fix, 
                                       addRandom = ran.ch,
                                       dropRandom = drop.ran, 
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

    #Find the optimal theta for rotating the penalty eigenvectors, fit model with the rotation
    if (rotateX && any(difforder == 2))
    {
      ran.rot.ch <- paste(paste0(sect.fac,  
                                 c(paste0("grp(TP.C.",1:difforder[1],"_frow)"), 
                                   paste0("grp(TP.R.",1:difforder[2],"_fcol)")), 
                                 collapse = " + "))
      #Fit the reduced random model
      rot.asrt <- do.call(changeTerms, 
                          args = list(init.asrt, 
                                      addFixed = fix.ch,
                                      dropFixed = drop.fix, 
                                      addRandom = ran.rot.ch,
                                      dropRandom = drop.ran, 
                                      group = grp,
                                      label = "Fit model for rotation gridding", 
                                      allow.unconverged = TRUE, 
                                      allow.fixedcorrelation = TRUE,
                                      checkboundaryonly = TRUE, 
                                      update = update, 
                                      maxit = maxit, 
                                      IClikelihood = IClikelihood, 
                                      which.IC = which.IC))
      rot.asr <- rot.asrt$asreml.obj
      dev <- deviance.asr(rot.asr)
      
      #Find the optimal thetas
      theta_opt <- rotate.penalty.U(rot.asr, data, sections = sections, ksect = ksect, 
                                    row.covar = row.covar, col.covar = col.covar,
                                    nsegs = nsegs, nestorder = nestorder,
                                    degree = degree, difforder = difforder,
                                    rotateX = rotateX, ngridangles = ngridangles, 
                                    which.rotacriterion = which.rotacriterion, 
                                    nrotacores = nrotacores, maxit = maxit, 
                                    stub = stub, mbf.env = sys.frame())
      theta.opt <- theta_opt$theta.opt
      cat("\n\n#### Optimal thetas:", paste(round(theta.opt,1), collapse = ","), "\n\n")

      #Fit the P-splines for the optimal theta
      mat <- makeTPPSplineMats(data, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar,
                               nsegs = nsegs, nestorder = nestorder,
                               degree = degree, difforder = difforder,
                               rotateX = rotateX, theta = theta.opt, 
                               asreml.opt = "grp")[[ksect]]
      grp <- mat$grp
      labrot <- gsub("tensor", 
                     paste0("rotated ", paste(theta.opt, collapse = ", "), 
                            " tensor"), lab)
      if (chooseOnIC)
      { 
        tspl.asrt <- updateOnIC.asrtests(tspl.asrt, data = mat$data.plus, 
                                         grp = grp, maxit = maxit, 
                                         label = labrot, IClikelihood = IClikelihood, 
                                         which.IC = which.IC)
        if (grepl("Unswapped", getTestEntry(tspl.asrt, label = labrot)$action))
          theta.opt <- c(0,0)
      } else
        tspl.asrt <- update.asrtests(tspl.asrt, data = mat$data.plus, 
                                     grp = grp, maxit = maxit, 
                                     label = labrot, IClikelihood = IClikelihood)
      #Check criteria
      # print(infoCriteria(list(old = tspl.asrt$asreml.obj, new = new.asr), IClikelihood = "full"))
      # new.dev <- deviance.asr(new.asr) 
      # print(c(dev, new.dev))
    }
  }
  
  #Prepare for fitting unstructured model to row and col marginal termsw
  vpars <- names(tspl.asrt$asreml.obj$vparameters)
  #repln <- as.data.frame(table(data[c(sections,row.covar,col.covar)]))
  repln <- 1
  
  if (!is.null(sections))
  {   
    klev <- levels(data[[sections]])[ksect]
    vpars <- vpars[grepl(klev, vpars)]
    vpars <- gsub(paste0("'",klev,"'"), ksect, vpars)
    repln <- length(levels(data[[sections]]))
  }
  
  #If more than one col variable in the marginal random row term in this section, try unstructured model
  rowmarg.vpar <- vpars[grepl("TP\\.C\\.", vpars)]
  nr <- length(rowmarg.vpar)
  if (nr > 1)
  {
    drop.ran <-paste(rowmarg.vpar, collapse = " + ") 
    add.ran <- paste0("str( ~ ", drop.ran, ", 
                      ~ ", ifelse(nr > 2, "corgh", "corh"),"(", length(rowmarg.vpar), 
                      "):id(", mat$dim['nbr']*repln,"))")
    if (asreml.opt == "mbf")
      tspl.asrt <- do.call(fitfunc, 
                           args = c(list(tspl.asrt, 
                                         addRandom = add.ran,
                                         dropRandom = drop.ran, 
                                         mbf = mbf.lis,
                                         label = "Try us variance for random row terms", 
                                         allow.unconverged = allow.unconverged, 
                                         allow.fixedcorrelation = allow.fixedcorrelation,
                                         checkboundaryonly = TRUE, 
                                         update = update, 
                                         maxit = maxit, 
                                         IClikelihood = IClikelihood, 
                                         which.IC = which.IC), 
                                    inargs))
    else
      tspl.asrt <- do.call(fitfunc, 
                           args = c(list(tspl.asrt, 
                                         addRandom = add.ran,
                                         dropRandom = drop.ran, 
                                         group = grp,
                                         label = "Try us variance for random row terms", 
                                         allow.unconverged = allow.unconverged, 
                                         allow.fixedcorrelation = allow.fixedcorrelation,
                                         checkboundaryonly = TRUE, 
                                         update = update, 
                                         maxit = maxit, 
                                         IClikelihood = IClikelihood, 
                                         which.IC = which.IC), 
                                    inargs))
  }
  
  #If more than one row variable in the marginal random col term in this section, try unstructured model
  colmarg.vpar <- vpars[grepl("TP\\.R\\.", vpars)]
  nc <- length(colmarg.vpar)
  if (length(colmarg.vpar) > 1)
  {
    drop.ran <-paste(colmarg.vpar, collapse = " + ") 
    add.ran <- paste0("str( ~ ", drop.ran, ", ~ ", ifelse(nr > 2, "corgh", "corh"), "(", nc, 
                      "):id(", mat$dim['nbc']*repln,"))")
    if (asreml.opt == "mbf")
      tspl.asrt <- do.call(fitfunc, 
                           args = c(list(tspl.asrt, 
                                         addRandom = add.ran,
                                         dropRandom = drop.ran, 
                                         mbf = mbf.lis,
                                         label = "Try us variance on random col terms", 
                                         allow.unconverged = allow.unconverged, 
                                         allow.fixedcorrelation = allow.fixedcorrelation,
                                         checkboundaryonly = TRUE, 
                                         update = FALSE, #to ensure clean refit
                                         maxit = maxit, 
                                         IClikelihood = IClikelihood, 
                                         which.IC = which.IC), 
                                    inargs))
    else
      tspl.asrt <- do.call(fitfunc, 
                           args = c(list(tspl.asrt, 
                                         addRandom = add.ran,
                                         dropRandom = drop.ran, 
                                         group = grp,
                                         label = "Try us variance on random col terms", 
                                         allow.unconverged = allow.unconverged, 
                                         allow.fixedcorrelation = allow.fixedcorrelation,
                                         checkboundaryonly = TRUE, 
                                         update = update, 
                                         maxit = maxit, 
                                         IClikelihood = IClikelihood, 
                                         which.IC = which.IC), 
                                    inargs))
  }
  tspl.asrt <- rmboundary(tspl.asrt, checkboundaryonly = checkboundaryonly, 
                          update = update, IClikelihood = IClikelihood)
  attr(tspl.asrt$asreml.obj, which = "theta.opt") <- theta.opt
  
  return(tspl.asrt)
}



#Fit a tensor-spline spatial model
fitTPPSMod <- function(asrtests.obj, sections = NULL, 
                       row.covar = "cRow", col.covar = "cCol", 
                       dropRowterm = NULL, dropColterm = NULL, 
                       nsegs = NULL, nestorder = c(1, 1), 
                       degree = c(3,3), difforder = c(2,2), 
                       rotateX = FALSE, ngridangles = c(18,18),
                       which.rotacriterion = "AIC", nrotacores = 1, 
                       asreml.opt = "mbf", 
                       tpps4mbf.obj = NULL, 
                       allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                       checkboundaryonly = FALSE, update = TRUE, 
                       maxit = 30, 
                       chooseOnIC = TRUE, 
                       IClikelihood = "full", which.IC = "AIC",
                       ...)
{ 

  inargs <- list(...)
  checkEllipsisArgs("makeTPPSplineMats.data.frame", inargs)
  checkEllipsisArgs("tpsmmb", inargs, pkg = "asremlPlus")
  
  #Check which.criterion options
  options <- c("deviance", "likelihood", "AIC", "BIC")
  which.rotacriterion <- options[check.arg.values(which.rotacriterion, options)]
  
  #Stop parallel processing for mbf
  if (asreml.opt == "mbf" && nrotacores > 1)
    stop(paste("Parallel processing has not been implemented for asreml.option set to mbf;",
               "nrotacores must be one"))

  #Check nsegs, nestorder, degree, difforder and ngridangles
  if (length(nsegs) != 2 && !is.null(nsegs))
    stop("nsegs must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(nestorder) != 2)
    stop("nestorder must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(degree) != 2)
    stop("degree must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(difforder) != 2)
    stop("difforder must specify exactly 2 values, one for each of the column and row dimensions")
  if (length(ngridangles) != 2)
    stop("ngridangles must specify exactly 2 values, one for each of the column and row dimensions")

  #Get the data from the original call and check that named columns are in the data
  dat.in <- asrtests.obj$asreml.obj$call$data
  if (is.symbol(dat.in))
    dat.in <- eval(dat.in)
  checkNamesInData(c(sections, row.covar, col.covar, dropRowterm, dropColterm), dat.in)
  #Check conformability of covars and factors
  if (!is.null(dropRowterm) && nlevels(dat.in[[dropRowterm]]) != length(unique(dat.in[[row.covar]])))
    stop(dropRowterm, " does not have the same number of levels as there are values of ", row.covar)
  if (!is.null(dropColterm) && nlevels(dat.in[[dropColterm]]) != length(unique(dat.in[[col.covar]])))
    stop(dropColterm, " does not have the same number of levels as there are values of ", col.covar)
  
  if (is.null(tpps4mbf.obj))
  {   
    #Create spline basis functions
    #do not set to NULL so that the mbf df will be assigned in the current environment
    # - needed for the rmboundary call at the end
    tps.XZmat <- makeTPPSplineMats(dat.in, sections = sections, 
                                   row.covar = row.covar, col.covar = col.covar,
                                   nsegs = nsegs, nestorder = nestorder,
                                   degree = degree, difforder = difforder, 
                                   rotateX = rotateX,
                                   asreml.opt = asreml.opt, 
                                   ...)
  }
  else #user supplied
    tps.XZmat <- tpps4mbf.obj
  dat <- tps.XZmat[[1]]$data.plus

  #Update the asreml.obj for the new data.frame
  asreml.obj  <- asrtests.obj$asreml.obj
  asreml.obj <- newfit(asreml.obj, data = dat)
  tspl.asrt <- as.asrtests(asreml.obj = asreml.obj, NULL, NULL, 
                           IClikelihood = "full", label = "Change to new data.frame with TPS bits")

  #Fit spatial TPPS to sections
  nsect <- calc.nsect(dat, sections)
  rotated <- rotateX & any(difforder > 1)
  if (rotated) theta.opt <- list()
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
    tspl.asrt <- fitTPSModSect(tspl.asrt, data = dat.in, mat = tps.XZmat[[i]], 
                               ksect = i, sect.fac = sect.fac, 
                               dropRowterm = dropRowterm, dropColterm = dropColterm, 
                               sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               nsegs = nsegs, nestorder = nestorder, 
                               degree = degree, difforder = difforder,
                               rotateX = rotateX, ngridangles = ngridangles, 
                               which.rotacriterion = which.rotacriterion, 
                               nrotacores = nrotacores, 
                               lab = lab, asreml.opt = asreml.opt, stub = stub, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, 
                               maxit = maxit, 
                               chooseOnIC = chooseOnIC, 
                               IClikelihood = IClikelihood, 
                               which.IC = which.IC, 
                               ...)
    if (rotated)
    { 
      theta.opt <- c(theta.opt, list(attr(tspl.asrt$asreml.obj, which = "theta.opt")))
    }
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
  return(tspl.asrt)
}

