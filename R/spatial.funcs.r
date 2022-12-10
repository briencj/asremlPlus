#The following functions have been developed for single and multiple section experiments

checkTrySpatial <- function(trySpatial)
{
  trySpat.opts <- c("none", "corr", "TPNCSS", "TPPS", "all")
  trySpatial <- trySpat.opts[unlist(lapply(trySpatial, check.arg.values, options=trySpat.opts))]
  
  if (length(intersect(trySpatial, trySpat.opts)) == 0)
    stop("trySpatial must be one of ", paste0(trySpat.opts, collapse = ", "))
  if ("all" %in% trySpatial)
    trySpatial <- c("corr", "TPNCSS", "TPPS")
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

addSpatialModelOnIC.asrtests <- function(asrtests.obj, spatial.model = "TPPS", 
                                         sections = NULL, 
                                         row.covar = "cRow", col.covar = "cCol", 
                                         row.factor = NULL, col.factor = NULL, 
                                         nsegs = NULL, asreml.option = "mbf", tpps4mbf.obj = NULL,  
                                         allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                                         checkboundaryonly = FALSE, update = FALSE, 
                                         IClikelihood = "full", which.IC = "AIC", 
                                         ...)
{    
  #Deal with arguments for tpsmmb and changeModelOnIC
  inargs <- list(...)
  checkEllipsisArgs(c("tpsmmb","changeModelOnIC", "asreml"), inargs)

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
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, 
                               IClikelihood = ic.lik, which.IC = ic.type, 
                               ...)
  #Fit a local spatial model involving TPNCSS
  if ("TPNCSS" %in% spatial.mod)
    spatial.asrt <- fitTPNCSSMod(asrtests.obj, sections = sections, 
                                 row.covar = row.covar, col.covar = col.covar, 
                                 row.factor = row.factor, col.factor = col.factor, 
                                 allow.unconverged = allow.unconverged, 
                                 allow.fixedcorrelation = allow.fixedcorrelation,
                                 checkboundaryonly = checkboundaryonly, 
                                 update = update, 
                                 IClikelihood = ic.lik, which.IC = ic.type, 
                                 ...)
  
  #Fit a residual spatial model involving TPPS
  if ("TPPS" %in% spatial.mod)
    spatial.asrt <- fitTPPSMod(asrtests.obj, sections = sections, 
                               row.covar = row.covar, col.covar = col.covar, 
                               row.factor = row.factor, col.factor = col.factor, 
                               nsegs = nsegs, asreml.opt = asreml.opt, 
                               tpps4mbf.obj = tpps4mbf.obj,
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, 
                               IClikelihood = ic.lik, which.IC = ic.type, 
                               ...)
  
  return(spatial.asrt)  
}

chooseSpatialModelOnIC.asrtests <- function(asrtests.obj, trySpatial = "all", 
                                            sections = NULL, 
                                            row.covar = "cRow", col.covar = "cCol", 
                                            row.factor = NULL, col.factor = NULL, 
                                            nsegs = NULL, asreml.option = "mbf", tpps4mbf.obj = NULL, 
                                            allow.unconverged = FALSE, allow.fixedcorrelation = FALSE,
                                            checkboundaryonly = FALSE, update = FALSE, 
                                            IClikelihood = "full", which.IC = "AIC", 
                                            return.asrts = "best", ...)
{    
  #Deal with arguments for tpsmmb and changeModelOnIC
  inargs <- list(...)
  checkEllipsisArgs(c("tpsmmb","changeModelOnIC", "asreml"), inargs)
  
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
    spatial.IC <- infoCriteria(asrtests.obj)
  } else #fit a spatial model
  {
    asreml::asreml.options(extra = 5, ai.sing = TRUE, fail = "soft")
    spatial.asrts <- list()
    #Fit a local spatial model involving correlated effects
    if ("corr" %in% trySpatial)
      spatial.asrts[["corr"]] <- fitCorrMod(asrtests.obj, sections = sections, 
                                            row.covar = row.covar, col.covar = col.covar, 
                                            allow.unconverged = allow.unconverged, 
                                            allow.fixedcorrelation = allow.fixedcorrelation,
                                            checkboundaryonly = checkboundaryonly, 
                                            update = update, 
                                            IClikelihood = ic.lik, which.IC = ic.type, 
                                            ...)
    #Fit a local spatial model involving TPNCSS
    if ("TPNCSS" %in% trySpatial)
      spatial.asrts[["TPNCSS"]] <- fitTPNCSSMod(asrtests.obj, sections = sections, 
                                                row.covar = row.covar, col.covar = col.covar, 
                                                row.factor = row.factor, col.factor = col.factor, 
                                                allow.unconverged = allow.unconverged, 
                                                allow.fixedcorrelation = allow.fixedcorrelation,
                                                checkboundaryonly = checkboundaryonly, 
                                                update = update, 
                                                IClikelihood = ic.lik, which.IC = ic.type, 
                                                ...)
    
    #Fit a residual spatial model involving TPPS
    if ("TPPS" %in% trySpatial)
      spatial.asrts[["TPPS"]] <- fitTPPSMod(asrtests.obj, sections = sections, 
                                            row.covar = row.covar, col.covar = col.covar, 
                                            row.factor = row.factor, col.factor = col.factor, 
                                            nsegs = nsegs, asreml.opt = asreml.opt, 
                                            tpps4mbf.obj = tpps4mbf.obj,
                                            allow.unconverged = allow.unconverged, 
                                            allow.fixedcorrelation = allow.fixedcorrelation,
                                            checkboundaryonly = checkboundaryonly, 
                                            update = update, 
                                            IClikelihood = ic.lik, which.IC = ic.type, 
                                            ...)
    
    asrt.names <- c("nonspatial", names(spatial.asrts))
    tmp.asrts <- c(list(nonspatial = asrtests.obj), spatial.asrts)
    spatial.IC <- infoCriteria(lapply(tmp.asrts, function(asrt) asrt$asreml.obj))
    # spatial.IC <- lapply(spatial.asrts, function(asrt) infoCriteria(asrt$asreml.obj))
    # spatial.IC <- do.call(rbind, spatial.IC)
    
    spatial.comp <- round(spatial.IC[[which.IC]], digits = 3)
    names(spatial.comp) <- rownames(spatial.IC)
    min.asrt <- which.min(spatial.comp)
    if (length(min.asrt) > 1)
    {
      #pick one in the order given below
      if ("nonspatial" %in% names(min.asrt)) min.asrt <- min.asrt["nonspatial"]
      if ("TPNCSS" %in% names(min.asrt)) min.asrt <- min.asrt["TPNCSS"]
      if ("corr" %in% names(min.asrt)) min.asrt <- min.asrt["corr"]
      if ("TPPS" %in% names(min.asrt)) min.asrt <- min.asrt["TPPS"]
    }
    #If only best, get the best astests.obj
    if (return.opt == "best")
      spatial.asrts <- tmp.asrts[names(min.asrt)]
    
  } 
  return(list(asrts = spatial.asrts, spatial.IC = spatial.IC, 
              best = names(min.asrt), best.AIC = spatial.comp[min.asrt]))
}
  
fitCorrMod <- function(asrtests.obj, sections = NULL,
                       row.covar = "cRow", col.covar = "cCol", 
                       allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                       checkboundaryonly = FALSE, update = TRUE, 
                       IClikelihood = "full", which.IC = "AIC", 
                       ...)
{
  #Check that named columns are in the data
  dat.in <- asrtests.obj$asreml.obj$call$data
  if (is.symbol(dat.in))
    dat.in <- eval(dat.in)
  checkNamesInData(c(sections, row.covar, col.covar), dat.in)
  
  nsect <- calc.nsect(dat.in, sections)

  #Loop over the sections
  corr.asrt <- asrtests.obj
  for (i in 1:nsect)
  {
    corr.term <- FALSE
    #Test for column exp in current section
    ran.term <- paste0(row.covar, ":exp(", col.covar, ")")
    lab <- "Try column exp"
    if (nsect > 1)
    {  
      ran.term <- paste0("at(", sections, ", ",i, "):", ran.term)
      lab <- paste0(lab, " for ", sections, " ",i)
    }
    corr.asrt <- changeModelOnIC(corr.asrt, 
                                 addRandom = ran.term, label = lab, 
                                 allow.unconverged = allow.unconverged, 
                                 allow.fixedcorrelation = allow.fixedcorrelation,
                                 checkboundaryonly = checkboundaryonly, 
                                 update = update, 
                                 IClikelihood = IClikelihood, 
                                 which.IC = which.IC, 
                                 ...)
    result <- getTestEntry(corr.asrt, label = lab)$action
    
    #Try row exp in current section
    lab <- "Try row exp"
    if (nsect > 1)
      lab <- paste0(lab, " for ", sections, " ",i)
    if (!grepl("Unswapped", result))
    { 
      corr.term <- TRUE
      last.term <- ran.term
      if (nsect > 1)
        last.term <- paste0(row.covar, ":at(", sections, ", ",i, "):exp(", col.covar, ")") 
      ran.term <- paste0("exp(", row.covar, "):exp(", col.covar, ")")
      if (nsect > 1)
        ran.term <- paste0("at(", sections, ", ",i, "):", ran.term)
      
      corr.asrt <- changeModelOnIC(corr.asrt, 
                                   addRandom = ran.term, 
                                   dropRandom = last.term, label = lab, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = checkboundaryonly, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC, 
                                   ...)
      if (!(grepl("Unswapped", getTestEntry(corr.asrt, label = lab)$action, 
                  fixed = TRUE)))
        last.term <- ran.term
    } else
    { 
      ran.term <- paste0("exp(", row.covar, "):", col.covar)
      if (nsect > 1)
        ran.term <- paste0("at(", sections, ", ",i, "):", ran.term)
      corr.asrt <- changeModelOnIC(corr.asrt, 
                                   addRandom = ran.term, label = lab, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = checkboundaryonly, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC, 
                                   ...)
      if (!(grepl("Unswapped", getTestEntry(corr.asrt, label = lab)$action, 
                  fixed = TRUE)))
      { 
        corr.term <- TRUE
        last.term <- ran.term
      }
    }
    #Try for units term
    if (corr.term)
    {
      lab <- "Try units term"
      if (nsect > 1)
        lab <- paste0(lab, " for ", sections, " ",i)
      ran.term <- "units"
      if (nsect > 1)
        ran.term <- paste0("at(", sections, ", ",i, "):", ran.term)
      corr.asrt <- changeModelOnIC(corr.asrt,
                                   addRandom = ran.term, label = lab,
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = checkboundaryonly, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC, 
                                   ...)
    }
  }
  
  return(corr.asrt)
}

fitTPNCSSMod <- function(asrtests.obj, sections = NULL, 
                         row.covar = "cRow", col.covar = "cCol", 
                         row.factor = NULL, col.factor = NULL, 
                         nsegs = NULL, 
                         allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                         checkboundaryonly = FALSE, update = TRUE, 
                         IClikelihood = "full", which.IC = "AIC", 
                         ...)
{ 
  #Check that named columns are in the data
  dat.in <- asrtests.obj$asreml.obj$call$data
  if (is.symbol(dat.in))
    dat.in <- eval(dat.in)
  checkNamesInData(c(sections, row.covar, col.covar, row.factor, col.factor), dat.in)
  
  #Check conformability of covars and factors
  if (!is.null(row.factor) && nlevels(dat.in[[row.factor]]) != length(unique(dat.in[[row.covar]])))
    stop(row.factor, " does not have the same number of levels as there are values of ", row.covar)
  if (!is.null(col.factor) && nlevels(dat.in[[col.factor]]) != length(unique(dat.in[[col.covar]])))
    stop(col.factor, " does not have the same number of levels as there are values of ", col.covar)
  
  #Are row.factor and col.factor already in the model?
  facs <- c(row.factor, col.factor)
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
  
  #Construct terms with sections
  spl.row <- paste0("spl(", row.covar, ")")
  spl.col <- paste0("spl(", col.covar, ")")
  spl.terms <- c(spl.row, spl.col, 
                 paste0("dev(", row.covar, ")"), paste0("dev(", col.covar, ")"), 
                 paste0(spl.row, ":", col.covar), paste0(row.covar, ":", spl.col), 
                 paste0(spl.row, ":", spl.col))
  
  #Add the fixed terms
  tspl.asrt <- changeTerms(asrtests.obj, 
                           addFixed = fix.terms, 
                           dropFixed = drop.fix, 
                           dropRandom = drop.ran, 
                           IClikelihood = "full", 
                           label = paste0("Add linear ",  row.covar, 
                                          " & ", col.covar, "terms")) 
  
  for (i in 1:nsect)
  {
    if (nsect > 1) 
      sect.fac <- paste0("at(", sections, ", ", i, "):")
    lab <- paste0("Try tensor NCS splines for ", sections, " ",i)
    #Fit TPNCSS to a section 
    tspl.asrt <- changeModelOnIC(tspl.asrt, 
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
    tspl.asrt <- rmboundary(tspl.asrt)
  }
  
  return(tspl.asrt)
}


#'### Function to create spline basis matrices and data for TPS splines
makeTPSPlineXZMats.data.frame <- function(data, sections = NULL, 
                                          row.covar, col.covar, 
                                          nsegs = NULL, asreml.option = "mbf", 
                                         ...)
{
  stub = "xx"

  #Check that named columns are in the data
  checkNamesInData(c(sections, row.covar, col.covar), data)
  
  #Make sure that row covar is not centred so that tpsmmb does not create a faulty index
  row.min <- min(data[[row.covar]], na.rm = TRUE)
  if (row.min < 0)
    data[row.covar] <- data[[row.covar]] - row.min + 1
  
  col.min <- min(data[[col.covar]], na.rm = TRUE)
  if (col.min < 0)
    data[col.covar] <- data[[col.covar]] - col.min + 1
  
  #Deal with arguments for tpsmmb and changeModelOnIC
  inargs <- list(...)
  checkEllipsisArgs(c("tpsmmb"), inargs)
  
  #Check nsegs
  if (length(nsegs) > 2)
    stop("nsegs must specify no more than 2 values")
  
  #Check asreml.option
  options <- "mbf"
  asreml.opt <- options[check.arg.values(asreml.option, options)]
  
  nsect <- calc.nsect(data, sections)
  
  rc.cols <- c(sections, row.covar, col.covar)
  dat.rc <- data[rc.cols]
  dat.rc <- unique(dat.rc)
  if (nsect != 1)
  {  
    tmp <- split(data, f = data[[sections]])
    sect.levs <- levels(data[[sections]])
    tps.XZmat  <- lapply(tmp, 
                         function(dat, columncoordinates, rowcoordinates, 
                                  sects, nsegments, asreml.opt)
                         { 
                           stub <- dat[[sects]][1]
                           XZ.mat <- tpsmmb(columncoordinates = columncoordinates, 
                                            rowcoordinates = rowcoordinates, 
                                            data = dat.rc, 
                                            stub = stub, nsegments = nsegments, 
                                            asreml = asreml.opt, 
                                            ...)
                           Zmat.names <- paste0(paste0(c("BcZ", "BrZ", "BcrZ"),stub), ".df")
                           assign(Zmat.names[1], XZ.mat$BcZ.df, envir = parent.frame(1))
                           assign(Zmat.names[2], XZ.mat$BrZ.df, envir = parent.frame(1))
                           assign(Zmat.names[3], XZ.mat$BcrZ.df, envir = parent.frame(1))
                           return(XZ.mat)
                         }, columncoordinates = col.covar, rowcoordinates = row.covar, 
                         sects = sect.levs, nsegments = nsegs, 
                         asreml.opt = asreml.opt)
  } else
  {
    tps.XZmat <- list(tpsmmb(columncoordinates = col.covar, 
                             rowcoordinates = row.covar, 
                             data = dat.rc, 
                             stub = stub, nsegments = nsegs, 
                             asreml = asreml.opt, 
                             ...))
    Zmat.names <- paste0(paste0(c("BcZ", "BrZ", "BcrZ"),stub), ".df")
    if (any(sapply(Zmat.names, exists, envir = parent.frame(1))))
      warning("THe following objects are being overwritten: ", 
              paste(Zmat.names[sapply(Zmat.names, exists, envir = parent.frame(2))], 
                    collapse = ", "))
    assign(Zmat.names[1], tps.XZmat[[1]]$BcZ.df, envir = parent.frame(1))
    assign(Zmat.names[2], tps.XZmat[[1]]$BrZ.df, envir = parent.frame(1))
    assign(Zmat.names[3], tps.XZmat[[1]]$BcrZ.df, envir = parent.frame(1))
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
  
  #Remove any previous Tensor Spline basis columns from the data
  if (any(grepl("TP\\.", names(data))) || any(grepl("TP\\_", names(data))))
    data <- data[,-c(grep("TP\\.", names(data)), grep("TP\\_", names(data)))]
  dat <- merge(data, dat, all.x = TRUE, by = rc.cols, sort = FALSE)
  
  #Add to returned object
  tps.XZmat <- lapply(tps.XZmat, function(mat, data = dat) c(mat, list(data.plus = dat)))
  
  return(tps.XZmat)
}

addPSdesign.mat <- function(dat, sections = NULL, nsect = 1, 
                            row.coords, col.coords, 
                            nsegs = NULL, asreml.opt = "grp", stub = "xx", 
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
                             asreml = asreml.opt, 
                             ...))
  }
  
  return(tps.XZmat)
}

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

fitTPSModSect <- function(tspl.asrt, mat, sect.fac, row.factor, col.factor, 
                          row.covar, col.covar, lab, 
                          asreml.opt = "mbf", stub = stub, 
                          allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                          checkboundaryonly = FALSE, update = TRUE, 
                          IClikelihood = "full", which.IC = "AIC", ...)
{
  inargs <- list(...)

  #Are row.factor and col.factor already in the model?
  facs <- c(row.factor, col.factor)
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

  if (asreml.opt == "mbf")
  {
    mbf.lis <- mat$mbflist
    tspl.asrt <- do.call(changeModelOnIC, 
                         args = c(list(tspl.asrt,
                                       addFixed = paste(paste0(
                                         sect.fac,
                                         c("TP.CR.2", "TP.CR.3", "TP.CR.4")),
                                         collapse = " + "),
                                       dropFixed = drop.fix, 
                                       addRandom = paste(paste0(
                                         sect.fac,
                                         c("TP.C.1:mbf(TP.row)", "TP.C.2:mbf(TP.row)",
                                           "TP.R.1:mbf(TP.col)", "TP.R.2:mbf(TP.col)",
                                           "mbf(TP.CxR)", 
                                           paste0("dev(",row.covar,")"), 
                                           paste0("dev(",col.covar,")"))),
                                         collapse = " + "),
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
    tspl.asrt <- do.call(changeModelOnIC, 
                         args = c(list(tspl.asrt, 
                                       addFixed = paste(paste0(
                                         sect.fac, 
                                         c("TP.CR.2", "TP.CR.3", "TP.CR.4")), 
                                         collapse = " + "),
                                       dropFixed = drop.fix, 
                                       addRandom = paste(paste0(
                                         sect.fac, 
                                         c("grp(TP.C.1_frow)", "grp(TP.C.2_frow)", 
                                           "grp(TP.R.1_fcol)", "grp(TP.R.2_fcol)",
                                           "grp(TP_fcol_frow)", 
                                           paste0("dev(",row.covar,")"), 
                                           paste0("dev(",col.covar,")"))), 
                                         collapse = " + "),
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
                       row.factor = NULL, col.factor = NULL, 
                       nsegs = NULL, asreml.opt = "mbf", 
                       tpps4mbf.obj = NULL, 
                       allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                       checkboundaryonly = FALSE, update = TRUE, 
                       IClikelihood = "full", which.IC = "AIC",
                       ...)
{ 
  stub = "xx"
  
  #Get data for mbf 
  if (asreml.opt == "mbf")
  {
    #stop('Sorry, but the mbf setting of asreml.opt is not functioning yet - use asreml.opt = "grp".')
    assign("tps.XZmat", tpps4mbf.obj)
    #Build the data.frame to be used in the analysis
    dat <- lapply(tpps4mbf.obj, function(mat) mat$data.plus)
    if (length(tpps4mbf.obj) > 1)
      dat <- do.call(rbind, dat)
    else
      dat <- dat[[1]]
    checkNamesInData(c(sections, row.covar, col.covar, row.factor, col.factor), dat)
    nsect <- calc.nsect(dat, sections)
  } else #doing grp
  {  
    #Check that named columns are in the data
    dat.in <- asrtests.obj$asreml.obj$call$data
    if (is.symbol(dat.in))
      dat.in <- eval(dat.in)
    checkNamesInData(c(sections, row.covar, col.covar, row.factor, col.factor), dat.in)
    
    nsect <- calc.nsect(dat.in, sections)
    
    #Make sure that row covar is not centred so that tpsmmb does not create a faulty index
    row.min <- min(dat.in[[row.covar]], na.rm = TRUE)
    if (row.min < 0)
      dat.in[row.covar] <- dat.in[[row.covar]] - row.min + 1
    
    col.min <- min(dat.in[[col.covar]], na.rm = TRUE)
    if (col.min < 0)
      dat.in[col.covar] <- dat.in[[col.covar]] - col.min + 1
    
    #Check conformability of covars and factors
    if (!is.null(row.factor) && nlevels(dat.in[[row.factor]]) != length(unique(dat.in[[row.covar]])))
      stop(row.factor, " does not have the same number of levels as there are values of ", row.covar)
    if (!is.null(col.factor) && nlevels(dat.in[[col.factor]]) != length(unique(dat.in[[col.covar]])))
      stop(col.factor, " does not have the same number of levels as there are values of ", col.covar)
    
    #Spatial local spatial model using tensor P-splines with the grp option 
    rc.cols <- c(sections, row.covar, col.covar)
    dat.rc <- dat.in[rc.cols]
    dat.rc <- unique(dat.rc)
    tps.XZmat <- addPSdesign.mat(dat.rc, sections = sections, nsect = nsect, 
                                 row.coords = row.covar, col.coords = col.covar, 
                                 nsegs = nsegs, asreml.opt = "grp", 
                                 stub = "xx", ...)
    
    #Build the data.frame to be used in the analysis
    kextra <- ncol(dat.in) - length(rc.cols)
    if (nsect == 1)
      dat <- tps.XZmat[[1]]$data
    else
    { 
      #Check, and if necessary, make data from different sections conformable
      tps.XZmat <- conformTPSSections(tps.XZmat)
      dat <- lapply(tps.XZmat, function(mat) mat$data)
      dat <- do.call(rbind, dat)
    }
    
    #Remove any previous Tensor Spline basis columns from the dat.in
    if (any(grepl("TP\\.", names(dat.in))) || any(grepl("TP\\_", names(dat.in))))
      dat.in <- dat.in[,-c(grep("TP\\.", names(dat.in)), grep("TP\\_", names(dat.in)))]
    dat <- merge(dat.in, dat, all.x = TRUE, by = rc.cols, sort = FALSE)
    
    #Adjust the grp columns for merging
    if (nsect == 1)
    { 
      tps.XZmat[[1]]$grp <- lapply(tps.XZmat[[1]]$grp, function(grp, kextra) grp <- grp + kextra, kextra = kextra)
    }
    else
      tps.XZmat <- lapply(tps.XZmat, 
                          function(mat, kextra) 
                          {
                            mat$grp <- lapply(mat$grp, 
                                              function(grp, kextra) grp <- grp + kextra, kextra = kextra)
                            return(mat)
                          }, kextra = kextra)
  }
  
  
  #Update the asreml.obj for the new data.frame
  asreml.obj  <- asrtests.obj$asreml.obj
  asreml.obj <- asreml::update.asreml(asreml.obj, data = dat)
  tspl.asrt <- as.asrtests(asreml.obj = asreml.obj, NULL, NULL, 
                           IClikelihood = "full", label = "Change to new data.frame with TPS bits")

  #Fit spatial TPPS to sections
  for (i in 1:nsect)
  {
    if (nsect == 1)
    { 
      sect.fac <- NULL
      lab <- paste0("Try tensor P-splines")
    }
    else
    { 
      sect.fac <- paste0("at(", sections, ",  ", i, "):")
      lab <- paste0("Try tensor P-splines for ", sections, " ",i)
    }
    tspl.asrt <- fitTPSModSect(tspl.asrt, mat = tps.XZmat[[i]], sect.fac = sect.fac, 
                               row.factor = row.factor, col.factor = col.factor, 
                               row.covar = row.covar, col.covar = col.covar, 
                               lab = lab, asreml.opt = asreml.opt, stub = stub, 
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation,
                               checkboundaryonly = checkboundaryonly, 
                               update = update, 
                               IClikelihood = IClikelihood, 
                               which.IC = which.IC, 
                               ...)
  }
  
  return(tspl.asrt)
}

