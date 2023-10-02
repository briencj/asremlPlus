"get.residuals" <- function(asreml.obj, units="ignore")
{ 
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)

  options <- c("addtoresiduals", "ignore")
  unit.opt <- options[check.arg.values(units, options)]
  if (length(asreml.obj$aom) == 0)
  { 
    res <- as.vector(asreml.obj$residuals)
    attr(res, which="restype") <- "Residuals"
    stdres <- FALSE
  } else
  { 
    res <- as.vector(asreml.obj$aom$R[,2])
    attr(res, which="restype") <- "Standardized conditional residuals"
    stdres <- TRUE
  }
  if (unit.opt == "addtoresiduals")
  { 
    #Check whether units is in random model
#    if (asr4)
      termno <- findterm("units", rownames(summary(asreml.obj)$varcomp))
#    else
#      termno <- findterm("units", rownames(summary(asreml.obj)$varcomp))
    if (termno > 0)
    { 
      if (asr4)
      {
        if (stdres)
          uBLUP <- asreml.obj$aom$G[,2] 
        else
          uBLUP <- asreml.obj$coefficients$random
      } else
      {
        if (!stdres)
          ucoeff <- 1
        else
        { 
#          ucoeff <- summary(asreml.obj)$varcomp[termno, "component"]
          ucoeff <- summary(asreml.obj)$varcomp[termno, "component"]
          ucoeff <- sqrt(ucoeff/(ucoeff^2 -1))
        }  
        uBLUP <- ucoeff*asreml.obj$coefficients$random
      }
      uBLUP <- uBLUP[grep("units_", substr(names(uBLUP), 1, 6), fixed=TRUE)]
      names(uBLUP) <- NULL
      res <- res + uBLUP
    }
  }
  res[is.nan(res)] <- NA
  return(res)
}


"variofaces.asreml" <- function(asreml.obj, means=NULL, V = NULL, nsim=100, seed = NULL, 
                                extra.matrix = NULL, ignore.terms = NULL, 
                                fixed.spline.terms = NULL, 
                                bound.exclusions = c("F","B","S","C"), 
                                tolerance = 1E-10, units = "ignore", 
                                update = TRUE, trace = FALSE, 
                                graphics.device=NULL, 
                                ncores = 2, ...)
#function to do the face variogram plots, including envelopes, described by 
#Stefanova et al (2010)
#asreml.obj is an asreml object from a call to asreml in which the data argument 
#   must have been set.
#means is a vector of predictions for fixed terms - should include spline terms
#   if null, fitted values from fitted model are used.
#V is the fitted variance matrix i.e. having the appropriate pattern and values 
#   given the model fitted and the estimates of the parameters obtained
#nsim is the number of data sets to be simulated in obtaining the envelopes
#   Note that only the results from data sets that converge in the fitting are 
#        included in computing the envelope
#tolerance is the value such that eigenvalues less than it are consdered to be zero
#... parameters to supply to asreml calls within variofaces.asreml
{ 
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2)
  #Check that have a valid object of class asreml
  validasr <- validAsreml(asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  options <- c("addtoresiduals", "ignore")
  unit.opt <- options[check.arg.values(units, options)] 
  #Check V
  if (is.null(V))
  {
    if (!asr4)
      stop("It appears that asreml 4.x is not loaded and so V must be supplied to variofaces.asreml.")
    asreml.obj <- asreml.obj
    V <- estimateV(asreml.obj, extra.matrix = extra.matrix, ignore.terms = ignore.terms, 
                   fixed.spline.terms = fixed.spline.terms, 
                   bound.exclusions = bound.exclusions)
    n <- nrow(V)
  }
  else #check supplied V
  {
    if (!isSymmetric(V))
      stop("Variance matrix must be symmetric")
    n <- length(asreml.obj$residuals)
    if (!is.null(means) & length(means) != n)
      stop("The lengths of means  and the response variable are not the same")
    if (!all(dim(V) == c(n, n))) 
      stop("V is not a square matrix whose order equals the length of the response variable")
  }

  #use eigenvalue decomposition to establish transformation matrix
  eigdecomp <- eigen(V, symmetric = TRUE)
  eigenval <- eigdecomp$values
  if (!all(eigenval >= -tolerance * abs(max(eigenval))))
     stop("Variance matrix is not nonnegative definite")
  R <- eigdecomp$vectors %*% diag(sqrt(pmax(eigenval, 0)), n)

 
  #check and get info in supplied call
  call <- asreml.obj$call
  if (!("data" %in% names(call)))
    stop("variofaces.asreml assumes that data has been set in call to asreml")
  env.dat <- eval(call$data)
 
  #set up call
  object.sim <- asreml.obj
  call <- object.sim$call
  elno <- grep("data", names(call))
  languageEl(call, which=elno) <- quote(env.dat)
  languageEl(call$fixed, which=2) <- quote(y.sim)
  call$trace <- trace

 #Deal with the R.param and G.param arguments
 if (update)
 { 
   #If update, set R.param and G.param
   languageEl(call, which = "R.param") <- asreml.obj$R.param
   languageEl(call, which = "G.param") <- asreml.obj$G.param
 } else
 { 
   #If R.param and G.param already set, make them NULL
   if (!is.null(languageEl(call, which = "R.param")))
     languageEl(call, which = "R.param") <- NULL
   if (!is.null(languageEl(call, which = "G.param")))
     languageEl(call, which = "G.param") <- NULL
 }
  
 #deal with args coming via ...
 #- this will overwite previously set values, except data and models are protectd
 tempcall <- list(...)
 if (length(tempcall)) 
 { 
   for (z in names(tempcall))
     if (z == "data")
       env.dat <- data
     else
       if (z %in% c("fixed", "random", "residual", "sparse"))
         stop("attempt to change model to be fitted")
     else
       languageEl(call, which = z) <- tempcall[[z]]
 }
 
 #investigate residual term to see if it is a two-factor term or a 
 #three-factor term with one factor indexing sections
 if (asr4)
   rterm <- languageEl(call, which = "residual")
 else
   rterm <- languageEl(call, which = "rcov")
 #Check for multiple terms
 rterm.form <- as.formula(rterm)
 res.specials <- c("dsum")
 common.specials <-  c("id", "diag", "us", 
                       "ar1", "ar2", "ar3", "sar","sar2",
                       "ma1", "ma2", "arma", "exp", "gau", 
                       "cor", "corb", "corg") 
 other.specials <- c("at", "sph", "chol", "ante", "sfa", "facv")
 rterm.obj <- as.terms.object(rterm, asreml.obj, 
                              specials=c(res.specials, common.specials, 
                                         paste(common.specials[4:length(common.specials)],"v",sep=""),
                                         paste(common.specials[4:length(common.specials)],"h",sep=""),
                                         other.specials))
 if (length(labels(rterm.obj)) != 1)
 {
   if (asr4)
     stop("In analysing ",asreml.obj$formulae$fixed[[2]],
          ", the residual model must involve a single term")
   else
     stop("In analysing ",asreml.obj$fixed.formula[[2]],
          ", the residual model must involve a single term")
 }
 if (asr4)
   kspecial <- attr(rterm.obj, which = "specials")$dsum
 else
 {
   if (length(labels(rterm.obj)) != 1)
     stop("In analysing ",object.sim$fixed.formula[[2]],
          ", the rcov model must involve a single term")
   kspecial <- attr(rterm.obj, which = "specials")$at
 }
 
 if (is.null(kspecial))
 {
   grid.facs <- rownames(attr(rterm.obj, which = "factors"))
   grid.facs <- unlist(lapply(grid.facs, rmFunction))
   sections <- 1
   fac.sec <- NULL
 } else
 {
   if (asr4)
   {
     rterm <- labels(rterm.obj)
     rterm <- substr(rterm, start=6, stop=nchar(rterm))
     rterm <- substr(rterm, start=1, stop =(nchar(rterm)-1))
     trim <- function (x) gsub("^\\s+|\\s+$", "", x)
     if (grepl(",", rterm, fixed=TRUE))
     {
       rterm <- trim(strsplit(rterm, ",", fixed=TRUE)[[1]])
     }
     kmodel <- grep("~", rterm, fixed=TRUE)
     if (length(kmodel) != 1)
       stop("residual model appears to lack a formula")
     rterm <- rterm[kmodel]
     rterm <- strsplit(rterm, "~", fixed=TRUE)[[1]][2]
     var <- trim(strsplit(rterm, "|", fixed=TRUE)[[1]])
     var <- trim(var)
     if (length(var) != 2)
       stop("residual model appears not to include `|' operator")
     fac.sec <- var[2]  
     sections <- levels(env.dat[[fac.sec]])
     rterm <- as.formula(paste("~",var[1],sep=""))
     rterm.obj <- as.terms.object(rterm, asreml.obj, 
                                  specials=c("dsum", "at", "ar1", "ar2", "ar3", "sar","sar2",
                                             "ma1", "ma2", "arma", "exp", "gau", 
                                             "cor", "corb", "corg", "diag", "us", 
                                             "chol", "ante", "sfa", "facv", "fa", "rr"))
     if (length(labels(rterm.obj)) != 1)
       stop("In analysing ",asreml.obj$formulae$fixed[[2]],
            ", the residual model must involve a single term")
     grid.facs <- rownames(attr(rterm.obj, which = "factors"))
     grid.facs <- unlist(lapply(grid.facs, rmFunction))
   } else
   { 
     grid.facs <- rownames(attr(rterm.obj, which = "factors"))
     grid.facs <- unlist(lapply(grid.facs, rmFunction))
     if (length(kspecial) != 1)
       stop("Can only have a single factor defining sections of the data (using at)")
     fac.sec <- grid.facs[kspecial[1]]  
     sections <- levels(env.dat[[fac.sec]])
     grid.facs <- grid.facs[-kspecial[1]]
   }
 }
 if (length(grid.facs) != 2)
   stop("Must have two dimensions in the residual model")

 #Set up for simulation
 #There will be two lists, one for each face
 #The list for each face  has number of sections components
 #Each component will be an nsim x nfaci matrix, where i is the face 
 fac1 <- grid.facs[1]
 fac2 <- grid.facs[2]
 if (any(table(env.dat[[fac1]])==0) | any(table(env.dat[[fac2]])==0))
   stop("Some levels in the the supplied factors are not observed") 
 n1 <- length(levels(env.dat[[fac1]]))
 n2 <- length(levels(env.dat[[fac2]]))
 ns1 <- nsim*n1
 ns2 <- nsim*n2
 env.var <- vector(mode = "list", length = 2)
 names(env.var) <- grid.facs
 env.var[[1]] <- vector(mode = "list", length = length(sections))
 names(env.var[[1]]) <- sections
 env.var[[2]] <- vector(mode = "list", length = length(sections))
 names(env.var[[2]]) <- sections

 
 #generate nsim data sets and save variogram face data
  conv <- FALSE
  res.dat <- data.frame(env.dat[[fac1]],env.dat[[fac2]])
  names(res.dat)[1:2] <- grid.facs
  res.dat[grid.facs] <- lapply(res.dat[grid.facs], as.numfac)
  if (!is.null(fac.sec))
  { 
    res.dat <- data.frame(env.dat[[fac.sec]], res.dat)
    names(res.dat)[1] <- fac.sec
  } 
  #Set up expectation
  if (is.null(means))
  {
    if (asr4)
      mu <- fitted.values(asreml.obj)
    else
      mu <- fitted(asreml.obj) #asreml::fitted.asreml(asreml.obj)
  }
  else
    mu <- means
  
  #Get observed residuals
  res.dat$res <- get.residuals(asreml.obj, units = units)
  restype <- attr(res.dat$res, which = "restype")
  r <- res.dat$res
  
  #Set OMP_NUM_THREADS to 1 on Windows systems
  if (grepl("Windows", Sys.getenv("OS")))
  {
    kTHREADS <- Sys.getenv("OMP_NUM_THREADS")
    Sys.setenv("OMP_NUM_THREADS" = 1)
  }
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  #process seed argument
  if (!is.null(seed))
  {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    ## start ncores workers
    s <- .Random.seed
    for (i in 1:ncores) {
      s <- nextRNGStream(s)
      # send s to worker i as .Random.seed
    }
  }
  
  if (asr4)
  {
    if (asr4.2)
    {
      kthreads <-   get("asr_options", envir = getFromNamespace(".asremlEnv", "asreml"))$threads
      asreml::asreml.options(threads = 1)
    }
      env.var <- foreach(i = 1:nsim, .packages = c("asreml","asremlPlus"))  %dopar%
      { 
        asreml::asreml.options(fail = "soft")
        while (!conv)
        { 
          env.dat <- within(env.dat, 
                            {y.sim <- as.vector(mu + R %*% rnorm(n))})
          sim.asreml <- eval(call)
          conv <- sim.asreml$converge
        }
        conv=FALSE
        #Get residuals
        res.dat <- within(res.dat, 
                          {res <- get.residuals(sim.asreml, units = units)})
        if (restype != attr(res.dat$res, which = "restype"))
          warning("Observed and simulated residuals are not of the same type")
        
        
        #Get variogram faces for this simulated data set
        for (k in sections)
        { 
          if (is.null(fac.sec))
            sect.dat <- res.dat
          else
            sect.dat <- res.dat[res.dat[[fac.sec]]==k, -1]
          #      sect.dat <- as.matrix(sect.dat)
          
          sim.var <- asreml::asr_varioGram(sect.dat)
          names(sim.var)[1:2] <- grid.facs
          env.var[[1]][[k]] <- sim.var[sim.var[[fac2]]==0,]$gamma
          env.var[[2]][[k]] <- sim.var[sim.var[[fac1]]==0,]$gamma
        }
        env.var
      }
    if (asr4.2)
      asreml::asreml.options(threads = kthreads)
  } else
  {
    env.var <- foreach (i = 1:nsim, .packages = c("asreml","asremlPlus"))  %dopar%
      { 
                          loadASRemlVersion(3, lib.loc = "D:\\Analyses\\R oldpkg")
                          while (!conv)
                          { 
                            env.dat <- within(env.dat, 
                                              {y.sim <- as.vector(mu + R %*% rnorm(n))})
                            sim.asreml <- eval(call)
                            conv <- sim.asreml$converge
                          }
                          conv=FALSE
                          #Get residuals
                          res.dat$res <- get.residuals(sim.asreml, units = units)
                          if (restype != attr(res.dat$res, which = "restype"))
                            warning("Observed and simulated residuals are not of the same type")
                          #Get variogram faces for this simulated data set
                          for (k in sections)
                          { 
                            if (is.null(fac.sec))
                              sect.dat <- res.dat
                            else
                              sect.dat <- res.dat[res.dat[[fac.sec]]==k, -1]
                            #      sect.dat <- as.matrix(sect.dat)
                            sim.var <- asreml::asreml.variogram(sect.dat)
                            names(sim.var)[1:2] <- grid.facs
                            env.var[[1]][[k]] <- sim.var[sim.var[[fac2]]==0,]$gamma
                            env.var[[2]][[k]] <- sim.var[sim.var[[fac1]]==0,]$gamma
                          }
                          env.var
                        }
  }
  stopCluster(cl)
  if (grepl("Windows", Sys.getenv("OS")))
    Sys.setenv("OMP_NUM_THREADS" = kTHREADS)
  if (abs(mean(r -res.dat$res, na.rm = TRUE)) > 1e-6)
    stop("res.dat$res has changed")
  
  #Do face plots for fac1 and fac2 in each section
  face.1 <- data.frame(matrix(nrow = 0, ncol=5))
  colnames(face.1) <- c("factor","observed","X2.5.","X50.","X97.5.")
  face.2 <- data.frame(matrix(nrow = 0, ncol=5))
  colnames(face.2) <- c("factor","observed","X2.5.","X50.","X97.5.")
  
  #Get the observed variogram and form the elements of the variofaces plot
  for (k in sections)
  { 
    if (is.null(fac.sec))
      sect.dat <- res.dat
    else
      sect.dat <- res.dat[res.dat[[fac.sec]]==k, -1]
    #    sect.dat <- as.matrix(sect.dat)
    if (asr4)
      object.var <- asreml::asr_varioGram(sect.dat)
    else
      object.var <- asreml::asreml.variogram(sect.dat)
    names(object.var)[1:2] <- grid.facs
    
    #Form variogram for current section
    sect.var <- vector(mode = "list", length = 2)
    names(sect.var) <- grid.facs
    sect.var[[1]] <- as.matrix(as.data.frame(lapply(env.var, 
                                                    function(x, k) {x[[1]][[k]]}, 
                                                    k = k)))
    sect.var[[2]] <- as.matrix(as.data.frame(lapply(env.var, 
                                                    function(x, k) {x[[2]][[k]]}, 
                                                    k = k)))

    #Form data.frame for current section
    face.1 <- rbind(face.1,
                    data.frame(factor = object.var[object.var[[fac2]]==0,][[fac1]],
                               observed = object.var[object.var[[fac2]]==0,]$gamma,   
                               t(apply(sect.var[[1]], 1, quantile, 
                                       probs=c(0.025, 0.50, 0.975),
                                       na.rm = TRUE)),
                               stringsAsFactors = FALSE))
    face.2 <- rbind(face.2,
                    data.frame(factor = object.var[object.var[[fac1]]==0,][[fac2]],   
                               observed = object.var[object.var[[fac1]]==0,]$gamma,   
                               t(apply(sect.var[[2]], 1, quantile, 
                                       probs=c(0.025, 0.50, 0.975),
                                       na.rm = TRUE)), 
                               stringsAsFactors = FALSE))
  }
  names(face.1)[1] <- fac1
  names(face.2)[1] <- fac2
  
  if (!is.null(fac.sec))
  { 
    face.1[[fac.sec]] <- factor(rep(1:length(sections), each=(n1)), labels=sections)
    face.2[[fac.sec]] <- factor(rep(1:length(sections), each=(n2)), labels=sections)
  }
  
  #Do plots
  if (!is.null(graphics.device) )
    do.call(graphics.device, list(record = FALSE))
  
  p <- ggplot(data=face.1) +
    theme_bw() +
    geom_line(aes(x = .data[[!!fac1]], y = .data[["observed"]]), linewidth=1) +
    geom_point(aes(x = .data[[!!fac1]], y = .data[["observed"]]), size=3) +
    geom_line(aes(x = .data[[!!fac1]], y = .data[["X2.5."]]), colour = "red") +
    geom_line(aes(x = .data[[!!fac1]], y = .data[["X97.5."]]), colour = "red") +
    geom_line(aes(x = .data[[!!fac1]], y = .data[["X50."]]), colour = "blue") +
    labs(title = paste("Variogram face of",restype,"for",fac1,sep=" "), 
         x = paste(fac1,"differences", sep=" "),y=NULL)
  if (!is.null(fac.sec))
    p <- p + facet_wrap(as.formula(paste("~",fac.sec)), ncol=2, scales="free")
  print(p)
  
  p <- ggplot(data=face.2) +
    theme_bw() +
    geom_line(aes(x = .data[[!!fac2]], y = .data[["observed"]]), linewidth=1) +
    geom_point(aes(x = .data[[!!fac2]], y = .data[["observed"]]), size=3) +
    geom_line(aes(x = .data[[!!fac2]], y = .data[["X2.5."]]), colour = "red") +
    geom_line(aes(x = .data[[!!fac2]], y = .data[["X97.5."]]), colour = "red") +
    geom_line(aes(x = .data[[!!fac2]], y = .data[["X50."]]), colour = "blue") +
    labs(title = paste("Variogram face of",restype,"for",fac2,sep=" "), 
         x = paste(fac2,"differences", sep=" "), y=NULL)
  if (!is.null(fac.sec))
    p <- p + facet_wrap(as.formula(paste("~",fac.sec)), ncol=2, scales="free")
  print(p)
  
  #return data frames containing the variogram values on which the plots are based
  invisible(list(face1 = face.1, face2 = face.2))
}

"simulate.asreml" <- function(object, nsim=100, seed = NULL, means=NULL, V, 
                              tolerance = 1E-10, update = TRUE, trace = FALSE, 
                              which = "data", units = "ignore", ncores = 2, 
                              ...)
  #function to obtain simulated data corresponding to a fitted asreml model
  #object is an asreml object from a call to asreml in which the data argument 
  #   must have been set.
  #means is a vector of predictions for fixed terms - should include spline terms
  #   if null, fitted values from fitted model are used.
  #V is the fitted variance matrix i.e. having the appropriate pattern and values 
  #   given the model fitted and the estimates of the parameters obtained
  #nsim is the number of data sets to be simulated in obtaining the envelopes
  #   Note that only the results from data sets that converge in the fitting are 
  #        included in computing the envelope
  #tolerance is the value such that eigenvalues less than it are consdered to be zero
#... parameters to supply to plot functions called within variofaces.asreml
{ 
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  #Check that have a valid object of class asreml
  validasr <- validAsreml(object)  
  if (is.character(validasr))
    stop(validasr)
  
  options <- c("data", "fitted", "residuals", "all")
  opt <- options[unlist(lapply(which, check.arg.values, options=options))]
  if ("all" %in% opt)
    opt <-  c("data", "fitted", "residuals")
  options <- c("addtoresiduals", "ignore")
  unit.opt <- options[check.arg.values(units, options)]
  n <- length(object$residuals)
  if (!is.null(means) & length(means) != n)
    stop("The lengths of means  and the response variable are not the same")
  #Check V
  if (!isSymmetric(V))
    stop("Variance matrix must be symmetric")
  if (!all(dim(V) == c(n, n))) 
    stop("V is not a square matrix whose order equals the length of the response variable")
  #use eigenvalue decomposition to establish transformation matrix
  eigdecomp <- eigen(V, symmetric = TRUE)
  eigenval <- eigdecomp$values
  if (!all(eigenval >= -tolerance * abs(max(eigenval))))
    stop("Variance matrix is not nonnegative definite")
  R <- eigdecomp$vectors %*% diag(sqrt(pmax(eigenval, 0)), n)
  
  
  #check and get info in supplied call
  call <- object$call
  if (!("data" %in% names(call)))
    stop("simulate.asreml assumes that data has been set in call to asreml")
  env.dat <- eval(call$data)
  
  #set up call
  object.sim <- object
  call <- object.sim$call
  elno <- grep("data", names(call))
  languageEl(call, which=elno) <- quote(env.dat)
  languageEl(call$fixed, which=2) <- quote(y.sim)
  call$trace <- trace
  
  #Deal with the R.param and G.param arguments
  if (update)
  { 
    #If update, set R.param and G.param
    languageEl(call, which = "R.param") <- object$R.param
    languageEl(call, which = "G.param") <- object$G.param
  } else
  { 
    #If R.param and G.param already set, make them NULL
    if (!is.null(languageEl(call, which = "R.param")))
      languageEl(call, which = "R.param") <- NULL
    if (!is.null(languageEl(call, which = "G.param")))
      languageEl(call, which = "G.param") <- NULL
  }
  
  #deal with args coming via ...
  #- this will overwite previously set values, except data and models are protectd
  tempcall <- list(...)
  if (length(tempcall)) 
  { 
    for (z in names(tempcall))
      if (z == "data")
        env.dat <- data
      else
      {
        if (z %in% c("fixed", "random", "residual", "rcov", "sparse"))
          stop("attempt to change model to be fitted")
        else
          languageEl(call, which = z) <- tempcall[[z]]
      }
  }

  
  #Set up for simulation and get observed residuals and fitted values
  #There will be a data.fame for each quantity saved, as well as the supplied 
  #data.frame to which the requested quantites have been added
  if ("fitted"%in% opt)
  { 
    if (asr4)
      env.dat$fitted <- fitted.values(object)
    else
      env.dat$fitted <- fitted(object)
#      env.dat$fitted <- asreml::fitted.asreml(object)
  } 
  if ("residuals"%in% opt)
  { 
    #Get observed residuals
    env.dat$residuals <- get.residuals(object, units = units)
    restype <- attr(env.dat$residuals, which = "restype")
  }
  
  #Set up expectation
  if (is.null(means))
  {
    if (asr4)
      mu <- fitted.values(object)
    else
      mu <- fitted(object)
#      mu <- asreml::fitted.asreml(object)
  }
  else
    mu <- means
  
  #Check type of residuals in simulation
  env.dat$y.sim <- as.vector(mu + R %*% rnorm(n))
  sim.asreml <- eval(call)
  if ("residuals" %in% opt)
    if (((length(sim.asreml$aom) == 0) & (restype == "Standardized conditional residuals")) | 
          ((length(sim.asreml$aom) != 0) & (restype == "Residuals")))
      warning("Observed and simulated residuals are not of the same type\n",
              "- check setting of aom")
  
  #generate nsim data sets and save variogram face data
  conv <- FALSE
  #Set OMP_NUM_THREADS to 1 on Windows systems
  if (grepl("Windows", Sys.getenv("OS")))
  {
    kTHREADS <- Sys.getenv("OMP_NUM_THREADS")
    Sys.setenv("OMP_NUM_THREADS" = 1)
  }
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  #process seed argument
  if (!is.null(seed))
  {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    ## start ncores workers
    s <- .Random.seed
    for (i in 1:ncores) {
      s <- nextRNGStream(s)
      # send s to worker i as .Random.seed
    }
  }
  if (asr4)
  {
    if (asr4.2)
    {
      kthreads <-   get("asr_options", envir = getFromNamespace(".asremlEnv", "asreml"))$threads
      asreml::asreml.options(threads = 1)
    }
    sim <- foreach(i = 1:nsim, .combine=rbind,
                   .packages = c("asreml","asremlPlus"))  %dopar%
      { 
        asreml::asreml.options(fail = "soft")
        while (!conv)
        { 
          env.dat <- within(env.dat, 
                            {y.sim <- as.vector(mu + R %*% rnorm(n))})
          sim.asreml <- eval(call)
          conv <- sim.asreml$converge
        }
        conv=FALSE
        sim <- vector(mode = "list", length = length(opt))
        names(sim) <- opt
        if ("data"%in% opt)
          sim[["data"]] <- env.dat$y.sim
        #Get residuals
        if ("residuals"%in% opt)
          sim[["residuals"]] <- get.residuals(sim.asreml, units=units)
        if ("fitted"%in% opt)
          sim[["fitted"]] <- fitted.values(sim.asreml)
        sim
      }
    if (asr4.2)
      asreml::asreml.options(threads = kthreads)
  } else
  {
    sim <- foreach(i = 1:nsim, .combine=rbind,
                   .packages = c("asreml","asremlPlus"))  %dopar%
      { 
        loadASRemlVersion(3, lib.loc = "D:\\Analyses\\R oldpkg")
        while (!conv)
        { 
          env.dat <- within(env.dat, 
                            {y.sim <- as.vector(mu + R %*% rnorm(n))})
          sim.asreml <- eval(call)
          conv <- sim.asreml$converge
        }
        conv=FALSE
        sim <- vector(mode = "list", length = length(opt))
        names(sim) <- opt
        if ("data"%in% opt)
          sim[["data"]] <- env.dat$y.sim
        #Get residuals
        if ("residuals"%in% opt)
          sim[["residuals"]] <- get.residuals(sim.asreml, units=units)
        if ("fitted"%in% opt)
          sim[["fitted"]] <- fitted(sim.asreml)
        #                        sim[["fitted"]] <- asreml::fitted.asreml(sim.asreml)
        sim
      }
  }
  stopCluster(cl)
  if (grepl("Windows", Sys.getenv("OS")))
    Sys.setenv("OMP_NUM_THREADS" = kTHREADS)
  
  out <- vector("list", length = 0)
  if ("residuals" %in% opt | "fitted" %in% opt)
  { 
    env.dat <- env.dat[,-match("y.sim",names(env.dat))]
    out[["observed"]] <- env.dat
  }
  if ("data" %in% opt)
  {
    out[["data"]] <- as.data.frame(lapply(1:nsim, 
                                          function(k, sim){sim[k, "data"]}, 
                                          sim = sim))
    names(out[["data"]]) <- paste("data",1:nsim,sep=".")
  }
  
  if ("fitted" %in% opt)
  {
    out[["fitted"]] <- as.data.frame(lapply(1:nsim, 
                                            function(k, sim){sim[k, "fitted"]}, 
                                            sim = sim))
    names(out[["fitted"]]) <- paste("fitted",1:nsim,sep=".")
  }
  if ("residuals" %in% opt)
  {
    out[["residuals"]] <- as.data.frame(lapply(1:nsim, 
                                               function(k, sim){sim[k, "residuals"]}, 
                                               sim = sim))
    names(out[["residuals"]]) <- paste("residuals",1:nsim,sep=".")
  }
  #return data frames containing the variogram values on which the plots are based
  invisible(out)
}


"plotVariofaces.data.frame" <- function(data, residuals, restype="Residuals", ...)
  #function to do the face variogram plots, including envelopes, described by 
  #Stefanova et al (2010)
  #data is a data.frame with 3 or 4 columns. If there are 4, the first is taken to index
  #     sections of the data, the second and the third to index the x and y dimensions 
  #     underlying the observations and the fourth the observed residuals. If there are 3, 
  #     then the same as for 4 columns, except that the factor index the sections has 
  #     been omitted.
  #residuals is a data frame with the same initial columns as data, the observed residuals 
  #     having been ommitted. Following the initial columns are nsim sets of residuals 
  #     derived from simulated data. 
  #restype is a character string describing the type of residuals supplied.
{ 
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2)
  
  if (!is.data.frame(data) | !is.data.frame(residuals))
    stop("Both data and residuals should be data frames")
  n <- nrow(data)
  nvars <- ncol(data)
  if (nvars < 3 & nvars > 4)
    stop("The number of columns in data is not 3 or 4")
  nsim <- ncol(residuals) - nvars + 1
  if (nsim < 1)
    stop("Could not find columns with simulated residuals")
  #Check compatibility
  if (n != nrow(residuals))
    stop("The number of rows in data and residuals are not equal")
  
  #Check that initial columns in data and residuals match
  if (!setequal(names(data)[1:(nvars-1)], names(residuals)[1:(nvars-1)]))
    stop("The names of the initial ",nvars-1," columns do not match")
  if (!all(unlist(lapply(data[1:(nvars-1)], is.factor))))
    stop("Some of the initial ",nvars-1," columns of data are not factors")
  if (!all(unlist(lapply(residuals[1:(nvars-1)], is.factor))))
    stop("Some of the initial ",nvars-1," columns of residuals are not factors")
  for (i in 1:(nvars-1))
    if (!all(data[i] == residuals[i]))
      stop("The values in column ",i," of data and residuals do not match")
  
  #Identify columns in data
  if (nvars == 3)
  { 
    fac.sec <- NULL
    sections <- 1
    grid.facs <- colnames(data)[1:2]
    resid <- colnames(data)[3]
  } 
  else
  { 
    fac.sec <- colnames(data)[1]
    sections <- levels(data[[fac.sec]])
    grid.facs <- colnames(data)[2:3]
    resid <- colnames(data)[4]
  }
  
  #Form variogram faces for each set of residuals
  #There will be two lists, one for each face
  #The list for each face  has number of sections components
  #Each component will be an nsim x nfaci matrix, where i is the face 
  fac1 <- grid.facs[1]
  fac2 <- grid.facs[2]
  if (any(table(data[[fac1]])==0) | any(table(data[[fac2]])==0))
    stop("Some levels in the the supplied factors are not observed") 
  n1 <- length(levels(data[[fac1]]))
  n2 <- length(levels(data[[fac2]]))
  ns1 <- nsim*n1
  ns2 <- nsim*n2
  variog <- vector(mode = "list", length = 2)
  variog[[1]] <- vector(mode = "list", length = length(sections))
  names(variog[[1]]) <- sections
  variog[[1]] <- lapply(variog[[1]][sections], 
                        function(variog, nrow, ncol)
                          variog <- matrix(rep(0, nrow*ncol), nrow=nrow, ncol=ncol),
                        nrow=nsim, ncol=n1)
  variog[[2]] <- vector(mode = "list", length = length(sections))
  names(variog[[2]]) <- sections
  variog[[2]] <- lapply(variog[[2]][sections], 
                        function(variog, nrow, ncol)
                          variog <- matrix(rep(0, nrow*ncol), nrow=nrow, ncol=ncol),
                        nrow=nsim, ncol=n2)
  data[grid.facs] <- lapply(data[grid.facs], as.numfac) 
  residuals[grid.facs] <- lapply(residuals[grid.facs], as.numfac) 
  
  for (k in sections)
  { 
    sect.dat <- data.frame(data[fac1], data[fac2])
    if (!is.null(fac.sec))
      sect.dat <- sect.dat[data[[fac.sec]] == k]
    names(sect.dat)[1:2] <- grid.facs
    for (i in 1:nsim)
    { 
      if (is.null(fac.sec))
        sect.dat$res <- as.vector(residuals[[(nvars+i-1)]])
      else
        sect.dat$res <- residuals[data[[fac.sec]] == k, (nvars+i-1)]
      if (asr4)
        sim.var <- asreml::asr_varioGram(sect.dat)
      else
        sim.var <- asreml::asreml.variogram(sect.dat)
      
      names(sim.var)[1:2] <- grid.facs
      variog[[1]][[k]][i,] <- sim.var[sim.var[[fac2]]==0,]$gamma
      variog[[2]][[k]][i,] <- sim.var[sim.var[[fac1]]==0,]$gamma
    }
  }
  
  #Do face plots for fac1 and fac2 in each section
  face.1 <- data.frame(matrix(nrow = 0, ncol=5))
  colnames(face.1) <- c("factor","observed","X2.5.","X50.","X97.5.")
  face.2 <- data.frame(matrix(nrow = 0, ncol=5))
  colnames(face.2) <- c("factor","observed","X2.5.","X50.","X97.5.")
  
  #Get the observed variogram and form the elements of the variofaces plot
  for (k in sections)
  { 
    if (is.null(fac.sec))
      sect.dat <- data
    else
      sect.dat <- data[data[[fac.sec]]==k, -1]
    if (asr4)
      object.var <- asreml::asr_varioGram(sect.dat)
    else
      object.var <- asreml::asreml.variogram(sect.dat)
    names(object.var)[1:2] <- grid.facs
    
    #Form data.frame for current section
    face.1 <- rbind(face.1,
                    data.frame(factor = object.var[object.var[[fac2]]==0,][[fac1]],
                               observed = object.var[object.var[[fac2]]==0,]$gamma,   
                               t(apply(variog[[1]][[k]], 2, quantile, 
                                       probs=c(0.025, 0.50, 0.975),
                                       na.rm = TRUE)),
                               stringsAsFactors = FALSE))
    face.2 <- rbind(face.2,
                    data.frame(factor = object.var[object.var[[fac1]]==0,][[fac2]],   
                               observed = object.var[object.var[[fac1]]==0,]$gamma,   
                               t(apply(variog[[2]][[k]], 2, quantile, 
                                       probs=c(0.025, 0.50, 0.975),
                                       na.rm = TRUE)), 
                               stringsAsFactors = FALSE))
  }
  names(face.1)[1] <- fac1
  names(face.2)[1] <- fac2
  
  if (!is.null(fac.sec))
  { 
    face.1[[fac.sec]] <- factor(rep(1:length(sections), each=(n1)), labels=sections)
    face.2[[fac.sec]] <- factor(rep(1:length(sections), each=(n2)), labels=sections)
  }
  
  #Do plots
  p <- ggplot(data=face.1) +
    theme_bw() +
    geom_line(aes(x = .data[[!!fac1]], y = .data[["observed"]]), linewidth=1) +
    geom_point(aes(x = .data[[!!fac1]], y = .data[["observed"]]), size=3) +
    geom_line(aes(x = .data[[!!fac1]], y = .data[["X2.5."]]), colour = "red") +
    geom_line(aes(x = .data[[!!fac1]], y = .data[["X97.5."]]), colour = "red") +
    geom_line(aes(x = .data[[!!fac1]], y = .data[["X50."]]), colour = "blue") +
    labs(title = paste("Variogram face of",restype,"for",fac1,sep=" "), 
         x = paste(fac1,"differences", sep=" "),y=NULL)
  if (!is.null(fac.sec))
    p <- p + facet_wrap(as.formula(paste("~",fac.sec)), ncol=2, scales="free")
  print(p)
  
  p <- ggplot(data=face.2) +
    theme_bw() +
    geom_line(aes(x = .data[[!!fac2]], y = .data[["observed"]]), linewidth=1) +
    geom_point(aes(x = .data[[!!fac2]], y = .data[["observed"]]), size=3) +
    geom_line(aes(x = .data[[!!fac2]], y = .data[["X2.5."]]), colour = "red") +
    geom_line(aes(x = .data[[!!fac2]], y = .data[["X97.5."]]), colour = "red") +
    geom_line(aes(x = .data[[!!fac2]], y = .data[["X50."]]), colour = "blue") +
    labs(title = paste("Variogram face of",restype,"for",fac2,sep=" "), 
         x = paste(fac2,"differences", sep=" "), y=NULL)
  if (!is.null(fac.sec))
    p <- p + facet_wrap(as.formula(paste("~",fac.sec)), ncol=2, scales="free")
  print(p)
  
  #return data frames containing the variogram values on which the plots are based
  invisible(list(face1 = face.1, face2 = face.2))
}

