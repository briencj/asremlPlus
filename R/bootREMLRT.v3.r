bootREMLRT.asreml <- function(h0.asreml.obj, h1.asreml.obj, 
                              nboot = 100, max.retries = 5, seed = NULL, 
                              means=NULL, V = NULL, extra.matrix = NULL, ignore.terms = NULL, 
                              fixed.spline.terms = NULL, 
                              bound.exclusions = c("F","B","S","C"), 
                              tolerance = 1E-10, update = TRUE, trace = FALSE, 
                              ncores = detectCores(), ...)
{
  #asreml codes:
  #  B - fixed at a boundary (!GP)
  #  F - fixed by user
  #  ? - liable to change from P to B    
  #  P - positive definite
  #  C - Constrained by user (!VCC)      
  #  U - unbounded
  #  S - Singular Information matrix
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  
  #Check that have a valid objects of class asreml
  validasr <- validAsreml(h0.asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  #Check that have a valid object of class asreml
  validasr <- validAsreml(h1.asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  #Check that fixed and sparse models are the same
  if (asr4)
  {
    fixed.labels <- lapply(list(h1.asreml.obj,h0.asreml.obj), 
                           function(x) {attr(terms(x$formulae$fixed), "term.labels")})
    sparse.labels <- lapply(list(h1.asreml.obj,h0.asreml.obj), 
                            function(x) {attr(terms(x$formulae$sparse), "term.labels")})
    mu <- sapply(list(h1.asreml.obj,h0.asreml.obj), 
                 function(x) {attr(terms(x$formulae$fixed), "intercept")})
  } else
  {
    fixed.labels <- lapply(list(h1.asreml.obj,h0.asreml.obj), 
                           function(x) {attr(terms(x$fixed.formula), "term.labels")})
    sparse.labels <- lapply(list(h1.asreml.obj,h0.asreml.obj), 
                            function(x) {attr(terms(x$sparse), "term.labels")})
    mu <- sapply(list(h1.asreml.obj,h0.asreml.obj), 
                 function(x) {attr(terms(x$fixed.formula), "intercept")})
  }
  
  if (!all(mu == mu[1])) 
    stop("Fixed models must be identical with respect to the intercept")
  if (!all(diff(sapply(fixed.labels, function(x) length(x))) == 0)) 
    stop("Fixed models differ in length")
  if (all(sapply(fixed.labels, function(x) length(x)) > 0))
  { 
    for (i in 2:length(fixed.labels)) 
    { 
      if (!all(is.element(fixed.labels[[1]], fixed.labels[[i]]))) 
        stop("Fixed models differ")
    }
  }
  if (!all(diff(sapply(sparse.labels, function(x) length(x))) == 0)) 
    stop("sparse models differ in length")
  if (all(sapply(sparse.labels, function(x) length(x)) > 0)) 
  { 
    for (i in 2:length(sparse.labels)) 
    { 
      if (!all(is.element(sparse.labels[[1]], sparse.labels[[i]]))) 
        stop("sparse models differ")
    }
  }
  
  #Set up expectation
  if (is.null(means))
  {
    if (!asr4 & !is.null(fixed.spline.terms))
      stop("Cannot incorporate fixed spline terms in asreml ver. 3 - must supply means yourself")
    if (asr4)
    {
      mu <- fitted.values(h0.asreml.obj)
      if (!is.null(fixed.spline.terms)) #add fixed spline terms
      {
        #Check whether have the design matrix
        if (is.null(h0.asreml.obj$design))
        {
          asreml::asreml.options(design = TRUE)
          h0.asreml.obj <- eval(h0.asreml.obj$call)
        }
        for (term in fixed.spline.terms)
        {
          if (!grepl("spl", term, fixed = TRUE))
            stop("A term in fixed.spline.terms does not include the spline function")
          Z <- as.matrix(getTermDesignMatrix(term, h0.asreml.obj))
          G.param <- h0.asreml.obj$G.param
          termvars <- names(G.param[[term]])[-1]
          cond.fac <- ""
          cols <- ""
          for (var in termvars)
          {
            if (cols == "")
              cols <- paste(var, "_", G.param[[term]][[var]]$levels, sep ="")
            else
              cols <- paste(cols, ":", var, "_", G.param[[term]][[var]]$levels, sep ="")
          }
          BLUPs <- h0.asreml.obj$coefficients$random
          BLUPs <- BLUPs[match(cols, rownames(BLUPs)), 1]
          mu <- mu + Z %*% BLUPs
        }
      }
    }
    else
      mu <- fitted(h0.asreml.obj) #asreml::fitted.asreml(h0.asreml.obj)
  }
  else
  {
    mu <- means
  }
  
  if (is.null(V))
  {
    if (!asr4)
      stop("It appears that asreml 4.x is not loaded and so V must be supplied to bootREMLRT.asreml.")
    V <- estimateV(h0.asreml.obj, extra.matrix = extra.matrix, ignore.terms = ignore.terms, 
                   fixed.spline.terms = fixed.spline.terms, 
                   bound.exclusions = bound.exclusions)
    n <- nrow(V)
  }
  else #check supplied V
  {
    if (!isSymmetric(V))
      stop("Variance matrix must be symmetric")
    n <- length(h0.asreml.obj$residuals)
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
  
  #check and get info in calls
  h1.call <- h1.asreml.obj$call
  h0.call <- h0.asreml.obj$call
  if (!("data" %in% names(h1.call) & "data" %in% names(h0.call)))
    stop("bootREMLRT.asreml assumes that data has been set in the two asreml objects")
  sim.dat <- eval(h1.call$data)
  h0.dat <- eval(h0.call$data)
  if (ncol(sim.dat) != ncol(h0.dat) || 
      nrow(suppressMessages(dplyr::anti_join(sim.dat,h0.dat)))!=0)
    stop("The data argument is not the same for the two asreml calls")
  #check that the y var is in both asreml objects and has the same values
  h1.y <- languageEl(h1.call$fixed, which = 2)
  if (1 - (abs(var(sim.dat[[as.character(h1.y)]], na.rm = TRUE) / 
               var(h0.dat[[as.character(h1.y)]], na.rm = TRUE))) > 1e-03)
    stop("The values of y in data for the two asreml objects appear not to be the same")
  
  #set up calls
  h1.object.sim <- h1.asreml.obj
  h1.call <- h1.object.sim$call
  elno <- grep("data", names(h1.call))
  languageEl(h1.call, which=elno) <- quote(sim.dat)
  languageEl(h1.call$fixed, which=2) <- quote(y.sim)
  h0.call$trace <- trace
  h0.object.sim <- h0.asreml.obj
  h0.call <- h0.object.sim$call
  elno <- grep("data", names(h0.call))
  languageEl(h0.call, which=elno) <- quote(sim.dat)
  languageEl(h0.call$fixed, which=2) <- quote(y.sim)
  h0.call$trace <- trace
  
  #Deal with the R.param and G.param arguments
  if (update)
  { 
    #If update, set R.param and G.param
    languageEl(h1.call, which = "R.param") <- h1.asreml.obj$R.param
    languageEl(h1.call, which = "G.param") <- h1.asreml.obj$G.param
    languageEl(h0.call, which = "R.param") <- h0.asreml.obj$R.param
    languageEl(h0.call, which = "G.param") <- h0.asreml.obj$G.param
  } else
  { 
    #If R.param and G.param already set, make them NULL
    if (!is.null(languageEl(h1.call, which = "R.param")))
      languageEl(h1.call, which = "R.param") <- NULL
    if (!is.null(languageEl(h1.call, which = "G.param")))
      languageEl(h1.call, which = "G.param") <- NULL
    if (!is.null(languageEl(h0.call, which = "R.param")))
      languageEl(h0.call, which = "R.param") <- NULL
    if (!is.null(languageEl(h0.call, which = "G.param")))
      languageEl(h0.call, which = "G.param") <- NULL
  }
  
  #deal with args coming via ...
  #- this will overwite previously set values, except data and models are protectd
  tempcall <- list(...)
  if (length(tempcall)) 
  { 
    for (z in names(tempcall))
      if (z == "data")
        sim.dat <- data
    else
    {
      if (z %in% c("fixed", "random", "residual", "rcov", "sparse"))
        stop("attempt to change model to be fitted")
      else
      {
        languageEl(h0.call, which = z) <- tempcall[[z]]
        languageEl(h1.call, which = z) <- tempcall[[z]]
      }
    }
  }
  
  #Get bound values
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  if (asr4)
  {
    if (asr4.2)
    {
      bound.h0 <- h0.asreml.obj$vparameters.con
      bound.h1 <- h1.asreml.obj$vparameters.con
    } else
    { 
      bound.h0 <- vpc.char(h0.asreml.obj)
      bound.h1 <- vpc.char(h1.asreml.obj)
    }
  }
  else
  {
    bound.h0 <- names(h0.asreml.obj$gammas.con)
    names(bound.h0) <- names(h0.asreml.obj$gammas)
    bound.h1 <- names(h1.asreml.obj$gammas.con)
    names(bound.h1) <- names(h1.asreml.obj$gammas)
  }
  #Calculate observed p-value
  DF.diff <- DFdiff(bound.h1, bound.h0, bound.exclusions = bound.exclusions)	
  DF <- DF.diff$DF
  NBound.h1 <- DF.diff$NBound.h1
  NBound.h0 <- DF.diff$NBound.h0
  
  REMLRT <- 2*(h1.asreml.obj$loglik-h0.asreml.obj$loglik)
  
  #Set up for simulation
  conv <- FALSE
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
  #generate nboot data sets and analyse them
  if (asr4)
  {
    if (asr4.2)
    {
      kthreads <-   get("asr_options", envir = getFromNamespace(".asremlEnv", "asreml"))$threads
      asreml::asreml.options(threads = 1)
    }
    REMLRT.out <- foreach (i = 1:nboot, .combine = rbind, .inorder=FALSE,
                           .packages = c("asreml","asremlPlus"))  %dopar%
      { 
        asreml::asreml.options(fail = "soft")
        while (!conv)
        { 
          nnonconv <- 0
          sim.dat <- within(sim.dat, 
                            {y.sim <- as.vector(mu + R %*% rnorm(n))})
          h1.asr <- eval(h1.call)
          conv <- h1.asr$converge
          if (conv)
          {
            h0.asr <- eval(h0.call)
            conv <- h0.asr$converge
          }
          if (!conv)
          {
            nnonconv <- nnonconv + 1
            REMLRT.sim <- NA
            if (max.retries > nnonconv)
              break()
          } else
          {
            #Perform the test
            REMLRT.sim <- 2*(h1.asr$loglik-h0.asr$loglik)
          }
        }
        conv=FALSE
        cbind(REMLRT.sim, nnonconv)
      }
    if (asr4.2)
      asreml::asreml.options(threads = kthreads)
  } else
  {
    REMLRT.out <- foreach (i = 1:nboot, .combine = rbind, .inorder=FALSE,
                           .packages = c("asreml","asremlPlus"))  %dopar%
      { 
        loadASRemlVersion(3, lib.loc = "D:\\Analyses\\R oldpkg")
        while (!conv)
        { 
          sim.dat <- within(sim.dat, 
                            {y.sim <- as.vector(mu + R %*% rnorm(n))})
          h1.asr <- eval(h1.call)
          conv <- h1.asr$converge
          if (conv)
          {
            h0.asr <- eval(h0.call)
            conv <- h0.asr$converge
          }
          if (!conv)
          {
            nnonconv <- nnonconv + 1
            REMLRT.sim <- NA
            if (max.retries > nnonconv)
              break()
          } else
          {
            #Perform the test
            REMLRT.sim <- 2*(h1.asr$loglik-h0.asr$loglik)
          }
        }
        conv=FALSE
        cbind(REMLRT.sim, nnonconv)
      }
  }
  stopCluster(cl)
  
  REMLRT.sim <- na.omit(REMLRT.out[,1])
  nunconverged <- REMLRT.out[,2]
  totalunconverged <- sum(nunconverged)
  ksim <- length(REMLRT.sim)
  
  #Calculate bootstrap p-value
  p <- (sum(REMLRT.sim >= REMLRT) + 1) / (ksim + 1)
  list(REMLRT = REMLRT, p = p, DF = DF, totalunconverged = totalunconverged, 
       REMLRT.sim = REMLRT.sim, nunconverged = nunconverged)
}

