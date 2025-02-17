DFdiff <- function(bound.h1, bound.h0, DF = NULL, bound.exclusions = c("F","B","S","C"))
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  
  NBound.h1 <- NBound.h0 <- NA
  if (asr4)
  {
    if (!(length(bound.exclusions) > 0) & !all(bound.exclusions %in% c("F","B","S","C")))
      warning("A code other than F, B S or C has been specified in bound.exclusions")
    Bound.h1 <- bound.h1  %in% bound.exclusions
    Bound.h0 <- bound.h0  %in% bound.exclusions
  } else #asr3
  {
    #Check bound.exclusions
    if (!all(bound.exclusions %in% c("F","B","S","C")))
    {
      stop("At least one bound.type is not one of those allowed with ASReml-R version 3")
    }
    else
    {
      bound.exclusions3 <- c("Fixed","Boundary","Singular","Constrained")
      bound.exclusions3 <- bound.exclusions3[c("F","B","S","C") %in% bound.exclusions]
      Bound.h1 <- bound.h1  %in% bound.exclusions3
      Bound.h0 <- bound.h0  %in% bound.exclusions3
    }
  }
  NBound.h1 <- sum(Bound.h1)
  Bound.h1 <- names(bound.h1)[Bound.h1] 
  NBound.h0 <- sum(Bound.h0)
  Bound.h0 <- names(bound.h0)[Bound.h0]
  Bound.diff <- c(Bound.h1[!(Bound.h1 %in% Bound.h0)], 
                  Bound.h0[!(Bound.h0 %in% Bound.h1)])
  NBound <- (NBound.h1 + NBound.h0 + length(Bound.diff))/2
  if (NBound > 0)
  {
    mess <- paste("There were a total of", NBound, "bound terms.", sep = " ")
    if (length(Bound.diff) > 0)
    {
      if (is.null(DF))
        warning(paste(mess, "\n  The following bound terms occur in only one of the models",
                      "compared and so were discounted:\n  ",
                      paste(Bound.diff, collapse = ", "),"\n"))
      else
        warning(paste(mess, "\n  The following bound terms occur in only one of the models",
                      "compared, but were not discounted:\n  ",
                      paste(Bound.diff, collapse = ", "),"\n"))
    } else
      if (length(Bound.diff) == 0)
        warning(paste(mess, "These bound terms occur in both models\n"))
  }
  #Calculate degrees of freedom
  DF.calc <- length(bound.h1) - length(bound.h0)
  DF.calc <- DF.calc - NBound.h1 + NBound.h0
  if (DF.calc == 0)
    warning("Calculated DF is zero indicating no difference between models in the number of parameters")
  else
    if (DF.calc <= 0)
      warning("Negative calculated degrees of freedom indicating the second model is more complex")
  if (is.null(DF))
    DF <- DF.calc
  return(list(DF = DF, NBound.h1 = NBound.h1, NBound.h0 = NBound.h0)) 
}

REMLRT.asreml <- function(h0.asreml.obj, h1.asreml.obj, 
                          positive.zero = FALSE, bound.test.parameters = "none", 
                          DF = NULL, bound.exclusions = c("F","B","S","C"), ...)
{
  #asreml codes (vparameters.con code) from ASReml4-SA User Guide Functional Sec. 16.4:
  # (1) P - positive definite
  # (2) ? - liable to change from P to B    
  # (3) U - unbounded
  # (4) F - fixed by user
  # (5) C - Constrained by user (!VCC)      
  # (6) S - Singular Information matrix
  # (7) B - fixed at a boundary (!GP)
  
  #Deal with deprecated constraints parameter
  tempcall <- list(...)
  if (length(tempcall)) 
  {
    if ("full.asreml.obj" %in% names(tempcall))
      stop("full.asreml.obj has been deprecated in REMLRT.asreml - use h1.asreml.obj")
    if ("reduced.asreml.obj" %in% names(tempcall))
      stop("reduced.asreml.obj has been deprecated in REMLRT.asreml - use h0.asreml.obj")
  }
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = FALSE)
  
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
  { for (i in 2:length(fixed.labels)) 
  { if (!all(is.element(fixed.labels[[1]], fixed.labels[[i]]))) 
    stop("Fixed models differ")
  }
  }
  if (!all(diff(sapply(sparse.labels, function(x) length(x))) == 0)) 
    stop("sparse models differ in length")
  if (all(sapply(sparse.labels, function(x) length(x)) > 0)) 
  { for (i in 2:length(sparse.labels)) 
  { if (!all(is.element(sparse.labels[[1]], sparse.labels[[i]]))) 
    stop("sparse models differ")
  }
  }
  
  #Check bound.test.parameters argument
  par.options <- c("none", "onlybound", "one-and-one")
  par.opt <- par.options[check.arg.values(bound.test.parameters, par.options)]
  if (positive.zero & par.opt == "none")
    par.opt <- "onlybound"
  
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
  DF.diff <- DFdiff(bound.h1, bound.h0, DF = DF, bound.exclusions = bound.exclusions)
  NBound.h1 <- DF.diff$NBound.h1
  NBound.h0 <- DF.diff$NBound.h0
  #Perform the test
  REMLRT <- 2*(h1.asreml.obj$loglik-h0.asreml.obj$loglik)
  if (is.null(DF))
    DF <- DF.diff$DF
  
  if (is.na(DF))	
  {
    stop("DF for REML test is NA")
  } else 
  {
    if (DF == 0)
      warning("DF for REML test is zero indicating no difference between models in the number of parameters")
    else
      if (DF <= 0)
        warning("Negative degrees of freedom for REML test indicating the second model is more complex")
  }
  
  if (par.opt != "none")
  { 
    if (REMLRT < 1e-10)
      p <- 1
    else
      if (par.opt == "onlybound")
      {
        if (DF == 1)
          #this option is used when testing a positively-contrained variance component is zero
          #it adjusts for testing on the boundary of the contraint by computing 0.5(ch1_0 + chi_2)
          #comes from Self & Liang (1987) Case 5
          p <- 0.5*(1-pchisq(REMLRT, DF))
        else
        {
          warning(paste("REMLRT on more than one positive variance parameters is being attempted,", 
                        "but is not available"))
          p <- NA
        }
        #following is for DF positively-contrained components equal to zero 
        #computes sum_i=0 to DF mix_i * chi_i where mix_i = choose(DF,  i)*2^(_DF)
        #comes from Self & Liang (1987) Case 9 with s = 0. However, assumes 
        #independence of information matrices
        # {  
        #   df <- seq(DF)
        #   p <- 1-2^(-DF)-sum(choose(DF, df)*2^(-DF)*pchisq(REMLRT, df))
        #   
        # }
      } else #one-and-one Case 6
      {
        if (DF == 2)
          p <- 1 - 0.5*(pchisq(REMLRT, 1) + pchisq(REMLRT, 2))
        else
          stop(paste("For one-and-one, DF must equal two",
                     "because there must be just two parameters tested"), sep = " ")
      }
  }
  else
    p <- 1-pchisq(REMLRT, DF)
  data.frame(REMLRT, DF, p, NBound.h0, NBound.h1)
}

#The code for the full likelihood was adapted from Verbyla (2019) ANZJS, File S1 (doi: 10.1111/anzs.12254) 
infoCriteria.asreml <- function(object, DF = NULL, 
                                bound.exclusions = c("F","B","S","C"), 
                                IClikelihood = "REML", fixedDF = NULL, varDF = NULL, ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  
  #Check that have a valid object of class asreml
  validasr <- validAsreml(object)  
  if (is.character(validasr))
    stop(validasr)
  
  asreml.obj <- object
  call.args <- rlang::call_match(asreml.obj$call, fn = asreml::asreml,  defaults = TRUE)
  #Check which of a Gaussian, GLM or other distribution/model
  if ("asr_gaussian" == rlang::call_name(call.args$family))
  {
    #Check IClikelihood option
    options <- c("REML", "full")
    loglik.opt <- options[check.arg.values(IClikelihood, options)] 
    
    if (loglik.opt == "full" & !asr4)
      stop("The full likelihood option has not been implemented for asreml-R version 3")
    
    #Get bound values
    if (asr4)
    {
      if (asr4.2)
        bound <- asreml.obj$vparameters.con
      else
        bound <- vpc.char(asreml.obj)
    } else
    {
      bound <- names(asreml.obj$gammas.con)
      names(bound) <- names(asreml.obj$gammas)
    }
    NBound <- NA
    if (asr4)
    {
      if (!(length(bound.exclusions) > 0) & !all(bound.exclusions %in% c("F","B","S","C")))
        warning("A code other than F, B S or C has been specified in bound.exclusions")
      Bound <- bound  %in% bound.exclusions
    } else #asr3
    {
      #Check bound.exclusions
      if (!(length(bound.exclusions) > 0) & !all(bound.exclusions %in% c("F","B","S","C")))
      {
        stop("At least one bound.type is not one of those allowed with ASReml-R version 3")
      } else
      {
        bound.exclusions3 <- c("Fixed","Boundary","Singular","Constrained")
        bound.exclusions3 <- bound.exclusions3[bound.exclusions %in% c("F","B","S","C")]
        Bound <- bound  %in% bound.exclusions3
      }
    }
    NBound <- sum(Bound)
    Bound <- names(bound)[Bound]
    #Calculate the varDF
    if (is.null(DF) && is.null(varDF))
    {
      varDF <- length(bound)
      varDF <- varDF - NBound
      if (NBound > 0)
        warning(paste("The following bound terms were discounted:\n", 
                      paste(Bound, collapse = ", ")))
    } else
    {
      if (is.null(varDF))
        varDF <- DF
      if (NBound > 0)
        warning(paste("The following bound terms were not discounted:\n", 
                      paste(Bound, collapse = ", ")))
    }
    #If full likelihood, calculate logdetC and fixedDF
    if (loglik.opt == "full")
    {
      if (asr4)
      {
        asreml::asreml.options(Cfixed = TRUE, gammaPar=FALSE)
        if (is.null(asreml.obj$Cfixed)) 
          asreml.obj <- asreml::update.asreml(asreml.obj, maxit=1)
        coefF <- summary(asreml.obj, coef=TRUE)$coef.fixed
        which.cF <- rownames(coefF)[!is.na(coefF[, "z.ratio"])]
        #      logdetC <- log(prod(svd(as.matrix(asreml.obj$Cfixed[which.cF, which.cF]))$d))
        if (length(asreml.obj$Cfixed) ==  1)
        { 
          logdetC <- log(asreml.obj$Cfixed)
          if (inherits(logdetC, what = "dspMatrix")) logdetC <- logdetC[1,1]
        } else
          logdetC <- sum(log(svd(as.matrix(asreml.obj$Cfixed[which.cF, which.cF]))$d))
      } else #asr3
      {
        if (is.null(asreml.obj$Cfixed)) 
          asreml.obj <- asreml::update.asreml(asreml.obj, maxit=1, Cfixed = TRUE)
        coefF <- summary(asreml.obj, all=TRUE)$coef.fixed
        which.cF <- !is.na(coefF[, "z ratio"])
        #asreml.obj$Cfixed is not a matrix and so this does not work 
        #       logdetC <- log(prod(svd(as.matrix(asreml.obj$Cfixed[which.cF, which.cF]))$d))
        if (any(which.cF)) 
          logdetC <- sum(log(svd(as.matrix(asreml.obj$Cfixed[which.cF, which.cF]))$d))
        else #all z-ratios are NA
        {
          fixedDF <- 0
          logdetC <- 0
          warning("The fixed effects variances are not estimable - reverting to REML likelihood")
        } 
      }
      if (is.null(fixedDF))
      {  
        if (is.logical(which.cF))
          fixedDF <- sum(which.cF)
        else
          fixedDF <- length(which.cF)
        }
    } else #REML
    {
      fixedDF <- 0
      logdetC <- 0
    }
    logREML <- object$loglik
    loglik <- logREML - logdetC/2
    #calculate AIC and BIC
    AIC <- -2 * loglik + 2 * (fixedDF + varDF)
    BIC <- -2 * loglik + (fixedDF + varDF) * log(object$nedf+fixedDF)
    info <- data.frame(fixedDF, varDF, NBound, AIC, BIC, loglik)
  } else #family that is not gaussian 
  {
    info <- data.frame(NA, NA, NA, NA, NA, NA)
    call.fam <- rlang::call_match(call.args$family, 
                                  fn = rlang::call_fn(call.args$family),  
                                  defaults = TRUE)
    distn <- rlang::call_name(call.fam)
    #use deviance for asr_binomial and asr_poison - suggested by Damian Collins
    if (distn %in% c("asr_binomial", "asr_poisson") && call.fam$dispersion == 1)
    {
      if (asr4)
        coefF <- summary(object, coef=TRUE)$coef.fixed
      else #asr3
        coefF <- summary(object, all=TRUE)$coef.fixed
      which.cF <- !is.na(coefF[, "z.ratio"])
      fixedDF <- sum(which.cF)
      #The deviance is, by definition +ve and is -2 * loglik
      if (asr4.2)
      { 
        AIC <- object$deviance + (2 * fixedDF)
        BIC <- object$deviance + fixedDF * log(object$nedf+fixedDF)
        info <- data.frame(fixedDF, 0, 0, AIC, BIC, object$deviance/(-2))
      } else
      { 
        #asreml < 4.2 stores - deviance so that it is 2 * loglik 
        AIC <- - object$deviance + (2 * fixedDF)
        BIC <- - object$deviance + fixedDF * log(object$nedf+fixedDF)
        info <- data.frame(fixedDF, 0, 0, AIC, BIC, (-object$deviance)/(-2))
      }

    } else
      warning("infoCriteria are not available for ", distn, " with dispersion equal to ", call.fam$dispersion, "; they have been set to NA.")
    names(info) <- c("fixedDF", "varDF", "NBound", "AIC", "BIC", "loglik")
  }
  return(info)
}


"infoCriteria.list" <- function(object, bound.exclusions = c("F","B","S","C"), 
                                IClikelihood = "REML", fixedDF = NULL, varDF = NULL, ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that is a list of asreml objects
  if (!all(unlist(lapply(object, function(m) "asreml" %in% class(m)))))
    stop("all elements of the list must be of class asreml")
  
  if (!is.null(fixedDF) & length(fixedDF) != 1)
    stop("Multiple fixedDF has not been implemented")
  if (!is.null(varDF) & length(varDF) != 1)
    stop("Multiple fixedDF has not been implemented")
  
  #Get information criteria
  if (asr4)
    asreml::asreml.options(trace = FALSE)
  ic <- lapply(object, 
               function(m, bound.exclusions, IClikelihood, fixedDF, varDF) 
                 suppressWarnings(suppressMessages(
                   infoCriteria(m, bound.exclusions = bound.exclusions, 
                                IClikelihood = IClikelihood, fixedDF = fixedDF, varDF = varDF, ...))), 
               bound.exclusions = bound.exclusions, 
               IClikelihood = IClikelihood, fixedDF = fixedDF, varDF = varDF)
  ic <- as.data.frame(do.call(rbind, ic))
  if (!is.null(names(object)))
    rownames(ic) <- names(object)
  return(ic)
}

"R2adj.asreml" <- function(asreml.obj, 
                           include.which.fixed = ~ ., orthogonalize = "hybrid", 
                           include.which.random = NULL, 
                           bound.exclusions = c("F","B","S","C"), ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = FALSE)
  if (!asr4)
    stop("R2adj.asreml requires asreml 4.x or later.")
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = FALSE)

  #Check that have a valid object of class asreml
  validasr <- validAsreml(asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  #Check whether Gaussian or other family
  call.args <- rlang::call_match(asreml.obj$call, fn = asreml::asreml,  defaults = TRUE)
  #Check which of a Gaussian, GLM or other distribution/model
  if ("asr_gaussian" != rlang::call_name(call.args$family))
    stop("R2adj.asreml has only been implemented for linear mixed models that assume normality")
  
  #Check that neither include.which.random nor include.which.fixed is NULL
  if (is.null(include.which.random) && is.null(include.which.fixed))
    stop("Cannot have both include.which.random and include.which.fixed set to NULL")
  
  options <- c("eigenmethods", "hybrid")
  orthogonalize <- options[check.arg.values(orthogonalize, options)]

  #Check whether have the design matrix
  if (is.null(asreml.obj$design))
  {
    asreml::asreml.options(design = TRUE)
    asreml.obj <- asreml::update.asreml(asreml.obj)
  }
  
  #Estimate the V for the fitted model
  missing.termmatrix <- NULL
  V <- estimateV(asreml.obj)
  if (all(!is.na(V)) && is.null(attr(V, which = "missing.termmatrix")))
  {  
    n <- nrow(V)
    ASV <- sum(diag((V - (matrix(rep(1, n*n), nrow = n)/n) %*% V))) / (n-1)
    
    #Get the design matrix and effects for the fitted fixed model
    beta <- asreml.obj$coefficients$fixed
    colnames <- rownames(beta)
    X <- as.matrix(asreml.obj$design[, colnames])
    fit <- X %*% beta
    Pfit <- fit - mean(fit)
    Vinv <- dae::mat.ginv(V)
    if (any(is.na(Vinv)))
    {  
      if (!all(is.null(include.which.fixed))) #var.beta not available to compute ASSBfix
      { 
        kresp <- as.character(asreml.obj$call$fixed)[2]
        warning("Unable to obtain the generalized inverse of the estimated variance matrix for ", kresp, 
                " that is required to compute the fixed contribution")
      }
      ASSB <- NA
    } else
    {
      var.beta <- t(X) %*% Vinv %*% X
      var.beta <- dae::mat.ginv(var.beta)
      PX <- X - (matrix(rep(1, n*n), nrow = n) /n) %*% X
      ASSB <- (t(Pfit) %*% Pfit - sum(diag(t(PX) %*% PX %*% var.beta))) / (n-1)
    }
    
    #Get the contribution of the specified fixed terms
    ASSBfix <- 0
    if (!all(is.null(include.which.fixed)) && !is.na(ASSB))
    {
      if (!inherits(include.which.fixed, what = "formula"))
        stop("include.which.fixed must be a formula")
      if (include.which.fixed == ~ .)
        ASSBfix <- ASSB
      else
      {
        #This converts the levels in at terms in both formulae to be double, not single, quoted
        fixterms <- attr(asreml.obj$coefficients$fixed, which = "terms")$tname[-1]
        include.which.fixed <- update.formula(as.formula(paste0('~',
                                                                paste(fixterms, 
                                                                      collapse = ' + '))), 
                                              include.which.fixed)
        #This reinstates single quotes around levels in at function
        incl.fixterms <- getTerms.formula(include.which.fixed)

        if (!all(is.null(incl.fixterms))) #include.which.fixed is not null
          incl.fixterms <- sapply(incl.fixterms, 
                                  function(term, fixterms)
                                  {
                                    term <- fixterms[findterm(term, fixterms)]
                                    if (length(term) == 0) term <- NA
                                    return(term)
                                  }, fixterms = fixterms)
        if (any(is.na(match(incl.fixterms, fixterms))))
        { 
          warning(paste("The following terms are not amongst the fitted fixed terms and will be ignored: ", 
                        paste(incl.fixterms[is.na(match(incl.fixterms, fixterms))], 
                              collapse = ", "), sep = ""))
          incl.fixterms <- incl.fixterms[incl.fixterms %in% fixterms]
        }
        if (length(setdiff(fixterms, incl.fixterms)) == 0)
          ASSBfix <- ASSB
        else
        {
          #use dae:pstructure to get the set of the fixed-term orthogonal projection operators 
          dat <- asreml.obj$call$data
          if (is.symbol(dat))
            dat <- eval(dat)
          Xs <- lapply(fixterms, getTermDesignMatrix, asreml.obj = asreml.obj)
          names(Xs) <- fixterms
          Q.G = projector(matrix(1, nrow=n, ncol=n)/n)
          Qs <- lapply(Xs, function(X, Q.G)
          {
            X <- as.matrix(X)
            X <- cbind(matrix(1, nrow = nrow(X), ncol = 1), X)
            Q <- X %*% dae::mat.ginv(t(X) %*% X) %*% t(X)
            Q <- projector(Q - Q.G)
            return(Q)
          }, Q.G = Q.G)
          #Othogonalize the Qs
          suppressWarnings(
            Qs <- dae::porthogonalize(projectors = Qs, formula = NULL, grandMean = FALSE, 
                                      orthogonalize = orthogonalize, labels = "terms", 
                                      marginality = NULL, check.marginality = FALSE, 
                                      omit.projectors = FALSE, 
                                      which.criteria = "none", 
                                      aliasing.print = FALSE)$Q)
          Q <- matrix(0, nrow = n, ncol = n)
          for (term in incl.fixterms)
            Q <- Q + Qs[[term]]
          QX <- Q %*% X
          Qfit <- Q %*% fit
          ASSBfix <- (t(Qfit) %*% Qfit - sum(diag(t(QX) %*% QX %*% var.beta))) / (n-1)
        }
      }
    }
    
    #Get the contribution of the specified random terms
    ASVran <- 0
    if (!is.null(include.which.random))
    {
      if (!inherits(include.which.random, what = "formula"))
        stop("include.which.random must be a formula")
      if (include.which.random == ~ .)
        G <- estimateV(asreml.obj, which.matrix = "G")
      else
      {
        #This converts the levels in at terms in both formulae to be double, not single, quoted
        ranterms <- attr(asreml.obj$coefficients$random, which = "terms")$tname
        include.which.random <- update.formula(as.formula(paste0('~', 
                                                                 paste(ranterms, 
                                                                       collapse = ' + '))), 
                                               include.which.random)
        #This reinstates single quotes around levels in at function
        incl.ranterms <- getTerms.formula(include.which.random)
        incl.ranterms <- convTerms2Vparnames(incl.ranterms)
        #This removes terms that are not in ranterms and results in terms named as in ranterms
        incl.ranterms <- lapply(incl.ranterms, findterm, termlist = ranterms)
        incl.ranterms <- sapply(incl.ranterms, function(term, ranterms) ranterms[term], 
                                ranterms = ranterms)
        if (!is.allnull(incl.ranterms))
          if (any(is.na(match(incl.ranterms, ranterms))))
            warning(paste("The following terms are not amongst the variance parameters ",
                          "and will be ignored: ", 
                          paste(incl.ranterms[is.na(match(incl.ranterms, ranterms))], 
                                collapse = ", "), sep = ""))
        
        #Work out ranterms to be ignored and get the G matrix
        ignore.terms <- setdiff(ranterms, incl.ranterms)
        if (length(ignore.terms) == 0)
          ignore.terms = NULL
        G <- estimateV(asreml.obj, which.matrix = "G", ignore.terms = ignore.terms)
      }
      
      #Calculate ASVran
      if (all(!is.na(G)) && is.null(attr(G, which = "missing.termmatrix")))
        ASVran <- sum(diag(G -  mean(G))) / (n-1)
      else
        missing.termmatrix <- attr(V, which = "missing.termmatrix")
    }
  } else
    missing.termmatrix <- attr(V, which = "missing.termmatrix")

  #Calculate R2adj
  if (is.null(missing.termmatrix))
  { 
    #Deal with case where ginv for V cannot be obtained and it is not needed for non-null include.which.fixed
    if (all(is.null(include.which.fixed)) && is.na(ASSB))
    {  
      if (is.symbol(dat))
        dat <- eval(dat)
      denom <- var(dat[as.character(asreml.obj$call$fixed)[2]], na.rm = TRUE)
    } else
      denom = ASSB + ASV
    R2.adj <-  as.vector((ASSBfix + ASVran) / denom * 100)
    attr(R2.adj, which = "fixed") <- include.which.fixed
    attr(R2.adj, which = "random") <- include.which.random
  } else
    R2.adj <- NA
  attr(R2.adj, which = "missing.termmatrix") <- missing.termmatrix
  
  return(R2.adj)
}

