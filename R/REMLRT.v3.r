DFdiff <- function(bound.h1, bound.h0, bound.exclusions = c("F","B","S","C"))
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  
  NBound.h1 <- NBound.h0 <- NA
  DF <- length(bound.h1) - length(bound.h0)
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
  DF <- DF - NBound.h1
  Bound.h1 <- names(bound.h1)[Bound.h1] 
  NBound.h0 <- sum(Bound.h0)
  DF <- DF + NBound.h0
  Bound.h0 <- names(bound.h0)[Bound.h0]
  Bound.diff <- c(Bound.h1[!(Bound.h1 %in% Bound.h0)], 
                  Bound.h0[!(Bound.h0 %in% Bound.h1)])
  NBound <- (NBound.h1 + NBound.h0 + length(Bound.diff))/2
  if (NBound > 0)
    warning(paste("There were a total of", NBound, "bound terms.", sep = " "))
  if (length(Bound.diff) > 0)
    warning(paste("The following bound terms occur in only one of the models",
                  "compared and so were discounted:\n",
                  paste(Bound.diff, collapse = ", "),"\n"))
  if (DF == 0)
    warning("DF is zero indicating no difference between models in the number of parameters")
  else
    if (DF <= 0)
      warning("Negative degrees of freedom indicating the second model is more complex")
 return(list(DF = DF, NBound.h1 = NBound.h1, NBound.h0 = NBound.h0)) 
}

REMLRT.asreml <- function(h0.asreml.obj, h1.asreml.obj, 
                          positive.zero = FALSE, bound.test.parameters = "none", 
                          DF = NULL, bound.exclusions = c("F","B","S","C"), ...)
{
#asreml codes:
#  B - fixed at a boundary (!GP)
#  F - fixed by user
#  ? - liable to change from P to B    
#  P - positive definite
#  C - Constrained by user (!VCC)      
#  U - unbounded
#  S - Singular Information matrix

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
    bound.h0 <- asreml::vpc.char(h0.asreml.obj)
    bound.h1 <- asreml::vpc.char(h1.asreml.obj)
  }
  else
  {
    bound.h0 <- names(h0.asreml.obj$gammas.con)
    names(bound.h0) <- names(h0.asreml.obj$gammas)
    bound.h1 <- names(h1.asreml.obj$gammas.con)
    names(bound.h1) <- names(h1.asreml.obj$gammas)
  }
  #Perform the test
	REMLRT <- 2*(h1.asreml.obj$loglik-h0.asreml.obj$loglik)
	if (is.null(DF))
	{
	  DF.diff <- DFdiff(bound.h1, bound.h0, bound.exclusions = bound.exclusions)
	  DF <- DF.diff$DF
	  NBound.h1 <- DF.diff$NBound.h1
	  NBound.h0 <- DF.diff$NBound.h0
	} else
	{
	  if (DF == 0)
	    warning("DF is zero indicating no difference between models in the number of parameters")
	  else
	    if (DF <= 0)
	      warning("Negative degrees of freedom indicating the second model is more complex")
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


infoCriteria.asreml <- function(asreml.obj, DF = NULL, 
                                bound.exclusions = c("F","B","S","C"), ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asreml
  validasr <- validAsreml(asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  #Get bound values
  if (asr4)
  {
    bound <- asreml::vpc.char(asreml.obj)
  }
  else
  {
    bound <- names(asreml.obj$gammas.con)
    names(bound) <- names(asreml.obj$gammas)
  }
  if (is.null(DF))
  {
    NBound <- NA
    DF <- length(bound)
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
      }
      else
      {
        bound.exclusions3 <- c("Fixed","Boundary","Singular","Constrained")
        bound.exclusions3 <- bound.exclusions3[bound.exclusions %in% c("F","B","S","C")]
        Bound <- bound  %in% bound.exclusions3
      }
    }
    NBound <- sum(Bound)
    DF <- DF - NBound
    Bound <- names(bound)[Bound]
    if (NBound > 0)
      warning(paste("The following bound terms were discounted:\n", 
                    paste(Bound, collapse = ", ")))
  }
	logREML <- asreml.obj$loglik
#calculate AIC and BIC
	AIC <- -2 * logREML + 2 * DF
	BIC <- -2 * logREML + DF * log(asreml.obj$nedf)
	data.frame(DF, NBound, AIC, BIC, logREML)
}
