addtoChooseSummary <- function(choose.summary, term, DF = NA, denDF = NA, p = NA, 
                               action, omit.DF)
{
  if (omit.DF)
    choose.summary <- addtoTestSummary(test.summary = choose.summary, terms = term, 
                                       p = p, action = action)
  else
    choose.summary <- addtoTestSummary(test.summary = choose.summary, terms = term, 
                                       DF = DF, denDF = denDF, 
                                       p = p, action = action)
  return(choose.summary)
}

#Function to process a df argument
"setdf" <- function(object, df, which.df = "DF")
{
  attr(object, which = which.df) <- which.df
  kattr <- attributes(object)
  
  #Process DF argument
  if (length(df) == 1)
  {
#    if (is.na(df))
#      attr(object, which = "omit.df") <- TRUE
#    else
#    {
      if (is.na(df) |is.numeric(df))
      {
        if (which.df %in% names(object))
          object[which.df] <- rep(df, nrow(object))
        else
        {
          kcols <- ncol(object) 
          object <- cbind(object, rep(df, nrow(object)))
          names(object)[kcols+1] <- which.df
        }
      } else
      {
        if (!is.character(df))
          stop(which.df," must be one of NA, character or numeric")
        if (!(df %in% names(object)))
          stop(df, " is not in ", deparse(substitute(object)))
        attr(object, which = which.df) <- df
      }
#    }
  } else #more than one element
  {
    if (all(is.na(df)) | is.numeric(df))
    {
      if (length(df) != nrow(object))
        stop("Length of ", which.df, "must 1 or equal to the number rows in ", deparse(substitute(object)))
      if (which.df %in% names(object))
        object[which.df] <- df
      else
      {
        kcols <- ncol(object) 
        object <- cbind(object, df)
        names(object)[kcols+1] <- which.df
      }
    } else  
      stop(which.df, " of length > 1 must be numeric")
  }
  newattr <- attributes(object)
  #Find missing atTributes in new alldiffs.obj and add them back in 
  kattr <- kattr[names(kattr)[!(names(kattr) %in% names(newattr))]]
  if (length(kattr) > 0)
  {
    newattr <- c(newattr,kattr)
    attributes(object) <- newattr
  }
  return(object)
}

"chooseModel.data.frame" <- function(object, terms = NULL, p.values = "Pr", 
                                     DF = "Df", denDF = "denDF", omit.DF = FALSE, 
                                     terms.marginality=NULL, alpha = 0.05, ...)
  #function to determine the set of significant terms taking into account marginality relations
  #terms.marginality should be a square matrix of ones and zeroes with row and column names 
  #   being the names of the terms. The diagonal elements should be one, indicating 
  #   that a term is marginal to itself. Elements should be one if the row term is 
  #   marginal to the column term. All other elements should be zero. 
{ 
  #check matrix  
  if (!is.matrix(terms.marginality) || 
      nrow(terms.marginality) != ncol(terms.marginality))
  {
    stop("Must supply a valid marginality matrix")
  } else
    if (is.null(rownames(terms.marginality)) || is.null(colnames(terms.marginality)))
    {
      stop("terms.marginality must have row and column names that are the terms to be tested")
    } else
      if (det(terms.marginality) == 0)
      {
        warning("Suspect marginalities of terms not properly specified - check")
      }
  #make sure have a terms.marginality matrix with lower triangle all zero
  terms.marginality <- permute.to.zero.lowertri(terms.marginality)
  
  # Get and Check terms supplied
  if (is.null(terms))
  {
    object <- cbind(terms = rownames(object), object)
    terms <- "terms"
  }
  if (!all(rownames(terms.marginality %in% object[terms])))
    stop("Supplied data.frame does not contain all terms in terms.marginality")

  #Process DF argument
  object <- setdf(object, df = DF, which.df = "DF")

  #Process denDF argument
  object <- setdf(object, df = denDF, which.df = "denDF")
  DF <- attr(object, which = "DF")
  denDF <- attr(object, which = "denDF")

  #perform tests
  sig.terms <- vector("list", length = 0)
  noterms <- dim(terms.marginality)[1]
  if (omit.DF)
  {
    choose.summary <- makeTestSummary(which.cols = c("terms", "p", "action"))
  } else
  {
    choose.summary <- makeTestSummary(which.cols = c("terms","DF","denDF","p","action"))
  }

  #Determine the significant terms
  if (omit.DF)
  {
    ndf <- NA; den.df <- NA
  } 
  j <- noterms
  #traverse the columns of terms.marginality
  while (j > 0)
  { 
    #get p-value for term for column j
    term <- (rownames(terms.marginality))[j]
    termno <- match(term, object[[terms]])
    p <- object[p.values][[1]][termno]
    if (!omit.DF)
    {
      ndf <- object[[DF]][termno]
      den.df <- object[[denDF]][termno]
    }
    
    #if significant, add to sig term list and work out which can be tested next
    if (!is.na(p)) 
    { 
      if (p <= alpha)
      { 
        action <- "Significant"
        sig.terms <- c(sig.terms, term)
        nonnest <- which(terms.marginality[1:(j-1), j] == 0)
        noterms <- length(nonnest)
        if (noterms == 0)
          j = 0
        else
        { 
          if (noterms == 1)
          { 
            nonnest.name <- rownames(terms.marginality)[nonnest]
            terms.marginality <- terms.marginality[nonnest, nonnest]
            dim(terms.marginality) <- c(1,1)
            rownames(terms.marginality) <- nonnest.name
          }
          else
          { 
            terms.marginality <- terms.marginality[nonnest, nonnest]
          }
          j <- noterms
        }
      }
      #if not significant, proceed to previous column
      else
      {
        action <- "Nonsignificant"
        j <- j - 1
      }
      choose.summary <- addtoChooseSummary(choose.summary, term = term, 
                                           DF = ndf, denDF = den.df, 
                                           p = p, action = action, omit.DF = omit.DF)
    }
    else
      j <- j - 1
  }
  class(choose.summary) <- c("test.summary", "data.frame")
  invisible(list(choose.summary = choose.summary, sig.terms = sig.terms))  
}

"chooseModel.asrtests" <- function(object, terms.marginality=NULL, alpha = 0.05, 
                                   allow.unconverged = TRUE, allow.fixedcorrelation = TRUE, 
                                   checkboundaryonly = FALSE, 
                                   drop.ran.ns=TRUE, positive.zero = FALSE, 
                                   bound.test.parameters = "none", 
                                   drop.fix.ns=FALSE, denDF = "numeric",  dDF.na = "none", 
                                   dDF.values = NULL, trace = FALSE, update = TRUE, 
                                   set.terms = NULL, ignore.suffices = TRUE, 
                                   bounds = "P", initial.values = NA, 
                                   IClikelihood = "none", ...)
  #function to determine the set of significant terms taking into account marginality relations
  #terms.marginality should be a square matrix of ones and zeroes with row and column names 
  #   being the names of the terms. The diagonal elements should be one, indicating 
  #   that a term is marginal to itself. Elements should be one if the row term is 
  #   marginal to the column term. All other elements should be zero. 
{ 
  #Deal with deprecated arguments
  tempcall <- list(...)
  if (length(tempcall)) 
  {
    if ("constraints" %in% names(tempcall))
      stop("constraints has been deprecated in chooseModel.asrtests - use bounds")
    if ("asrtests.obj" %in% names(tempcall))
      stop("asrtests.obj has been deprecated in chooseModel.asrtests - use object")
  }
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(object)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #check matrix  
  if (asr4)
    kresp <- object$asreml.obj$formulae$fixed[[2]]
  else
    kresp <- object$asreml.obj$fixed.formula[[2]]
  if (!is.matrix(terms.marginality) || 
      nrow(terms.marginality) != ncol(terms.marginality))
  {
    stop("In analysing ",kresp,
         ", must supply a valid marginality matrix")
  } else
    if (is.null(rownames(terms.marginality)) || is.null(colnames(terms.marginality)))
    {
      stop("In analysing ",kresp,
           ", terms.marginality must have row and column names that are the terms to be tested")
    } else
      if (det(terms.marginality) == 0)
      {
        warning("In analysing ",kresp,
                ", Suspect marginalities of terms not properly specified - check")
      }
  #make sure have a terms.marginality matrix with lower triangle all zero
  terms.marginality <- permute.to.zero.lowertri(terms.marginality)
  #perform tests
  sig.terms <- vector("list", length = 0)
  noterms <- dim(terms.marginality)[1]
  current.asrt <- object
  j <- noterms
  #traverse the columns of terms.marginality
  while (j > 0)
  { 
    #get p-value for term for column j and, if random, drop if ns and drop.ran.ns=TRUE
    term <- (rownames(terms.marginality))[j]
    current.asrt <- testranfix.asrtests(asrtests.obj = current.asrt, term, 
                                        alpha=alpha, allow.unconverged = allow.unconverged, 
                                        allow.fixedcorrelation = allow.fixedcorrelation, 
                                        checkboundaryonly = checkboundaryonly,  
                                        drop.ran.ns = drop.ran.ns, 
                                        positive.zero = positive.zero, 
                                        bound.test.parameters = bound.test.parameters, 
                                        drop.fix.ns = drop.fix.ns, 
                                        denDF = denDF, dDF.na = dDF.na, 
                                        dDF.values = dDF.values, trace = trace, 
                                        update = update, set.terms = set.terms, 
                                        ignore.suffices = ignore.suffices, 
                                        bounds = bounds, IClikelihood = IClikelihood, 
                                        initial.values = initial.values, ...)
    p <- getTestPvalue(current.asrt, label = term)
    #if significant, add to sig term list and work out which can be tested next
    if (!is.na(p)) 
    { 
      if (p <= alpha)
      { 
        sig.terms <- c(sig.terms, term)
        nonnest <- which(terms.marginality[1:(j-1), j] == 0)
        noterms <- length(nonnest)
        if (noterms == 0)
          j = 0
        else
        { 
          if (noterms == 1)
          { 
            nonnest.name <- rownames(terms.marginality)[nonnest]
            terms.marginality <- terms.marginality[nonnest, nonnest]
            dim(terms.marginality) <- c(1,1)
            rownames(terms.marginality) <- nonnest.name
          }
          else
          { 
            terms.marginality <- terms.marginality[nonnest, nonnest]
          }
          j <- noterms
        }
      }
      #if not significant, proceed to previous column
      else
        j <- j - 1
    }
    else
      j <- j - 1
  }
  invisible(list(asrtests.obj = current.asrt, sig.terms = sig.terms))  
}

"changeModelOnIC.asrtests" <- function(asrtests.obj, 
                                       dropFixed = NULL, addFixed = NULL, 
                                       dropRandom = NULL,  addRandom = NULL, 
                                       newResidual = NULL, 
                                       allow.absentDropTerms = FALSE, label = "Changed terms", 
                                       allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                                       checkboundaryonly = FALSE, 
                                       trace = FALSE, update = TRUE, denDF = "numeric", 
                                       set.terms = NULL, ignore.suffices = TRUE, 
                                       bounds = "P", initial.values = NA, 
                                       which.IC = "AIC", IClikelihood = "REML", 
                                       fixedDF = NULL, varDF = NULL, 
                                       bound.exclusions = c("F","B","S","C"),  
                                       ...)
  #Uses information criteria to select the best model after comparing that fitted in the asrtests.obj 
  #and the one obtained after modifying the model using a combination of adding and removing sets of 
  #terms from one or both of the fixed or random asreml models and replacing the residual model.
  #The function changeTerms is used to change the model.
{ 
  #Check IClikelihood options, because "none" is not allowed here
  options <- c("REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  options <- c("AIC", "BIC") #, "both")
  ic.type <- options[check.arg.values(which.IC, options)]
  
  #Check for fixed correlations in supplied asrtests.obj
  if (!isFixedCorrelOK(asrtests.obj$asreml.obj, allow.fixedcorrelation = allow.fixedcorrelation))
  {
    if ("formulae" %in% names(asrtests.obj$asreml.obj))
      kresp <- asrtests.obj$asreml.obj$formulae$fixed[[2]]
    else
    {
      if ("fixed.formula" %in% names(asrtests.obj$asreml.obj))
        kresp <- asrtests.obj$asreml.obj$fixed.formula[[2]]
      else
        kresp <- NULL
    }
    warning(paste("The estimated value of one or more correlations in the supplied asreml fit for", kresp,
                  "is fixed and allow.fixedcorrelation is FALSE"))
  }
  
  #Calculate the IC for the incoming fit
  old.IC <- infoCriteria(asrtests.obj$asreml.obj, IClikelihood = ic.lik, 
                         bound.exclusions = bound.exclusions, 
                         fixedDF = fixedDF, varDF = varDF, ...)
  old.IC <- as.vector(old.IC[c("fixedDF", "varDF", "AIC", "BIC")], mode = "numeric")
  names(old.IC) <- c("DF", "denDF", "AIC", "BIC")
  nlines.test <- nrow(asrtests.obj$test.summary)
  
  #Check that terms to drop are present - otherwise do not change the model, unless 
  #     allow.absentDropTermsis TRUE
  absent <- FALSE
  if (!allow.absentDropTerms&& (!is.null(dropFixed) || !is.null(dropRandom)))
  {
    if (!is.null(dropFixed))
    {
      #Test whether any terms in dropFixed are absent
      fix.form <- as.formula(paste("~ (", dropFixed, ")", sep = ""))
      fix.terms <- getTerms.formula(fix.form)
      fixterms.obj <- as.terms.object(languageEl(asrtests.obj$asreml.obj$call, which="fixed"), 
                                      asrtests.obj$asreml.obj)
      if (any(unlist(lapply(fix.terms, 
                            function (term, terms.obj) findterm(term, labels(terms.obj)) == 0, 
                            terms.obj = fixterms.obj))))
        absent <- TRUE
    }
    if (!is.null(dropRandom))
    {
      #Test whether any terms in dropRandom are absent
      ran.form <- as.formula(paste("~ (", dropRandom, ")", sep = ""))
      ran.terms <- getTerms.formula(ran.form)
      ranterms.obj <- as.terms.object(languageEl(asrtests.obj$asreml.obj$call, which="random"), 
                                      asrtests.obj$asreml.obj)
      if (any(unlist(lapply(ran.terms, 
                            function (term, terms.obj) findterm(term, labels(terms.obj)) == 0, 
                            terms.obj = ranterms.obj))))
        absent <- TRUE
    }
  }
 
  if (absent)
  {
    ic.NA <- data.frame(fixedDF = NA, varDF = NA, AIC = NA, BIC = NA)
    action <- "Absent"
    ic <- ic.NA
    test.summary <- asrtests.obj$test.summary
    test.summary <- addtoTestSummary(test.summary, terms = label, 
                                     DF=NA, denDF = NA, p = NA, AIC = NA, BIC = NA, 
                                     action = "Absent")
  } else
  {
    #Use changeTerms to change the model
    new.asrtests.obj <- changeTerms(asrtests.obj, 
                                    dropFixed = dropFixed, addFixed = addFixed, 
                                    dropRandom = dropRandom,  addRandom = addRandom, 
                                    newResidual = newResidual, label = label, 
                                    allow.unconverged = allow.unconverged, 
                                    allow.fixedcorrelation = allow.fixedcorrelation, 
                                    checkboundaryonly = checkboundaryonly, 
                                    trace = trace, update = update, denDF = denDF, 
                                    set.terms = set.terms, ignore.suffices = ignore.suffices, 
                                    bounds = bounds, initial.values = initial.values, 
                                    IClikelihood = IClikelihood, 
                                    bound.exclusions = bound.exclusions,  
                                    ...)
    #Obtain IC for new model
    new.IC <- infoCriteria(new.asrtests.obj$asreml.obj, IClikelihood = ic.lik, 
                           bound.exclusions = bound.exclusions, 
                           fixedDF = fixedDF, varDF = varDF, ...)
    new.IC <- as.vector(new.IC[c("fixedDF", "varDF", "AIC", "BIC")], mode = "numeric")
    names(new.IC) <- c("DF", "denDF", "AIC", "BIC")
    
    #Extract asreml.objects
    asreml.obj <- asrtests.obj$asreml.obj
    if (!is.null(asreml.obj$mf) && !is.null(attr(asreml.obj$mf, which = "mbf.env")))
      attr(new.asrtests.obj$asreml.obj$mf, 
           which = "mbf.env") <- attr(asreml.obj$mf, which = "mbf.env")
    new.asreml.obj <- new.asrtests.obj$asreml.obj

    change <- FALSE
    action <- getTestEntry(new.asrtests.obj, label = label)$action 
    diff.IC <- new.IC - old.IC
    #check convergence
    # if (!allow.unconverged && (!asreml.obj$converge | !new.asreml.obj$converge))
    # {
    #   if (!asreml.obj$converge)
    #   {
    #     if (!new.asreml.obj$converge)
    #     {
    #       action <- "Unchanged - both unconverged"
    #       change <- FALSE
    #     } else
    #     {
    #       action <- "Swappped - old unconverged"
    #       change <- TRUE
    #     }
    #   } else
    #   {
    #     action <- "Unchanged - new unconverged"
    #     change <- FALSE
    #   }
    # } else
    # {
    #Check fixed correlation
      if (grepl("Unchanged - fixed correlation", action))
      {
        #new.asrtests.obj <- asrtests.obj
        change <- FALSE
      } else
      {
        if (grepl("Unchanged", action))
        {
          #new.asrtests.obj <- asrtests.obj
          change <- FALSE
        } else
        {
          if ((ic.type == "AIC" & diff.IC["AIC"] < 0) || 
              (ic.type == "BIC" & diff.IC["BIC"] < 0))
          {
            change <- TRUE
            action <- "Swapped"
          } else
            action <- "Unswapped"
          
          #check convergence, when it is allowed
          if (allow.unconverged)
          {
            if (!asreml.obj$converge && !new.asreml.obj$converge)
            {
              action <- paste(action, " - both unconverged", sep="")
            } else
            {
              if (!asreml.obj$converge)
                action <- paste(action, " - old unconverged", sep="")
              else
              {
                if (!new.asreml.obj$converge)
                  action <- paste(action, " - new unconverged", sep="")
              }
            }
          }
        }
      }
 #   }

    #It is not checked that removing a boundary term will result in convergence here as is done
    #in, for example, testresidual.asrtests
    if (change)
      asrtests.obj <- new.asrtests.obj
    
    #Modify action to reflect model selection, except for fixed correlation
    test.summary <- new.asrtests.obj$test.summary
    krows <- (nlines.test+1):nrow(test.summary)
    k <- krows[test.summary$terms[krows] == label & 
                 test.summary$action[krows] != "Boundary"]
    
    test.summary[k, "DF"] <- diff.IC["DF"]
    test.summary[k, "denDF"] <- diff.IC["denDF"]
    test.summary[k, "AIC"] <- diff.IC["AIC"]
    test.summary[k, "BIC"] <- diff.IC["BIC"]
    test.summary[k, "action"] <- action
  }
  asrtests.obj$test.summary <- test.summary

  return(asrtests.obj)
}
