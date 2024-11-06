#### The functions in this file are documented utility functions for manipulating asremlPlus objects  

"as.asrtests" <- function(asreml.obj, wald.tab = NULL, test.summary = NULL, 
                          denDF = "numeric", label = NULL, 
                          IClikelihood = "none", bound.exclusions = c("F","B","S","C"), 
                          ...)
{ 
  
  # Get original data name.
  orig.name <- asreml.obj$call$data                                       ## add VSNi 14/03/2024
  
  # Expose dataset to avoid conflicts.
  asreml.obj$call$data <- str2lang("asreml.obj$mf")                       ## add VSNi 14/03/2024
  
  #Check that have a valid object of class asreml
  validasr <- validAsreml(asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  if (asr4)
    asreml::asreml.options(trace = FALSE)
  
  #Check IClikelihood options
  options <- c("none", "REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  ic.NA <- data.frame(fixedDF = NA, varDF = NA, AIC = NA, BIC = NA)
  
  #Process test.summary and label arguments
  if (is.null(test.summary))
    test.summary <- makeTestSummary(which.cols = c("terms","DF","denDF",
                                                   "p","AIC","BIC","action"))
  else
  {
    validtests <- validTestSummary(test.summary)
    if (is.character(validtests))
      stop(validtests)
  }
  if (!is.null(label))
  {
    if (ic.lik != "none")
      ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                         bound.exclusions = bound.exclusions)
    else
    {
      ic <- ic.NA
      ic$varDF <- infoCriteria(asreml.obj, IClikelihood = "REML", 
                               bound.exclusions = bound.exclusions)$varDF
    }
    test.summary <- addtoTestSummary(test.summary, terms = label, 
                                     DF=ic$fixedDF, denDF = ic$varDF, p = NA, 
                                     AIC = ic$AIC, BIC = ic$BIC, 
                                     action = "Starting model")
  }
  
  #Deal with wald.tab
  if (!is.null(wald.tab))
  {  
    #Check that have a valid wald.tab object
    validwald <- validWaldTab(wald.tab)  
    if (is.character(validwald))
      stop(validwald)
  } else #form wald.tab
  { 
    wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = FALSE, ...)
    wald.tab <- chkWald(wald.tab)
  }
  
  # Add original name back to asreml.obj.
  asreml.obj$call$data <- orig.name                                     ## add VSNi 14/03/2024
  
  #Put together the asrtests.obj, with updated wald tab
  test <- list(asreml.obj = asreml.obj, wald.tab=wald.tab, test.summary = test.summary)
  class(test) <- c("asrtests", "list")
  test$wald.tab <- recalcWaldTab(test, ...)
  
  #Reset trace to default
  if (asr4)
    asreml::asreml.options(trace = TRUE)
  
  
  # Return.
  return(test)
}

"is.asrtests" <- function(object)
{
  inherits(object, "asrtests")
}

"validAsrtests" <- function(object)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  isasrtests <- TRUE 
  
  #Check that have a valid object of class asrtests
  if (is.null(object) || !is.asrtests(object))
  {
    if (asr4)
      msg <- paste0("\n  In analysing ",object$asreml.obj$formulae$fixed[[2]],
                    ", must supply an object of class asrtests", sep = "")
    else
      msg <- paste0("\n  In analysing ",object$asreml.obj$fixed.formula[[2]],
                    ", must supply an object of class asrtests", sep = "")
    isasrtests[1] <- FALSE
    isasrtests <- c(isasrtests, msg)
  }    
  
  #Check have appropriate components
  if (!all(c("asreml.obj", "wald.tab", "test.summary") %in% names(object)))
  {
    isasrtests[1] <- FALSE
    isasrtests <- c(isasrtests, 
                    paste("\n ", deparse(substitute(object)), 
                          "is not a list with named components",
                          "asreml.obj, wald.tab and test.summary"))
  }    
  if (length(isasrtests) > 1)
    isasrtests[1] <- "Error in validAsrtests : "
  return(isasrtests)
}

setOldClass("asrtests")

"validAsreml" <- function(object)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  isasr <- TRUE 
  #Check class
  if (!inherits(object, "asreml") || is.null(object))
  {
    isasr[1] <- FALSE 
    isasr <- c(isasr, "\n ",deparse(substitute(object))," is not of class 'asreml'")
  }
  #Check have corresponding asreml version
  if (asr4)
  {
    if (!("vparameters" %in% names(object)))
    { 
      isasr[1] <- FALSE 
      isasr <- c(isasr, 
                 paste0("\n ",deparse(substitute(object)), 
                        " is not an asreml object compatible ASReml-R version 4; ",
                        "\nuse convASRemlobjVersion.asreml to create a compatible asreml object"))
    } else
      if (asr4.2 && ("errtxt" %in% names(object)))
      {
        isasr[1] <- FALSE 
        isasr <- c(isasr, 
                   paste0("\n ",deparse(substitute(object)), 
                          " is not an asreml object compatible with ASReml-R version 4.2; ",
                          "\nuse convASRemlobjVersion.asreml to create a compatible asreml object"))
      }
  } else #not asr4
    if ("vparameters" %in% names(object))
    {
      isasr[1] <- FALSE 
      isasr <- c(isasr, 
                 paste0("\n ",deparse(substitute(object)), " is not compatible with ",
                        "ASReml-R ",packageVersion("asreml"), ", the currently loaded version; ",
                        "\nuse convASRemlobjVersion.asreml to create a compatible asreml object"))
    }
  if (length(isasr) > 1)
    isasr[1] <- "Error(s) in validAsreml(object) : "
  return(isasr)
}

"convAsremlobj.asreml" <- function(asreml.obj, ...)
{
  y <- runif(10); x <- runif(10); a <- asreml::asreml(y~x)
  call <- asreml.obj$call
  call[[1]] <- a$call[[1]]  
  asreml.obj <- eval(call)
  validAsreml(asreml.obj)
  return(asreml.obj)
}

"validWaldTab" <- function(object)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = FALSE)
  iswald <- TRUE 
  if (!is.null(object) && (!is.data.frame(object) || all(ncol(object) != c(4,6))))
  {
    iswald[1] <- FALSE
    iswald <- c(iswald, 
                "wald.tab should be a 4- or 6-column data.frame -- perhaps extract Wald component from list")
  }
  if (length(iswald) > 1)
    iswald[1] <- "Error(s) in validWaldTab(object) : "
  return(iswald)
}

"validTestSummary" <- function(object)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = FALSE)
  isTSumm <- TRUE 
  if (!is.null(object) && (!is.data.frame(object) || all(ncol(object) != c(3,5,7))))
  {
    isTSumm[1] <- FALSE
    isTSumm <- c(isTSumm, 
                 "test.summary should be a 3-, 5- or 7-column data.frame")
  }
  which.names <- names(object) %in% c("terms","DF","denDF","p", "AIC", "BIC","action")
  if(!all(which.names))
  {
    isTSumm <- c(isTSumm, 
                 paste0("test.summary contains the illegal column(s) ", 
                        paste(names(object)[!which.names], collapse = ",")))
  }
  if (length(isTSumm) > 1)
    isTSumm[1] <- "Error(s) in validWaldTab(object) : "
  return(isTSumm)
}

"makeTestSummary" <- function(which.cols = c("terms","DF","denDF","p","AIC","BIC","action"))
{
  test.summary <- as.data.frame(matrix(nrow = 0, ncol = length(which.cols)))
  names(test.summary) <- which.cols
  class(test.summary) <- c("test.summary", "data.frame")
  return(test.summary)
}

"addto.test.summary" <- function(test.summary, terms, DF = 1, denDF = NA, 
                                 p = NA, AIC = NA, BIC = NA,
                                 action = "Boundary")
{
  test.summary <- addtoTestSummary(test.summary = test.summary, terms = terms, 
                                   DF = DF, denDF = denDF, p = p, AIC = AIC, BIC = BIC,
                                   action = action)
  return(test.summary)    
}

"addtoTestSummary" <- function(test.summary, terms, DF = 1, denDF = NA, 
                               p = NA, AIC = NA, BIC = NA,
                               action = "Boundary")
{
  which.cols <- names(test.summary)
  new.row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(test.summary)))
  names(new.row) <- which.cols
  if ("terms" %in% which.cols) new.row$terms <- terms
  if ("DF" %in% which.cols) new.row$DF <- DF
  if ("denDF" %in% which.cols) new.row$denDF <- denDF
  if ("p" %in% which.cols) new.row$p <- p
  if ("AIC" %in% which.cols) new.row$AIC <- AIC
  if ("BIC" %in% which.cols) new.row$BIC <- BIC
  if ("action" %in% which.cols) new.row$action <- action
  
  test.summary <- rbind(test.summary, new.row, stringsAsFactors = FALSE)
  class(test.summary) <- c("test.summary", "data.frame")
  return(test.summary)
}

#Function to revert to a previous model fit, transferring the test.summary to the last 
#asrtests object, with an extra line added to the test summary using terms and the action 
#arguments.
#
#Additional arguments to addtoTestSummary can be supplied to revert2previousFit using the ellipsis.
#If asrtests.prev = NULL, then only a line is added to asrtests.last (same as directly using 
#addtoTestSummary)
revert2previousFit <- function(asrtests.last, asrtests.prev = NULL, terms, action, DF = NA, ...)
{
  test.summary <- addtoTestSummary(asrtests.last$test.summary, terms = terms, 
                                   DF = DF, action = action)
  if (!is.null(asrtests.prev)) asrtests.last <- asrtests.prev
  asrtests.last$test.summary <- test.summary
  return(asrtests.last)
}

"chkWald" <- function(wald.tab)
{
  if (!is.null(wald.tab)) 
  {
    if (is.list(wald.tab))
      wald.tab <- wald.tab$Wald
    else
      if (is.matrix(wald.tab))
      {
        hd <- attr(wald.tab, which = "heading")
        wald.tab <- as.data.frame(wald.tab, stringsAsFactors = FALSE)
        nofixed <- dim(wald.tab)[1]
        wald.tab <- wald.tab[1:(nofixed-1), c(1,3,4)] #Remove Residual line
        wald.tab$F.inc <- wald.tab$'Wald statistic'/wald.tab$Df
        hd[1] <- "Conservative Wald F tests for fixed effects \n"
        attr(wald.tab, which = "heading") <- hd
      }
    
    #Check that have a valid wald.tab object
    validwald <- validWaldTab(wald.tab)  
    if (is.character(validwald))
      stop(validwald)
  }
  
  if (!is.null(wald.tab))
    class(wald.tab) <- c("wald.tab", "data.frame")
  
  return(wald.tab)  
}

"print.test.summary" <- function(x,  which.print = c("title", "table"), 
                                 omit.columns = NULL, response = NULL, ...)
{
  options <- c("title", "table", "all")
  opt <- options[unlist(lapply(which.print, check.arg.values, options=options))]
  
  #make change to control printing
  class(x) <- c("test.summary", "data.frame")
  x$p <- round(x$p, digits=4)
  
  #remove unwanted columns
  if (!is.null(omit.columns))
  {
    which.cols <- names(x)
    which.cols <- which.cols[!(which.cols %in% omit.columns)]
    x <- x[which.cols]
  }
  
  if (any(c("title", "all") %in% opt))
  {
    if (is.null(response)) 
      cat("\n\n####  Sequence of model investigations \n\n")
    else
      cat("\n\n####  Sequence of model investigations for", response, "\n\n")
    if (any(c("AIC", "BIC") %in% names(x)))
      cat("(If a row has NA for p but not denDF, DF and denDF relate to fixed and variance parameter numbers)\n\n")
  }
  
  print.data.frame(x, ...)
  invisible()
}

"print.wald.tab" <- function(x, which.wald = c("title", "heading", "table"), 
                             colourise = FALSE, ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = FALSE)
  
  options <- c("title", "heading", "table", "all")
  opt <- options[unlist(lapply(which.wald, check.arg.values, options=options))]
  
  #make change to control printing
  class(x) <- c("wald", "data.frame")
  x$Pr <- round(x$Pr, digits=4)
  
  if (any(c("title", "all") %in% opt))
    cat("\n\n####  Pseudo-anova table for fixed terms \n\n")
  if  (!any(c("table", "all") %in% opt)) #no table to be printed
  {
    if ("heading" %in% opt)
    {
      hd <- attr(x, which = "heading")
      for (i in 1:length(hd))
        cat(hd[i],"\n")
    }
  }
  else #print table, possibly with heading
  {
    if (any(c("heading", "all") %in% opt) && !is.null(asr4) && asr4)
    {
      if (asr4)
      {
        asr.col <- asreml::asreml.options()$colourise
        if (length(asr.col) == 0)
        {
          if (colourise) 
            asreml::asreml.options(colourise = colourise)
        } else 
          if (xor(colourise,asr.col))
            asreml::asreml.options(colourise = colourise)
        print(x, ...)
        asreml::asreml.options(colourise = asr.col)
      } else
        print(x, ...)
    } else
    {
      if (any(c("heading", "all") %in% opt))
      {
        hd <- attr(x, which = "heading")
        for (i in 1:length(hd))
          cat(hd[i],"\n")
      } else
      {
        cat("\n")      
      }
      print.data.frame(x, ...)
    }
  }
  
  invisible()
}

"print.asrtests" <- function(x, which = "key", colourise = FALSE, ...)
{ 
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(x)  
  if (is.character(validasrt))
    stop(validasrt)
  
  options <- c("asremlsummary", "vparametersummary", "pseudoanova", "wald.tab", 
               "testsummary", "key", "all")
  opt <- options[unlist(lapply(which, check.arg.values, options=options))]
  if (all(c("key", "all") %in% opt))
    stop("Can only specify one of key and all for which argument")
  if ("wald.tab" %in% opt)
  {
    opt[match("wald.tab", opt)] <- "pseudoanova"
    opt <- unique(opt)
  }
  
  kresp <- x$asreml.obj$formulae$fixed[[2]]
  
  #print summary of asreml.obj
  if (any(c("asremlsummary", "all") %in% opt))
    print(summary(x$asreml.obj), ...)
  
  #print vparameter summary of asreml.obj
  if (any(c("vparametersummary", "key") %in% opt))
  {
    cat("\n\n####  Summary of the fitted variance parameters for", kresp, "\n\n")
    print(summary(x$asreml.obj)$varcomp, ...)
  }
  
  #print wald.tab
  if (any(c("pseudoanova", "key", "all") %in% opt))
    print.wald.tab(x$wald.tab, colourise = colourise, ...)
  
  #print test.summary
  if (any(c("testsummary", "key", "all") %in% opt))
    print.test.summary(x$test.summary, which.print = "all", 
                       response = kresp, ...)
  
  invisible()
}

#This function identifies the names of the rows or columns that correspond to the effects for term.
#The argument `use` identifies which component of an asreml.obj is to be used to obtain the effect names.
#Possible values are: fixed.coeffs, random.coeffs, G.aom or design. 
"getTermEffectNames" <- function(term, asreml.obj, use = "design.matrix", sep = ":")
{
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  if (asr4.2)
  { 
    vars <- fac.getinTerm(term, asr4.2 = asr4.2)
    # vars <- gsub("\\(", "\\\\(", vars)
    # vars <- gsub("\\)", "\\\\)", vars)
    nvars <- length(vars)
    
    #Get the effects names from the object in use
    if (use == "fixed.coeffs")
      effnames <- rownames(asreml.obj$coefficients$fixed)
    else
    {
      if (use == "random.coeffs")
        effnames <- rownames(asreml.obj$coefficients$random)
      else
      {
        if (use == "G.aom")
          effnames <- rownames(asreml.obj$aom$G)
        else if (use == "design.matrix")
        {
          effnames <- colnames(asreml.obj$design)
        } else
          stop("Cannot use ", use)
      }
    }
    if (is.null(effnames))
      stop("No effects names found")
    col.vars <- strsplit(effnames, split = sep, fixed = TRUE)
    #reduce to columns of same length as term
    len.nvars <- unlist(lapply(col.vars, length)) == nvars
    effnames <- effnames[len.nvars]
    if (length(effnames) > 0)
    {
      col.vars <- col.vars[len.nvars]
      col.vars <- as.data.frame(do.call(rbind, col.vars))
      #Determine which columns have the same vars, in any order, as the term; select effnames that do
      which.cols <- apply(col.vars, MARGIN = 1, 
                          FUN = function(krow, vars)
                          {
                            all(sapply(vars, function(v, krow) 
                            {
                              if (any(startsWith(krow, paste0(v,"_")))) #term with a level
                                TRUE
                              else #must match exactly if does not include a level
                                any(krow == v)
                            }, krow = krow))
                          }, vars = vars)
      effnames <- effnames[which.cols]
    }
  } else
  { 
    if (use == "fixed.coeffs")
      effnames <- rownames(asreml.obj$coefficients$fixed)[
        startsWith(rownames(asreml.obj$coefficients$fixed), term)]
    else
    {
      if (use == "random.coeffs")
        effnames <- rownames(asreml.obj$coefficients$random)[
          startsWith(rownames(asreml.obj$coefficients$random), term)]
      else
      {
        if (use == "G.aom")
          effnames <- rownames(asreml.obj$aom$G)[startsWith(rownames(asreml.obj$aom$G), term)]
        else if (use == "design")
        {
          effnames <- colnames(asreml.obj$design)[startsWith(colnames(asreml.obj$design), term)]
        } else
          stop("Cannot use", use)
      }
      restnams <- substring(effnames, first = nchar(term)+1)
      effnames <- effnames[grepl("^[0-9]*$", restnams)]
    }
  }
  if (length(effnames) == 0)
    effnames <- NULL
  # else
  # {
  #   #grp terms do not have a known number of columns
  #   if (term %in% names(asreml.obj$noeff) && length(effnames) != asreml.obj$noeff[term])
  #     stop(paste("Error in finding the columns in ", use,  " for ", term, sep=""))
  # }
  return(effnames)
}

#This function produces a data frame of the factors involved a set of effnames
#Was as.data.frame.effnames 
#The argument `use` identifies which component of an asreml.obj is to be used to obtain the effect names.
#Possible values are: fixed.coeffs, random.coeffs, G.aom or design. 
convEffectNames2DataFrame.asreml <- function(asreml.obj, term, use = "design.matrix", sep = ":", ...)
{
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  if (!asr4.2)
    stop("This funtion requires asreml v4.2 to run")
  
  #Check that have a valid object of class asreml
  validasr <- validAsreml(asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  #Check use options
  options <- c("fixed.coeffs", "random.coeffs", "G.aom", "design.matrix")
  use <- options[check.arg.values(use, options)]
  
  #get effect names for term and determine the factors in the term
  effnames <- getTermEffectNames(term = term, asreml.obj = asreml.obj, use = use, sep = sep)
  if (is.null(effnames))
    facs.vars <- NULL
  else
  {
    col.vars <- strsplit(effnames, split = sep, fixed = TRUE)
    col.vars <- as.data.frame(do.call(rbind, col.vars))
    
    #Generate the data.frame
    facs.vars <- as.data.frame(lapply(col.vars, function(var)
    {
      if (startsWith(var[1], "at(")) #Deal with at, assuming only one level
      {
        if (stringr::str_count(var[1], "at\\(") > 1 || grepl("\\+", var[1]))
          stop("term is too complex for convEffectNames2DataFrame.asreml to process")
        var <- stringr::str_sub(var, 4, nchar(var)-1)
        var <- gsub("\\'", "", var)
        var <- strsplit(var, ", ")
      } else
        var <- strsplit(var, split = "\\_")
      var.nam <- var[[1]][1]
      if (!all(sapply(var, function(effnames, var.nam) effnames[1] == var.nam, var.nam = var.nam)))
        stop("The effects names are not based on the same variables")
      var.levs <- as.data.frame(do.call(rbind, lapply(var, function(v) v[2])))
      names(var.levs) <- rmFunction(var.nam)
      var.levs[,1] <- factor(var.levs[,1], levels = unique(var.levs[,1]))
      return(var.levs)
    }))
    
    #Add the effect names as an attribute
    attr(facs.vars, which = "effect.names") <- effnames
  }
  return(facs.vars)
}

"getTermDesignMatrix" <- function(term, asreml.obj, sep = ":")
{
  colnams <- getTermEffectNames(term, asreml.obj, use = "design.matrix", sep = sep)
  terms <- rbind(attr(asreml.obj$coefficients$fixed, which = "terms"),
                 attr(asreml.obj$coefficients$random, which = "terms"))
  #grp terms do not have a known number of columns
  #Does the no. cols in the  design matrix equal the number of effects in the coefficents? 
  if (term %in% names(asreml.obj$noeff) && !is.null(colnams) && 
      length(colnams) != terms$n[terms$tname == term])
  {
    #A fixed term?
    if (term %in% attr(asreml.obj$coefficients$fixed, which = "terms"))
      coeffnams <- getTermEffectNames(term, asreml.obj, use = "fixed.coeffs", sep = sep)
    else #Random term
    {  
      if (term %in% attr(asreml.obj$coefficients$random, which = "terms")$tname)
        coeffnams <- getTermEffectNames(term, asreml.obj, use = "random.coeffs", sep = sep)
      else #Can't find term
        stop("Could not find term on fixed and random coeficients")
    }
    
    #Have both sets of names, so add zero cols to the design matrix to make them the same
    if (all(colnams %in% coeffnams))
    { 
      Z <- matrix(0, nrow = nrow(asreml.obj$design), ncol = terms$n[terms$tname == term])
      colnames(Z) <- coeffnams
      for (nam in colnams)
        Z[,nam] <- asreml.obj$design[ ,nam]
    } else
      stop(paste0("For ", term, ", cannot match the ", length(colnams), 
                  " columns in the design matrix, with the ", terms$n[terms$tname == term], 
                  " elements in the coefficients component"))
  } else
  { 
    if (is.null(colnams))
      Z <- NULL
    else
      Z <- asreml.obj$design[, colnams]
  } 
  return(Z)
}

