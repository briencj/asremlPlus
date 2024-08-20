#### The functions in this file are documented utility functions for manipulating asremlPlus objects  

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

