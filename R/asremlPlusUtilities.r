"check.arg.values" <- function(arg.val, options)
  #Function to check that arg.val is one of the allowed values
  #and to return the position of the argument in the set of values
  #that is stored in options
{ 
  kopt <- pmatch(arg.val, options)
  if (any(is.na(kopt)))
    stop("Value ",paste(arg.val, collapse = ","), " is either not unique or is not an allowed option for its argument")
  if (length(kopt) > 1)
  {
    warning(paste("Only one value allowed for argument where", 
                  paste(arg.val, collapse = ","), "have been supplied", 
                  sep = " "))
    kopt <- kopt[1]
  }
  return(kopt)
}

#Get the loaded version of asreml
"getASRemlVersionLoaded" <- function(nchar = NULL, notloaded.fault = FALSE)
{
  if (isNamespaceLoaded("asreml"))
  {
    sInf <- sessionInfo(package = "asreml")
    vers <- sInf$otherPkgs$asreml$Version
    if (!is.null(nchar))
      vers <- substring(vers, first = 1, last = nchar)
  } else
  {
    if (notloaded.fault)
      stop("No version of asreml loaded")
    else
      vers <- NULL
  }
  return(vers)
}

#Checks whether the loaded version is a subset of version, excluding incompatible versions 
"isASRemlVersionLoaded" <- function(version = 4, notloaded.fault = FALSE)
{
  vers <- getASRemlVersionLoaded(nchar = NULL, notloaded.fault = notloaded.fault)
  if (!is.null(vers))
  {
    if (substring(vers, first = 1, last = 3) == "4.0")
      stop("This version of asremlPlus does not support asreml 4.0.x - asremlPlus 4.0.x does")
    version <- as.character(version)
    vers <- substring(vers, first = 1, last = nchar(version)) == version
  }
  return(vers)
}

#Checks if loaded version of asreml begins with version; if not unloads asreml and loads the version in lib.loc 
"loadASRemlVersion" <- function(version = "4", ...)
{
  loadASReml <- FALSE
  vers <- isASRemlVersionLoaded(version)
  if (is.null(vers))
    loadASReml <- TRUE
  else
  {
    if (!vers)
    {
      unloadNamespace("asreml")
      loadASReml <- TRUE
    }
  }
  if (loadASReml)
  {
    gotASReml <- requireNamespace("asreml", ...)
    if (gotASReml)
      attachNamespace("asreml")
  }
  if (!isASRemlVersionLoaded(version))
    stop("Unable to load asreml, version ", version)
  vers <- getASRemlVersionLoaded()
  invisible(vers)
}

"getFormulae.asreml" <- function(asreml.obj, which = c("fixed", "random", "residual"), 
                                 expanded = FALSE, envir = parent.frame(), ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)

  #Process which argument
  which.options <- c("fixed", "random", "residual", "sparse", "all")
  forms.opt <- which.options[unlist(lapply(which, check.arg.values, options=which.options))]
  if ("all" %in% forms.opt)
    forms.opt <- c("fixed", "random", "residual", "sparse")
  if (!asr4)
    forms.opt <- gsub("residual", "rcov", forms.opt, fixed = TRUE)
  
  #Get formula(e)
  modlist <- lapply(forms.opt, 
                    function(form, asreml.obj) 
                    {
                      modl <- asreml.obj$call[[form]]
                      #evaluate mod in the frame of the function call that led to here 
                      modl <- eval(modl, envir = envir) 
                      return(modl)
                    }, 
                    asreml.obj = asreml.obj)
  names(modlist) <- forms.opt
  
  #Expand if required
  if (expanded)
    modlist <- lapply(modlist, 
                      function(x) 
                      {  
                        if (!is.null(x)) x <- update.formula(x, ~., ...) 
                        else x <- x
                      })
  return(modlist)
}

"printFormulae.asreml" <- function(asreml.obj, which = c("fixed", "random", "residual"), 
                                   expanded = FALSE, envir = parent.frame(), ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  mod <- getFormulae.asreml(asreml.obj, which = which, expanded = expanded, envir = envir,  ...)
  mod.ch <- lapply(mod, 
                   function(m) 
                   {
                     attributes(m) <- NULL
                     m <- capture.output(m)
#                     m <- m[-length(m)]
                     if (length(m) > 1)
                     {
                       m <- unlist(lapply(m, function(m) m <- stringr::str_trim(m, side = "left")))
                       m <- paste0(m, collapse = "")
                     }
                     return(m)
                   })
  if ("random" %in% names(mod.ch))
    mod.ch$random <- gsub("~", "~ ", mod.ch$random)
  if (asr4)
  {
    if ("residual" %in% names(mod.ch))
      mod.ch$residual <- gsub("~", "~ ", mod.ch$residual)
  } else
  {
    if ("rcov" %in% names(mod.ch))
      mod.ch$rcov <- gsub("~", "~ ", mod.ch$rcov)
  }
  if ("sparse" %in% names(mod.ch))
    mod.ch$sparse <- gsub("~", "~ ", mod.ch$sparse)
  m.ch <- unlist(lapply(names(mod.ch), 
                        function(mname, mod.ch) 
                          paste0(mname,": ", mod.ch[[mname]]), mod.ch))
  
  cat("\n\n#### Formulae from asreml object\n\n")
  cat(paste0(m.ch, collapse = "\n"), "\n\n\n")
  invisible(m.ch)
}



"as.terms.object" <- function(terms = NULL, asreml.obj = NULL, ...)
{ 
  if (is.character(terms))
  { 
    terms <- as.formula(paste("~ ",terms, sep=""))
  }
  else
  { 
    if (!is.null(terms))  
      terms <- as.formula(terms)
  }
  if (is.null(terms))
    terms.obj <- NULL
  else
  {
    terms.obj <- terms(terms, 
                       keep.order = T, 
                       data = eval(languageEl(asreml.obj$call, which="data"), 
                                   envir = .GlobalEnv), ...)
  }
  return(terms.obj)
}

"rmTermDescription" <- function(term)
{ 
  #Remove description, if there is one, from term in an asreml termlist
  if (length(grep("!", term, fixed=TRUE))!=0) 
    term <- (strsplit(term, "!", fixed=TRUE) )[[1]][1]
  return(term)
}

"getTerms.formula" <- function(form)
{
  terms <- as.character(update.formula(form, ~.))
  terms <- stringr::str_trim(unlist(strsplit(terms[length(terms)], 
                                             split = "+", fixed = TRUE)))
  return(terms)
}

"findterm" <- function(term, termlist, rmDescription=TRUE)
  #This function finds the position of a term in an asreml termlist 
  #It strips off stuff to the right of an !, provided it is not to the right of  
  #a term starting with R!
{ 
  if (length(termlist) == 0 | is.null(termlist))
  { 
    k <- 0
  }
  else
  { 
    if (rmDescription)
    { 
      if(substr(term, 1, 2) != "R!")  term <- rmTermDescription(term)
      k <- which(sapply(termlist, 
                        FUN=function(kterm, term)
                        { 
                          if (substr(kterm, 1, 2) != "R!")  
                            kterm <- rmTermDescription(kterm)
                          haveterm <- setequal(fac.getinTerm(term), 
                                               fac.getinTerm(kterm))
                        }, 
                        term=term))
    }
    else
    { 
      k <- which(sapply(termlist, 
                        FUN=function(kterm, term)
                        { haveterm <- setequal(fac.getinTerm(term), 
                                               fac.getinTerm(kterm))
                        }, 
                        term=term))
    }
    if (length(k) == 0) k <- 0
  }
  return(k)
}

"fac.getinTerm" <- function(term, rmfunction=FALSE)
  #function to return the set of factors/variables in a term separated by ':"
{ 
  if (length(term) != 1)
    stop("Multiple terms supplied where only one allowed")
  vars <- unlist(strsplit(term, ":", fixed=TRUE))
  if (rmfunction)
    vars <- unlist(lapply(vars, rmFunction))
  return(vars)
}

"getTermVars" <- function(term)
  #This function gets the vars in each random term from an asreml termlist 
  #It strips off stuff to the right of an !, provided it is not to the right of  
  #a term starting with R!
{ 
  if(substr(term, 1, 2) != "R!")  term <- rmTermDescription(term)
  vars <- fac.getinTerm(term)
  return(vars)
}

"fac.formTerm" <- function(factors)
  #function to form a term from a set of factors/variables in a structure
{ 
  term <- paste(factors,collapse=":")
  return(term)
}

"separateFunction" <- function(var)
  #A function to separate the name of a function and the argument to the function
{ 
  #Remove description, if there is one, from term in an asreml termlist
  if (length(grep("(", var, fixed=TRUE))!=0) 
  { 
    var <- (strsplit(var, "(", fixed=TRUE) )[[1]]
    var[2] <- (strsplit(var[2], ")", fixed=TRUE) )[[1]][1]
  }
  return(var)
}

"rmFunction" <- function(var, asreml.obj)
  #A function that returns the variable without any function
{ 
  var <- separateFunction(var)
  if (length(var)==2)
  { 
    var <- var[2]
    #Check for further arguments and strip, if found
    if (length(grep(",", var, fixed=TRUE))!=0) 
    { 
      var <- (strsplit(var, ",", fixed=TRUE) )[[1]]  
      var <- var[1]
    } 
  }  
  return(var)
}

"getVarCode" <- function(var, asreml.obj)
  #A function that returns a code for a variable
  # 0 for an integer or numeric with no function
  # 5 for any other variable with no function
  # Added to these codes is 1,2,3 or 4 for lin, pol, spl and dev, respectively
{ 
  sys.funcs <- c("lin","pol","spl","dev")
  #Does it have a function
  if (length(grep("(", var, fixed=TRUE))!=0)
  { 
    var.comp <- separateFunction(var)
    k <- match(var.comp[1], sys.funcs, nomatch = 0)
    var <- var.comp[2]
  } else
    k <- 0
  type <- class(eval(languageEl(asreml.obj$call, which="data"))[[match(var, 
                                                                       colnames(eval(languageEl(asreml.obj$call, which="data"))))]])
  if (type == "integer" | type == "numeric")
    code <- k
  else
    code <- k + 5
  return(code)
}

"getVarsCodes" <- function(term, asreml.obj)
  #A function to call getVarCode for each of the vars in a term that have been stored in character vector
{ 
  codes <- unlist(lapply(term[1:length(term)], getVarCode, asreml.obj = asreml.obj))
}

"num.recode" <- function(x, new.values)
  #function to form a new variate by changing its unique values to new.values
{ 
  x.values <- sort(unique(x))
  nval <- length(x.values)
  if (nval != length(new.values))
    stop("Must supply a new value for every unique value in the supplied vector")
  new.x <- x
  for (i in 1:nval)
    new.x[x == x.values[i]] <- new.values[i]
  return(new.x)
}

"trend.terms.types" <- function(terms, trend.num = NULL, devn.fac = NULL)
  #Examines a set of terms to see if, out of devn.fac and trend.num, whether it 
  # some terms with devn.fac and/or sum with trend.num.
{ 
  has.devn.fac <- has.trend.num <- FALSE
  if (length(terms) > 0)
  { 
    for (term in terms)
    { 
      factors <- fac.getinTerm(term)
      if (!is.null(devn.fac) && any(devn.fac %in%  factors))
        has.devn.fac <- TRUE
      else
        if (!is.null(trend.num) && any(grepl(trend.num, factors, fixed = TRUE)))
          has.trend.num <- TRUE
    }
  }
  term.types <- list(has.devn.fac = has.devn.fac, has.trend.num = has.trend.num)
  return(term.types)
}

"getTermDesignMatrix" <- function(term, asreml.obj)
{
  colnams <- colnames(asreml.obj$design)[startsWith(colnames(asreml.obj$design), term)]
  restnams <- substring(colnams, first = nchar(term)+1)
  colnams <- colnams[grepl("^[0-9]*$", restnams)]
  if (length(colnams) != asreml.obj$noeff[term])
    stop(paste("Error in finding the columns in the design matrix for ",term,sep=""))
  Z <- asreml.obj$design[, colnams]
  return(Z)
}

"makeTestSummary" <- function(which.cols = c("terms","DF","denDF","p","AIC","BIC","action"))
{
  test.summary <- as.data.frame(matrix(nrow = 0, ncol = length(which.cols)))
  names(test.summary) <- which.cols
  class(test.summary) <- c("test.summary", "data.frame")
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

"ginv" <- function(x, tol = .Machine$double.eps ^ 0.5)
{ 
  # computes Moore-Penrose inverse of a matrix
  if (!is.matrix(x) | length(dim(x)) != 2 )
    stop("x must be a matrix")
  svd.x <- svd(x)
  nonzero.x <- (svd.x$d > svd.x$d[1] * tol)
  rank.x <- sum(nonzero.x)
  geninv.x <- matrix(0, dim(x)[1], dim(x)[2])
  if (rank.x)
  { 
    i <- matrix((1:length(nonzero.x))[nonzero.x], rank.x, 2)
    geninv.x[i] <- 1/svd.x$d[nonzero.x]
    if (all(nonzero.x))
      geninv.x <- svd.x$v %*% geninv.x %*% t(svd.x$u)
    else 
      geninv.x <- svd.x$v[, nonzero.x] %*% geninv.x[nonzero.x, nonzero.x] %*% 
      t(svd.x$u[, nonzero.x])
  }
  geninv.x
}
