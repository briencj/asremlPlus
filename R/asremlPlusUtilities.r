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

#Function to allow testing for NULL when objects of length > 1 are possible
is.allnull <- function(x) all(is.null(x))

#Test for a large change in at least one variance parameter
largeVparChange <- function(asreml.obj, threshold = 0.75)
{
  largeChg <- FALSE
  chg <- summary(asreml.obj)$varcomp$"%ch"
  if (!all(is.na(chg)))
  {
    chg <- chg[!is.na(chg)]
    largeChg <- any(chg > threshold)
  }
  return(largeChg)
}

#Check for fixed correlations
isFixedCorrelOK.asreml <- function(asreml.obj, allow.fixedcorrelation = TRUE, ...)  
{
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  correlOK <- TRUE
  
  if (!allow.fixedcorrelation)
  {
    if (asr4.2)
    {
      corrs <- names(asreml.obj$vparameters.type)[asreml.obj$vparameters.type 
                                                  %in% c("R", "P")]
      if (length(corrs) > 0) #have a correlation
      {
        corrs <- asreml.obj$vparameters.con[names(asreml.obj$vparameters.con) 
                                            %in% corrs]
      } 
    }else #not 4.2
    {
      corrs <- names(vpt.char(asreml.obj))[vpt.char(asreml.obj) 
                                           %in% c("R", "P")]
      if (length(corrs) > 0) #have a correlation
      {
        corrs <- vpc.char(asreml.obj)[names(vpc.char(asreml.obj)) 
                                      %in% corrs]
      }
    }
    if (any(corrs %in% c("F", "B", "S")))
      correlOK <- FALSE
  }
  
  return(correlOK)
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
# "isASRemlVersionLoaded" <- function(version = 4, notloaded.fault = FALSE)
# {
#   vers <- getASRemlVersionLoaded(nchar = NULL, notloaded.fault = notloaded.fault)
#   if (!is.null(vers)) vers <- "4.2"
#   if (!is.null(vers))
#   {
#     if (substring(vers, first = 1, last = 3) == "4.0")
#       stop("This version of asremlPlus does not support asreml 4.0.x - asremlPlus 4.0.x does")
#     version <- as.character(version)
#     vers <- substring(vers, first = 1, last = nchar(version)) == version
#   }
#   return(vers)
# }

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

#Checks whether the loaded version is greater than version
"isASReml4_2Loaded" <- function(version = 4.2, notloaded.fault = FALSE)
{
  vers <- getASRemlVersionLoaded(nchar = 3, notloaded.fault = notloaded.fault)
  if (!is.null(vers)) 
    vers <- as.numeric_version(vers) >= 4.2
  return(vers)
}

#Checks if loaded version of asreml begins with version; if not unloads asreml and 
#loads the version in lib.loc 
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

#Function to check that arguments that come in via the ellipsis are all arguments to the allowed functions
#Usually funcs should be the current function and any to which ... is allowed to pass arguments
checkEllipsisArgs <- function(funcs, inargs, pkg = NULL)
{
  inargs <- names(inargs)
  if (length(inargs))
  {
    if (is.null(pkg))
      args <- unique(unlist(lapply(funcs, formalArgs)))
    else
      args <- unique(unlist(lapply(funcs, 
                                   function(func) 
                                     formalArgs(getFunction(func, 
                                                            where = rlang::ns_env(pkg))))))
    args <- args[-match("...", args)]
    foreignArgs <- inargs[!(inargs %in% args)]
    if (length(foreignArgs))
      stop("the argument(s) ", paste0(foreignArgs, collapse = ", "), " are not legal arguments for ", 
           paste0(paste0("'",funcs, collapse = "', "), "'"))
  }
  invisible()
}

checkNamesInData <- function(Names, data)
{
  if (!all(Names %in% names(data)))
  { 
    mess <- paste0("The following required columns are not in data: ", #deparse(substitute(data)))
                   paste0(Names[!(Names %in% names(data))], collapse = ", "), ".\n")
    stop(mess)
  }
  invisible()
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
    if (is.null(asreml.obj))
      terms.obj <- terms(terms, 
                         keep.order = T)
    else
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

"termsplit" <- function(term)
{
  #If a term has neither of these separators, the return vector will consist of the term and ""
  sep <- ifelse(grepl("!R", term), "!R", "!")
  if (grepl(sep, term)) #have a compound name
  {
    term.parts <- strsplit(term, sep)[[1]]
    term.source <- term.parts[1]
    term.end <- gsub(term.source, "", term)
    
  } else #not a compound term
  { 
    term.source <- term
    term.end <- ""
  }
  return(c(term.source, term.end))
}

"findterm" <- function(term, termlist, rmDescription=TRUE)
  #This function finds the position of a term in an asreml termlist 
  #It strips off stuff to the right of an !, provided it is not to the right of  
  #a term starting with R!
{ 
  if (length(termlist) == 0 || is.null(termlist))
  { 
    k <- 0
  }
  else
  { 
    if (rmDescription)
    { 
      #Remove description if not an ASR3 residual
      if (substr(term, 1, 2) != "R!")  term <- rmTermDescription(term)
      k <- which(sapply(termlist, 
                        FUN=function(kterm, term)
                        { 
                          if (substr(kterm, 1, 2) != "R!")  
                            kterm <- rmTermDescription(kterm)
                          if (grepl("\\:",kterm) && grepl("\\:",term))
                            haveterm <- setequal(fac.getinTerm(term), 
                                                 fac.getinTerm(kterm))
                          else
                            haveterm <- term == kterm
                        }, 
                        term=term))
    }
    else
    { 
      #ASR4 Residual term - allow for description after !R - not sure if this occurs
      #Perhaps !R is indicative of a residual variance
      term.parts <- termsplit(term)
      k <- which(sapply(termlist, 
                        FUN=function(kterm, term.parts)
                        { 
                          kterm.parts <- termsplit(kterm)
                          if (grepl("\\:",kterm.parts[1]) && grepl("\\:",term.parts[1]))
                            haveterm.source <- setequal(fac.getinTerm(term.parts[1]), 
                                                        fac.getinTerm(kterm.parts[1]))
                          else
                            haveterm.source <- term.parts[1] == kterm.parts[1]
                          if (haveterm.source && kterm.parts[2] == term.parts[2])
                            haveterm <- TRUE
                          else
                            haveterm <- FALSE
                          return(haveterm)
                        }, term.parts=term.parts))
    }
    if (length(k) == 0) k <- 0
  }
  return(k)
}

"fac.getinTerm" <- function(term, asreml.obj = NULL, rmfunction=FALSE, asr4.2 = FALSE)
  #function to return the set of factors/variables in a term separated by ':"
{ 
  if (length(term) != 1)
    stop("Multiple terms supplied where only one allowed")
  t.obj <- as.terms.object(term, asreml.obj)
  vars <- as.character(parse(text = attr(t.obj, which = "variables")))[-1]
  if (asr4.2) vars <- gsub('\\\"', "'", vars)
  #  vars <- unlist(strsplit(term, ":", fixed=TRUE))
  if (rmfunction)
    vars <- unlist(lapply(vars, rmFunction))
  return(vars)
}

chk4TermInFormula <- function(form, term, asreml.obj)
{
  terms <- getTerms.formula(form)
  facs.terms <- lapply(terms, fac.getinTerm, asreml.obj = asreml.obj)
  names(facs.terms) <- terms
  which.term <- sapply(facs.terms, 
                       function(tfacs, last.facs){all(last.facs %in% tfacs)},
                       last.facs = fac.getinTerm(term))
  if (sum(which.term) == 1) #there is a single term with the same factors with  functions
    term <- terms[which.term]
  return(term)
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
  # 0 for an integer or numeric with no sys function
  # 5 for any other variable with no sys function
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
  type <- class(eval(languageEl(asreml.obj$call, 
                                which="data"))[[match(var, 
                                                      colnames(eval(languageEl(asreml.obj$call, 
                                                                               which="data"))))]])
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

"getTermDesignMatrix" <- function(term, asreml.obj, split = ":")
{
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  if (asr4.2)
  { 
    vars <- fac.getinTerm(term, asr4.2 = asr4.2)
    nvars <- length(vars)
    colnams <- colnames(asreml.obj$design)
    col.vars <- strsplit(colnams, split = split, fixed = TRUE)
    #reduce to columns of same length as term
    len.nvars <- unlist(lapply(col.vars, length)) == nvars
    colnams <- colnams[len.nvars]
    col.vars <- col.vars[len.nvars]
    col.vars <- as.data.frame(do.call(rbind, col.vars))
    #Determine which columns have the same vars, in any order, as the term; select colnams that do
    which.cols <- apply(col.vars, MARGIN = 1, 
                        FUN = function(krow, vars)
                        {
                          all(sapply(vars, function(v, krow) any(startsWith(krow, v)), krow = krow))
                        }, vars = vars)
    colnams <- colnams[which.cols]
  }
  else
  { 
    colnams <- colnames(asreml.obj$design)[startsWith(colnames(asreml.obj$design), term)]
    restnams <- substring(colnams, first = nchar(term)+1)
    colnams <- colnams[grepl("^[0-9]*$", restnams)]
  }
  #grp terms do not have a known number of columns
  if (term %in% names(asreml.obj$noeff) && length(colnams) != asreml.obj$noeff[term])
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

#These are asreml functions, included here to save having to use asreml::
guzpfx <- function (i) 
{
  con <- c(" ", "P", "?", "U", "F", "C", "S", "B")
  pc <- apply(matrix((i + 1), ncol = 1), 1, function(x) {
    if (is.na(x)) 
      return(3)
    y <- max(x, 1)
    if (y > 8 && y%%10 == 2) 
      y <- 3
    if (y > 8) 
      y <- 8
    y
  })
  con[pc]
}

vpc.char <- function (object) 
{
  if (!inherits(object, "asreml")) 
    stop("'object' must be of class 'asreml'")
  nm <- names(object$vparameters.con)
  
  ch <- guzpfx(object$vparameters.con)
  names(ch) <- nm
  return(ch)
}

vpt.char <- function (object) 
{
  if (!inherits(object, "asreml")) 
    stop("'object' must be of class 'asreml'")
  nm <- names(object$vparameters.type)
  typeCode <- c("V", "G", "R", "C", "P", "L")
  ch <- typeCode[object$vparameters.type]
  names(ch) <- nm
  attr(ch, "legend") <- data.frame(Code = typeCode, 
                                   Type = c("variance", "variance ratio", "correlation", 
                                            "covariance", "positive correlation", "loading"), 
                                   stringsAsFactors = FALSE)
  return(ch)
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

"setToZero" <- function(x, zero.tolerance = .Machine$double.eps ^ 0.5)
{
  zeroes <- abs(x) < zero.tolerance
  if (any(zeroes))
    x[zeroes] <- 0
  return(x)
}
subset.list <- function(x, select = 1:length(x), ...)
{
  selx <- select
  if (inherits(select, what = "logical"))
    selx <- c(1:length(x))[select]
  y <- lapply(selx, function(k, x){x[[k]]}, x = x)
  namx <- names(x)
  if (!is.null(namx))
  {
    if (inherits(select, what = "character"))
      selx <- namx %in% select
    names(y) <- names(x)[selx]
  }
  return(y)
}
