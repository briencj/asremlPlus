.asremlPlusEnv <- new.env()
#Set up allowed specials
corr.specs <- c("id","ar1", "ar2", "ar3", "sar","sar2",
                "ma1", "ma2", "arma", "cor", "corb", "corg",
                "exp", "gau", "lvr", "iexp", "igau", "ieuc", 
                "sph", "circ", "aexp", "agau")
var.specs <- c("diag", "us", "ante", "chol", "sfa", 
               "fa", "facv", "rr", "vm", "dsum")
mod.specs <- c("and","at", "C", "con", "dev", "gpf", "grp", 
               "leg", "lin", "ma", "mbf", "pol", "pow", 
               "sbs", "spl", "uni")
all.specs <- c(corr.specs, paste0(corr.specs,"v"), 
                  paste0(corr.specs,"h"), var.specs,
                  mod.specs)
assign("corr.specials",  corr.specs, envir=.asremlPlusEnv)
assign("var.specials",  var.specs, envir=.asremlPlusEnv)
assign("mod.specials",  mod.specs, envir=.asremlPlusEnv)
assign("all.specials",  all.specs, envir=.asremlPlusEnv)

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

#Function to test for compound symmetry in a matrix
isCompoundSymmetric.matrix <- function(object, tol = 100 * .Machine$double.eps, ...)
{
  cs <- FALSE
  if (isSymmetric(object))
  { 
    d <- diag(object)
    offd <- object[upper.tri(object)]
    cs <- all(abs(d[-1] - d[1]) < 1e-08) && all(abs(offd[-1] - offd[1]) < 1e-08)
  }
  return(cs)
}



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

#Checks whether the loaded version is greater than or equal "4.2" 
# - version is ignored, but is retained because it is throughout the calling code
"isASReml4_2Loaded" <- function(version = "4.2", notloaded.fault = FALSE)
{
  vers <- getASRemlVersionLoaded(nchar = 3, notloaded.fault = notloaded.fault)
  if (!is.null(vers)) 
    vers <- as.numeric_version(vers) >= "4.2"
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
      args <- unique(unlist(lapply(
        funcs, 
        function(f) 
        {
          args <- { if (f == "asreml.options")
            names(asreml::asreml.options()) #not being used, but left in case
          else
            formalArgs(f)}
        })))
    else
      args <- unique(unlist(lapply(funcs, 
                                   function(f) 
                                     args <- 
                                     { 
                                       if (f == "asreml.options")
                                         names(asreml::asreml.options()) #not being used, but left in case
                                       else
                                         formalArgs(getFunction(f, 
                                                            where = rlang::ns_env(pkg)))})))
    if ("..." %in% args)
      args <- args[-match("...", args)]
    foreignArgs <- inargs[!(inargs %in% args)]
    if (length(foreignArgs))
      stop("the argument(s) ", paste0(foreignArgs, collapse = ", "), " are not legal arguments for ", 
           paste0(paste0("'",funcs, collapse = "', "), "'"))
  }
  invisible()
}

getTpsmmb.args <- function(inargs)
{   
  if (length(names(inargs)))
  {
    tpsmmb.args <- formalArgs(getFunction("tpsmmb", where = rlang::ns_env("asremlPlus")))
    tpsmmb.args <- inargs[intersect(tpsmmb.args, names(inargs))]
  } else
    tpsmmb.args <- list()
  return(tpsmmb.args)
}


#Function to separate checking of tpsmmb args from other arguments
# - needed because tpsmmb is not an asremlPlus exported function
checkEllipsisArgs_tpsmmb <- function(funcs, inargs, pkg = NULL)
{
  if (length(names(inargs)))
  {
    tpsmmb.args <- getTpsmmb.args(inargs)
    other.inargs <- inargs[setdiff(names(inargs), names(tpsmmb.args))]
    checkEllipsisArgs(funcs, other.inargs, pkg = NULL)
  } else
    tpsmmb.args <- list()
  invisible(tpsmmb.args)
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

"getTerms.formula" <- function(form, ...)
{
  # terms <- as.character(update.formula(form, ~.))
  # terms <- stringr::str_trim(unlist(strsplit(terms[length(terms)], 
  #                                            split = "+", fixed = TRUE)))
  form <- update.formula(form, ~.)
  terms <- attr(as.terms.object(form, ...), "term.labels")
  terms <- gsub('\\\"', "'", terms)
  
  return(terms)
}

#So far, this function only needs to deal with str terms
"convTerms2Vparnames" <- function(terms)
{
  
  #Deal with any str function in the formula
  if (any(grepl("str\\(", terms)))
  {
    terms <- sapply(terms, function(term)
    {
      if (grepl("str\\(", term))
      {
        start <- stringr::str_locate(term, "~")[1,1]
        end <- stringr::str_locate(term, ",")[1,1]
        term <- stringr::str_sub(term, start = start, end = end-1)
        term <- as.formula(paste("~",term))
        term <- (update.formula(term, ~.))
        term <- gsub(" ", "", term)[2]
      }
      return(term)
    })
  }
  
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

getVpars <- function(asreml.obj, asr4.2)
{
  if (asr4.2)
  { 
    vpc <- asreml.obj$vparameters.con
    vpt <- asreml.obj$vparameters.type
  } else
  {
    vpc <- vpc.char(asreml.obj)
    vpt <- vpt.char(asreml.obj)
  }
  return(list(vpc = vpc, vpt = vpt))
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

#Function that replaces startsWith, startsWith not being able to deal with special regular expressiion characters
"beginsWith" <- function(x, prefix)
  grepl(paste0("^",prefix), x)


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

#This version uses svd::propack.svd which doesn't crash as much as svd
"ginv" <- function(x, tol = .Machine$double.eps ^ 0.5)
{ 
  # computes Moore-Penrose inverse of a matrix
  if (!inherits(x, what = "matrix") || length(dim(x)) != 2)
    stop("x must be a matrix")
  #propack.svd return only the nonzero eigenvalues and their eigenvectors
  suppressWarnings(
    svd.x <- svd::propack.svd(x, opts = list(tol = tol)))
  nonzero.x <- (svd.x$d > svd.x$d[1] * tol)
  if (!all(nonzero.x))
  {
    svd.x$d <- svd.x$d[nonzero.x]
    svd.x$u <- svd.x$u[,nonzero.x]
    svd.x$v <- svd.x$v[,nonzero.x]
  }
  rank.x <- length(svd.x$d)
  #set diagonal elements of geninv.x, which may not be square
  geninv.x <- matrix(0, rank.x, rank.x)
  diag(geninv.x) <- 1/svd.x$d
  #form inverse
  if (isSymmetric(x))
    geninv.x <- svd.x$v %*% geninv.x %*% t(svd.x$v)
  else
    geninv.x <- svd.x$v %*% geninv.x %*% t(svd.x$u)
  
  return(geninv.x)
}

"ginv" <- function(x, tol = .Machine$double.eps ^ 0.5)
{ 
  # computes Moore-Penrose inverse of a matrix
  if (!inherits(x, what = "matrix") || length(dim(x)) != 2)
    stop("x must be a matrix")
  svd.x <- tryCatchLog(svd(x),
                       error = function(e) 
                       {print("Computing a generalized inverse using svd has failed; NA returned"); NA}, 
                       include.full.call.stack = FALSE, include.compact.call.stack = TRUE)
  if (all(!is.na(svd.x)))
  { 
    nonzero.x <- (svd.x$d > svd.x$d[1] * tol)
    rank.x <- sum(nonzero.x)
    geninv.x <- matrix(0, dim(x)[1], dim(x)[2])
    if (rank.x)
    { 
      #set diagonal elements of geninv.x
      i <- matrix((1:length(nonzero.x))[nonzero.x], rank.x, 2)
      geninv.x[i] <- 1/svd.x$d[nonzero.x]
      #form inverse
      if (all(nonzero.x))
        geninv.x <- svd.x$v %*% geninv.x %*% t(svd.x$u)
      else 
      {  
        v <- svd.x$v[, nonzero.x]
        if (rank.x == 1)
          v <- matrix(v, ncol = 1)
        geninv.x <- v %*% geninv.x[nonzero.x, nonzero.x] %*% t(v)
      }
    }
    attr(geninv.x, which = "rank") <- rank.x
  } else
    geninv.x <- NA
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
