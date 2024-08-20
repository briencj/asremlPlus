checkTrySpatial <- function(trySpatial)
{
  trySpat.opts <- c("none", "corr", "TPNCSS", "TPPCS", "TPPSC2", "TPP1LS", "TPPSL1", "all")
  trySpatial <- trySpat.opts[unlist(lapply(trySpatial, check.arg.values, options=trySpat.opts))]
  if ("TPPCS" %in% trySpatial)
  { 
    trySpatial <- gsub("TPPCS", "TPPSC2", trySpatial)
    warning("TPPCS in trySpatial has been changed to the new abbreviation TPPSC2")
  }
  if ("TPP1LS" %in% trySpatial)
  { 
    trySpatial <- gsub("TPP1LS", "TPPSL1", trySpatial)
    warning("TPP1LS in trySpatial has been changed to the new abbreviation TPPSL1")
  }
  trySpatial <- unique(trySpatial)
  
  if (length(intersect(trySpatial, trySpat.opts)) == 0)
    stop("trySpatial must be one of ", paste0(trySpat.opts, collapse = ", "))
  if ("all" %in% trySpatial)
    trySpatial <- c("corr", "TPNCSS", "TPPSC2", "TPPSL1")
  if ("none" %in% trySpatial && length(trySpatial) > 1)
    trySpatial <= "none"
  return(trySpatial)
}

calc.nsect <- function(dat, sections)
{
  if (!is.null(sections))
  {    
    if (!is.factor(dat[[sections]]))
      stop("sections, if not NULL, must be a factor")
    else
    {
      nsect <- nlevels(dat[[sections]]) 
    }
  } else #Sections is NULL
    nsect <- 1
  
  return(nsect)
}

checkSections4RanResTerms <- function(asrtests.obj, sections = NULL, asr4)
{
  #Check that separate residual terms have been fitted for multiple sections 
  if (!is.null(sections))
  {
    ran <- getFormulae(asrtests.obj$asreml.obj)$random
    if (!is.null(ran))
    { 
      ran <- getTerms.formula(ran)
      ran <- ran[ran != sections] #Remove sections term, if present
      if (!any((grepl("idh\\(", ran) | grepl("at\\(", ran)) & grepl(sections, ran)))
        warning(paste0("There are no 'idh' or 'at' terms for ", sections, 
                       " in the supplied random model; the fitting may be improved ",
                       "if separate random block terms are specified for each section"))
    }
    if (asr4)
    {  
      res <- getFormulae(asrtests.obj$asreml.obj)$residual
      if (!is.null(res))
      { 
        res <- getTerms.formula(res)
        if (!any((grepl("idh\\(", res) | grepl("dsum\\(", res)) & grepl(sections, res)))
          warning(paste0("There are no 'idh' or 'dsum' terms involving ", sections, 
                         " in the supplied residual model; the fitting may be improved ",
                         "if separate units terms are specified for each section"))
      }
    }  
    else
    {
      res <- getFormulae(asrtests.obj$asreml.obj)$rcov
      if (!is.null(res))
      { 
        res <- getTerms.formula(res)
        if (!any((grepl("idh\\(", res) | grepl("at\\(", res)) & grepl(sections, res)))
          warning(paste0("There are no 'idh' or 'at' terms for ", sections, 
                         " in the supplied rcov model; the fitting may be improved ",
                         "if separate units terms are specified for each section"))
      }
    }
  }
  
  invisible()
}

#Function to identify residual and random correlation model terms (i.e. variance & correlations) 
#   that are currently fitted
#Returns NULL if the residual model cannot accommodate nugget variances
# when there are multiple residual terms, but none include sections
getSectionVpars <- function(asreml.obj, sections, stub, corr.facs, which = c("res", "ran"), 
                            asr4.2)
{
  vpar <- getVpars(asreml.obj, asr4.2 = asr4.2)
  vpc <- vpar$vpc
  vpt <- vpar$vpt
  
  #Get the residual term
  if ("res" %in% which)
  {
    #Get residual variance term using vpc and vpt
    vpc.res <- vpc[grepl("!R$", names(vpc))] #no sections or dsum used
    if (!is.null(sections) && any(grepl(paste0("!", sections), names(vpc))))
      vpc.res <- c(vpc.res, vpc[grepl(paste0("!", sections), names(vpc))]) #idh used
    vpt.res <- vpt[names(vpc.res)]
    vpt.res <- vpt.res[vpt.res == "V"]
    vpc.res <- vpc.res[names(vpt.res)]
    
    if (length(vpc.res) > 1) 
    { 
      if (!is.null(sections))
      { 
        #Determine the terms that include sections
        vpc.facs <- lapply(names(vpc.res), 
                          function(x) 
                          { 
                            x <- stringr::str_split_1(x,"!")[1]
                            x <- fac.getinTerm(x, rmfunction = TRUE, asr4.2 = TRUE)
                            x <- sapply(x, function(fac) stringr::str_split_1(fac, "\\_")[1])
                          })
        names(vpc.facs) <- names(vpc.res)
        which.sectres <- sapply(vpc.facs, 
                                function(x, sections) any(x == sections), 
                                sections = sections)
        #Got any residual terms with sections?
        if (!any(which.sectres))
          vpc.res <- NULL
        else
        {  
          vpc.sect <- vpc.res[which.sectres]
          if (any(grepl(paste0("!",sections), names(vpc.res)))) #! at start implies idh or?
          {
            term <- rmTermDescription(names(vpc.sect[1]))
            nt <- length(startsWith(vpc.res, term))
          } else
            nt <- length(vpc.sect)
          #Get the sections term for the current section
          if (length(vpc.res) == nt)
          { 
            kres.term <- grepl(sections, names(vpc.res)) & grepl(stub, names(vpc.res))
            vpc.res <- vpc.res[kres.term]
          } else
          { 
            warning("There are multiple residual terms, but at least some of these do not involve ",
                    sections, "; consequently the model cannot accommodate nugget variances.")
            vpc.res <- NULL
          }
        }
      }
    }
    if (length(vpc.res) != 1)
      warning("Could not find a residual term for ", sections, " ", stub)
    
  } else #res not required
    vpc.res <- NULL
  
  #Get the random correlation terms (if sections, from all sections)
  if ("ran" %in%  which)
  { 
    vpc.ran <- vpc[!grepl("!R$", names(vpc))]
    if (!is.null(sections) && any(grepl(paste0("!", sections), names(vpc))))
      vpc.ran <- vpc.ran[!grepl(paste0("!", sections), names(vpc.ran))] #idh used
    
    if (length(vpc.ran) > 0)
      vpc.ran <- vpc.ran[sapply(names(vpc.ran),
                                function(term, corr.facs)
                                {
                                  facs <- fac.getinTerm(rmTermDescription(term), asr4.2 = asr4.2)
                                  ran.corr <- 
                                    {
                                      if (length(facs) == length(corr.facs))
                                        all(sapply(corr.facs,
                                                   {
                                                     function(fac, term)
                                                       any(grepl(fac, term))
                                                   }, term = term))
                                      else
                                        FALSE
                                    }
                                  return(ran.corr)
                                }, corr.facs = c(sections, corr.facs))]
    if (length(vpc.ran) == 0)
      vpc.ran <- NULL
  } else
    vpc.ran <- NULL
  return(list(ran = vpc.ran, res = vpc.res))
}

#add set.terms arguments to inargs
#setterms must be a list of single-valued argument settings
#inargs can be a list of multivalued argument settings, but all must be of the same length
addSetterms2inargs <- function(setterms, inargs) 
{
  if ("set.terms" %in% names(inargs)) #already have set.terms in inargs
  {
    if (!all(sapply(setterms, function(x) length(x) == 1)))
      stop("setterms must be a list of single-valued argument settings")
    
    if (setterms$set.terms %in% inargs$set.terms) #same term in set.terms of inargs
    {
      if (length(inargs$set.terms) > 1) #multiple terms in inargs?
      {
        kset.term <- match(setterms$set.terms, inargs$set.terms)
        inargs$ignore.suffices[kset.term] <- setterms$ignore.suffices
        inargs$bounds[kset.term] <- setterms$bounds
        inargs$initial.values[kset.term] <- setterms$initial.values
      } else
      {
        inargs$ignore.suffices <- setterms$ignore.suffices
        inargs$bounds <- setterms$bounds
        inargs$initial.values <- setterms$initial.values
      } 
    } else
    {
      inargs$set.terms <- c(inargs$set.terms, setterms$set.terms)
      inargs$ignore.suffices <- c(inargs$ignore.suffices, setterms$ignore.suffices)
      inargs$bounds <- c(inargs$bounds, setterms$bounds)
      inargs$initial.values <- c(inargs$initial.values, setterms$initial.values)
    }
  } else #no set.terms in inargs
    inargs <- c(inargs, setterms)
  
  return(inargs)
}

#Function to check that a resdual term that is P, does not have is.na(std.error)
isValidResTerm <- function(asreml.obj, vpc.res)
{
  vpc.res <- vpc.res
  vcomp <- summary(asreml.obj)$varcomp
  vcomp <- vcomp[rownames(vcomp) == names(vpc.res), "std.error"]
  res.ok <- !(vpc.res == "P" & is.na(vcomp))
  return(res.ok)
}

#Function to fix a single residual term or the residual variance for current level of section 
chgResTermBound <- function(corr.asrt, sections, stub, asr4, asr4.2, 
                            fitfunc = "changeTerms", 
                            corr.facs, vpc.res, maxit = 30, 
                            allow.unconverged = FALSE, 
                            allow.fixedcorrelation = FALSE,
                            checkboundaryonly = TRUE, 
                            update = TRUE, 
                            IClikelihood = "full", 
                            which.IC = "AIC", ...)
{
  inargs <- list(...)
  #Exclude used fitfunc args from the ellipsis
  # fitfunc.args <- c("asrtests.obj", "newResidual", "label", 
  #                   "set.terms", "initial.values", "bounds", "ignore.suffices", 
  #                   "maxit", "allow.unconverged", "allow.fixedcorrelation",
  #                   "checkboundaryonly", "update", "IClikelihood", "which.IC")
  # other.args <- inargs[setdiff(names(inargs), fitfunc.args)]
  
  bounds.excl <- c("S","B")
  
  #If already fixed, the try P
  if (vpc.res == "F")
  {  bound <- "P"; initial.values = 0.1}
  else 
  {  bound <- "F"; initial.values = 1}
  bound.lab <- ifelse(bound == "F", "fixed", "positive")
  #Add set.terms args to inargs
  inargs <- addSetterms2inargs(setterms = list(set.terms = names(vpc.res), 
                                               ignore.suffices = FALSE, 
                                               bounds = bound, 
                                               initial.values = initial.values),
                               inargs)
  #A single residual variance
  if (is.null(sections))
  { 
    resmod <- as.character(getFormulae(corr.asrt$asreml.obj)$residual)
    if (length(resmod) == 0)
      resmod <- "units"
    else
      resmod <- resmod[2]
    
    if (vpc.res %in% bounds.excl) #will try to fix residual
      fitfunc <- "changeTerms"
    lab <- paste("Try", bound.lab, "nugget (residual) variance")
    if (fitfunc == "changeTerms")
      lab <- gsub("Try", "Force", lab)
    tmp.asrt <- tryCatchLog(
      do.call(fitfunc, 
              c(list(corr.asrt, label = lab, 
                     newResidual = resmod, 
                     maxit = maxit, 
                     allow.unconverged = allow.unconverged, 
                     allow.fixedcorrelation = allow.fixedcorrelation,
                     checkboundaryonly = TRUE, 
                     update = update, 
                     IClikelihood = IClikelihood, 
                     which.IC = which.IC), 
                inargs)),
      error = function(e) {print(paste("Failed attempting to fit correlations to both dimensions;",
                                       "continued analysis without them")); NULL}, 
      include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
    if (!is.asrtests(tmp.asrt))
    {
      test.summary <- addtoTestSummary(corr.asrt$test.summary, terms = lab, 
                                       DF = NA, denDF = NA, p = NA, 
                                       AIC = NA, BIC = NA, 
                                       action = "Unchanged - Singular")
      tmp.asrt <- corr.asrt
      tmp.asrt$test.summary <- test.summary
    }
  } else #sections with multiple residuals - fix the one for this section
  {
    if (vpc.res %in% bounds.excl)
      fitfunc <- "changeTerms"
    kres.term <- names(vpc.res)
    lab <- paste("Try", bound.lab, "nugget with", kres.term)
    #      fitfunc <- "changeModelOnIC"
    if (fitfunc == "changeTerms")
      lab <- gsub("Try", "Force", lab)
    #newResidual not needed here because the residual must be named i.e. cannot be units
    #- only need to set the residual
    tmp.asrt <- tryCatchLog(
      do.call(fitfunc, 
              c(list(corr.asrt, label = lab, 
                     maxit = maxit, 
                     allow.unconverged = allow.unconverged, 
                     allow.fixedcorrelation = allow.fixedcorrelation,
                     checkboundaryonly = TRUE, 
                     update = update, 
                     IClikelihood = IClikelihood, 
                     which.IC = which.IC), 
                inargs)),
      error = function(e) {print(paste("Failed attempting to fit correlations to both dimensions;",
                                       "continued analysis without them")); NULL}, 
      include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
  }
  if (!is.asrtests(tmp.asrt))
  {
    test.summary <- addtoTestSummary(corr.asrt$test.summary, terms = lab, 
                                     DF = NA, denDF = NA, p = NA, 
                                     AIC = NA, BIC = NA, 
                                     action = "Unchanged - Singular")
    tmp.asrt <- corr.asrt
    tmp.asrt$test.summary <- test.summary
  } else
  { 
    if (largeVparChange(tmp.asrt$asreml.obj, 0.75))
      tmp.asrt <- iterate(tmp.asrt)
    #Get bound for Res in tmp.asrt
    vpc.res <- getSectionVpars(tmp.asrt$asreml.obj, 
                               sections = sections, stub = stub,
                               corr.facs = corr.facs,
                               asr4.2 = asr4.2)$res
    if (allow.unconverged || 
        (tmp.asrt$asreml.obj$converge && isValidResTerm(corr.asrt$asreml.obj, vpc.res = vpc.res)))
      corr.asrt <- tmp.asrt
    else
    { 
      lastline <- tail(tmp.asrt$test.summary, n = 1)
      lastline$action <- gsub("Changed", "Unchanged", lastline$action)
      lastline$action <- gsub("Swapped", "Unswapped", lastline$action)
      test.summary <- rbind(corr.asrt$test.summary, lastline)
      corr.asrt$test.summary <- test.summary
    }
  }
  return(corr.asrt)
}

rmRanTerm <- function(corr.asrt, vpbound, 
                      maxit = 30, 
                      allow.unconverged, allow.fixedcorrelation,
                      checkboundaryonly, update,
                      IClikelihood, which.IC, 
                      inargs)
{ 
  for (bound in vpbound)
  { 
    #Find term in random model formula to delete
    ran.terms <- getFormulae(corr.asrt$asreml.obj)$random
    ran.terms <- getTerms.formula(ran.terms)
    facs.bound <- fac.getinTerm(rmTermDescription(bound))
    facs.bound <- gsub('\\\"', "'", facs.bound)
    facs.bound <- gsub('\\(', "\\\\(", facs.bound)
    facs.bound <- gsub('\\)', "\\\\)", facs.bound)
    terms.bound <-  sapply(ran.terms,
                           function(term, facs.bound)
                           {
                             all(sapply(facs.bound,
                                        {
                                          function(fac, term)
                                            grepl(fac, term)
                                        }, term = term))
                           }, facs.bound = facs.bound)
    terms.bound <- names(terms.bound)[terms.bound]
    terms.bound <- paste(terms.bound, collapse = " + ")
    
    #Remove a bound term
    if (!is.null(terms.bound) && terms.bound != "")
      corr.asrt <- do.call(changeTerms,
                           c(list(corr.asrt,
                                  dropRandom = terms.bound,
                                  label = paste("Drop bound",terms.bound),
                                  maxit = maxit, 
                                  allow.unconverged = allow.unconverged,
                                  allow.fixedcorrelation = allow.fixedcorrelation,
                                  checkboundaryonly = FALSE,
                                  update = update,
                                  IClikelihood = IClikelihood,
                                  which.IC = which.IC),
                             inargs))
  }
  return(corr.asrt)
}

#This function assumes that row.factor and col.factor are in the data
makeCorrSpec1D <- function(corr.funcs, corr.orders, dimension, 
                           row.covar, col.covar, row.factor, col.factor, 
                           met.funcs, unimpl.funcs)
{  
  if (corr.funcs[dimension] == "")
    corr1D <- ifelse(dimension == 1, row.factor, col.factor)
  else
  {
    if (corr.funcs[dimension] %in% met.funcs)
      corr1D <- paste0(corr.funcs[dimension],"(",
                       ifelse(dimension == 1, row.covar, col.covar), ")")
    else
    { 
      if (corr.funcs[dimension] != "")
        corr1D <- paste0(corr.funcs[dimension],"(",
                         ifelse(dimension == 1, row.factor, col.factor),")")
    }
    if (grepl("^corb", corr.funcs[dimension])) #at the moment, corb is the only function with an order
    {
      if (is.null(corr.orders) || corr.orders[dimension] == 0)
        corr.orders[dimension] <- 1
      corr1D <- gsub("\\)", paste0(", b = ", corr.orders[dimension], ")"), corr1D)
    }
  }
  return(corr1D)
}  

chk4SingularCorrTerms <- function(asrtests.obj, corr.asrt, label, 
                                  sections, stub, corr.facs, asr4, asr4.2)
{
  if (!is.allnull(asrtests.obj))
  { 
    vpc.corr <- getSectionVpars(asrtests.obj$asreml.obj, 
                                sections = sections, stub = stub, 
                                corr.facs = corr.facs, 
                                asr4.2 = asr4.2)
    
    
    #Determine the correlation terms, if any
    vpt.corr <- getVpars(asrtests.obj$asreml.obj, asr4.2)$vpt
    vpt.ran <- vpt.corr[names(vpc.corr$ran)]
    vpt.r <- vpt.ran[vpt.ran %in% c("R", "P", "C")]
    vpc.r <- vpc.corr$ran[names(vpt.r)]
    #Are there singular r terms
    if (length(vpc.r) > 0 && any(unlist(vpc.r) %in% "S"))
    {
      entry <- getTestEntry(asrtests.obj, label = label)
      entry$action <- "Unchanged - singular term(s)"
      corr.asrt$test.summary <- rbind(corr.asrt$test.summary, entry) 
    } else #no S terms
      corr.asrt <- asrtests.obj
  }
  return(corr.asrt)
}

#ran.term must involve either 2 or 3 variables, at least one of which involves a corb function
#rorder must be zero for the dimension to be changed
#lab must be the llabel for the current fit and be based on ran.term
#dimension can be 0, 1 or 2
# 0 indicates there is only one corb in ran.term so just change it.
# 1 indicates change the part of the ran.term for the first grid dimension. 
# 2 indicates change the part of the ranterm for the second grid dimension. 
"fitCorbPlus1" <- function(corr.asrt, ran.term, rorder, lab, result, dimension = 0, 
                           IClikelihood = "full", trace = FALSE,  ...)
{ 
  inargs <- list(...)
  b <- 1
  last.term <- ran.term
  last.lab <- lab
  if (grepl("corb", ran.term) && (rorder == 0) && #corb with b == 0
      (!grepl("Unswapped", result) && !grepl("Unchanged", result))) #corb is fitted
  {
    if (trace) cat("\n#### Try fitting additional bands to corb\n\n")
    if (!(dimension %in% 0:2))
      stop("dimension must be  eiother 0, 1 or 2")
    if (dimension > 0)
    { 
      ran.parts <- stringr::str_split_1(ran.term, ":")
      if (!(length(ran.parts)%in% 2:3))
        stop("Can only deal with 2 or 3 variables in a random spatial term")
      if (length(ran.parts) == 3)
        dimension <- dimension + 1
    }
    for (b in 2:10)
    {
      last.lab <- lab
      last.term <- ran.term
      if (dimension == 0)
      {     
        if (stringr::str_count(ran.term, "corb") != 1)
          stop("There is not just one corb in ", ran.term, " and dimension is 0")
        lab <- gsub(as.character(b-1), as.character(b), lab)
        ran.term <- gsub(as.character(b-1), as.character(b), ran.term)
      } else
      { 
        #change rand parts for the current dimension
        old.ran.part <- ran.parts[dimension]
        ran.parts[dimension] <-  gsub(as.character(b-1), as.character(b), 
                                      ran.parts[dimension])
        ran.term <- stringr::str_c(ran.parts, collapse = ":")
        lab <- gsub(old.ran.part, ran.parts[dimension], last.lab, fixed = TRUE)
      }
      corr.asrt <- tryCatchLog(
        do.call(changeModelOnIC,
                c(list(corr.asrt, 
                       dropRandom = last.term, 
                       addRandom = ran.term, 
                       label = lab, 
                       allow.fixedcorrelation = FALSE, allow.unconverged = FALSE, 
                       IClikelihood = IClikelihood),
                  inargs)),
        error = function(e) {print("Analysis continued"); NULL}, 
        include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
      if (!grepl("Swapped", getTestEntry(corr.asrt, label = lab)$action, fixed = TRUE)) 
        break
    }
    result <- getTestEntry(corr.asrt, label = last.lab)$action
  }
  return(list(asrt = corr.asrt, last.term = last.term, last.b = (b-1), 
              last.lab = last.lab, result = result))
}

