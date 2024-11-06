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
getSectionVpars <- function(asreml.obj, sections, stub, sterm.facs, which = c("res", "ran"), 
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
  
  #Get the random spatial terms (if sections, from all sections)
  if ("ran" %in%  which)
  { 
    vpc.ran <- vpc[!grepl("!R$", names(vpc))]
    if (!is.null(sections) && any(grepl(paste0("!", sections), names(vpc))))
      vpc.ran <- vpc.ran[!grepl(paste0("!", sections), names(vpc.ran))] #idh used
    
    if (length(vpc.ran) > 0)
      vpc.ran <- vpc.ran[unlist(sapply(names(vpc.ran),
                                       function(term, sterm.facs)
                                       {
                                         facs <- fac.getinTerm(rmTermDescription(term), 
                                                               rmfunction = TRUE, 
                                                               asr4.2 = asr4.2)
                                         ran.corr <- setequal(sterm.facs, facs)
                                         return(ran.corr)
                                       }, sterm.facs = c(sections, sterm.facs)))]
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

#Function to determine if a correlation has been fitted by looking for correlation terms 
#of the form spat.term!.
setCorrTerm <- function(asreml.obj, spat.var, asr4.2)
{
  corr.term <- FALSE
  # vpc.ran <- getSectionVpars(asreml.obj, 
  #                            sections = sections, stub = stub, 
  #                            which = "ran", asr.4.2 = asr.4.2)$ran
  # 
  # if (any((sapply(paste0(spat.term, "!", sterm.facs), grepl, names(vpc.ran)))))
  #   corr.term <- TRUE
  
  #Get correlation terms
  vpt.corr <- getVpars(asreml.obj, asr4.2)$vpt
  spat.term <- names(findterm(spat.var, names(vpt.corr),  #allows for changed order
                              rmDescription = FALSE))
  vpt.corr <- vpt.corr[vpt.corr %in% c("R", "P")]
  
  #Any involve "spat.term!"?
  spat.term <-   gsub("\\\'", "\\\\'", spat.term)
  spat.term <-   gsub("\\(", "\\\\(", spat.term)
  spat.term <-   gsub("\\)", "\\\\)", spat.term)
  corr.term <-  (length(vpt.corr) > 0 && length(spat.term) > 0  && 
                   any(grepl(paste0(spat.term, "!"), names(vpt.corr))))
  return(corr.term)
}

#Function to check that a residual term that is P, does not have is.na(std.error)
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
                            sterm.facs, vpc.res, maxit = 30, 
                            allow.unconverged = FALSE, 
                            allow.fixedcorrelation = FALSE,
                            checkboundaryonly = TRUE, 
                            update = TRUE, 
                            IClikelihood = "full", 
                            which.IC = "AIC", bounds.excl, ...)
{
  inargs <- list(...)
  #Exclude used fitfunc args from the ellipsis
  # fitfunc.args <- c("asrtests.obj", "newResidual", "label", 
  #                   "set.terms", "initial.values", "bounds", "ignore.suffices", 
  #                   "maxit", "allow.unconverged", "allow.fixedcorrelation",
  #                   "checkboundaryonly", "update", "IClikelihood", "which.IC")
  # other.args <- inargs[setdiff(names(inargs), fitfunc.args)]
  
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
      error = function(e) {print(paste(
        "Failed attempting to change bound on nugget (residual) variance;",
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
                               sterm.facs = sterm.facs,
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
    if (!is.allnull(ran.terms))
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

#Function to remove individual bound correlation terms in fitting spatial models
rmboundCorrVpar <- function(corr.asrt, vpcbound,  maxit = 30, 
                            allow.unconverged, allow.fixedcorrelation,
                            checkboundaryonly, update,
                            IClikelihood, 
                            inargs)
{ 
  if (!allow.fixedcorrelation)
  { 
    asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
    asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
    rfuncs <- c(get.specials("tseries.specials"), get.specials("met.specials"))
    ignored.rfuncs <- c("corb", "corg")
    ignored.rfuncs <- c(ignored.rfuncs, 
                        sapply(ignored.rfuncs, function(f) paste0(f, c("v","h"))))
    rfuncs <- rfuncs[-match(ignored.rfuncs, rfuncs)]
    
    #Need to first reduce to unique correlation functions
    vpc.bC <- sapply(vpcbound, 
                     function(x)
                       paste(stringr::str_split_1(
                         x, "\\!")[1:2], collapse = "!"))
    vpc.bC <- unique(vpc.bC)
    
    #Get the factors in the bound variance parameters
    facs.vpc.bC <- lapply(vpc.bC, function(bC) 
    { 
      bC <- fac.getinTerm(rmTermDescription(bC), asr4.2 = asr4.2)
    })
    names(facs.vpc.bC) <- vpc.bC
    
    #Reduce random terms to those that include a bound correlation
    asreml.obj <- corr.asrt$asreml.obj
    terms <- getTerms.formula(asreml.obj$call$random)
    terms <- unlist(
      lapply(terms, 
             function(kterm, facs.vpc.bC, asreml.obj)
             {
               term <- NULL
               facs.kterm <- convTerm2VparFacs(kterm, asr4.2 = asr4.2)
               haveterm <- sapply(facs.vpc.bC, 
                                  function(facs.bC, facs.kterm) 
                                    setequal(facs.bC, facs.kterm), 
                                  facs.kterm = facs.kterm)
               if (any(haveterm)) term <- kterm
               return(term)
             }, facs.vpc.bC = facs.vpc.bC, asreml.obj = asreml.obj))
    
    #Cannot find a term to remove
    if (is.allnull(terms)) 
    { 
      if (asr4)
        kresp <- asreml.obj$formulae$fixed[[2]]
      else
        kresp <- asreml.obj$fixed.formula[[2]]    
      
      warning("\nIn analysing ", kresp,
              ", cannot remove the following fixed/boundary/singular parameters: ", 
              vpcbound, "\n\n")
    }
    else
    {
      #Modify the terms in the random formula to remove the correlation functions corresponding 
      #   to the bound correlations 
      new.terms <- terms
      for (bC in vpc.bC)
      {              
        #find term in random formula corresponding to current bound correlation
        which.term  <- sapply(
          new.terms, 
          function(kterm, facs.term, asreml.obj) 
            setequal(facs.term, convTerm2VparFacs(kterm, asr4.2 = asr4.2)), 
          facs.term = facs.vpc.bC[[bC]], asreml.obj = asreml.obj)
        if (!any(which.term))
        { 
          if (asr4)
            kresp <- asreml.obj$formulae$fixed[[2]]
          else
            kresp <- asreml.obj$fixed.formula[[2]]    
          warning("\nIn analysing ", kresp,
                  ", cannot remove the following fixed/boundary/singular parameter: ", 
                  bC, "\n\n")
        } else 
        {
          #Find if the dimension of the current bound correlation has a (correlation) function
          bC.dim <- stringr::str_split_i(bC, "\\!", 2)
          ran.corr <- new.terms[which.term]
          comps.corr <- fac.getinTerm(ran.corr, 
                                      asreml.obj = asreml.obj, asr4.2 = asr4.2)
          facs.corr <- fac.getinTerm(ran.corr, 
                                     asreml.obj = asreml.obj, 
                                     rmfunction = TRUE, asr4.2 = asr4.2)
          which.ran.dim <- grep(bC.dim, facs.corr)
          #Does ran,corr have a function in the dimension
          if (stringr::str_detect(comps.corr[which.ran.dim], "(?<=[:alnum:])\\(")) 
            #(separateFunction(comps.corr[which.ran.dim])[1] %in% rfuncs)
          {
            func <- separateFunction(comps.corr[which.ran.dim])[1]
            #Remove (correlation) function for a bound dimension
            if (func %in% rfuncs)
            { 
              comps.corr[which.ran.dim] <- facs.corr[which.ran.dim]
              #if the correlation func ends in "h", then add the function "idh" to the dimension
              if (grepl("h$", func))
                comps.corr[which.ran.dim] <- paste0("idh(", comps.corr[which.ran.dim], ")")
              #Make the new term
              new.terms[which.term] <- fac.formTerm(comps.corr)
            } else
              new.terms[which.term] <- new.terms[!which.term]
          }
        }
      }
      
      #If have any terms to use in replacing old terms
      if (length(new.terms) > 0)
      { 
        #Find terms to be replaced by new.terms
        old.terms <- unlist(
          sapply(new.terms, 
                 function(knew.term, terms)
                 {
                   term <- unlist(
                     lapply(terms, 
                            function(kterm, knew.term)
                            {
                              term <-NULL
                              if (all(convTerm2VparFacs(kterm, asr4.2 = asr4.2) == 
                                      convTerm2VparFacs(knew.term, asr4.2 = asr4.2)))
                                term <- kterm
                              return(term)
                            }, knew.term = knew.term))
                 }, terms = terms))
        if (length(old.terms) == 0)
          old.terms <- NULL
        
        #Check if any replacement terms no longer involve a correlation function
        is.corrfunc <- unlist(
          sapply(new.terms, 
                 function(kterm, asreml.obj)
                 {
                   comps.kterm <- fac.getinTerm(kterm, asreml.obj = asreml.obj, 
                                                asr4.2 = asr4.2)
                   havecorr <- sapply(comps.kterm, 
                                      function(kcomp) 
                                      {
                                        havecorr <- FALSE
                                        if(stringr::str_detect(kcomp, "(?<=[:alnum:])\\(") && 
                                           separateFunction(kcomp)[1] %in% rfuncs)
                                          havecorr <- TRUE
                                        return(havecorr)
                                      })
                   havecorr <- any(havecorr)
                   return(havecorr)
                 }, asreml.obj = asreml.obj))
        new.terms <- new.terms[is.corrfunc]
        if (length(new.terms) == 0)
          new.terms <- NULL
        else
          new.terms <- paste(new.terms, collapse = " + ")
        old.terms <- paste(old.terms, collapse = " + ")
        
        #Remove correlations that are bound from terms
        if (!is.null(old.terms) || !is.null(new.terms))
        { 
          if (is.allnull(new.terms))
            lab <- "Drop bound spatial term(s)"
          else
            lab <- "Drop bound correlation from spatial term(s)"
          corr.asrt <- do.call(changeTerms,
                               c(list(corr.asrt,
                                      dropRandom = old.terms,
                                      addRandom = new.terms,
                                      label = lab,
                                      maxit = maxit, 
                                      allow.unconverged = allow.unconverged,
                                      allow.fixedcorrelation = allow.fixedcorrelation,
                                      checkboundaryonly = FALSE,
                                      update = update,
                                      IClikelihood = IClikelihood),
                                 inargs))
        }
      }
    }
  }
  return(corr.asrt)
}

#This function assumes that row.factor and col.factor are in the data
makeCorrSpec1D <- function(corr.funcs, corr.orders, dimension, 
                           row.covar, col.covar, row.factor, col.factor)
{  
  if (corr.funcs[dimension] == "")
    corr1D <- ifelse(dimension == 1, row.factor, col.factor)
  else
  {
    if (corr.funcs[dimension] %in% get.specials("met.specials"))
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

# (i) Checks for singular spatial variance and attempts to set it to F
# (ii) Checks for singular R, P, C terms: 
#      asrtests.obj should be for the model with the correlation term fitted.
#         if random is S, returns asrtest.obj with an entry in test summary; 
#         if residual is S, attempts to change it to F using chgResTermBound
# corr.asrt should be for the model before a correlation term was fitted;
chk4SingularSpatTerms <- function(asrtests.obj, corr.asrt, label, 
                                  sections, stub, sterm.facs, asr4, asr4.2, 
                                  maxit = 30, 
                                  allow.unconverged = FALSE, 
                                  allow.fixedcorrelation = FALSE,
                                  checkboundaryonly = TRUE, 
                                  update = TRUE, 
                                  IClikelihood = "full", 
                                  which.IC = "AIC", 
                                  bounds.excl, sing.excl)
{
  if (!is.allnull(asrtests.obj))
  { 
    vpc.corr <- getSectionVpars(asrtests.obj$asreml.obj, 
                                sections = sections, stub = stub, 
                                sterm.facs = sterm.facs, 
                                asr4.2 = asr4.2)
    vpt.corr <- getVpars(asrtests.obj$asreml.obj, asr4.2)$vpt

    #Is spatial variance singular? Try to set fixed
    vpc.var <- vpc.corr$ran
    vpc.var <- vpc.var[names(vpc.var) %in% intersect(names(vpt.corr)[vpt.corr == "G"], 
                                                     names(vpc.var))]
    if (length(vpc.var) == 1 && vpc.var %in% sing.excl)
    {
      tmp.asrt <- changeTerms(asrtests.obj, set.terms = names(vpc.var), bounds = "F",
                              initial.values = 1, ignore.suffices = FALSE, 
                              update = FALSE, 
                              allow.unconverged = allow.unconverged, 
                              allow.fixedcorrelation = TRUE,
                              checkboundaryonly = checkboundaryonly)
      #Only if the spatial variance is F in tmp.asrt, make it the current fit
      if (tmp.asrt$asreml.obj$vparameters[names(vpc.var)] == "F")
      { 
        asrtests.obj <- tmp.asrt
        vpc.corr <- getSectionVpars(asrtests.obj$asreml.obj, 
                                    sections = sections, stub = stub, 
                                    sterm.facs = sterm.facs, 
                                    asr4.2 = asr4.2)
        vpt.corr <- getVpars(asrtests.obj$asreml.obj, asr4.2)$vpt
      }
    }

    #Determine the correlation terms, if any
    vpt.ran <- vpt.corr[names(vpc.corr$ran)]
    vpt.r <- vpt.ran[vpt.ran %in% c("R", "P", "C")]
    vpc.r <- vpc.corr$ran[names(vpt.r)]
    #Are there singular r terms
    if (length(vpc.r) > 0 && any(unlist(vpc.r) %in% sing.excl))
    {
      entry <- getTestEntry(asrtests.obj, label = label)
      entry$action <- "Unchanged - singular term(s)"
      corr.asrt$test.summary <- rbind(corr.asrt$test.summary, entry) 
    } else #no Singular correlation terms
      corr.asrt <- asrtests.obj
    if (length(vpc.corr$res) > 0 && (vpc.corr$res %in% c("B",sing.excl)))
    {  
      corr.asrt <- chgResTermBound(corr.asrt, sections = sections, stub = stub, 
                                   asr4 = asr4, asr4.2 = asr4.2, 
                                   fitfunc = "changeTerms", 
                                   sterm.facs = sterm.facs, vpc.res = vpc.corr$res, 
                                   maxit = 30, 
                                   allow.unconverged = allow.unconverged, 
                                   allow.fixedcorrelation = allow.fixedcorrelation,
                                   checkboundaryonly = checkboundaryonly, 
                                   update = update, 
                                   IClikelihood = IClikelihood, 
                                   which.IC = which.IC, bounds.excl = bounds.excl)
    }
  }
  return(corr.asrt)
}

#A function that checks whether all spatial vpars for a correlation model for a section 
#  are in all.bounds.excl. If they are, this is considered irretrievable and the 
#  initial model is retained.
allBoundSectionVpars <- function(asrtests.last, asrtests.prev, 
                                 lab, action = "Swapped - all bound", 
                                 sections = NULL, stub = NULL, 
                                 sterm.facs, all.bounds.excl, asr4.2)
{
  #If all correlation model terms are bound, reinstate initial model
  vpars <- unlist(getSectionVpars(asrtests.last$asreml.obj, sections = sections, stub = stub, 
                                  sterm.facs = sterm.facs, asr4.2 = asr4.2))
  if (all(vpars %in% all.bounds.excl))
    asrtests.last <- revert2previousFit(asrtests.last, asrtests.prev, terms = lab,
                                    action = action)
  return(asrtests.last)
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

