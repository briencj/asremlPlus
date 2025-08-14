
"findEntry.asrtests" <- function(asrtests.obj, label, ...)
{
  tests <- asrtests.obj$test.summary
  #Assume if label has spaces then cannot be just a term that is factors separated by colons. 
  #Thus labels that are a combination of text and terms or terms with spaces (e.g. str and 
  #  at) must match the label exactly. 
  #  - the problem is to identify terms and separate them from ordinary text when there are 
  #    spaces in label.
  if (grepl(" ", label)) 
    k <- tail(which(as.character(tests$terms)==label),1)
  else #no spaces in label
  {  
    k <- 0
    nospace.terms <- !grepl(" ",as.character(asrtests.obj$test.summary$terms))
    if (any(nospace.terms))
    {  
      tests <- tests[nospace.terms,]
      #findterm allows for differences in factor order between label and tests$terms
      k <- tail(findterm(label, as.character(tests$terms)),1) 
    }
  }
  if (length(k) == 0 || k == 0)
    entry <- NULL
  else
    entry <- tests[k,]
  return(entry)
}

"getTestPvalue.asrtests" <- function(asrtests.obj, label, ...)
{
  entry <- findEntry.asrtests(asrtests.obj, label, ...)
  if (is.allnull(entry))
      stop(label, " not found in test.summary of supplied asrtests.obj")
  p <- entry$p
  return(p)
}

"getTestEntry.asrtests" <- function(asrtests.obj, label, error.absent = TRUE, ...)
{
  entry <- findEntry.asrtests(asrtests.obj, label, ...)
  if (is.allnull(entry))
  {
    if (error.absent)
      stop(label, " not found in test.summary of supplied asrtests.obj")
    else 
      entry <- NULL
  } else
  { 
    class(entry) <- "data.frame"
  }
  return(entry)
}

"recalcWaldTab.asrtests" <- function(asrtests.obj, recalc.wald = FALSE, 
                                     denDF="numeric", dDF.fault = "none", 
                                     dDF.values = NULL, trace = FALSE, ...)
{ 
  tempcall <- list(...)
  if (length(tempcall)) 
  {
    #Check for deprecated dDF.na argument
    if ("dDF.na" %in% names(tempcall))
      stop("The argument dDF.na has been deprecated; use dDF.fault instead.")
  }

  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  if (is.null(asrtests.obj))
    stop("Must supply an asrtests object")
  
  #Initialize
  asreml.obj <- asrtests.obj$asreml.obj
  #Call wald.asreml if recalc.wald is TRUE
  if (recalc.wald)
  {
    wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = FALSE, ...)
    wald.tab <- chkWald(wald.tab)  
  }
  else #extract wald.tab from the asrtests object
    wald.tab <- asrtests.obj$wald.tab
  nofixed <- dim(wald.tab)[1]
  hd <- attr(wald.tab, which = "heading")
  options <- c("none", "residual", "maximum", "supplied")
  opt <- options[check.arg.values(dDF.fault, options)]
  if (opt == "supplied" & is.null(dDF.values))
    stop('Need to set dDF.values because have set dDF.fault = \"supplied\"')
  if (opt == "supplied")
    if (length(dDF.values) != nofixed)
      stop("Number of  dDF.values must be the same as the number of terms in wald.tab")
  den.df <- NA
  if (is.null(wald.tab) || nrow(wald.tab) == 0)
  {
    wald.tab <- as.data.frame(matrix(nrow = 0, ncol = 4))
    names(wald.tab) <- c("Df" , "denDF","F.inc","Pr")
    warning("Wald.tab is empty - probably the calculations have failed")
  } else
  {
    if (!("denDF" %in% colnames(wald.tab))) #no denDF
    { 
      #Get denDF
      if (opt == "supplied")
      { wald.tab$denDF <- dDF.values
      warning("Supplied dDF.values used for denDF")
      }
      else
        if (opt == "maximum" | opt == "residual") 
        { 
          wald.tab$denDF <- asreml.obj$nedf
          warning("Residual df used for denDF")
        }
    }
    else #have denom. df
    { 
      den.df.fault <- is.na(wald.tab$denDF) | is.infinite(wald.tab$denDF)  | wald.tab$denDF < 0.01
      if (any(den.df.fault)) #some denDf are NA or < 0.01 and so need to use approximate DF
      { 
        if (opt == "supplied")
        { 
          wald.tab$denDF[den.df.fault] <- dDF.values[den.df.fault]
          warning("At least some supplied dDF.values used for denDF that are NA, Inf or < 0.01")
        }
        else
        { 
          if (opt == "maximum") 
          { 
            if (!all(den.df.fault))
            { 
              den.df <- max(wald.tab$denDF[-1], na.rm=TRUE) 
              warning("Maximum denDF used for for denDF that are NA, Inf or < 0.01")
            }
            else
            { 
              den.df <- asreml.obj$nedf
              warning("Residual df used for denDF for denDF that are NA, Inf or < 0.01")
            }
          }
          else
          { 
            if (opt == "residual")
            { 
              den.df <- asreml.obj$nedf
              warning("Residual df used for denDF for at least some terms")
            }
          }
          wald.tab$denDF[den.df.fault] <- den.df
        }
      }
    }
    attr(wald.tab, which = "heading") <- hd
    if (nrow(wald.tab) == 0)
      warning("Wald.tab is empty - probably the calculations have failed")
    else
    {
      #Calc Pr
      if ("denDF" %in% colnames(wald.tab))
      {
        if ("F.con" %in% colnames(wald.tab))
          wald.tab$Pr <- 1 - pf(wald.tab$F.con, wald.tab$Df, wald.tab$denDF)
        else
          wald.tab$Pr <- 1 - pf(wald.tab$F.inc, wald.tab$Df, wald.tab$denDF)
      }
    }
  }
  return(wald.tab)
}


"powerTransform" <- function(var.name, power = 1, offset = 0, scale = 1, 
                             titles = NULL, data)
  #Function to perform a power transformation on a variable whose name is given as 
  #a character string in var.name. The transformed variable is stored in data
{ 
  k <- match(var.name, names(data))
  if (!is.null(titles) & !is.na(match(var.name, names(titles))))
    title <- titles[[var.name]]
  else
    title <- names(data)[k]
  #Get current variable and its name and title
  tvar.name <- var.name
  ttitle <- title
  y <- data[[k]]
  #Apply scale transformation
  if (scale != 1)
  { 
    y <- scale * y
    if (scale == -1)
    { 
      tvar.name <- paste("neg.",tvar.name,sep="")
      ttitle <- paste("Negative of ", ttitle,sep="")
    } else
    { 
      tvar.name <- paste("scaled.",tvar.name,sep="")
      ttitle <- paste("Scaled ", ttitle,sep="")
    }
  }
  #Apply offset
  if (offset != 0)
  { 
    y <- offset + y
    tvar.name <- paste("offset.",tvar.name,sep="")
    if (scale != 1)
      ttitle <- paste("and ", ttitle,sep="")
    ttitle <- paste("Offset  ", ttitle,sep="")
  }
  #Apply power transformation
  if (power != 1)
  { 
    if (power == 0)
    { 
      tvar.name <- paste("log.",tvar.name,sep="")
      if (min(y, na.rm = "TRUE") < 1e-04)
        stop("Negative values for log transform - could use offset/scale")
      else
        y <- log(y)
      ttitle <- paste("Logarithm of ", ttitle,sep="")
    } else
    { 
      y <- y^power
      if (power == 0.5)
      { 
        tvar.name <- paste("sqrt.",tvar.name,sep="")
        ttitle <- paste("Square root of ", ttitle,sep="")
      } else
      { 
        if (power == -1)
        { 
          tvar.name <- paste("recip.",tvar.name,sep="")
          ttitle <- paste("Reciprocal of ", ttitle,sep="")
        } else
        { 
          tvar.name <- paste("power.",tvar.name,sep="")
          ttitle <- paste("Power ",power," of ", ttitle, sep="")
        }
      }
    }

    #Add transformed variable to data
    tk <- match(tvar.name, names(data))
    if (is.na(tk))
      tk <- ncol(data) + 1
    data[[tk]]  <- y
    names(data)[tk] <- tvar.name
    #Add transformed title to titles
    nt <- length(titles)
    titles[[nt+1]] <- ttitle
    names(titles)[nt+1] <- tvar.name
  }
  return(list(data = data, tvar.name = tvar.name, titles = titles))
}

#Function to construct a vparameters table that can be supplied to asreml with new bounds and initial values
# - it is useful for having setvarianceterms.call use previous estimates of the variance parameters
#Need to decide if it is worth providing for users as a way of doing the equivalent of update = TRUE when 
# using setvarianceterms.call.
#IF it is to be implemented, should all singular terms be reset when terms is not NULL?
makeVpar.table.asreml <- function(asreml.obj, terms = NULL, 
                                  ignore.suffices = TRUE, bounds = "P", 
                                  initial.values = NA)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  if (!asr4)
    stop("makeVparameters.table.asreml is only available for asreml versions greater than 4.1")
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = FALSE)
  
  vpars.table <- unique(data.frame(Component = names(asreml.obj$vparameters), 
                                   Value = asreml.obj$vparameters, 
                                   Constraint = asreml.obj$vparameters.con))
  
  if (!asr4.2)
    vpars.table$Component <- vpc.char(asreml.obj)
  
  #Change singular terms to appropriate constraint and near-zero values
  # - might not need this, although might be best to reset in case changes to some terms affect other terms 
  if (any("S" %in% vpars.table$Constraint))
  {
    Sterms <- vpars.table$Component[vpars.table$Constraint == "S"]
    vpars.table$Constraint[vpars.table$Component %in% Sterms & 
                             asreml.obj$vparameters.type %in% c("V", "G", "P")] <- "P"
    vpars.table$Value[vpars.table$Component %in% Sterms & 
                        asreml.obj$vparameters.type %in% c("V", "G", "P")] <- 1e-04
    vpars.table$Constraint[vpars.table$Component %in% Sterms & 
                             asreml.obj$vparameters.type %in% c("R", "C", "L")] <- "U"
    vpars.table$Value[vpars.table$Component %in% Sterms & 
                        asreml.obj$vparameters.type %in% c("R", "C", "L")] <- 1e-04
  }
  
  if (!is.null(terms))
  { 
    #test for compatibility of arguments
    nt <- length(terms)
    if (length(ignore.suffices) == 1 & nt != 1)
      ignore.suffices <- rep(ignore.suffices, nt)
    if (length(ignore.suffices) != nt)
      stop("ignore.suffices specification is not consistent with terms")
    if (length(bounds) == 1 & nt != 1)
      bounds <- rep(bounds, nt)
    if (length(bounds) != nt)
      stop("bounds specification is not consistent with terms")
    if (length(initial.values) == 1 & nt != 1)
      initial.values <- rep(initial.values, nt)
    if (length(initial.values) != nt)
      stop("initial.values specification is not consistent with terms")
    if (!is.allnull(bounds) && any(!(bounds %in% c("B", "F", "P", "C", "U", "?"))))
      stop("the bound(s) ", 
           paste0(bounds[!(bounds %in% c("B", "F", "P", "C", "U", "?"))], collapse = ", "), 
           " is/are not used by asreml")
    
    #Change settings for terms  
    k <- unlist(lapply(1:nt, 
                       FUN=function(i, terms, termslist, ignore.suffices=TRUE)
                       { 
                         k <- which(termslist == terms[i])
                         return(k)
                       }, 
                       terms=terms,
                       termslist=vpars.table$Component, 
                       ignore.suffices=ignore.suffices))
    if (any(length(k)==0))
      stop(paste("Could not find", paste(terms[k==0], collapse=", ")))
    else
    { 
      if (!all(is.na(bounds)))
      { 
        kk <- k[!is.na(bounds)]
        vpars.table$Constraint[kk] <- bounds[!is.na(bounds)]
      }
      if (!all(is.na(initial.values)))
      { 
        kk <- k[!is.na(initial.values)]
        vpars.table$Value[kk] <- initial.values[!is.na(initial.values)]
      }
    }
  }
  return(vpars.table)
}

#Function that finds set.terms in a lsit of terms
# - allows for different factor/variable order
findsetterms <- function(setvparameters, termlist)
{
  sterms <- mapply(function(term, ignore.suffices, termlist) 
  {
    findterm(term, termlist, rmDescription = ignore.suffices)
  }, 
  term = setvparameters$set.terms, 
  ignore.suffices = setvparameters$ignore.suffices, 
  MoreArgs = list(termlist = termlist), SIMPLIFY = FALSE)
  sterms <- unlist(lapply(sterms, function(x) names(x)))
  return(sterms)
}

"setvarianceterms.call" <- function(call, terms, ignore.suffices = TRUE, bounds = "P", 
                                    initial.values = NA, ...)
  # call is an unevaluated call to asreml (can create using the call function)
  # terms is a vector of characters with the names of the gammas to be set
  # bounds specifies the bounds to be applied to the terms 
  #asreml codes:
  #  B - fixed at a boundary (!GP)
  #  F - fixed by user
  #  ? - liable to change from P to B    
  #  P - positive definite
  #  C - Constrained by user (!VCC)      
  #  U - unbounded
  #  S - Singular Information matrix
# initial.values specifies the initial values for the terms
# - ignore.suffices, bounds and initial.values must be of length 1 or 
#   the same length as terms
# - if any of bounds or initial.values is set to NA, then they are 
#   left unchanged for those terms
{ 
  #Deal with deprecated constraints parameter
  tempcall <- list(...)
  if (length(tempcall)) 
  { 
    if ("constraints" %in% names(tempcall))
      stop("constraints has been deprecated in setvarianceterms.call - use bounds")
    for (z in names(tempcall))
      languageEl(call, which = z) <- tempcall[[z]]
  }

  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  if (!inherits(call, "call"))
    stop("Need to supply an argument of class call")
  
  #Check for setvparameters in call
  if (is.null(call$setvparameters))
    setvparameters <- NULL
  else
    setvparameters <- call$setvparameters

  #Identify any variance parameters whose bounds have been deliberately changed
  # asreml.obj <- eval(call, sys.parent())
  # delib.bound <- getTermsDeliberateBound (asreml.obj, update = TRUE,
  #                                         asr4 = asr4, asr4.2 = asr4.2)

  #test for compatibility of arguments
  nt <- length(terms)
  if (length(ignore.suffices) == 1 & nt != 1)
    ignore.suffices <- rep(ignore.suffices, nt)
  if (length(ignore.suffices) != nt)
    stop("ignore.suffices specification is not consistent with terms")
  if (length(bounds) == 1 & nt != 1)
    bounds <- rep(bounds, nt)
  if (length(bounds) != nt)
    stop("bounds specification is not consistent with terms")
  if (length(initial.values) == 1 & nt != 1)
    initial.values <- rep(initial.values, nt)
  if (length(initial.values) != nt)
    stop("initial.values specification is not consistent with terms")
  if (!is.allnull(bounds) && any(!(bounds %in% c("B", "F", "P", "C", "U", "?"))))
    stop("the bound(s) ", 
         paste0(bounds[!(bounds %in% c("B", "F", "P", "C", "U", "?"))], collapse = ", "), 
         " is/are not used by asreml")
  
  #add start.values to call so can get vparamaeters in the current call
  start.call <- call
  languageEl(start.call, which = "start.values") <- TRUE
  gamma.start <- eval(start.call, sys.parent())
  #Get the gammas in the current model
  if (asr4)
  {
    gamma.table <- gamma.start$vparameters.table
    gammas <- gamma.table$Component
  }
  else
  {
    gamma.table <- gamma.start$gammas.table
    gammas <- gamma.table$Gamma
  }
  #Check if setvparameters terms are still in fitted vparameters - if not, remove them 
  if (!is.allnull(setvparameters) && !is.allnull(findsetterms(setvparameters, gammas)))
   {
    setvpars.terms <- findsetterms(setvparameters, gammas) #get setvpar.terms in gammas
    if (length(setvpars.terms) > 0)
    {
      #Reduce setvparameters to set.terms not in terms and make set.terms equal to those found
      setvparameters <- setvparameters[setvparameters$set.terms %in% names(setvpars.terms),]
      setvparameters$set.terms <- setvpars.terms
    }
    else
      setvparameters <- NULL
  } else
    setvparameters <- NULL
  
  # Add the setvparameters terms not in terms to terms
  if (!is.null(setvparameters))
  {
    if (!is.null(terms))
    {
      setvpars.terms <- setvparameters$set.terms
      setvpars.terms <- setdiff(setvpars.terms, terms) #get setvpar.terms not in terms
      if (length(setvpars.terms) > 0)
      { 
        #Reduce setvparameters to set.terms not in terms
        setvparameters <- setvparameters[setvparameters$set.terms %in% setvpars.terms,]
        setvparameters <- rbind(setvparameters, 
                                data.frame(set.terms = terms,
                                           ignore.suffices = ignore.suffices,
                                           bounds = bounds,
                                           initial.values = initial.values))
      }
    }
    terms <- setvparameters$set.terms
    ignore.suffices <- setvparameters$ignore.suffices
    initial.values <- setvparameters$initial.values
    bounds <- setvparameters$bounds
  }
  nt <- length(terms)
    
  #Apply bounds to the gammas specified by set.terms
  call <- setGRparam.call(call = call, set.terms = terms, 
                          ignore.suffices = ignore.suffices, bounds = bounds, 
                          initial.values = initial.values, 
                          asr4 = asr4, asr4.2 = asr4.2)
  #Set vparameters in the call
  call[["setvparameters"]] <- data.frame(set.terms = terms, 
                                         ignore.suffices = ignore.suffices, 
                                         bounds = bounds, 
                                         initial.values = initial.values)
  
  #deal with args coming via ...
  if (length(tempcall)) 
  { 
    for (z in names(tempcall))
      languageEl(call, which = z) <- tempcall[[z]]
  }
  
  #Evaluate the call
  new.reml <- eval(call, sys.parent())
  new.reml$call$setvparameters <- data.frame(set.terms = terms, 
                                             ignore.suffices = ignore.suffices, 
                                             bounds = bounds, 
                                             initial.values = initial.values)
  if (!new.reml$converge || largeVparChange(new.reml))
    new.reml <- newfit(new.reml)
  
   #Check that setvpars bounds have stuck- if not, warn that they have not stuck
  if (!is.null(setvparameters))
  {
    vpars <- getVpars(new.reml, asr4.2)
    #Find the location of setvparameters$set.terms in vpars$vpc
    ndb <- nrow(setvparameters)
    k <- unlist(lapply(1:ndb, 
                       FUN=function(i, set.terms, termslist, ignore.suffices=TRUE)
                       { 
                         k <- findterm(set.terms[i], termslist, 
                                       rmDescription=ignore.suffices[i])
                         return(k)
                       }, 
                       set.terms=setvparameters$set.terms,
                       termslist=names(vpars$vpc), 
                       ignore.suffices=setvparameters$ignore.suffices))
    if (any(k==0))
    { 
      warning(paste("Could not find", paste(setvparameters$set.terms[k==0], collapse=", ")))
      setvparameters <- setvparameters[k != 0, ]
      k <- k[k != 0]
    }
    bounds <- vpars$vpc[k]
    
    #If bounds have changed from setvparameters$bounds, then output a warning
    if (any(bounds != setvparameters$bounds))
      warning("The bounds have changed for the following terms: ", 
              paste0(bounds[bounds != setvparameters$bounds], collapse = ", "))
  }
  
  invisible(new.reml)
}

get.atargs <- function(at.term, dd, always.levels = FALSE)
{
  kargs <- stringr::str_replace(at.term, "at\\(", "")
  #kargs <- strsplit(kargs, "\\:")[[1]][1]
  kargs <- stringr::str_sub(kargs, 1, stringr::str_length(kargs)-1)
  #  kargs <- substr(kargs, 1, nchar(kargs)-1)
  kargs <- stringr::str_split(kargs, ",", n = 2)[[1]]
  obj <- kargs[1]
  obj.levs <- levels(dd[[obj]])
  lvls <- stringr::str_trim(kargs[2])
  
  if (grepl("\\\\", lvls) || grepl("\\\\'", lvls) ||
      grepl("c\\(", lvls) || all(!is.na(suppressWarnings(as.numeric(lvls)))))
    lvls <- eval(parse(text = lvls))
  if (is.character(lvls))
  { 
    lvls <- stringr::str_replace_all(lvls, "\"", "")
    lvls <- stringr::str_trim(lvls)
  }
  
  if (is.numeric(lvls))
  {
    if (always.levels)
    {
      levs.indx <- lvls
      lvls <- obj.levs[levs.indx]
    } else
    {
      if (any(is.na(suppressWarnings(as.numeric(obj.levs)))))
        levs.indx <- lvls
      else
      {
        if (all(as.character(lvls) %in% obj.levs))
          levs.indx <- which(obj.levs %in% as.character(lvls))
        else
          levs.indx <- lvls
      }
    }
  } else
    levs.indx <- which(obj.levs %in% as.character(lvls))
  
  #Check that levs.indx is legal  
  if (length(levs.indx) != 0 && (min(levs.indx) < 0 || max(levs.indx) > length(obj.levs)))
    stop('at has numeric values that are more than the number of levels')
  
  return(list(obj = obj, lvls = lvls, obj.levs = obj.levs, levs.indx = levs.indx))
}

#Fit model with quotes around AMF_plus 
# NB must have quotes for character levels, 
# as from v4.2 must also have in testranfix term 
#   because wald.tab, vparameters and varcomp rownames now always have single quotes.
# However, if fit with a character level call$fixed/random has double quotes and 
#     and the terms.object in formulae$fixed/random has \".
#If fit model with numeric specifying the position in the levels list, 
#   both call$fixed/random and the terms.object in formulae$fixed/random have the numeric.
#   This latter fact means that when updating a formula, needs to give index if was given in the 
#   formula to be updated.


#The purpose of this function is to make sure that any new "at" term being changed in a formula update
# matches that in the model that is being updated.
# new is the term to be added/dropped
# old is the formula to be updated
#It assumes that the new "at" term has the actual levels, rather than an index 1:no.levels.
#
#Testing of formulae with at functions is in test42Selection
atLevelsMatch <- function(new, old, call, single.new.term = FALSE, always.levels = TRUE)
{
  new.ch <- deparse(new)
  new.ch <- paste0(stringr::str_trim(new.ch, side = "left"), collapse = "")
  if (any(grepl("at\\(", new.ch))) #only process if new involves an at
  {
    ##Split new into pieces based on whether they are separated by "+"s
    dd <- eval(languageEl(call, which = "data")) #needed for levels
    if (single.new.term)
      new.split <- unlist(strsplit(new.ch, "[-]")) #introduced on 1/11/2023
    else
      new.split <- unlist(strsplit(new.ch, "[-~/+]")) #removed "*" on 17/9/2022
    at.parts <- stringr::str_trim(new.split[unlist(lapply(new.split, grepl, 
                                                          pattern = "at"))])
    ##Look for stray end comma if str in formula
    if (any(sapply(new.split, function(x) grepl("str\\(", x))))
      at.parts <- sapply(at.parts, function(x) gsub(",$", "", x))
    names(at.parts) <- at.parts

    #remove beginning and ending parentheses for multiple random terms
    at.parts <- sapply(at.parts, function(part)
    {
      part.len <- stringr::str_length(part)
      if (stringr::str_count(part,"\\(") > stringr::str_count(part,"\\)"))
      { 
        if (stringr::str_locate(part, "\\(")[1] == 1) #"(" at start?
          stringr::str_sub(part,1,1) <- ""
        else
          stop("Parentheses do not balance in ",part)
      }
      if (stringr::str_count(part,"\\(") < stringr::str_count(part,"\\)")) 
      {
        if  (any(stringr::str_locate_all(part, "\\)")[[1]] == part.len)) #")" at end?
          stringr::str_sub(part,part.len,part.len) <- ""
        else
          stop("Parentheses do not balance in ",part)
      }
      #Also remove parentheses from around a single term
      part.len <- stringr::str_length(part)
      if (stringr::str_sub(part,1,1) ==  "(" && stringr::str_sub(part,part.len,part.len) ==  ")")
        part <- stringr::str_sub(part,2,part.len-1)
      return(part)
    })
    at.parts <- at.parts[grepl("at\\(", at.parts)]
    
    ##Find the sets of factors association with terms in old that involve one or more at functions
    at.old.terms <- getTerms.formula(old)
    at.old.terms <- at.old.terms[grepl("at\\(", at.old.terms)]
    #    at.old.terms.vars <- strsplit(at.old.terms, split = ":")
    at.old.terms.vars <- lapply(at.old.terms, fac.getinTerm)
    at.old.terms.vars <- lapply(at.old.terms.vars, 
                                function(vars) vars[!grepl("at\\(", vars)])
    names(at.old.terms.vars) <- at.old.terms
    
    #Loop over the pieces of new
    for (piece in at.parts)
    {
      #Find the sets of factors association with terms in new that involve one or more at functions
      at.new.terms <- stringr::str_trim(unlist(strsplit(piece[length(piece)], split = "\\+")))
      at.new.terms <- at.new.terms[grepl("at\\(", at.new.terms)]
      at.new.terms.vars <- lapply(at.new.terms, fac.getinTerm)
      at.new.terms.vars <- lapply(at.new.terms.vars, 
                                  function(vars) vars[!grepl("at\\(", vars)])
      names(at.new.terms.vars) <- at.new.terms
      
      #check if any new terms in this piece have the same non-at variables as an old term
      for (kterm in at.new.terms)
      {
        matches <- unlist(lapply(at.old.terms.vars, 
                                 function(old.term, knew.term, at.new.terms.vars){
                                   same <- setequal(at.new.terms.vars[[knew.term]], old.term)
                                   return(same)
                                 }, knew.term = kterm, at.new.terms.vars = at.new.terms.vars))
        if (any(matches)) #have old term(s) whose non-at vars match kterm; do the at variables match?
        {
          matches <- names(at.old.terms.vars)[matches]
          at.new.term <- fac.getinTerm(kterm)
          at.new.term <- at.new.term[grepl("at\\(", at.new.term)]
          for (kmatch in matches)
          {
            at.kmatch <- fac.getinTerm(kmatch)
            at.kmatch <- at.kmatch[grepl("at\\(", at.kmatch)]
            if (at.new.term == at.kmatch) break
            at.new.args <- get.atargs(at.new.term, dd, always.levels = always.levels)
            at.kmatch.args <- get.atargs(at.kmatch, dd)
            #if same var and both have levels, try matching levels
            if (at.new.args$obj == at.kmatch.args$obj && !(all(is.na(at.new.args$lvls)) || all(is.na(at.kmatch.args$lvls)))) 
            {
              if (all(at.new.args$lvls == at.kmatch.args$lvls))
              {
                new.ch <- gsub(at.new.term, at.kmatch, new.ch, fixed = TRUE) #now substitute the old at term
              } else
              {
                #Check if kmatch is a set of levels indexes & new at is a set of levels
                if (is.numeric(at.kmatch.args$lvls) && 
                    all(as.character(at.new.args$lvls) %in% at.new.args$obj.levs))
                {
                  #Check that new at levels are those corresponding to the kmatch levels indexes 
                  if (setequal(as.character(at.new.args$lvls), 
                               at.kmatch.args$obj.levs[at.kmatch.args$levs.indx]))
                    new.ch <- gsub(at.new.term, at.kmatch, new.ch, fixed = TRUE) #now substitute the old at term
                } else
                {
                  #Check if kmatch is a set of levels & new at is a set of levels indices
                  if (!always.levels && all(as.character(at.kmatch.args$lvls) %in% at.kmatch.args$obj.levs) && 
                      is.numeric(at.new.args$lvls))
                  {
                    #Check if new at indices correspond to the same levels as kmatch
                    if (setequal(at.new.args$obj.levs[at.new.args$levs.indx], at.kmatch.args$lvls))
                      new.ch <- gsub(at.new.term, at.kmatch, new.ch, fixed = TRUE) #now substitute the old at term
                  }
                }
              }
            }
          }
        }
      }
    }
    new <- as.formula(new.ch)
  }
  return(new)
}

"my.update.formula" <- function(old, new, call, keep.order = TRUE, ...) 
  #function to update a formula
{ 
  old <- as.formula(old)
  old <- update.formula(old, ~ .)
  env <- environment(old)
  if (missing(new)) new <- old
  new <- atLevelsMatch(new, old, call, single.new.term = FALSE)
  #update formula expands the formula using keep.order = FALSE (cannot be changed)
  if (old == formula(~1))
    out <- old
  else
  { 
    tmp <- update.formula(old, new, keep.order = keep.order) 
    out <- formula(terms.formula(tmp, keep.order = keep.order))
 }
  environment(out) <- env
  return(out)
}

setGRparam.call <- function(call, set.terms = NULL, ignore.suffices = TRUE, 
                            bounds = "P", initial.values = NA, asr4, asr4.2)
{
  nt <- length(set.terms)
  #add start.values to call and unconstrain the gammas specified by set.terms
  start.call <- call
  languageEl(start.call, which = "start.values") <- TRUE
  gamma.start <- eval(start.call, sys.parent())
  #get the vparameters in the current call
  if (asr4)
  {
    gamma.table <- gamma.start$vparameters.table
    gamma.table <- gamma.table[!grepl("<NotEstimated>", gamma.table$Component),]
    rownames(gamma.table) <- NULL
    gammas <- gamma.table$Component
  } else
  {
    gamma.table <- gamma.start$gammas.table
    gammas <- gamma.table$Gamma
  }
  #Check that there are no set.terms absent from the current fit
  k <- unlist(lapply(1:nt, 
                     FUN=function(i, set.terms, termslist, ignore.suffices=TRUE)
                     { 
                       k <- findterm(set.terms[i], termslist, 
                                     rmDescription=ignore.suffices[i])
                       return(k)
                     }, 
                     set.terms=set.terms,
                     termslist=gammas, 
                     ignore.suffices=ignore.suffices))
  if (any(k==0))
  { 
    warning(paste("Could not find", paste(set.terms[k==0], collapse=", ")))
    bounds <- bounds[k != 0]
    initial.values <- initial.values[k != 0]
    k <- k[k != 0]
  }
  #Apply bounds to set.terms
  if (!all(is.na(bounds)))
  { 
    kk <- k[!is.na(bounds)]
    gamma.table$Constraint[kk] <- bounds[!is.na(bounds)]
  }
  if (!all(is.na(initial.values)))
  { 
    kk <- k[!is.na(initial.values)]
    gamma.table$Value[kk] <- initial.values[!is.na(initial.values)]
  }
  #modify call to set variance parameters
  languageEl(call, which = "G.param") <- gamma.table
  languageEl(call, which = "R.param") <- gamma.table
  
  return(call)
}

"newfit.asreml" <- function(asreml.obj, fixed., random., sparse., 
                            residual., rcov., update = TRUE, trace = FALSE, 
                            allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                            keep.order = TRUE, 
                            set.terms = NULL, ignore.suffices = TRUE, 
                            bounds = "P", initial.values = NA, ...)
  #a function to refit an asreml model with modified model formula
  #using either update.asreml or a direct call to asreml
  #- the principal difference is that the latter does not enforce the 
  #  use of previous values of the variance parameters as initial values.
  #- ... is used to pass arguments to asreml.
  # set.terms is a vector of characters with the names of the gammas to be set
  # bounds specifies the bounds to be applied to the terms 
  # initial.values specifies the initial values for the terms
  # - ignore.suffices, bounds and initial.values must be of length 1 or 
  #   the same length as terms
  # - if any of bounds or initial.values is set to NA, then they are 
#   left unchanged for those terms
{ 
  #Deal with deprecated constraints parameter
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("constraints" %in% names(tempcall))
      stop("constraints has been deprecated in setvarianceterms.call - use bounds")
  
  #For asr4, need to set trace using asreml.options in here and as.asrtests
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  if (asr4 & !trace)
    asreml::asreml.options(trace = trace)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  
  #Check that have a valid object of class asreml
  validasr <- validAsreml(asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  #Check for fixed correlations in supplied asrtests.obj
  if (!isFixedCorrelOK.asreml(asreml.obj, allow.fixedcorrelation = allow.fixedcorrelation))
  {
    if (asr4)
      kresp <- asreml.obj$formulae$fixed[[2]]
    else
      kresp <- asreml.obj$fixed.formula[[2]]
    warning(paste("The estimated value of one or more correlations in the supplied asreml fit for", 
                  kresp, "is bound, singular or fixed and allow.fixedcorrelation is FALSE"))
  }
  
  #Check for setvparameters in call
  if (is.null(asreml.obj$call$setvparameters))
    setvparameters <- NULL
  else
    setvparameters <- asreml.obj$call$setvparameters
  
  #Identify any variance parameters whose bounds have been deliberatley changed
  # delib.bound <- getTermsDeliberateBound(asreml.obj, update = update,
  #                                        asr4 = asr4, asr4.2 = asr4.2)
  
  if (is.null(call <- asreml.obj$call) && 
      is.null(call <- attr(asreml.obj, "call"))) 
    stop("Need an object with call component or attribute")
  #Evaluate formulae in case they are stored in a user-supplied object
  if (!is.null(languageEl(call, which = "fixed")))
    languageEl(call, which = "fixed") <- eval(languageEl(call, which = "fixed"))
  if (!is.null(languageEl(call, which = "random")))
    languageEl(call, which = "random") <- eval(languageEl(call, which = "random"))
  if (!is.null(languageEl(call, which = "sparse")))
    languageEl(call, which = "sparse") <- eval(languageEl(call, which = "sparse"))
  if (asr4)
  {
    if (!is.null(languageEl(call, which = "residual")))
      languageEl(call, which = "residual") <- eval(languageEl(call, which = "residual"))
  } else
  {
    if (!is.null(languageEl(call, which = "rcov")))
      languageEl(call, which = "rcov") <- eval(languageEl(call, which = "rcov"))
  }
  
  #Make sure that call has current value of update
  languageEl(call, which = "update") <- update
  
  #If data argument in call to newfit.asreml, use it to replace data in call
  if (length(tempcall)) 
    if ("data" %in% names(tempcall))
    { 
      dat <- tempcall$data
      if (is.symbol(data))
        dat <- eval(data)
      call$data <- dat
    }
  
  #Now update formulae
  if (!missing(fixed.)) 
  {
    if (fixed. != ". ~ . ")
      languageEl(call, which = "fixed") <- 
        my.update.formula(as.formula(languageEl(call, which = "fixed")), 
                          fixed., call = call, keep.order = keep.order)
  }
  
  if (!missing(random.)) 
  {
    if (is.null(random.))
      languageEl(call, which = "random") <- NULL
    else
    { 
      languageEl(call, which = "random") <- 
        { 
          if (!is.null(languageEl(call, which = "random"))) 
          { 
            my.update.formula(as.formula(languageEl(call, which = "random")), 
                              random., call = call, keep.order = keep.order)
          }
          else 
            random.
        }
      if (call$random == formula(~1))
        call$random <- NULL
    }
  }
  if (!missing(sparse.)) 
    languageEl(call, which = "sparse") <- 
    { if (!is.null(languageEl(call, which = "sparse"))) 
      my.update.formula(as.formula(languageEl(call, which = "sparse")), 
                        sparse., call = call, keep.order = keep.order)
      else 
        sparse.
    }
  if (asr4)
  {
    if (!missing(rcov.)) 
      stop("ASReml verson 4 is loaded - residual rather than deprecated rcov should be set")
    if (!missing(residual.)) 
    {
      if (is.null(residual.))
        languageEl(call, which = "residual") <- NULL
      else
        languageEl(call, which = "residual") <- 
          { if (!is.null(languageEl(call, which = "residual"))) 
            my.update.formula(as.formula(languageEl(call, which = "residual")), 
                              residual., call = call, keep.order = keep.order)
            else 
              residual.
          }
    }
  } else
  {
    if (!missing(residual.)) 
      stop("ASReml verson 3 is loaded - rcov should be set, rather than residual which requires version 4")
    if (!missing(rcov.)) 
    {
      if (is.null(rcov.))
        languageEl(call, which = "rcov") <- NULL
      else
        languageEl(call, which = "rcov") <- 
          { if (!is.null(languageEl(call, which = "rcov"))) 
            my.update.formula(as.formula(languageEl(call, which = "rcov")), 
                              rcov., call = call, keep.order = keep.order)
            else 
              rcov.
          }
    }
  }
  
  if (is.null(set.terms) && is.null(setvparameters))
  { 
    if (update)
    {
      #If R.param and G.param set in asreml.obj, copy them
      if (!is.null(asreml.obj$R.param))
        languageEl(call, which = "R.param") <- asreml.obj$R.param
      if (!is.null(asreml.obj$G.param))
        languageEl(call, which = "G.param") <- asreml.obj$G.param
    } else
    {
      #If R.param and G.param already set, make them NULL
      if (!is.null(languageEl(call, which = "R.param")))
        languageEl(call, which = "R.param") <- NULL
      if (!is.null(languageEl(call, which = "G.param")))
        languageEl(call, which = "G.param") <- NULL
    }
  } else
  { 
    #set variance terms
    if (!is.null(set.terms))
    { 
      #test for compatibility of arguments
      nt <- length(set.terms)
      if (length(ignore.suffices) == 1 & nt != 1)
        ignore.suffices <- rep(ignore.suffices, nt)
      if (length(ignore.suffices) != nt)
        stop("ignore.suffices specification is not consistent with set.terms")
      if (length(bounds) == 1 & nt != 1)
        bounds <- rep(bounds, nt)
      if (length(bounds) != nt)
        stop("bounds specification is not consistent with set.terms")
      if (length(initial.values) == 1 & nt != 1)
        initial.values <- rep(initial.values, nt)
      if (length(initial.values) != nt)
        stop("initial.values specification is not consistent with set.terms")
    }
    
    #add start.values to call so can get vparamaeters in the current call
    start.call <- call
    languageEl(start.call, which = "start.values") <- TRUE
    gamma.start <- eval(start.call, sys.parent())
    if (asr4)
    {
      gamma.table <- gamma.start$vparameters.table
      gammas <- gamma.table$Component
    }
    else
    {
      gamma.table <- gamma.start$gammas.table
      gammas <- gamma.table$Gamma
    }
    
    #Check if setvparameters terms are still in fitted vparameters - if not, remove them 
    if (!is.allnull(setvparameters) && !is.allnull(findsetterms(setvparameters, gammas)))
    {
      setvpars.terms <- findsetterms(setvparameters, gammas) #get setvpar.terms in gammas
      if (length(setvpars.terms) > 0)
      {
        #Reduce setvparameters to set.terms not in terms and make set.terms equal to those found
        setvparameters <- setvparameters[setvparameters$set.terms %in% names(setvpars.terms),]
        setvparameters$set.terms <- setvpars.terms
      }
      else
        setvparameters <- NULL
    } else
      setvparameters <- NULL
    
    # Add the setvparameters terms not in set.terms to set.terms
    if (!is.allnull(setvparameters))
    {
      if (!is.allnull(set.terms))
      {
        setvpars.terms <- setvparameters$set.terms
        setvpars.terms <- setdiff(setvpars.terms, set.terms) #get setvpar.terms not in set.terms
        if (length(setvpars.terms) > 0)
        { 
          #Reduce setvparameters to the terms not in set.terms
          setvparameters <- setvparameters[setvparameters$set.terms %in% setvpars.terms,]
          setvparameters <- rbind(setvparameters, 
                                  data.frame(set.terms = set.terms,
                                             ignore.suffices = ignore.suffices,
                                             bounds = bounds,
                                             initial.values = initial.values))
        }
      }
      set.terms <- setvparameters$set.terms
      ignore.suffices <- setvparameters$ignore.suffices
      initial.values <- setvparameters$initial.values
      bounds <- setvparameters$bounds
    }
    
    nt <- length(set.terms)
    if (nt > 0)
    {   
      #Apply bounds to the gammas specified by set.terms
      call <- setGRparam.call(call = call, set.terms = set.terms, 
                              ignore.suffices = ignore.suffices, bounds = bounds, 
                              initial.values = initial.values, 
                              asr4 = asr4, asr4.2 = asr4.2)
      #Set vparameters in the call
      call[["setvparameters"]] <- data.frame(set.terms = set.terms, 
                                             ignore.suffices = ignore.suffices, 
                                             bounds = bounds, 
                                             initial.values = initial.values)
    }
  }

  #deal with args coming via ...
  tempcall <- list(...)
  if (length(tempcall)) 
  { 
    for (z in names(tempcall))
      languageEl(call, which = z) <- tempcall[[z]]
  }

  #Check whether formulae have been reduced to no terms and, if have, update the formulae
  if (!is.null(languageEl(call, which = "random")) && 
      length(attr(terms(
        my.update.formula(as.formula(languageEl(call, which = "random")), 
                          random., call = call, keep.order = keep.order)), 
        "factors")) == 0) 
    languageEl(call, which = "random") <- NULL
  if (!is.null(languageEl(call, which = "sparse")) && 
      length(attr(terms(
        my.update.formula(as.formula(languageEl(call, which = "sparse")), 
                          sparse., call = call, keep.order = keep.order)), 
        "factors")) == 0) 
    languageEl(call, which = "sparse") <- NULL
  if (asr4)
  {
    if (!is.null(languageEl(call, which = "residual")) && 
        length(attr(terms(
          my.update.formula(as.formula(languageEl(call, which = "residual")), 
                            residual., call = call, keep.order = keep.order)), 
          "factors")) == 0) 
      languageEl(call, which = "residual") <- NULL
  } else
  {
    if (!is.null(languageEl(call, which = "rcov")) && 
        length(attr(terms(
          my.update.formula(as.formula(languageEl(call, which = "rcov")), 
                            rcov., call = call, keep.order = keep.order)), 
          "factors")) == 0) 
      languageEl(call, which = "rcov") <- NULL
  }

  #Evaluate the call
  if (asr4 & keep.order)
    asreml::asreml.options(keep.order = TRUE)
  class(asreml.obj) <- c("asreml", "list") #needed for within to work
  asreml.new.obj <- tryCatchLog(
    eval(call, sys.parent()),
    error = function(e) {print("Analysis continued"); NULL}, 
    include.full.call.stack = FALSE, include.compact.call.stack = FALSE,
    finally = within(asreml.obj, converge <- FALSE))
  asreml.new.obj$call <- call
  class(asreml.obj) <- "asreml" #return to single class
  
  #Check if updating and all Residual variance parameters are either F, B or S, 
  # and not equal to 1 (as in spatial models with no nugget variance)
  if (update)
  {
    Ranterms <- names(asreml.new.obj$vparameters)
    Resterms <- grepl("\\!R$", names(asreml.new.obj$vparameters))
    Resterms <- Ranterms[Resterms]
    Ranterms <- setdiff(Ranterms, Resterms)
    
    #Check Resterms
    Resterms <- unique(unlist(lapply(Resterms, rmTermDescription)))
    if (asr4.2)
    { 
      resterms.chk <- unlist(lapply(Resterms, grepl, names(asreml.new.obj$vparameters), fixed = TRUE)) & 
        asreml.new.obj$vparameters.type == "V"
      if (all(asreml.new.obj$vparameters.con[resterms.chk] %in% c("F","B","S")) && 
          all(abs(asreml.new.obj$vparameters[resterms.chk] - 1) > 1e-05))
        stop("There is no unbound residual variance - perhaps try update = FALSE")
    } else
    { 
      resterms.chk <- unlist(lapply(Resterms, grepl, names(asreml.new.obj$vparameters), fixed = TRUE)) & 
        vpt.char(asreml.new.obj) == "V"
      if (all(vpc.char(asreml.new.obj)[resterms.chk] %in% c("F","B","S")) && 
          all(abs(asreml.new.obj$vparameters[resterms.chk] - 1) > 1e-05))
        stop("There is no unbound residual variance - perhaps try update = FALSE")
    }
    
    #See if any bound or singular random terms, and if fit can be improved 
    #by fitting with update = FALSE
    #Have boundary terms?
    boundinfo <-  findboundary.asreml(asreml.new.obj, asr4 = asr4, asr4.2 = asr4.2)
    allvcomp <- boundinfo$allvcomp
    bound.terms <- boundinfo$bound.terms
    vcomp <- allvcomp[bound.terms, ]
    nbound <- nrow(vcomp)
    if (nbound > 0 && !is.null(asreml.new.obj$call$random) && is.null(set.terms))
    {   
      #Try a newfit with update = FALSE to remove the bound terms
      #If R.param and G.param already set, make them NULL
      if (!is.null(languageEl(call, which = "R.param")))
        languageEl(call, which = "R.param") <- NULL
      if (!is.null(languageEl(call, which = "G.param")))
        languageEl(call, which = "G.param") <- NULL
      new.asreml.obj <- tryCatchLog(
        eval(call, sys.parent()),
        error = function(e) {print("Analysis continued"); NULL}, 
        include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
      new.asreml.obj$call <- call
      new.boundinfo <-  findboundary.asreml(new.asreml.obj, asr4 = asr4, asr4.2 = asr4.2)
      new.allvcomp <- new.boundinfo$allvcomp
      new.bound.terms <- new.boundinfo$bound.terms
      new.nbound <- sum(new.bound.terms)
      if (new.nbound == 0 || 
          (all(new.allvcomp[new.bound.terms, "bound"] %in% allvcomp[bound.terms, "bound"]) && 
           new.nbound < nbound)) 
        asreml.new.obj <- new.asreml.obj
#      if (all(allvcomp[!bound.terms, "bound"] %in% new.allvcomp[!new.bound.terms, "bound"]))
    }    
  }
  
  #If not converged or large change in the variance parameters, try more iterations
  if (!asreml.new.obj$converge || largeVparChange(asreml.new.obj, threshold = 0.75))
  { 
    call <- asreml.new.obj$call
    asreml.new.obj <- tryCatchLog(
      eval(call, sys.parent()),
      error = function(e) {print("Analysis continued"); NULL}, 
      include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
  }
  #    asreml.new.obj <- asreml::update.asreml(asreml.new.obj)
  
  #Make sure that mbf.env is set in asreml.new.obj
  if (!is.null(asreml.obj$mf) && !is.null(attr(asreml.obj$mf, which = "mbf.env")))
    attr(asreml.new.obj$mf, 
         which = "mbf.env") <- attr(asreml.obj$mf, which = "mbf.env")

  #If still not converged, issue warning
  if (!asreml.new.obj$converge)
  {
    warning(asreml.new.obj$last.message)
    if (!allow.unconverged)
      asreml.new.obj <- asreml.obj
  } else
  {
    #Check for fixed correlation
    if (!isFixedCorrelOK.asreml(asreml.new.obj, allow.fixedcorrelation = allow.fixedcorrelation))
    {
      warning("At least one correlation's estimated value is bound, singular or fixed")
      asreml.new.obj <- asreml.obj
    }
  }
  
  #Check that setvpars bounds have stuck- if not, try to reset
  if (!is.null(setvparameters))
  {
    vpars <- getVpars(asreml.new.obj, asr4.2)
    #Find the location of setvparameters$set.terms in vpars$vpc 
    ndb <- nrow(setvparameters)
    k <- unlist(lapply(1:ndb, 
                       FUN=function(i, set.terms, termslist, ignore.suffices=TRUE)
                       { 
                         k <- findterm(set.terms[i], termslist, 
                                       rmDescription=ignore.suffices[i])
                         return(k)
                       }, 
                       set.terms=setvparameters$set.terms,
                       termslist=names(vpars$vpc), 
                       ignore.suffices=setvparameters$ignore.suffices))
    if (any(k==0))
    { 
      warning(paste("Could not find", paste(setvparameters$set.terms[k==0], collapse=", ")))
      setvparameters <- setvparameters[k != 0, ]
      k <- k[k != 0]
    }
    if (length(k) > 0)
      bounds <- vpars$vpc[k]
    else 
      bounds <- NULL
    
    #If bounds have changed from setvparameters$bounds, then try to reset to setvparameters$bounds
    if (!is.allnull(bounds) && any(bounds != setvparameters$bounds))
    {
      setvparameters$set.terms <- names(bounds)
      call <- asreml.new.obj$call
      #add start.values to call and unconstrain the gammas specified by set.terms
      call <- setGRparam.call(call = call, set.terms = setvparameters$set.terms, 
                              ignore.suffices = setvparameters$ignore.suffices, 
                              bounds = setvparameters$bounds, 
                              initial.values = setvparameters$initial.values, 
                              asr4 = asr4, asr4.2 = asr4.2)
      asreml.new.obj <- tryCatchLog(
        eval(call, sys.parent()),
        error = function(e) {print("Analysis continued"); NULL}, 
        include.full.call.stack = FALSE, include.compact.call.stack = FALSE)
    }
  }
  
  #Reset trace to default on the way out
  if (asr4)
    asreml::asreml.options(trace = TRUE)
  invisible(asreml.new.obj)
}

"iterate.asrtests" <- function(asrtests.obj, denDF = "numeric", trace = FALSE, ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  if (asr4 & !trace)
    asreml::asreml.options(trace = trace)
  
  #Check that have a valid asrtests object
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #initialize
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  
  #Update the asreml.obj
  asreml.obj <- newfit(asreml.obj)
  #  call <- asreml.obj$call
  #  asreml.obj <- eval(call, sys.parent())
  #If not converged, issue warning
  if (!asreml.obj$converge)
    warning(asreml.obj$last.message)
  
  #Update wald.tab
  wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
  wald.tab <- chkWald(wald.tab)
  
  #Update asrtests.object
  results <- as.asrtests(asreml.obj = asreml.obj, 
                         wald.tab = wald.tab, 
                         test.summary = test.summary,
                         denDF = denDF, trace = trace, ...)
  
  #Reset trace to default on the way out
  if (asr4)
    asreml::asreml.options(trace = TRUE)
  
  invisible(results)
}

findboundary.asreml <- function(asreml.obj, asr4, asr4.2)
{
  rancomps <- getResidRandomTerms(asreml.obj, which = "ran", asr4.2 = asr4.2)$ran
  if (asr4)
  {
    if (asr4.2)
    { 
      allvcomp <- data.frame(bound = asreml.obj$vparameters.con,  
                             component = asreml.obj$vparameters, 
                             stringsAsFactors = FALSE)
    } else
    {
      allvcomp <- data.frame(bound = vpc.char(asreml.obj),  
                             component = asreml.obj$vparameters, 
                             stringsAsFactors = FALSE)
    }
    #Reduce to random terms only
    allvcomp <- allvcomp[rownames(allvcomp) %in% names(rancomps),]
    bound.terms <- allvcomp$bound %in% c("B", "S") 
    #                 | bound == "F"
    #bound.terms[grep("?", bound, fixed=TRUE)] <- TRUE
  } else
  {
    allvcomp <- data.frame(bound = names(asreml.obj$gammas.con), 
                           component = asreml.obj$gammas,
                           stringsAsFactors = FALSE)
    allvcomp <- allvcomp[rownames(allvcomp) %in% names(rancomps),]
    bound.terms <- allvcomp$bound == "Boundary" | allvcomp$bound == "Singular"
    #                 | bound == "Fixed"
    #bound.terms[grep("?", bound, fixed=TRUE)] <- TRUE
  }
  
  return(list(allvcomp = allvcomp, bound.terms = bound.terms)) 
}

"rmboundary.asrtests" <- function(asrtests.obj, checkboundaryonly = FALSE, 
                                  IClikelihood = "none", trace = FALSE, update = TRUE, 
                                  set.terms = NULL, ignore.suffices = TRUE, 
                                  bounds = "P", initial.values = NA, ...)
  #Removes any boundary (B) or singular (S) terms from the fit stored in asreml.obj, 
  #one by one from largest to smallest
  #For a list of bounds codes see setvarianceterms.call (or ASREML entry in Ecco)
{ 
  
  #Deal with deprecated constraints parameter
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("constraints" %in% names(tempcall))
      stop("constraints has been deprecated in setvarianceterms.call - use bounds")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = TRUE)
  
  #Check that a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check IClikelihood options
  options <- c("none", "REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  ic.NA <- data.frame(AIC = NA, BIC = NA)
  
  #check input arguments
  if (asr4)
    kresp <- asrtests.obj$asreml.obj$formulae$fixed[[2]]
  else
    kresp <- asrtests.obj$asreml.obj$fixed.formula[[2]]    
  
  #Initialize
  asreml.obj <- asrtests.obj$asreml.obj
  summary(asreml.obj)$varcomp
  reml <- asreml.obj$loglik
  test.summary <- asrtests.obj$test.summary
  wald.tab <- asrtests.obj$wald.tab
  
  #Look for bound terms in random
  num.unremoveable <- 0
  #Loop  until have removed all removable boundary terms in random
  repeat
  { 
    #Find boundary terms
    boundinfo <-  findboundary.asreml(asreml.obj, asr4 = asr4, asr4.2 = asr4.2)
    allvcomp <- boundinfo$allvcomp
    bound.terms <- boundinfo$bound.terms
    vcomp <- allvcomp[bound.terms, ]
    nbound <- nrow(vcomp)
    #No boundary terms, so get out
    if (nbound <= 0 || checkboundaryonly || is.null(asreml.obj$call$random)) break

    #Reduce bound terms to set of bound terms that can be removed
    k <- 1
    while (k <= nbound)
    { 
      #Check if bound term is in the random formula
      term <- rownames(vcomp)[k]
      call <- asreml.obj$call
      if (substr(term, 1, 2) != "R!")  term <- rmTermDescription(term)
      termform <- atLevelsMatch(new = as.formula(paste("~", term)), 
                                old = as.formula(languageEl(call, which = "random")), 
                                call = call, single.new.term = TRUE)
      ranterms <- getTerms.formula(call$random)
      ranforms <- lapply(ranterms, 
                         function(term) 
                           ran <- as.formula(paste("~", 
                                                   paste0(term, collapse = " + "))))

      #Are there any nonremovable terms? i.e. not an R term or not a recognizable G term
      if (!any(sapply(ranforms, 
                      function(rform, termform) 
                        setequal(fac.getinTerm(getTerms.formula(rform), asr4.2 = asr4.2),
                                 fac.getinTerm(gsub('\\\"', "'", as.character(termform)[2]), 
                                               asr4.2 = asr4.2)),
                      termform = termform)))
      { 
        vcomp <- vcomp[-k, ]
        k <- k - 1
        nbound <- nbound - 1
        #Store any unremoveable terms
        if (num.unremoveable == 0)
        { 
          terms.unremoveable <- term
          num.unremoveable <- num.unremoveable + 1
        } else
        { 
          if (!(term %in% terms.unremoveable))
          { 
            terms.unremoveable <- c(terms.unremoveable, term)
            num.unremoveable <- num.unremoveable + 1
          }
        }
      }
      k <- k + 1
    }
    #No removeable boundary terms, so get out
    if (nbound <= 0) break
    #If single term, set it as the term to remove
    if (nbound == 1)
    { 
      term <- rmTermDescription(rownames(vcomp)[1])
    } else #Choose a term to remove
    { 
      #Classify terms as involving random or not because it involves a covariate
      vcomp.vars <- lapply(rownames(vcomp)[1:nbound], getTermVars)
      vcomp.codes <- lapply(vcomp.vars, getVarsCodes, asreml.obj = asreml.obj)
      vcomp <- within(vcomp, {
        terms.random <- unlist(lapply(vcomp.codes, function(term){all(term == 5)}))
        varnos <- unlist(lapply(vcomp.codes, length))
      })
      #Identify terms of the same type that are next to be removed
      max.no.factor <- TRUE
      #Check for terms that involve the dev function
      this.type <- unlist(lapply(vcomp.codes, function(term){any(term == 4 | term == 9)}))
      if (any(this.type))
      { 
        vcomp <- subset(vcomp, this.type)
      } else #Check for random terms (factors only)
      {   
        this.type <- unlist(lapply(vcomp.codes, function(term){all(term == 5)}))
        if (any(this.type))
        { 
          vcomp <- subset(vcomp, this.type)
          max.no.factor <- FALSE
        } else #Check for terms with spl
        { 
          this.type <- unlist(lapply(vcomp.codes, function(term){any(term == 3 | term == 8)}))
          if (any(this.type))
            vcomp <- subset(vcomp, this.type)
          else #Check for terms with pol
          { 
            this.type <- unlist(lapply(vcomp.codes, function(term){any(term == 2 | term == 7)}))
            if (any(this.type))
              vcomp <- subset(vcomp, this.type)
            else #Check for terms with covariate or lin
            { 
              this.type <- unlist(lapply(vcomp.codes, function(term){any(term < 2 | term == 6)}))
              if (any(this.type))
                vcomp <- subset(vcomp, this.type)
            }
          }
        }
      }
      #Reduce to subset of terms of the type to be removed and choose a term  
      if (max.no.factor)
      { 
        #Get smallest value term amongst those with the most vars
        vcomp <- subset(vcomp, vcomp$varnos == max(vcomp$varnos))
        vcomp <- with(vcomp, vcomp[order(-as.numeric(component)),])
        term <- rmTermDescription(tail(rownames(vcomp), 1))
      } else #Get smallest value term amongst those with the least vars
      { 
        vcomp <- subset(vcomp, vcomp$varnos == min(vcomp$varnos))
        vcomp <- with(vcomp, vcomp[order(-(component)),])
        term <- rmTermDescription(tail(rownames(vcomp), 1))
      }
    }
    #Remove chosen term
    if (ic.lik != "none")
    { 
      ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik)
    } else
    {  
      #if want to include bound.exclusions then have to make it an arg of rmboundary.asrtests
      #, bound.exclusions = bound.exclusions) 
      ic <- ic.NA
    }
    test.summary <- addtoTestSummary(test.summary, terms = term, 
                                     DF = 1, denDF = NA, p = NA, 
                                     AIC = ic$AIC, BIC = ic$BIC, 
                                     action = "Boundary")
    mod.ran <- as.formula(paste("~ . - ", term, sep=""))
    # mod.ran <- atLevelsMatch(as.formula(paste("~ . - ", term)), as.formula(languageEl(call, which = "random")), call)
    asreml.obj <- newfit.asreml(asreml.obj, random. = mod.ran, trace = trace, 
                                update = update, set.terms = set.terms, 
                                ignore.suffices = ignore.suffices, 
                                bounds = bounds, 
                                initial.values = initial.values, ...)
  }
  #Output warning if there are unremovable bound terms
  if (num.unremoveable > 0)
  {
    warning("\nIn analysing ", kresp,
            ", cannot remove the following boundary/singular term(s): ", 
            paste(terms.unremoveable, collapse = "; "), "\n\n")
  } else
  {
    if (checkboundaryonly & nbound > 0)
      warning("In analysing ", kresp,
              ", the following term(s) were boundary/singular, ","
              but not removed because checkboundaryonly = TRUE:\n",
              paste(rownames(vcomp), collapse = "; "), "\n")
  }
  #Check for variance terms that have been fixed and issue warning
  if (any(grepl("^F", allvcomp$bound))) #check for F or Fixed
  { 
    fix.vcomp <- allvcomp[grepl("^F", allvcomp$bound), ]
    #Leave if equal to one
    if (!any(abs(fix.vcomp$component - 1) < 1e-06))
    {  
      fix.terms <- rownames(allvcomp)[grepl("^F", allvcomp$bound)]
      #See if newfit removes any fixed terms
      new.asreml.obj <- newfit(asreml.obj, update = FALSE)
      new.boundinfo <-  findboundary.asreml(new.asreml.obj, asr4 = asr4, asr4.2 = asr4.2)
      new.allvcomp <- new.boundinfo$allvcomp
      
      fix.vcomp <- new.allvcomp[rownames(new.allvcomp) %in% fix.terms, ]
      if (!any(grepl("^S", fix.vcomp$bound)) && any(!grepl("^F", fix.vcomp$bound)))
      {
        asreml.obj <- new.asreml.obj
        allvcomp <- new.allvcomp
      } 
      if (any(grepl("^F", allvcomp$bound)))
        warning("In analysing ", kresp,
                ", estimates of the following parameter(s) are fixed:\n",
                rownames(allvcomp)[grepl("^F", allvcomp$bound)])
    }
  }
  #check for change in log likelihood and if changed issue warning
  change <- abs(reml -asreml.obj$loglik)/reml*100
  if (!is.na(change) & change > 1)
  { 
    warning(paste("Removing boundary terms has changed the log likelihood by "),
            change,"%")
    wald.tab <- NULL #as.asrtests will recalulcate the wald tab
  }  
  
  results <- as.asrtests(asreml.obj = asreml.obj, 
                         wald.tab = wald.tab, 
                         test.summary = test.summary, ...)
  invisible(results)
}

"changeTerms.asrtests" <- function(asrtests.obj, 
                                   dropFixed = NULL, addFixed = NULL, 
                                   dropRandom = NULL,  addRandom = NULL, 
                                   newResidual = NULL, label = "Changed terms", 
                                   allow.unconverged = TRUE, allow.fixedcorrelation = TRUE, 
                                   checkboundaryonly = FALSE, 
                                   trace = FALSE, update = TRUE, denDF = "numeric", 
                                   set.terms = NULL, ignore.suffices = TRUE, 
                                   bounds = "P", initial.values = NA, 
                                   IClikelihood = "none", 
                                   bound.exclusions = c("F","B","S","C"),  
                                   ...)
#Adds or removes sets of terms from one or both of the fixed or random asreml models
{ 
  
  #Deal with deprecated constraints parameter
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("constraints" %in% names(tempcall))
      stop("constraints has been deprecated in setvarianceterms.call - use bounds")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check IClikelihood options
  options <- c("none", "REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  ic.NA <- data.frame(fixedDF = NA, varDF = NA, AIC = NA, BIC = NA)
  
  #check input arguments
  if (asr4)
    kresp <- asrtests.obj$asreml.obj$formulae$fixed[[2]]
  else
    kresp <- asrtests.obj$asreml.obj$fixed.formula[[2]]
  #Check for fixed correlations in supplied asrtests.obj
  if (!isFixedCorrelOK.asreml(asrtests.obj$asreml.obj, 
                              allow.fixedcorrelation = allow.fixedcorrelation))
    warning(paste("The estimated value of one or more correlations in the supplied asreml fit for", 
                  kresp,
                  "is bound, singular or fixed and allow.fixedcorrelation is FALSE"))
  all.terms <- c(dropFixed, addFixed, dropRandom, addRandom, newResidual,set.terms)
  if (is.allnull(all.terms))
    stop("In analysing ", kresp, ", must supply terms to be removed/added")
  if (any(substr(trimws(all.terms), 1, 1) == "~"))
    stop("In analysing ", kresp, ", a leading tilde (~) has been included")
  if (!is.character(all.terms))
    stop("In analysing ", kresp, ", must supply terms as character")

  #initialize
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  
  #Form formula with terms to be added and removed
  fix.form <- ". ~ . "
  if (!is.null(dropFixed) || !is.null(addFixed))
  {
    if (!is.null(dropFixed))
      fix.form <- paste(fix.form, " - (", dropFixed, ")", sep = "")
    if (!is.null(addFixed))
      fix.form <- paste(fix.form, " + (", addFixed, ")", sep = "")
    fix.form <- as.formula(fix.form)
  }
  
  nullran <- is.null(languageEl(asreml.obj$call, which = "random"))
  ran.form <-" ~ . "
  if (((is.null(dropRandom) && is.null(addRandom)) ||
       (!is.null(dropRandom) && is.null(addRandom))) && nullran)
  {
    ran.form <- NULL
    if (!is.null(dropRandom))
      warning("No random terms to drop")
  } else
  {
    if (!is.null(dropRandom))
    {
      if (nullran)
        warning("No random terms to drop")
      else
        ran.form <- paste(ran.form, " - (", dropRandom, ")", sep = "")
    }
    if (!is.null(addRandom))
    { 
      if (nullran)
        ran.form <- paste("~ (", addRandom, ")", sep = "")
      else
        ran.form <- paste(ran.form , " + (", addRandom, ")", sep = "")
     }
    ran.form <- as.formula(ran.form)
  }
  
  res.form <- " ~ "
  if (!is.null(newResidual))
    res.form <- as.formula(paste(res.form, newResidual, sep=""))
  else 
  {
    if (asr4)
    {
      if (!is.null(languageEl(asreml.obj$call, which = "residual")))
        res.form <- " ~ . "
      else
        res.form <- NULL
    } else
    {
      if (!is.null(languageEl(asreml.obj$call, which = "rcov")))
        res.form <- " ~ . "
      else
        res.form <- NULL
    }
  }
  
  if (!is.allnull(all.terms) && length(setdiff(all.terms, set.terms)) == 0) #set.terms only
    action <- "Fix terms"
  else
  {  
    #Check models to update by determining if they have any terms in the update formula
    forms <- stringr::str_trim(as.character(c(fix.form, ran.form, res.form), side = "both"))
    lastch <- stringr::str_sub(forms, start = nchar(forms))
    lastch <- !(lastch == "." | lastch == "~")
    if (is.null(res.form))
      lastch <- c(lastch, FALSE)
    if (is.null(ran.form))
      lastch <- c(lastch[1], FALSE, lastch[2])
    if (!any(lastch))
      action <- "No changes"
    else
    {
      if (sum(lastch) == 1)
        action <- paste("Changed ", c("fixed", "random","residual")[lastch], sep = "")
      else
        action <- paste("Changed ", 
                        paste(c("fixed", "random","residual")[lastch], collapse = ", "), 
                        sep = "")
    }
  }
  
  #Update the models
  if (action == "No changes")
  {
    if (ic.lik != "none")
      ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                         bound.exclusions = bound.exclusions)
    else
      ic <- ic.NA
    test.summary <- addtoTestSummary(test.summary, terms = label, 
                                     DF=ic$fixedDF, denDF = ic$varDF, 
                                     p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                     action = action)
  } else
  {
    if (action == "Fix terms") #set.terms only
    {
      asreml.obj <- newfit(asreml.obj)
      call <- asreml.obj$call
      asreml.new.obj <- setvarianceterms(call, 
                                         terms = set.terms, 
                                         ignore.suffices = ignore.suffices, 
                                         bounds = bounds, 
                                         initial.values = initial.values)
    } else #have changes to model terms
    { 
      asreml.obj <- newfit(asreml.obj)
      
      if (asr4)
        asreml.new.obj <- newfit.asreml(asreml.obj, 
                                        fixed. = fix.form, random. = ran.form, 
                                        residual. = res.form, 
                                        trace = trace, update = update, 
                                        allow.unconverged = TRUE, 
                                        allow.fixedcorrelation = TRUE, 
                                        checkboundaryonly = checkboundaryonly,
                                        set.terms = set.terms, 
                                        ignore.suffices = ignore.suffices, 
                                        bounds = bounds, 
                                        initial.values = initial.values, ...)
      else
        asreml.new.obj <- do.call(newfit.asreml,
                                  args = list(asreml.obj, 
                                              fixed. = fix.form, random. = ran.form, 
                                              rcov. = res.form, 
                                              trace = trace, update = update, 
                                              allow.unconverged = TRUE, 
                                              allow.fixedcorrelation = TRUE, 
                                              set.terms = set.terms, 
                                              ignore.suffices = ignore.suffices, 
                                              bounds = bounds, 
                                              initial.values = initial.values, ...))
    }
    
    #Update results, checking for convergence
    if (asreml.new.obj$converge || allow.unconverged)
    {
      #Check fixed correlation
      if (!isFixedCorrelOK.asreml(asreml.new.obj, 
                                  allow.fixedcorrelation = allow.fixedcorrelation))
      {
        action <- "Unchanged - fixed correlation"
        if (ic.lik != "none")
          ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                             bound.exclusions = bound.exclusions)
        else
          ic <- ic.NA
        test.summary <- addtoTestSummary(test.summary, terms = label, 
                                         DF=ic$fixedDF, denDF = ic$varDF, 
                                         p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                         action = action)
      } else
      {
        asreml.obj <- asreml.new.obj
        #Update wald.tab
        wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
        wald.tab <- chkWald(wald.tab)
        if (!asreml.obj$converge)
          action <- paste(action, " - old uncoverged", sep="")
        if (ic.lik != "none")
        {
          ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                             bound.exclusions = bound.exclusions)
          test.summary <- addtoTestSummary(test.summary, terms = label, 
                                           DF=ic$fixedDF, denDF = ic$varDF, 
                                           p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                           action = action)
        } else
          test.summary <- addtoTestSummary(test.summary, terms = label, DF=NA, denDF = NA, 
                                           p = NA, action = action)
        #Check for boundary terms
        temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.obj, wald.tab, test.summary, ...), 
                                         checkboundaryonly = checkboundaryonly, 
                                         IClikelihood = IClikelihood, 
                                         trace = trace, update = update, 
                                         set.terms = set.terms, 
                                         ignore.suffices = ignore.suffices, 
                                         bounds = bounds, 
                                         initial.values = initial.values, ...)
        if (nrow(temp.asrt$test.summary) > nrow(test.summary))
        {
          if (asr4)
            warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                    ", boundary terms removed")
          else
            warning("In analysing ",asreml.obj$fixed.formula[[2]],
                    ", boundary terms removed")
        }
        asreml.obj <- temp.asrt$asreml.obj
        test.summary <- temp.asrt$test.summary
        #Update wald.tab
        wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
        wald.tab <- chkWald(wald.tab)
      }
    } else #unconverged and not allowed
    {
      #Check if get convergence with any boundary terms removed
      temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.new.obj, wald.tab, test.summary, ...), 
                                       checkboundaryonly = checkboundaryonly, 
                                       IClikelihood = IClikelihood, 
                                       trace = trace, update = update, 
                                       set.terms = set.terms, 
                                       ignore.suffices = ignore.suffices, 
                                       bounds = bounds, 
                                       initial.values = initial.values, ...)
      if (nrow(temp.asrt$test.summary) > nrow(test.summary))
      {
        if (asr4)
          warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                  ", boundary terms removed")
        else
          warning("In analysing ",asreml.obj$fixed.formula[[2]],
                  ", boundary terms removed")
      }
      if (temp.asrt$asreml.obj$converge)
      {
        if (!isFixedCorrelOK.asreml(temp.asrt$asreml.obj, 
                                    allow.fixedcorrelation = allow.fixedcorrelation))
        {
          action <- "Unchanged - fixed correlation"
          if (ic.lik != "none")
            ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                               bound.exclusions = bound.exclusions)
          else
            ic <- ic.NA
          test.summary <- addtoTestSummary(test.summary, terms = label, 
                                           DF=ic$fixedDF, denDF = ic$varDF, 
                                           p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                           action = action)
        } else
        { 
          asreml.obj <- temp.asrt$asreml.obj
          test.summary <- temp.asrt$test.summary
          #Update wald.tab
          wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
          wald.tab <- chkWald(wald.tab)
        }
      } else
      {
        p <- NA
        action <- "Unchanged - new unconverged"
      }
      if (ic.lik != "none")
        ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                           bound.exclusions = bound.exclusions)
      else
        ic <- ic.NA
      test.summary <- addtoTestSummary(test.summary, terms = label, 
                                       DF=ic$fixedDF, denDF = ic$varDF, 
                                       p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                       action = action)
    }
  }
  results <- as.asrtests(asreml.obj = asreml.obj, 
                         wald.tab = wald.tab, 
                         test.summary = test.summary,
                         denDF = denDF, trace = trace, ...)
  invisible(results)
}

"update.asrtests" <- function(asrtests.obj, label,  
                              allow.unconverged = TRUE, allow.fixedcorrelation = TRUE, 
                              checkboundaryonly = FALSE, 
                              trace = FALSE, update = TRUE, denDF = "numeric", 
                              #bounds = "P", initial.values = NA, 
                              bound.exclusions = c("F","B","S","C"),  
                              IClikelihood = "none", 
                              ...)
  #calls newfit.asreml to update the asreml.obj stored in asrtests.obj
{ 
  
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check whether any update.asreml args have been supplied
  tempcall <- list(...)
  upd.asr.args <- FALSE
  if (length(tempcall)) 
  {
    inargs <- names(tempcall)
    args <- unique(unlist(lapply(c("asreml", "update.asreml"), 
                                 function(func) 
                                   formalArgs(getFunction(func, 
                                                          where = rlang::ns_env("asreml"))))))
    if ("..." %in% args)
      args <- args[-match("...", args)]
    asremlArgs <- inargs[inargs %in% args]
    if (!length(asremlArgs))
      warning("No asreml or update.asreml args have been supplied to update.asrtests")
    else
      upd.asr.args <- TRUE
  }
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  
  #Check IClikelihood options
  options <- c("none", "REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  ic.NA <- data.frame(fixedDF = NA, varDF = NA, AIC = NA, BIC = NA)
  
  #initialize
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  
  #Update the asreml.obj
  asreml.new.obj <- newfit(asreml.obj, ...)
  
  #Update the asreml.obj
  #if (abs(asreml.new.obj$loglik - asreml.obj) <- 1.5e-08)  else
  if (!upd.asr.args)
  {
    action <- "No changes"
    if (ic.lik != "none")
      ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                         bound.exclusions = bound.exclusions)
    else
      ic <- ic.NA
    test.summary <- addtoTestSummary(test.summary, terms = label, 
                                     DF=ic$fixedDF, denDF = ic$varDF, 
                                     p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                     action = action)
  } else  #have changes to asreml.obj
  {
    #Update results, checking for convergence
    action <- "Changed"
    if (asreml.new.obj$converge || allow.unconverged)
    {
      #Check fixed correlation
      if (!isFixedCorrelOK.asreml(asreml.new.obj, allow.fixedcorrelation = allow.fixedcorrelation))
      {
        action <- "Unchanged - fixed correlation"
        if (ic.lik != "none")
          ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                             bound.exclusions = bound.exclusions)
        else
          ic <- ic.NA
        test.summary <- addtoTestSummary(test.summary, terms = label, 
                                         DF=ic$fixedDF, denDF = ic$varDF, 
                                         p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                         action = action)
      }
      else
      {
        asreml.obj <- asreml.new.obj
        #Update wald.tab
        wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
        wald.tab <- chkWald(wald.tab)
        if (!asreml.obj$converge)
          action <- paste(action, " - old uncoverged", sep="")
        if (ic.lik != "none")
        {
          ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                             bound.exclusions = bound.exclusions)
          test.summary <- addtoTestSummary(test.summary, terms = label, 
                                           DF=ic$fixedDF, denDF = ic$varDF, 
                                           p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                           action = action)
        } else
          test.summary <- addtoTestSummary(test.summary, terms = label, DF=NA, denDF = NA, 
                                           p = NA, action = action)
        #Check for boundary terms
        temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.obj, wald.tab, test.summary, 
                                                     IClikelihood = IClikelihood), 
                                         checkboundaryonly = checkboundaryonly, 
                                         IClikelihood = IClikelihood, 
                                         trace = trace, update = update, 
                                         ...)
        if (nrow(temp.asrt$test.summary) > nrow(test.summary))
        {
          if (asr4)
            warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                    ", boundary terms removed")
          else
            warning("In analysing ",asreml.obj$fixed.formula[[2]],
                    ", boundary terms removed")
        }
        asreml.obj <- temp.asrt$asreml.obj
        test.summary <- temp.asrt$test.summary
        #Update wald.tab
        wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
        wald.tab <- chkWald(wald.tab)
      }
    } else #unconverged and not allowed
    {
      #Check if get convergence with any boundary terms removed
      temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.new.obj, wald.tab, test.summary, ...), 
                                       checkboundaryonly = checkboundaryonly, 
                                       IClikelihood = IClikelihood, 
                                       trace = trace, update = update, 
                                       ...)
      if (nrow(temp.asrt$test.summary) > nrow(test.summary))
      {
        if (asr4)
          warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                  ", boundary terms removed")
        else
          warning("In analysing ",asreml.obj$fixed.formula[[2]],
                  ", boundary terms removed")
      }
      if (temp.asrt$asreml.obj$converge)
      {
        asreml.obj <- temp.asrt$asreml.obj
        test.summary <- temp.asrt$test.summary
        #Update wald.tab
        wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
        wald.tab <- chkWald(wald.tab)
      } else
      {
        p <- NA
        action <- "Unchanged - new unconverged"
      }
      if (ic.lik != "none")
        ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                           bound.exclusions = bound.exclusions)
      else
        ic <- ic.NA
      test.summary <- addtoTestSummary(test.summary, terms = label, 
                                       DF=ic$fixedDF, denDF = ic$varDF, 
                                       p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                       action = action)
    }
  }
  results <- as.asrtests(asreml.obj = asreml.obj, 
                         wald.tab = wald.tab, 
                         test.summary = test.summary,
                         denDF = denDF, trace = trace, ...)
  invisible(results)
}

"updateOnIC.asrtests" <- function(asrtests.obj, 
                                  label = "Changed terms", 
                                  allow.unconverged = TRUE, allow.fixedcorrelation = TRUE,
                                  checkboundaryonly = FALSE, 
                                  trace = FALSE, update = TRUE, denDF = "numeric", 
                                  which.IC = "AIC", IClikelihood = "REML", 
                                  fixedDF = NULL, varDF = NULL, #material.diff = NA, 
                                  bound.exclusions = c("F","B","S","C"),  
                                  ...)
  #Uses information criteria to select the best model after comparing that the model in the 
  #asreml.obj stored in asrtests.obj and the one obtained after updating this asreml.obj 
  #directly using asreml and newfit.asreml arguments.
  #The function update.asrtests is used to change the model.
{ 
  #Check IClikelihood options, because "none" is not allowed here
  options <- c("REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  options <- c("AIC", "BIC") #, "both")
  ic.type <- options[check.arg.values(which.IC, options)]
  
  asreml.obj <- asrtests.obj$asreml.obj
  
  #Check for fixed correlations in supplied asrtests.obj
  if (!isFixedCorrelOK.asreml(asreml.obj, allow.fixedcorrelation = allow.fixedcorrelation))
  {
    if ("formulae" %in% names(asreml.obj))
      kresp <- asreml.obj$formulae$fixed[[2]]
    else
    {
      if ("fixed.formula" %in% names(asreml.obj))
        kresp <- asreml.obj$fixed.formula[[2]]
      else
        kresp <- NULL
    }
    warning(paste("The estimated value of one or more correlations in the supplied asreml fit for", kresp,
                  "is fixed and allow.fixedcorrelation is FALSE"))
  }
  
  #Calculate the IC for the incoming fit
  old.IC <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                         bound.exclusions = bound.exclusions, 
                         fixedDF = fixedDF, varDF = varDF, ...)
  old.IC <- as.vector(old.IC[c("fixedDF", "varDF", "AIC", "BIC")], mode = "numeric")
  names(old.IC) <- c("DF", "denDF", "AIC", "BIC")
  nlines.test <- nrow(asrtests.obj$test.summary)
  
  #Check that update.asreml arguments have been supplied - otherwise do not change the the asreml.obj
  absent <- FALSE
  #Check whether any update.asreml args have been supplied
  tempcall <- list(...)
  upd.asr.args <- FALSE
  if (length(tempcall)) 
  {
    inargs <- names(tempcall)
    args <- unique(unlist(lapply(c("asreml", "update.asreml"), 
                                 function(func) 
                                   formalArgs(getFunction(func, 
                                                          where = rlang::ns_env("asreml"))))))
    if ("..." %in% args)
      args <- args[-match("...", args)]
    asremlArgs <- inargs[inargs %in% args]
    if (!length(asremlArgs))
      warning("No asreml or update.asreml args have been supplied to updateOnIC.asrtests")
    else
      upd.asr.args <- TRUE
  }
  
  if (!upd.asr.args)
  {
    ic.NA <- data.frame(fixedDF = NA, varDF = NA, AIC = NA, BIC = NA)
    action <- "Absent"
    ic <- ic.NA
    test.summary <- asrtests.obj$test.summary
    test.summary <- addtoTestSummary(test.summary, terms = label, 
                                     DF=NA, denDF = NA, p = NA, AIC = NA, BIC = NA, 
                                     action = "Unswapped")
  } else
  {
    #Use update.asrts to update the model in the asreml.obj
    new.asrtests.obj <- update(asrtests.obj, label = label,
                               allow.unconverged = allow.unconverged, 
                               allow.fixedcorrelation = allow.fixedcorrelation, 
                               checkboundaryonly = checkboundaryonly, 
                               trace = trace, update = update, denDF = denDF, 
                               IClikelihood = IClikelihood, 
                               bound.exclusions = bound.exclusions,  
                               ...)
    new.asreml.obj <- new.asrtests.obj$asreml.obj
    #Obtain IC for new model
    new.IC <- infoCriteria(new.asreml.obj, IClikelihood = ic.lik, 
                           bound.exclusions = bound.exclusions, 
                           fixedDF = fixedDF, varDF = varDF, ...)
    new.IC <- as.vector(new.IC[c("fixedDF", "varDF", "AIC", "BIC")], mode = "numeric")
    names(new.IC) <- c("DF", "denDF", "AIC", "BIC")
    
    #Extract asreml.objects
    if (!is.null(asreml.obj$mf) && !is.null(attr(asreml.obj$mf, which = "mbf.env")))
      attr(new.asreml.obj$mf, 
           which = "mbf.env") <- attr(asreml.obj$mf, which = "mbf.env")
    new.asreml.obj <- new.asreml.obj
    
    change <- FALSE
    action <- getTestEntry(new.asrtests.obj, label = label)$action 
    diff.IC <- new.IC - old.IC
    #Check fixed correlation
    if (grepl("Unchanged - fixed correlation", action))
    {
      change <- FALSE
    } else
    {
      if (grepl("Unchanged", action))
      {
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

"testranfix.asrtests" <- function(asrtests.obj, term=NULL, alpha = 0.05, 
                                  allow.unconverged = TRUE, allow.fixedcorrelation = TRUE, 
                                  checkboundaryonly = FALSE, 
                                  drop.ran.ns = TRUE, positive.zero = FALSE, 
                                  bound.test.parameters = "none", 
                                  bound.exclusions = c("F","B","S","C"), REMLDF = NULL, 
                                  drop.fix.ns = FALSE, denDF="numeric", dDF.fault = "none", 
                                  dDF.values = NULL, IClikelihood = "none", 
                                  trace = FALSE, update = TRUE, 
                                  set.terms = NULL, ignore.suffices = TRUE, 
                                  bounds = "P", initial.values = NA, ...)
  #function to test for a single term, using a REMLRT for a random term or based 
  #on Wald statistics for a fixed term. Note that fixed terms are, by default, not dropped.
{ 
  #Deal with deprecated constraints parameter
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("constraints" %in% names(tempcall))
      stop("constraints has been deprecated in setvarianceterms.call - use bounds")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check IClikelihood options
  options <- c("none", "REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  ic.NA <- data.frame(fixedDF = NA, varDF = NA, AIC = NA, BIC = NA)
  
  #Initialize
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  
  #Check for fixed correlations in supplied asrtests.obj
  if (!isFixedCorrelOK.asreml(asrtests.obj$asreml.obj, allow.fixedcorrelation = allow.fixedcorrelation))
  {
    if (asr4)
      kresp <- asrtests.obj$asreml.obj$formulae$fixed[[2]]
    else
      kresp <- asrtests.obj$asreml.obj$fixed.formula[[2]]
    warning(paste("The estimated value of one or more correlations in the supplied asreml fit for", kresp,
                  "is bound or fixed and allow.fixedcorrelation is FALSE"))
  }
  
  #Check for multiple terms
  term.form <- as.formula(paste("~ ",term, sep=""))
  term.obj <- as.terms.object(term, asreml.obj)
  if (length(labels(term.obj)) != 1)
  {
    if (asr4)
      stop("In analysing ",asreml.obj$formulae$fixed[[2]],
           ", multiple terms not allowed in testranfix.asrtests")
    else
      stop("In analysing ",asreml.obj$fixed.formula[[2]],
           ", multiple terms not allowed in testranfix.asrtests")
  } else
  {
    if (grepl("at(", term, fixed = TRUE)) #have an at term
    {
      at.facs <- rownames(attr(term.obj, which = "factors"))
      at.facs <- at.facs[grepl("at(", at.facs, fixed = TRUE)]
      lvls <- stringr::str_trim(stringr::str_split(at.facs, ", ", n = 2)[[1]])[2]
      lvls <- substr(lvls, 1, nchar(lvls)-1)
      if (substr(lvls, 1, 2) == "c(")
        lvls <- eval(parse(text = lvls))
      if (length(lvls) > 1)
      {
        if (asr4)
          stop("In analysing ",asreml.obj$formulae$fixed[[2]],
               ", an at term involving multiple levels will result in multiple terms and cannot be tested in testranfix.asrtests")
        
        else
          stop("In analysing ",asreml.obj$fixed.formula[[2]],
               ", an at term involving multiple levels will result in multiple terms and cannot be tested in testranfix.asrtests")
      }
    }
  }
  
  #Test whether term is in random model
  ranterms.obj <- as.terms.object(languageEl(asreml.obj$call, which="random"), asreml.obj)
  termno <- findterm(term, labels(ranterms.obj))
  if (termno == 0)
  { 
    #See if in fixed model by searching wald.tab - need to use wald.tab rather than the model
    #because will need the p-value from the wald.tab
    if (asr4)
      termno <- findterm(term, rownames(wald.tab))
    else
      termno <- findterm(term, rownames(wald.tab))
    #Term is not in either model
    if (termno == 0)
    {
      #Term is not in either model
      if (ic.lik != "none")
        ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                           bound.exclusions = bound.exclusions)
      else
        ic <- ic.NA
      test.summary <- addtoTestSummary(test.summary, terms = term, 
                                       DF=ic$fixedDF, denDF = ic$varDF, p = NA, 
                                       AIC = ic$AIC, BIC = ic$BIC, 
                                       action = "Absent")
    } else
      #Have a fixed term
    { 
      wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
      wald.tab <- chkWald(wald.tab)
      termno <- findterm(term, rownames(wald.tab)) #in case order has changed
      options <- c("none", "residual", "maximum", "supplied")
      opt <- options[check.arg.values(dDF.fault, options)]
      if (opt == "supplied" & is.null(dDF.values))
        stop('Need to set dDF.values because have set dDF.fault = \"supplied\"')
      #Compute p-value
      p <- wald.tab[termno, 
                    colnames(wald.tab)[grepl("Pr", colnames(wald.tab), fixed = TRUE)]]
      ndf <- wald.tab$Df[termno]
      den.df <- NA
      if ("denDF" %in% colnames(wald.tab) & !is.na(wald.tab$denDF[termno]))
        den.df <- wald.tab$denDF[termno]
      else
      { 
        if (opt == "supplied")
          den.df <- dDF.values[termno]
        else
        { 
          if (opt == "maximum") 
          { 
            if ("denDF" %in% colnames(wald.tab) & !all(is.na(wald.tab$denDF)))
              den.df <- max(wald.tab$denDF[-1], na.rm=TRUE) 
            else
              den.df <- asreml.obj$nedf
          } else
          { 
            if (opt == "residual")
              den.df <- asreml.obj$nedf
            else
              p <- NA
          }
        }
      }
      #Calc F, if necessary, and p
      if (!is.na(den.df))
      { 
        if ("denDF" %in% colnames(wald.tab))
        {
          if ("F.con" %in% colnames(wald.tab))
            test.stat <- wald.tab$F.con[termno]
          else
            test.stat <- wald.tab$F.inc[termno]
        }
        else
          test.stat <- wald.tab$'Wald statistic'[termno]/ndf
        p <- 1 - pf(test.stat, ndf, den.df)
      }
      
      #Add record for test to test.summary and, if drop.fix.ns is TRUE, remove term
      if (is.na(p))
      {
        if (ic.lik != "none")
          ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                             bound.exclusions = bound.exclusions)
        else
          ic <- ic.NA
        test.summary <- addtoTestSummary(test.summary, terms = rownames(wald.tab)[termno], 
                                         DF = ndf, denDF = NA, p = p, 
                                         AIC = ic$AIC, BIC = ic$BIC, 
                                         action = "Absent")
      } else
      {
        if (p <= alpha)
        {
          if (drop.fix.ns)
          {
            if (ic.lik != "none")
              ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                                 bound.exclusions = bound.exclusions)
            else
              ic <- ic.NA
            test.summary <- addtoTestSummary(test.summary, terms = rownames(wald.tab)[termno], 
                                             DF = ndf, denDF = den.df, p = p, 
                                             AIC = ic$AIC, BIC = ic$BIC, 
                                             action = "Retained")
          } else
          {
            if (ic.lik != "none")
              ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                                 bound.exclusions = bound.exclusions)
            else
              ic <- ic.NA
            test.summary <- addtoTestSummary(test.summary, terms = rownames(wald.tab)[termno], 
                                             DF = ndf, denDF = den.df, p = p, 
                                             AIC = ic$AIC, BIC = ic$BIC, 
                                             action = "Significant")
          }
        } else
        {
          if (drop.fix.ns)
          { 
            term.form <- as.formula(paste(". ~ . - ",term, sep=""))
            asreml.new.obj <- newfit.asreml(asreml.obj, fixed. = term.form, trace = trace, 
                                            update = update, 
                                            allow.unconverged = TRUE, 
                                            allow.fixedcorrelation = TRUE,
                                            set.terms = set.terms, 
                                            ignore.suffices = ignore.suffices, 
                                            bounds = bounds, 
                                            initial.values = initial.values, ...)
            if (asreml.new.obj$converge | allow.unconverged)
            {
              #Check for fixed correlation
              if (!isFixedCorrelOK.asreml(asreml.new.obj, allow.fixedcorrelation = allow.fixedcorrelation))
              {
                action <- "Unchanged - fixed correlation"
                asreml.new.obj <- asreml.obj
              } else
              {
                action <- "Dropped"
                asreml.obj <- asreml.new.obj
                #Update wald.tab
                wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, 
                                                 trace = trace, ...)
                wald.tab <- chkWald(wald.tab)
                if (!asreml.obj$converge)
                  action <- paste(action, " - unconverged", sep="")
                if (ic.lik != "none")
                  ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                                     bound.exclusions = bound.exclusions)
                else
                  ic <- ic.NA
                test.summary <- addtoTestSummary(test.summary, terms = term, 
                                                 DF = ndf, denDF = den.df, p = p, 
                                                 AIC = ic$AIC, BIC = ic$BIC, 
                                                 action = action)
                
                #Check for boundary terms
                temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.obj, wald.tab, 
                                                             test.summary, ...), 
                                                 checkboundaryonly = checkboundaryonly, 
                                                 IClikelihood = IClikelihood, 
                                                 trace = trace, update = update, 
                                                 set.terms = set.terms, 
                                                 ignore.suffices = ignore.suffices, 
                                                 bounds = bounds, 
                                                 initial.values = initial.values, ...)
                if (nrow(temp.asrt$test.summary) > nrow(test.summary))
                {
                  if (asr4)
                    warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                            ", Boundary terms removed")
                  else
                    warning("In analysing ",asreml.obj$fixed.formula[[2]],
                            ", Boundary terms removed")
                }
                asreml.obj <- temp.asrt$asreml.obj
                test.summary <- temp.asrt$test.summary
              }
            } else #unconverged and not allowed
            {
              #Check for boundary terms
              temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.new.obj, wald.tab, 
                                                           test.summary, ...), 
                                               checkboundaryonly = checkboundaryonly, 
                                               IClikelihood = IClikelihood, 
                                               trace = trace, update = update, 
                                               set.terms = set.terms, 
                                               ignore.suffices = ignore.suffices, 
                                               bounds = bounds, 
                                               initial.values = initial.values, ...)
              if (nrow(temp.asrt$test.summary) > nrow(test.summary))
              {
                if (asr4)
                  warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                          ", Boundary terms removed")
                else
                  warning("In analysing ",asreml.obj$fixed.formula[[2]],
                          ", Boundary terms removed")
              }
              if (temp.asrt$asreml.obj$converge) #have we now got convergence
              {
                asreml.obj <- temp.asrt$asreml.obj
                test.summary <- temp.asrt$test.summary
                action <- "Dropped"
              } else
              {
                p <- NA
                action = "Unchanged - unconverged"
              }
              if (ic.lik != "none")
                ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                                   bound.exclusions = bound.exclusions)
              else
                ic <- ic.NA
              test.summary <- addtoTestSummary(test.summary, terms = term, 
                                               DF = ndf, denDF = den.df, p = p, 
                                               AIC = ic$AIC, BIC = ic$BIC, 
                                               action = action)
            }
          } else
          {
            if (ic.lik != "none")
              ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                                 bound.exclusions = bound.exclusions)
            else
              ic <- ic.NA
            test.summary <- addtoTestSummary(test.summary, terms = rownames(wald.tab)[termno], 
                                             DF = ndf, denDF = den.df, p = p, 
                                             AIC = ic$AIC, BIC = ic$BIC, 
                                             action = "Nonsignificant")
          }
        }
      }
    }
  } else
    #Remove random term and test
  { 
    term.form <- as.formula(paste("~ . - ",term, sep=""))
    asreml.new.obj <- newfit.asreml(asreml.obj, random. = term.form, trace = trace, 
                                    update = update, 
                                    allow.unconverged = TRUE, allow.fixedcorrelation = TRUE, 
                                    set.terms = set.terms, 
                                    ignore.suffices = ignore.suffices, 
                                    bounds = bounds, 
                                    initial.values = initial.values, ...)
    
    
    
    
    
    if (!drop.ran.ns)
    {
      #Perform test
      test <- REMLRT.asreml(h1.asreml.obj = asreml.obj, h0.asreml.obj = asreml.new.obj, 
                            positive.zero = positive.zero, 
                            bound.exclusions = bound.exclusions, DF = REMLDF, 
                            bound.test.parameters = bound.test.parameters)
      if (!allow.unconverged && (!asreml.obj$converge | !asreml.new.obj$converge))
      {
        p <- NA
        if (!asreml.obj$converge)
        {
          if (!asreml.new.obj$converge)
          {
            action <- "Both unconverged"
          } else
          {
            action <- "Old unconverged"
          }
        } else
        {
          action <- "New unconverged"
        }
      } else
      { 
        #Check for fixed correlation
        if (!isFixedCorrelOK.asreml(asreml.new.obj, allow.fixedcorrelation = allow.fixedcorrelation))
        {
          p <- NA
          action <- "Fixed correlation"
        }
        else        
        {
          if (test$DF <= 0)
            p <- NA
          else
            p <- test$p
          if (is.na(p) || p <= alpha)
          {
            action <- "Significant"
          } else
          { 
            action <- "Nonsignificant"
          }
          if (!asreml.new.obj$converge)
            action <- paste(action, " - new unconverged", sep="")
        }
      }
    } else #drop.ran.ns
    {
      if (!allow.unconverged && !asreml.new.obj$converge)
      {
        #If new model not converged then see if removing boundary terms will result in convergence
        temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.new.obj, wald.tab, 
                                                     test.summary, ...), 
                                         checkboundaryonly = checkboundaryonly, 
                                         IClikelihood = IClikelihood, 
                                         trace = trace, update = update, 
                                         set.terms = set.terms, 
                                         ignore.suffices = ignore.suffices, 
                                         bounds = bounds, 
                                         initial.values = initial.values, ...)
        if (nrow(temp.asrt$test.summary) > nrow(test.summary))
        {
          if (asr4)
            warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                    ", Boundary terms removed")
          else
            warning("In analysing ",asreml.obj$fixed.formula[[2]],
                    ", Boundary terms removed")
        }
        if (temp.asrt$asreml.obj$converge)
        {
          asreml.obj <- temp.asrt$asreml.obj
          test.summary <- temp.asrt$test.summary
          term.form <- as.formula(paste("~ . + ",term, sep=""))
          asreml.new.obj <- newfit.asreml(asreml.obj, random. = term.form, trace = trace, 
                                          update = update, 
                                          allow.unconverged = TRUE, 
                                          allow.fixedcorrelation = TRUE,
                                          set.terms = set.terms, 
                                          ignore.suffices = ignore.suffices, 
                                          bounds = bounds, 
                                          initial.values = initial.values, ...)
        }
      }
      
      #Perform the test
      test <- REMLRT.asreml(h1.asreml.obj = asreml.obj, h0.asreml.obj = asreml.new.obj, 
                            positive.zero = positive.zero, 
                            bound.exclusions = bound.exclusions, DF = REMLDF, 
                            bound.test.parameters = bound.test.parameters)
      
      #Check convergence and force a converged model, if possible
      if (!allow.unconverged && (!asreml.obj$converge | !asreml.new.obj$converge))
      {
        p <- NA
        if (!asreml.new.obj$converge)
        {
          if (!asreml.obj$converge)
          {
            action <- "Retained - both unconverged"
            change <- FALSE
          } else
          {
            action <- "Retained - new unconverged"
            change <- FALSE
          }
        } else
        {
          action <- "Dropped - old unconverged"
          asreml.obj <- asreml.new.obj
          #Update wald.tab
          wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
          wald.tab <- chkWald(wald.tab)
        }
      } else #Evaluate test for drop.ran.ns
      {
        #Check for fixed correlation
        if (!isFixedCorrelOK.asreml(asreml.new.obj, allow.fixedcorrelation = allow.fixedcorrelation))
        {
          p <- NA
          action <- "Unchanged - fixed correlation"
        }
        else
        {
          if (test$DF <= 0)
            p <- NA
          else
            p <- test$p
          if (is.na(p) || p <= alpha)
          {
            if (drop.ran.ns)
            {
              action <- "Retained"
            } else
            {
              action <- "Significant"
            }
          } else
          { 
            action <- "Dropped"
            asreml.obj <- asreml.new.obj
            #Update wald.tab
            wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
            wald.tab <- chkWald(wald.tab)
          }
          if (!asreml.new.obj$converge)
            action <- paste(action, " - new unconverged", sep="")
        }
      }
    }
    
    #Update summary
    if (ic.lik != "none")
      ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                         bound.exclusions = bound.exclusions)
    else
      ic <- ic.NA
    test.summary <- addtoTestSummary(test.summary, terms = term, 
                                     DF=test$DF, denDF = NA, p = p, 
                                     AIC = ic$AIC, BIC = ic$BIC, 
                                     action = action)
    
    #Check for boundary terms
    temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.obj, wald.tab, 
                                                 test.summary, ...), 
                                     checkboundaryonly = checkboundaryonly, 
                                     trace = trace, update = update, 
                                     IClikelihood = IClikelihood, 
                                     set.terms = set.terms, 
                                     ignore.suffices = ignore.suffices, 
                                     bounds = bounds, 
                                     initial.values = initial.values, ...)
    if (nrow(temp.asrt$test.summary) > nrow(test.summary))
    {
      if (asr4)
        warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                ", Boundary terms removed")
      else
        warning("In analysing ",asreml.obj$fixed.formula[[2]],
                ", Boundary terms removed")
    }
    asreml.obj <- temp.asrt$asreml.obj
    test.summary <- temp.asrt$test.summary
  }
  
  results <- as.asrtests(asreml.obj = asreml.obj, 
                         wald.tab = wald.tab, 
                         test.summary = test.summary,
                         denDF = denDF, dDF.fault = dDF.fault, 
                         dDF.values = dDF.values, trace = trace, ...)
  invisible(results)
}

"testswapran.asrtests" <- function(asrtests.obj, oldterms = NULL, newterms = NULL, 
                                   label = "Swap in random model", simpler = FALSE, 
                                   alpha = 0.05, allow.unconverged = TRUE, 
                                   allow.fixedcorrelation = TRUE, checkboundaryonly = FALSE, 
                                   positive.zero = FALSE, bound.test.parameters = "none", 
                                   bound.exclusions = c("F","B","S","C"), REMLDF = NULL, 
                                   denDF="numeric", IClikelihood = "none", 
                                   trace = FALSE, update = TRUE, 
                                   set.terms = NULL, ignore.suffices = TRUE, 
                                   bounds = "P", initial.values = NA, ...)
  #function to test difference between current random model and one in which oldterms are dropped 
  #and newterms are added, using a REMLRT.
{ 
  #Deal with deprecated constraints parameter
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("constraints" %in% names(tempcall))
      stop("constraints has been deprecated in setvarianceterms.call - use bounds")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check IClikelihood options
  options <- c("none", "REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  ic.NA <- data.frame(AIC = NA, BIC = NA)
  
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  
  #Check for fixed correlations in supplied asrtests.obj
  if (!isFixedCorrelOK.asreml(asrtests.obj$asreml.obj, allow.fixedcorrelation = allow.fixedcorrelation))
  {
    if (asr4)
      kresp <- asrtests.obj$asreml.obj$formulae$fixed[[2]]
    else
      kresp <- asrtests.obj$asreml.obj$fixed.formula[[2]]
    warning(paste("The estimated value of one or more correlations in the supplied asreml fit for", kresp,
                  "is bound or fixed and allow.fixedcorrelation is FALSE"))
  }
  
  #Test whether oldterms are in random model
  oldterms.obj <- as.terms.object(oldterms, asreml.obj)
  ranterms.obj <- as.terms.object(languageEl(asreml.obj$call, which="random"), asreml.obj)
  termno <- findterm(oldterms, labels(ranterms.obj))
  if (any(lapply(labels(oldterms.obj), findterm, termlist=labels(ranterms.obj)) == 0))
  {
    if (asr4)
      stop("In analysing ",asrtests.obj$asreml.obj$formulae$fixed[[2]],
           ", some random terms in oldterms not in random model")
    else
      stop("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
           ", some random terms in oldterms not in random model")
  }
  
  #Remove random oldterms and add random newterms
  new.form <- as.formula(paste("~ . - ",oldterms," + ",newterms, sep=""))
  asreml.new.obj <- newfit.asreml(asreml.obj, random. = new.form, 
                                  trace = trace, update = update, 
                                  allow.unconverged = TRUE, 
                                  allow.fixedcorrelation = TRUE,
                                  set.terms = set.terms, 
                                  ignore.suffices = ignore.suffices, 
                                  bounds = bounds, 
                                  initial.values = initial.values, ...)
  change <- FALSE
  #Perform the test
  if (simpler)
  { 
    test <- REMLRT.asreml(h1.asreml.obj = asreml.obj, h0.asreml.obj = asreml.new.obj, 
                          positive.zero = positive.zero,
                          bound.exclusions = bound.exclusions, DF = REMLDF, 
                          bound.test.parameters = bound.test.parameters)
  } else
  {
    test <- REMLRT.asreml(h1.asreml.obj = asreml.new.obj, h0.asreml.obj = asreml.obj, 
                          positive.zero = positive.zero,
                          bound.exclusions = bound.exclusions, DF = REMLDF, 
                          bound.test.parameters = bound.test.parameters)
  }
  
  #check convergence
  if (!allow.unconverged && (!asreml.obj$converge | !asreml.new.obj$converge))
  {
    p <- NA
    if (!asreml.obj$converge)
    {
      if (!asreml.new.obj$converge)
      {
        action <- "Unchanged - both unconverged"
        change <- FALSE
      } else
      {
        action <- "Swapped - old unconverged"
        change <- TRUE
      }
    } else
    {
      action <- "Unchanged - new unconverged"
      change <- FALSE
    }
  } else
  {
    #Check fixed correlation
    if (!isFixedCorrelOK.asreml(asreml.new.obj, allow.fixedcorrelation = allow.fixedcorrelation))
    {
      p <- NA
      action <- "Unchanged - fixed correlation"
      change <- FALSE
    } else
    {
      #Evaluate the test
      if (simpler)
      { 
        if (test$DF <= 0)
          p <- NA
        else
          p <- test$p
        if (is.na(p) | p <= alpha)
          action <- "Unswapped"
        else
        { 
          action <- "Swapped"
          change <- TRUE
        }
      }
      else
      { 
        if (test$DF <= 0)
          p <- NA
        else
          p <- test$p
        if (!is.na(p) & p <= alpha)
        { 
          action = "Swapped"
          change <- TRUE
        }
        else
        { 
          action = "Rejected"
        }
      }
      #check convergence, when it is allowed
      if (allow.unconverged)
      {
        if (!asreml.obj$converge & !asreml.new.obj$converge)
        {
          action <- paste(action, " - both unconverged", sep="")
        } else
        {
          if (!asreml.obj$converge)
            action <- paste(action, " - old unconverged", sep="")
          else
          {
            if (!asreml.new.obj$converge)
              action <- paste(action, " - new unconverged", sep="")
          }
        }
      }
    }
  }
  if (ic.lik != "none")
    ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                       bound.exclusions = bound.exclusions)
  else
    ic <- ic.NA
  test.summary <- addtoTestSummary(test.summary, terms = label, 
                                   DF=test$DF, denDF = NA, p = p, 
                                   AIC = ic$AIC, BIC = ic$BIC, 
                                   action = action)
  
  #Update results
  if (change)
  { 
    #Check for boundary terms
    temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.new.obj, wald.tab, 
                                                 test.summary, ...), 
                                     checkboundaryonly = checkboundaryonly, 
                                     IClikelihood = IClikelihood, 
                                     trace = trace, update = update, 
                                     set.terms = set.terms, 
                                     ignore.suffices = ignore.suffices, 
                                     bounds = bounds, 
                                     initial.values = initial.values, ...)
    if (nrow(temp.asrt$test.summary) > nrow(test.summary))
    {
      if (asr4)
        warning("In swapping random terms for ",asreml.obj$formulae$fixed[[2]],
                ", Boundary terms removed")
      else
        warning("In swapping random terms for ",asreml.obj$fixed.formula[[2]],
                ", Boundary terms removed")
    }
    asreml.obj <- temp.asrt$asreml.obj
    test.summary <- temp.asrt$test.summary
    #Update wald.tab
    wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
    wald.tab <- chkWald(wald.tab)
  }
  results <- as.asrtests(asreml.obj = asreml.obj, 
                         wald.tab = wald.tab, 
                         test.summary = test.summary, 
                         denDF = denDF, trace = trace, ...)
  invisible(results)
}

"testresidual.asrtests" <- function(asrtests.obj, terms = NULL, label = "R model", 
                                    simpler = FALSE, alpha = 0.05, 
                                    allow.unconverged = TRUE, allow.fixedcorrelation = TRUE, 
                                    checkboundaryonly = FALSE, 
                                    positive.zero = FALSE, bound.test.parameters = "none", 
                                    bound.exclusions = c("F","B","S","C"), REMLDF = NULL, 
                                    denDF="numeric", IClikelihood = "none", 
                                    update = TRUE, trace = FALSE, 
                                    set.terms = NULL, ignore.suffices = TRUE, 
                                    bounds = "P", initial.values = NA, ...)
  #Fits new residual formula and tests whether the change is significant
{ 
  #Deal with deprecated constraints parameter
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("constraints" %in% names(tempcall))
      stop("constraints has been deprecated in testresidual.asreml - use bounds")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid asrtests object
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check IClikelihood options
  options <- c("none", "REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  ic.NA <- data.frame(AIC = NA, BIC = NA)
  
  #check input arguments
  if (asr4)
    kresp <- asrtests.obj$asreml.obj$formulae$fixed[[2]]
  else
    kresp <- asrtests.obj$asreml.obj$fixed.formula[[2]]    
  #Check for fixed correlations in supplied asrtests.obj
  if (!isFixedCorrelOK.asreml(asrtests.obj$asreml.obj, allow.fixedcorrelation = allow.fixedcorrelation))
    warning(paste("The estimated value of one or more correlations in the supplied asreml fit for", kresp,
                  "is bound or fixed and allow.fixedcorrelation is FALSE"))
  if (is.null(terms))
    stop("In analysing ", kresp, ", must supply terms to be tested")
  else
    if (!is.character(terms))
      stop("In analysing ", kresp, ", must supply terms as character")
  
  #initialize
  asreml.obj <- asrtests.obj$asreml.obj
  wald.tab <- asrtests.obj$wald.tab
  test.summary <- asrtests.obj$test.summary
  term.form <- as.formula(paste("~ ", terms, sep=""))
  #Update the R model
  if (asr4)
    asreml.new.obj <- newfit.asreml(asreml.obj, residual. = term.form, 
                                    trace = trace, update = update, 
                                    allow.unconverged = TRUE, 
                                    allow.fixedcorrelation = TRUE,
                                    set.terms = set.terms, 
                                    ignore.suffices = ignore.suffices, 
                                    bounds = bounds, 
                                    initial.values = initial.values, ...)
  else
    asreml.new.obj <- newfit.asreml(asreml.obj, rcov. = term.form, 
                                    trace = trace, update = update, 
                                    allow.unconverged = TRUE,
                                    allow.fixedcorrelation = TRUE,
                                    set.terms = set.terms, 
                                    ignore.suffices = ignore.suffices, 
                                    bounds = bounds, 
                                    initial.values = initial.values, ...)
  
  change <- FALSE
  #Perform the test
  if (simpler)
  { 
    test <- REMLRT.asreml(h1.asreml.obj = asreml.obj, h0.asreml.obj = asreml.new.obj, 
                          positive.zero = positive.zero, 
                          bound.exclusions = bound.exclusions, DF = REMLDF, 
                          bound.test.parameters = bound.test.parameters)
  } else
  {
    test <- REMLRT.asreml(h1.asreml.obj = asreml.new.obj, h0.asreml.obj = asreml.obj, 
                          positive.zero = positive.zero, 
                          bound.exclusions = bound.exclusions, DF = REMLDF, 
                          bound.test.parameters = bound.test.parameters)
  }
  
  #check convergence
  if (!allow.unconverged && (!asreml.obj$converge | !asreml.new.obj$converge))
  {
    p <- NA
    if (!asreml.obj$converge)
    {
      if (!asreml.new.obj$converge)
      {
        action <- "Unchanged - both unconverged"
        change <- FALSE
      } else
      {
        action <- "Swappped - old unconverged"
        change <- TRUE
      }
    } else
    {
      action <- "Unchanged - new unconverged"
      change <- FALSE
    }
  } else
  {
    #Check fixed correlation
    if (!isFixedCorrelOK.asreml(asreml.new.obj, allow.fixedcorrelation = allow.fixedcorrelation))
    {
      p <- NA
      action <- "Unchanged - fixed correlation"
      change <- FALSE
    } else
    {
      #Evaluate the test
      if (test$DF <= 0)
      {
        p <- NA
        action <- "Unswapped"
      } else
      {
        if (simpler)
        { 
          p <- test$p
          if (!is.na(p) & p <= alpha)
            action <- "Unswapped"
          else
          { 
            action = "Swapped"
            change <- TRUE
          }
        }
        else
        { 
          p <- test$p
          if (!is.na(p) & p <= alpha)
          { 
            action = "Swapped"
            change <- TRUE
          }
          else
            action = "Rejected"
        }
      }
      
      #check convergence, when it is allowed
      if (allow.unconverged)
      {
        if (!asreml.obj$converge && !asreml.new.obj$converge)
        {
          action <- paste(action, " - both unconverged", sep="")
        } else
        {
          if (!asreml.obj$converge)
            action <- paste(action, " - old unconverged", sep="")
          else
          {
            if (!asreml.new.obj$converge)
              action <- paste(action, " - new unconverged", sep="")
          }
        }
      }
    }
  }
  if (ic.lik != "none")
    ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                       bound.exclusions = bound.exclusions)
  else
    ic <- ic.NA
  test.summary <- addtoTestSummary(test.summary, terms = label, 
                                   DF=test$DF, denDF = NA, p = p, 
                                   AIC = ic$AIC, BIC = ic$BIC, 
                                   action = action)
  
  #Update results
  if (change)
  { 
    #Check for boundary terms
    temp.asrt <- rmboundary.asrtests(as.asrtests(asreml.new.obj, wald.tab, test.summary), 
                                     checkboundaryonly = checkboundaryonly, 
                                     IClikelihood = IClikelihood, 
                                     trace = trace, update = update, 
                                     set.terms = set.terms, 
                                     ignore.suffices = ignore.suffices, 
                                     bounds = bounds, 
                                     initial.values = initial.values, ...)
    if (nrow(temp.asrt$test.summary) > nrow(test.summary))
    {
      if (asr4)
        warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                ", boundary terms removed")
      else
        warning("In analysing ",asreml.obj$fixed.formula[[2]],
                ", boundary terms removed")
    }
    asreml.obj <- temp.asrt$asreml.obj
    test.summary <- temp.asrt$test.summary
    #Update wald.tab
    wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
    wald.tab <- chkWald(wald.tab)
  } 
  results <- as.asrtests(asreml.obj = asreml.obj, 
                         wald.tab = wald.tab, 
                         test.summary = test.summary,
                         denDF = denDF, trace = trace, ...)
  invisible(results)
}


"permute.square" <- function(x, permutation)
{ #function to permute the rows and coluns of a square matrix
  if (!is.matrix(x) || nrow(x) != ncol(x))
    stop("x must be a square matrix")
  permuted <- x[permutation, permutation]
  if (!is.null(rownames(x)))
    rownames(x) <- (rownames(x))[permutation]
  if (!is.null(colnames(x)))
    colnames(x) <- (colnames(x))[permutation]
  return(permuted)
}

"permute.to.zero.lowertri" <- function(x)
{ 
  #function to permute a square matrix until all the lower triangular elements are zero
  if (!is.matrix(x) || nrow(x) != ncol(x))
    stop("x must be a square matrix")
  m <- dim(x)[1]
  noperm <- m*(m-1)/2
  if (length(which(x == 0)) < noperm)
    stop("not enough zero elements to permute to all-zero, lower triangular form")
  x.lt <- matrix(0, nrow=m, ncol=m)
  x.lt[lower.tri(x.lt)] <- x [lower.tri(x)]
  i = 0
  #perform permutations until zero lower triangle 
  #or have done twice the no. of possible permutations
  while (any(x.lt == 1) && i <= noperm*2)
  { 
    #find the non-zero element in rightmost column (c) of the lowest row (r)
    nonzero <- which(x.lt == 1, arr.ind=TRUE)
    nonzero <- nonzero[nonzero[, 1] == max(nonzero[, 1]), ]
    if (!is.null(nrow(nonzero)))
      nonzero <- nonzero[nonzero[, 2] == max(nonzero[, 2]), ]
    #perform the permutation that in interchanges rows r and c and columns r and c
    permtn <- 1:m
    permtn[nonzero[1]] <- nonzero[2]
    permtn[nonzero[2]] <- nonzero[1]
    x <- permute.square(x, permtn)
    #get new lower triangle
    x.lt <- matrix(0, nrow=m, ncol=m)
    x.lt[lower.tri(x.lt)] <- x [lower.tri(x)]
  } 
  if (any(x.lt == 1))
    stop("Unable to form matrix with all zero lower triangular elements")
  return(x)
}

"reparamSigDevn.asrtests" <- function(asrtests.obj, terms = NULL, 
                                      trend.num = NULL, devn.fac = NULL, 
                                      allow.unconverged = TRUE, 
                                      allow.fixedcorrelation = TRUE,
                                      checkboundaryonly = FALSE, 
                                      denDF = "numeric", IClikelihood = "none", 
                                      trace = FALSE, update = TRUE, 
                                      set.terms = NULL, ignore.suffices = TRUE, 
                                      bounds = "P", initial.values = NA, ...)
  #reparamterizes a deviations term to a fixed term
  #It assumes that the deviations term are deviations from trend in trend.num and 
  #  that there is a random deviations term that involves devn.fac
  #Also assumes that the same term, but with the devn.fac substitued by trend.num
  #   and spl(trend.num), are in the fixed and random models respectively.
  #Retains the trend.num term if there are any significant terms involving it.
{ 
  #Deal with deprecated constraints parameter
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("constraints" %in% names(tempcall))
      stop("constraints has been deprecated in setvarianceterms.call - use bounds")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asrtests
  validasrt <- validAsrtests(asrtests.obj)  
  if (is.character(validasrt))
    stop(validasrt)
  
  #Check IClikelihood options
  options <- c("none", "REML", "full")
  ic.lik <- options[check.arg.values(IClikelihood, options)]
  ic.NA <- data.frame(fixedDF = NA, varDF = NA, AIC = NA, BIC = NA)
  
  #Check for fixed correlations in supplied asrtests.obj
  if (!isFixedCorrelOK.asreml(asrtests.obj$asreml.obj, allow.fixedcorrelation = allow.fixedcorrelation))
  {
    if (asr4)
      kresp <- asrtests.obj$asreml.obj$formulae$fixed[[2]]
    else
      kresp <- asrtests.obj$asreml.obj$fixed.formula[[2]]
    warning(paste("The estimated value of one or more correlations in the supplied asreml fit for", kresp,
                  "is bound or fixed and allow.fixedcorrelation is FALSE"))
  }
  
  #find out if any terms have either a devn.fac or a trend.num term 
  #  - (marginality implies cannot be both)
  if (is.null(terms))
    stop("The terms argument is NULL")
  ranterms.obj <- as.terms.object(languageEl(asrtests.obj$asreml.obj$call, which="random"), 
                                  asrtests.obj$asreml.obj)
  asrtests.old.obj <- asrtests.obj
  term.types <- trend.terms.types(terms = terms, devn.fac = devn.fac, trend.num = trend.num)
  #Check have terms involving devn.fac or trend.num
  if (term.types$has.devn.fac || term.types$has.trend.num)
    #Now reparamaterize deviations terms, the method depending on whether or not have a mixture
  { 
    for (term in terms)
    { 
      termno <- findterm(term, labels(ranterms.obj))
      if (termno == 0)
        stop(paste(term,"is not a random term"))
      factors <- fac.getinTerm(term)
      #Check this term involves a devn fac
      if (!is.null(devn.fac) && any(devn.fac %in%  factors))
      { 
        #Have mixture so convert random deviations to lin(trend.num) + fixed devn.fac deviations
        #form spl(time.num) and random deviation
        spl.term <- sub(devn.fac,paste("spl(",trend.num,")", sep=""),term)
        ran.term <- paste(spl.term, term, sep = " + " )
        if (term.types$has.devn.fac && term.types$has.trend.num)
        { 
          #Add time.num term to be sure - will do nothing if already there
          lin.term <- sub(devn.fac,trend.num,term)
          asrtests.obj <- changeTerms.asrtests(asrtests.obj, 
                                               addFixed = paste(lin.term, term, sep = " + " ), 
                                               dropRandom = ran.term, 
                                               trace = trace, allow.unconverged = TRUE, 
                                               allow.fixedcorrelation = TRUE,
                                               update = update, set.terms = set.terms, 
                                               ignore.suffices = ignore.suffices, 
                                               bounds = bounds, IClikelihood = IClikelihood, 
                                               initial.values = initial.values, ...)
        } else
        {
          #remove spl(time.num) and random deviation and add devn term to fixed model
          asrtests.obj <- changeTerms.asrtests(asrtests.obj, addFixed = term, 
                                               dropRandom = ran.term, 
                                               trace = trace, allow.unconverged = TRUE, 
                                               allow.fixedcorrelation = TRUE, 
                                               update = update, set.terms = set.terms, 
                                               ignore.suffices = ignore.suffices, 
                                               bounds = bounds, IClikelihood = IClikelihood, 
                                               initial.values = initial.values, ...)
          
        }
      }
    }
  }
  if (!allow.unconverged)
  {
    if (!asrtests.obj$asreml.obj$converge)
    {
      #Check for boundary terms
      temp.asrt <- rmboundary.asrtests(asrtests.obj, 
                                       checkboundaryonly = checkboundaryonly,  
                                       IClikelihood = IClikelihood, 
                                       trace = trace, update = update, 
                                       set.terms = set.terms, 
                                       ignore.suffices = ignore.suffices, 
                                       bounds = bounds, 
                                       initial.values = initial.values, ...)
      if (nrow(temp.asrt$test.summary) > nrow(asrtests.obj$test.summary))
      {
        if (asr4)
          warning("In analysing ",asrtests.obj$asreml.obj$formulae$fixed[[2]],
                  ", boundary terms removed")
        else
          warning("In analysing ",asrtests.obj$asreml.obj$fixed.formula[[2]],
                  ", boundary terms removed")
      }
      if (temp.asrt$asreml.obj$converge)
      {
        asrtests.obj <- temp.asrt
        asreml.obj <- asrtests.obj$asreml.obj
        #Update wald.tab
        wald.tab <- asreml::wald.asreml(asreml.obj, denDF = denDF, trace = trace, ...)
        wald.tab <- chkWald(wald.tab)
        asrtests.obj$wald.tab <- wald.tab
      } else
      {
        asrtests.obj <- asrtests.old.obj
      }
    }
  }
  #Check fixed correlation
  if (!isFixedCorrelOK.asreml(asrtests.obj$asreml.obj, allow.fixedcorrelation = allow.fixedcorrelation))
  {
    asrtests.obj <- asrtests.old.obj
    if (ic.lik != "none")
      ic <- infoCriteria(asreml.obj, IClikelihood = ic.lik, 
                         bound.exclusions = c("F","B","S","C"))
    else
      ic <- ic.NA
    test.summary <- addtoTestSummary(asrtests.old.obj$test.summary, terms = terms, 
                                     DF=ic$fixedDF, denDF = ic$varDF, 
                                     p = NA, AIC = ic$AIC, BIC = ic$BIC, 
                                     action = "Unchanged - fixed correlation")
    
    asrtests.obj$test.summary <- test.summary
  }
  else
    asrtests.obj$wald.tab <- recalcWaldTab(asrtests.obj, denDF = denDF, trace = trace, ...)
  invisible(asrtests.obj)
}

"predictASR4" <- function(asreml.obj, classify, levels = list(), 
                          sed = TRUE, vcov = FALSE, 
                          trace = FALSE, ...)
{
  if (sed & vcov)
  {
    pred <- predict(asreml.obj, classify=classify, levels=levels, 
                    sed = FALSE, vcov = TRUE, 
                    trace = trace, ...)
    if (!("vcov" %in% names(pred)))
      stop(paste0("predict.asreml has not returned the variance matrix of the predictions as requested\n",
                  "(possibly no estimable predicted values)"))
    vc <- as.matrix(pred$vcov)
    n <- nrow(vc)
    dvcov <- diag(vc)
    pred$sed <- matrix(rep(dvcov, each = n), nrow = n) + 
      matrix(rep(dvcov, times = n), nrow = n) - 2 * vc
    pred$sed <- sqrt(pred$sed)
    diag(pred$sed) <- NA_real_
  } else
  {
    pred <- predict(asreml.obj, classify=classify, levels=levels, 
                    sed = sed, vcov = vcov, 
                    trace = trace, ...)
    if (vcov && !("vcov" %in% names(pred)))
      stop(paste0("predict.asreml has not returned the variance matrix of the predictions as requested\n",
                  "(possibly no estimable predicted values)"))
  }
  
  if (sed && !("sed" %in% names(pred)))
    stop(paste0("predict.asreml has not returned the sed component for the predictions as requested\n",
                "(possibly no estimable predicted values)"))
  if (inherits(pred$pvals, "list"))
  {  
    attr.hd <- attr(pred$pvals, which = "heading")
    pred$pvals <- as.data.frame(pred$pvals)
    attr(pred$pvals, which = "heading") <- attr.hd
  }
  class(pred$pvals) <- c("predictions.frame", class(pred$pvals))
  return(pred)
}

"predictASR3" <- function(asreml.obj, classify, #levels = list(), 
                          sed = TRUE, vcov = FALSE, 
                          trace = FALSE, ...)
{
  pred <- predict(asreml.obj, classify=classify, #levels=levels, 
                  sed = sed, vcov = vcov, 
                  trace = trace, ...)$predictions
  class(pred$pvals) <- c("predictions.frame", "asreml.predict", "data.frame")
  return(pred)
}

"predictPlus.asreml" <- function(asreml.obj, classify, term = NULL, 
                                 inestimable.rm = TRUE, 
                                 linear.transformation = NULL, EGLS.linTransform = TRUE, 
                                 error.intervals = "Confidence", alpha = 0.05, 
                                 wald.tab = NULL, dDF.fault = "residual",  dDF.values = NULL, 
                                 pairwise = TRUE, Vmatrix = FALSE, 
                                 avsed.tolerance = 0.25, accuracy.threshold = NA, 
                                 LSDtype = "overall", LSDsupplied = NULL, LSDby = NULL, 
                                 LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
                                 x.num = NULL, x.fac = NULL,  
                                 x.pred.values = NULL, x.plot.values = NULL, 
                                 titles = NULL,  tables = "all" , level.length = NA, 
                                 transform.power = 1, offset = 0, scale = 1, 
                                 transform.function = "identity", 
                                 sortFactor = NULL, sortParallelToCombo = NULL, 
                                 sortNestingFactor = NULL, sortOrder = NULL, 
                                 decreasing = FALSE, trace = FALSE, ...)
#a function to get asreml predictions when there a parallel vector and factor are involved
{ 
  #Check for deprecated argument meanLSD.type and warn
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("meanLSD.type" %in% names(tempcall))
      stop("meanLSD.type has been deprecated - use LSDtype")
  
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asreml
  validasr <- validAsreml(asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction", "supplied")
  avLSD <- AvLSD.options[check.arg.values(LSDtype, AvLSD.options)]
  if (!is.null(LSDby) &&  !is.character(LSDby))
    stop("LSDby must be a character")
  
  #Get ... argument so can check for levels argument
  tempcall <- list(...)
  if ("levels.length" %in% names(tempcall))
    stop("levels.length has been deprecated - use level.length")
  if ("vcov" %in% names(tempcall))
    stop("Use Vmatrix to request that the vcov matrix be saved")
  
  #Check that have a valid wald.tab object
  validwald <- validWaldTab(wald.tab)  
  if (is.character(validwald))
    stop(validwald)
  if (!is.null(x.pred.values) && !is.null(x.plot.values))
    if (length(x.pred.values) != length(x.plot.values))
    {
      if (asr4)
        stop("In analysing ",asreml.obj$formulae$fixed[[2]],
             ", length of x.pred.values and x.plot.values should be equal")
      else
        stop("In analysing ",asreml.obj$fixed.formula[[2]],
             ", length of x.pred.values and x.plot.values should be equal")
    }
  if (!is.na(avsed.tolerance) & (avsed.tolerance <0 | avsed.tolerance > 1))
    stop("avsed.tolerance should be between 0 and 1")
  int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
  int.opt <- int.options[check.arg.values(error.intervals, int.options)]
  #Get table option and  check if must form pairwise differences
  tab.options <- c("none", "predictions", "vcov", "backtransforms", 
                   "differences", "p.differences", "sed", "LSD", "all")
  table.opt <- tab.options[unlist(lapply(tables, check.arg.values, options=tab.options))]
  if ("all" %in% table.opt || "differences" %in% table.opt || 
      int.opt == "halfLeastSignificant")
    pairwise <- TRUE
  if (classify == "(Intercept)")
  {
    pairwise <- FALSE
    get.vcov <- FALSE
    vars <- classify
  } else
  {
    #Make sure no functions in classify
    vars <- fac.getinTerm(classify, rmfunction=TRUE)
    classify <- fac.formTerm(vars)
    #Check if need vcov matrix
    if (!is.null(linear.transformation) || Vmatrix)
      get.vcov <- TRUE
    else
      get.vcov <- FALSE
  }
  
  #Get the predicted values when x.num is not involved in classify
  if (is.null(x.num) || !(x.num %in% vars))
  { 
    if (asr4)
      pred <- predictASR4(asreml.obj, classify=classify, 
                          sed=pairwise, vcov = get.vcov, 
                          trace = trace, ...)
    else
      pred <- predictASR3(asreml.obj, classify=classify, 
                          sed=pairwise, vcov = get.vcov,  
                          trace = trace, ...)
    if (!is.null(x.fac) && x.fac %in% vars)
    { 
      k <- match(x.fac, names(pred$pvals))
      #Set x values for plotting and table labels in x.fac
      if (is.numeric(pred$pvals[[k]]))
      { 
        if (!is.null(x.plot.values))
          pred$pvals[[k]] <- num.recode(pred$pvals[[k]], x.plot.values)
      }
      else
      { 
        if (is.character(pred$pvals[[k]]))
          pred$pvals[[k]] <- factor(pred$pvals[[k]])
        if (!is.null(x.plot.values))
        { 
          pred$pvals[[k]] <- fac.recode(pred$pvals[[k]], x.plot.values)
          pred$pvals[[k]] <- as.numfac(pred$pvals[[k]])
        }
        else
          if (all(is.na((as.numfac(pred$pvals[[k]])))))
            stop(paste("the levels of ",x.fac, "are not numeric in nature"))
        else
          pred$pvals[[k]] <- as.numfac(pred$pvals[[k]])
      }
    }
    else #make sure all variables are factors
    { 
      if (classify != "(Intercept)")
        pred$pvals[vars] <- lapply(1:length(vars), 
                                   function(k, vars, data){
                                     if (!is.factor(data[[vars[k]]]))
                                       data[vars[k]] <- factor(data[[vars[k]]], 
                                                               levels = unique(data[[vars[k]]]))
                                     return(data[[vars[k]]])}, 
                                   data=pred$pvals, vars=vars)
    }
  } else    #Get the predicted values when x.num is involved in classify
  { 
    #if levels in ... ignore x.pred.values
    if ("levels" %in% names(tempcall)) 
    {
      if (asr4)
        pred <- predictASR4(asreml.obj, classify=classify,
                            sed = pairwise, vcov = get.vcov, 
                            trace = trace, ...)
      else
        pred <- predictASR3(asreml.obj, classify=classify, 
                            sed = pairwise,  vcov = get.vcov, 
                            trace = trace, ...)
    } else
    {
      if (is.null(x.pred.values)) 
      {
        if (asr4)
          warning("In analysing ",asreml.obj$formulae$fixed[[2]],
                  ", predictions involve ",x.num,
                  " - you may want to specify x.pred.values or levels")
        else
          warning("In analysing ",asreml.obj$fixed.formula[[2]],
                  ", predictions involve ",x.num,
                  " - you may want to specify x.pred.values or levels")
      }
      x.list <- list(x.pred.values)
      names(x.list) <- x.num
      if (asr4)
        pred <- predictASR4(asreml.obj, classify=classify, levels=x.list, 
                            sed = pairwise, vcov = get.vcov, 
                            trace = trace, ...)
      else
        pred <- predictASR3(asreml.obj, classify=classify, levels=x.list, 
                            sed = pairwise, vcov = get.vcov, 
                            trace = trace, ...)
    }
    k <- match(x.num, names(pred$pvals))
    #Set x values for plotting and table labels in x.num
    if (!is.null(x.plot.values))
      pred$pvals[[k]] <- num.recode(pred$pvals[[k]], x.plot.values)
  }
  #Check that x.num and x.fac are conformable
  if ((!is.null(x.num) && x.num %in% names(pred$pvals)) && 
      (!is.null(x.fac) && x.fac %in% names(pred$pvals)))
  { 
    kn <- match(x.num, names(pred$pvals))
    kf <- match(x.fac, names(pred$pvals))  
    if (is.factor(pred$pvals[[kf]]))
      if (length(levels(pred$pvals[[kf]])) != length(unique(pred$pvals[[kn]])))
      {
        if (asr4)
          stop("In analysing ",asreml.obj$formulae$fixed[[2]],
               ", length of ",x.num," and the number of levels of ",x.fac," must be equal")
        else
          stop("In analysing ",asreml.obj$fixed.formula[[2]],
               ", length of ",x.num," and the number of levels of ",x.fac," must be equal")
      }
    else
      if (length(unique(pred$pvals[[kf]])) != length(unique(pred$pvals[[kn]])))
      {
        if (asr4)
          stop("In analysing ",asreml.obj$formulae$fixed[[2]],
               ", length of ",x.num," and the number of levels of ",x.fac," must be equal")
        else
          stop("In analysing ",asreml.obj$fixed.formula[[2]],
               ", length of ",x.num," and the number of levels of ",x.fac," must be equal")
      }
  }
  #Get denominator degrees of freedom
  if (is.null(term))
    term <- classify
  denom.df <- NA
  #Use wald.tab, if possible
  if (!is.null(wald.tab))
  {
    i <- findterm(term, rownames(wald.tab))
    if (i != 0 && "denDF" %in% colnames(wald.tab) && !is.na(wald.tab$denDF[i]))
      denom.df <- wald.tab$denDF[i]
    else  if (i == 0)
    { 
      {
        if (asr4)
          warning("For ",asreml.obj$formulae$fixed[[2]],
                  ", ",term," is not a fixed term that has been fitted")
        else
          warning("For ",asreml.obj$fixed.formula[[2]],
                  ", ",term," is not a fixed term that has been fitted")
      }
    }
  }
  if (is.na(denom.df))
  { 
    #Get options for dDF.fault
    options <- c("none", "residual", "maximum", "supplied")
    opt <- options[check.arg.values(dDF.fault, options)]
    if (opt == "supplied" && is.null(dDF.values))
      stop('Need to set dDF.values because have set dDF.fault = \"supplied\"')
    warn <- paste("Denominator degrees of freedom obtained using dDF.fault method", opt)
    if (is.null(wald.tab))
      warn <- paste(warn, "\n- no wald.tab supplied")
    warning(warn)
    #Compute denom.df
    denom.df <- NA
    if (opt == "supplied")
      denom.df <- dDF.values[i]
    else
    { 
      if (opt == "maximum") 
      { 
        if ("denDF" %in% colnames(wald.tab) & !all(is.na(wald.tab$denDF)))
          denom.df <- max(wald.tab$denDF[-1], na.rm=TRUE) 
        else
          denom.df <- asreml.obj$nedf
      }
      else
      { 
        if (opt == "residual")
          denom.df <- asreml.obj$nedf
      }
    }
  }
  
  #Initialize for setting up alldiffs object
  if (asr4)
    response <- as.character(asreml.obj$formulae$fixed[[2]])
  else
    response <- as.character(asreml.obj$fixed.formula[[2]])
  if (length(response) > 1)
    response <- paste(response, collapse = ".")
  if (!is.null(titles) && !is.na(match(response, names(titles))))
    response.title <- titles[[response]]
  else
  {
    response.title <- response
    if (!is.null(linear.transformation))
      response.title <- paste(response.title, "transform(s)", sep = " ")
  }
  
  #Form alldiffs object for predictions
  if (!get.vcov)
    lintrans.vcov <- NULL
  else
    lintrans.vcov <- pred$vcov
  diffs <- allDifferences(predictions = pred$pvals, vcov = lintrans.vcov, 
                          sed = pred$sed, LSDtype = LSDtype, LSDsupplied = LSDsupplied, 
                          LSDby = LSDby, LSDstatistic = LSDstatistic, LSDaccuracy = LSDaccuracy, 
                          response = response, response.title =  response.title, 
                          term = term, classify = classify, 
                          tdf = denom.df, 
                          x.num = x.num, x.fac = x.fac,
                          level.length = level.length, 
                          pairwise = pairwise,
                          transform.power = transform.power, 
                          transform.function = transform.function, 
                          offset = offset, scale = scale, 
                          inestimable.rm = inestimable.rm, 
                          alpha = alpha)
  if (is.null(linear.transformation))
  {
    #Add lower and upper uncertainty limits - send transform info to addBacktransforms.alldiffs 
    #so that backtransformed limits are updated
    diffs <- redoErrorIntervals.alldiffs(diffs, error.intervals = error.intervals, alpha = alpha, 
                                         avsed.tolerance = avsed.tolerance, 
                                         accuracy.threshold = accuracy.threshold, 
                                         LSDtype = LSDtype, LSDsupplied = LSDsupplied, LSDby = LSDby, 
                                         LSDstatistic = LSDstatistic, LSDaccuracy = LSDaccuracy) 
  } else
  {
    #Linear transformation required - send transform info to addBacktransforms.alldiffs
    diffs <- linTransform.alldiffs(alldiffs.obj = diffs, classify = classify,  term = term, 
                                   linear.transformation = linear.transformation, 
                                   EGLS.linTransform = EGLS.linTransform, 
                                   Vmatrix = Vmatrix, 
                                   error.intervals = error.intervals, 
                                   avsed.tolerance = avsed.tolerance, 
                                   accuracy.threshold = accuracy.threshold, 
                                   LSDtype = LSDtype, LSDsupplied = LSDsupplied, LSDby = LSDby, 
                                   LSDstatistic = LSDstatistic, LSDaccuracy = LSDaccuracy, 
                                   response = response, response.title = response.title, 
                                   x.num = x.num, x.fac = x.fac, 
                                   tables = "none", level.length = level.length, 
                                   pairwise = pairwise, alpha = alpha, 
                                   inestimable.rm = inestimable.rm)
  }
  
  #Sort alldiffs components, if required
  if (!is.null(sortFactor))
    diffs <- sort(diffs, decreasing = decreasing, sortFactor = sortFactor, 
                  sortParallelToCombo = sortParallelToCombo, 
                  sortNestingFactor = sortNestingFactor, sortOrder = sortOrder)
  
  #Outut tables according to table.opt and save alldiff object for current term
  if (!("none" %in% table.opt))
    print(diffs, which = table.opt)
  return(diffs)
}



"plotPredictions.data.frame" <- function(data, classify, y, x.num = NULL, x.fac = NULL, 
                                         nonx.fac.order = NULL, colour.scheme = "colour", 
                                         panels = "multiple", graphics.device = NULL,
                                         error.intervals = "Confidence", 
                                         interval.annotate = TRUE, 
                                         titles = NULL, y.title = NULL, 
                                         filestem = NULL, printPlot = TRUE, 
                                         ggplotFuncs = NULL, ...)
  #a function to plot asreml predictions and associated statistics
{ 
  
  #Get ... argument so can check for linear.transformation argument
  tempcall <- list(...)
  if ("linear.transformation" %in% names(tempcall))
    warning("linear.transformation is not an argument to plotPredictions - perhaps use linTransform.alldiffs")
  
  #Change asreml4 names to asreml3 names
  data <- as.predictions.frame(data, classify = classify, se = "std.error", est.status = "status")
  #Check that a valid object of class predictions.frame
  validpframe <- validPredictionsFrame(data)  
  if (is.character(validpframe))
    stop(validpframe)
  #check have some predictions
  if (nrow(data) == 0)
    stop("no rows in object supplied to data argument")
  
  #Check options
  scheme.options <- c("black", "colour")
  scheme.opt <- scheme.options[check.arg.values(colour.scheme, scheme.options)]
  panel.options <- c("single", "multiple")
  panel.opt <- panel.options[check.arg.values(panels, panel.options)]
  int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
  int.opt <- int.options[check.arg.values(error.intervals, int.options)]
  intervals <- !(int.opt == "none")
  #check that x.fac is numeric in nature when no x.num
  if (is.null(x.num))
    if (!is.null(x.fac) && all(is.na((as.numfac(data[[x.fac]])))))
      stop(paste("the levels of ",x.fac, "are not numeric in nature"))
  
  #work out columns for intervals and their names along with labels for them
  if (intervals)
  { 
    klow <- pmatch("lower", names(data))
    kupp <- pmatch("upper", names(data))
    if (klow == 0 || kupp == 0)
      stop(paste("Cannot find a column in data starting with lower ",
                 "and another with upper to plot intervals"))
    ylow <- names(data)[klow]
    yupp <- names(data)[kupp]
    low.parts <- (strsplit(ylow, ".", fixed=TRUE))[[1]]
    upp.parts <- (strsplit(yupp, ".", fixed=TRUE))[[1]]
    if (!setequal(low.parts[-1], upp.parts[-1]))
      stop("Names of columns for lower and upper limits are not consistent")
    if (low.parts[2] == "StandardError")
    { 
      labend <- "+/- standard errors"
      abbrev <- "SE"
    }
    else 
      if (low.parts[2] == "Confidence")
      { 
        labend <- "confidence intervals"    
        abbrev <- "CI"
      }
    else
      if (low.parts[2] == "halfLeastSignificant") 
      { 
        labend <- "+/- half mean LSD"
        abbrev <- "LSI"
      }
    else
      stop("Names for interval limit are not consistent with the types of interval allowed")
  }
  #Do plots
  #Make sure no functions in classify
  if (classify == "(Intercept)")
  {
    vars <- "Intercept"
    names(data)[match("(Intercept)", names(data))] <- "Intercept"
  } else
  {
    vars <- fac.getinTerm(classify, rmfunction = TRUE)
    classify <- fac.formTerm(vars)
  }
  non.x.terms <- setdiff(vars, c(x.fac, x.num))
  n.non.x <- length(non.x.terms)
  if ((length(vars) != n.non.x && n.non.x > 2) || (length(vars) == n.non.x && n.non.x > 3)) 
    stop("Sorry but plotting for prediction terms with more than 3 variables ", 
         "not implemented; \n", 
         "One possibility is to use facCombine.alldiffs to combine some factors.")
  #Determine plotting order of non.x.vars
  if (!is.null(nonx.fac.order)) 
  { 
    if (!setequal(nonx.fac.order, non.x.terms))
      stop("The factors nonx.fac.order are not the same as the set in the term being plotted")
    non.x.terms <- nonx.fac.order
    nos.lev <- sapply(1:n.non.x, 
                      FUN=function(k, non.x.terms, data){return(length(levels(data[[non.x.terms[k]]])))},
                      non.x.terms = non.x.terms, data = data)
  }
  else
    #Default order - sort non.x.vars according to their numbers of levels 
  { 
    nos.lev <- 1
    if(n.non.x >0)
    { nos.lev <- sapply(1:n.non.x, 
                        FUN=function(k, non.x.terms, data)
                        {return(length(levels(data[[non.x.terms[k]]])))},
                        non.x.terms = non.x.terms, data = data)
    non.x.terms <- non.x.terms[order(nos.lev, decreasing = TRUE)]
    nos.lev <- nos.lev[order(nos.lev, decreasing = TRUE)]
    }
  }
  #Append x.vars 
  if (length(vars) != n.non.x)
  { 
    x.terms <- setdiff(vars, non.x.terms)
    vars <- c(non.x.terms, x.terms)
  }
  else
    vars <- non.x.terms
  # meanLSD <- attr(data, which = "meanLSD")
  # if (!is.null(meanLSD) && !is.na(meanLSD))
  #   meanLSD <- sqrt(mean(meanLSD*meanLSD, na.rm = TRUE))
  # else
  #   meanLSD <- NA
  #Plot predicted values
  cbPalette <- rep(c("#CC79A7", "#56B4E9", "#009E73", "#E69F00", "#0072B2", "#D55E00", "#000000"), times=2)
  symb <- rep(c(18,17,15,3,13,8,21,9,3,2,11,1,7,5,10,0), times=10)
  if (!is.null(graphics.device) )
    do.call(graphics.device, list(record = FALSE))
  barplot <- FALSE
  if (length(vars) == length(non.x.terms))
    #Do bar chart
  { 
    barplot <- TRUE
    if (!is.null(titles) & !is.na(match(vars[1], names(titles))))
      x.title <- titles[[vars[1]]]
    else
      x.title <- vars[[1]]
    pred.plot <-  ggplot(data=data, ...) + theme_bw() + 
      scale_y_continuous(y.title)
    if (classify == "(Intercept)")
      pred.plot <- pred.plot + scale_x_discrete(x.title, labels = NULL)
    else
      pred.plot <- pred.plot + scale_x_discrete(x.title)
    if (scheme.opt == "black")
      pred.plot <-  pred.plot + theme_bw() +
      theme(panel.grid.major = element_line(colour = "grey95", 
                                            linewidth= 0.5), 
            panel.grid.minor = element_line(colour = "grey98", 
                                            linewidth = 0.5))
    if (n.non.x == 1)
    { 
      pred.plot <-  pred.plot + 
        aes(x = .data[[!!vars[1]]], y = .data[[!!y]]) + 
        scale_fill_manual(values=cbPalette) + 
        theme(axis.text.x=element_text(angle=90, hjust=1))
      if (scheme.opt == "black")
        pred.plot <-  pred.plot + geom_bar(stat="identity", fill="grey50") 
      else
        pred.plot <-  pred.plot + geom_bar(stat="identity", fill=cbPalette[1])
      if (intervals)
      { 
        if (scheme.opt == "black")
          pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                     ymax=.data[[!!yupp]]),   
                                                 linetype = "solid", colour = "black") 
        else
          pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                     ymax=.data[[!!yupp]]),   
                                                 linetype = "solid", colour = cbPalette[5]) 
        if (interval.annotate)
        {
          pred.plot <- pred.plot + 
            annotate("text", x=Inf, y=-Inf,  hjust=1, vjust=-0.3, size = 2, 
                     label = paste("Error bars are ", labend, sep=""))
          # if (low.parts[2] == "halfLeastSignificant" && !is.na(meanLSD))
          #   pred.plot <- pred.plot + 
          #     annotate("text", x=-Inf, y=-Inf,  hjust=-0.01, vjust=-0.3, size = 2, 
          #              label = paste("mean LSD = ",signif(meanLSD, digits=3)))
        }
      }
    } else
      if (n.non.x == 2)
      { 
        pred.plot <-  pred.plot + 
          aes(x = .data[[!!vars[1]]], y = .data[[!!y]]) + 
          scale_fill_manual(values=cbPalette) + 
          facet_grid(as.formula(paste("~ ",vars[2],sep=""))) + 
          theme(axis.text.x=element_text(angle=90, hjust=1))
        if (scheme.opt == "black")
          pred.plot <-  pred.plot + geom_bar(stat="identity", fill="grey50") 
        else
          pred.plot <-  pred.plot + geom_bar(stat="identity", fill=cbPalette[1])
        if (intervals)
        { 
          non.x.lev <- levels(data[[vars[2]]])
          annot <- data.frame(Inf, -Inf, 
                              factor(non.x.lev[length(non.x.lev)], levels = non.x.lev))
          names(annot) <- c(vars[1], y, vars[2])
          if (scheme.opt == "black")
            pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                       ymax=.data[[!!yupp]]),   
                                                   linetype = "solid", colour = "black") 
          else
            pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                       ymax=.data[[!!yupp]]),   
                                                   linetype = "solid", colour = cbPalette[5]) 
          if (nos.lev[2] > 6)
            pred.plot <- pred.plot +
            theme(strip.text.x = element_text(size = 8))
          
          if (interval.annotate)
          {
            pred.plot <- pred.plot + 
              geom_text(data = annot, label = "Error bars are", 
                        hjust=1, vjust=-1.3, size = 2) +
              geom_text(data = annot, label = labend, 
                        hjust=1, vjust=-0.3, size = 2)
            # if (low.parts[2] == "halfLeastSignificant" && !is.na(meanLSD)) 
            # { 
            #   annot <- data.frame(-Inf, -Inf, 
            #                       factor(non.x.lev[1], levels = non.x.lev))
            #   names(annot) <- c(vars[1], y, vars[2])
            #   pred.plot <- pred.plot + 
            #     geom_text(data = annot, hjust=-0.01, vjust=-0.3, size = 2, 
            #               label = paste("mean LSD = ",signif(meanLSD, digits=3)))
            # }
          }
        }
      } else
        if (n.non.x == 3)
        { 
          pred.plot <-  pred.plot + 
            aes(x = .data[[!!vars[1]]], y = .data[[!!y]]) + 
            scale_fill_manual(values=cbPalette) + 
            facet_grid(as.formula(paste(vars[3]," ~ ",vars[2],sep=""))) + 
            theme(axis.text.x=element_text(angle=90, hjust=1))
          if (scheme.opt == "black")
            pred.plot <-  pred.plot + geom_bar(stat="identity", fill="grey50") 
          else
            pred.plot <-  pred.plot + geom_bar(stat="identity", fill=cbPalette[1])
          if (nos.lev[2] > 6 || nos.lev[3] > 6)
            pred.plot <- pred.plot +
              theme(strip.text.x = element_text(size = 8),
                    strip.text.y = element_text(size = 8, angle = -90))
          if (intervals)
          { 
            non.x.lev1 <- levels(data[[vars[3]]])
            non.x.lev2 <- levels(data[[vars[2]]])
            annot <- data.frame(Inf, -Inf, 
                                factor(non.x.lev1[length(non.x.lev1)], levels = non.x.lev1),
                                factor(non.x.lev2[length(non.x.lev2)], levels = non.x.lev2))
            names(annot) <- c(vars[1], y, vars[c(3,2)])
            if (scheme.opt == "black")
              pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                         ymax=.data[[!!yupp]]),   
                                                     linetype = "solid", colour = "black") 
            else
              pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                         ymax=.data[[!!yupp]]),   
                                                     linetype = "solid", colour = cbPalette[5]) 
            
            if (interval.annotate)
            {
              pred.plot <- pred.plot + 
                geom_text(data = annot, label = "Error bars are", 
                          hjust=1, vjust=-2.1, size = 2) +
                geom_text(data = annot, label = labend, 
                          hjust=1, vjust=-1.1, size = 2)
              # if (low.parts[2] == "halfLeastSignificant" && !is.na(meanLSD)) 
              # { 
              #   annot <- data.frame(-Inf, -Inf, 
              #                       factor(non.x.lev1[length(non.x.lev1)], levels = non.x.lev1),
              #                       factor(non.x.lev2[1], levels = non.x.lev2))
              #   names(annot) <- c(vars[1], y, vars[c(3,2)])
              #   pred.plot <- pred.plot + 
              #     geom_text(data = annot, hjust=-0.01, vjust=-1.1, size = 2, 
              #               label = paste("mean LSD = ",signif(meanLSD, digits=3)))
              # }
            }
          }
        }
  }
  else
    #Do line plot
  { 
    if (intervals)
      if (!(is.null(x.num)) && x.num %in% vars)
        int.width <- (max(data[[x.num]]) - min(data[[x.num]]))*0.025
    else
      int.width <- length(levels(data[[x.fac]]))*0.025
    if (!(is.null(x.num)) && x.num %in% vars)
      x.var <- x.num
    else
      x.var <- x.fac
    if (!is.null(titles) & !is.na(match(x.var,names(titles))))
      x.title <- titles[[x.var]]
    else
      x.title <- x.var
    pred.plot <-  ggplot(data=data, ...) + theme_bw() +
      scale_x_continuous(x.title) + scale_y_continuous(y.title)
    if (scheme.opt == "black")
      pred.plot <-  pred.plot + theme_bw() +
      theme(panel.grid.major = element_line(colour = "grey95", 
                                            linewidth = 0.5), 
            panel.grid.minor = element_line(colour = "grey98", 
                                            linewidth = 0.5))
    #If no non.x.terms ignore single & multiple
    if (n.non.x == 0)
    { 
      pred.plot <-  pred.plot + 
        aes(x = .data[[!!x.var]], y = y) + 
        scale_colour_manual(values=cbPalette) + 
        scale_shape_manual(values=symb) + 
        geom_point(shape = symb[7])
      if (scheme.opt == "black")
        pred.plot <-  pred.plot + geom_line(colour = "black") 
      else
        pred.plot <-  pred.plot + geom_line(colour = cbPalette[1])
      if (intervals)
      { 
        if (scheme.opt == "black")
          pred.plot <-  pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                      ymax=.data[[!!yupp]]),  
                                                  linetype = "solid", colour = "black", width=int.width) 
        else
          pred.plot <-  pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                      ymax=.data[[!!yupp]]),  
                                                  linetype = "solid", colour = cbPalette[5], width=int.width) 
        
        if (interval.annotate)
        {
          pred.plot <- pred.plot + 
            annotate("text", x=Inf, y=-Inf,  hjust=1, vjust=-0.3, size = 2, 
                     label =  paste("Error bars are ", labend, sep=""))
          # if (low.parts[2] == "halfLeastSignificant" && !is.na(meanLSD)) 
          #   pred.plot <- pred.plot + 
          #     annotate("text", x=-Inf, y=-Inf,  hjust=-0.01, vjust=-0.3, size = 2, 
          #              label = paste("mean LSD = ",signif(meanLSD, digits=3)))
        }
      }
    } else
      if (panel.opt == "single")
        #Single with non.x.terms
      { 
        if (n.non.x == 1)
        { 
          if (scheme.opt == "black")
            pred.plot <-  pred.plot + aes(x = .data[[!!x.var]], y = .data[[!!y]], 
                                          linetype=non.x.terms, shape=non.x.terms)
          else
            pred.plot <-  pred.plot + aes(x = .data[[!!x.var]], y = .data[[!!y]], 
                                          colour = non.x.terms, 
                                          linetype=non.x.terms, shape=non.x.terms) 
          pred.plot <-  pred.plot + 
            labs(colour=non.x.terms, linetype=non.x.terms, shape=non.x.terms) +
            # scale_colour_manual(values=cbPalette) + 
            scale_shape_manual(values=symb) + 
            geom_line() + 
            geom_point()
          if (intervals)
          { 
            pred.plot <- pred.plot + 
              geom_errorbar(aes(ymin=.data[[!!ylow]], ymax=.data[[!!yupp]]),   
                            linetype = "solid", position = position_dodge(1))
            if (interval.annotate)
            {
              pred.plot <- pred.plot + 
                annotate("text", x=Inf, y=-Inf,  hjust=1, vjust=-0.3, size = 2, 
                         label =  paste("Error bars are ", labend, sep=""))
              # if (low.parts[2] == "halfLeastSignificant" && !is.na(meanLSD)) 
              # { pred.plot <- pred.plot + 
              #   annotate("text", x=-Inf, y=-Inf,  hjust=-0.01, vjust=-0.3, size = 2, 
              #            label = paste("mean LSD = ",signif(meanLSD, digits=3)))
              # }
            }
          }
        }
      }
    else #"multiple"
    { if (n.non.x == 1)
    { 
      pred.plot <- pred.plot + 
        aes(x = .data[[!!x.var]], y = .data[[!!y]]) +
        scale_colour_manual(values=cbPalette) + 
        scale_shape_manual(values=symb) + 
        geom_point(shape = symb[7])  + 
        facet_wrap(as.formula(paste("~ ",non.x.terms,sep="")))
      if (scheme.opt == "black")
        pred.plot <- pred.plot + geom_line(colour = "black")
      else
        pred.plot <- pred.plot + geom_line(colour = cbPalette[1])
      if (intervals)
      { 
        non.x.lev <- levels(data[[non.x.terms]])
        annot <- data.frame(Inf, -Inf, 
                            factor(non.x.lev[length(non.x.lev)], levels = non.x.lev))
        names(annot) <- c(x.var, y, non.x.terms)
        if (scheme.opt == "black")
          pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                     ymax=.data[[!!yupp]]),   
                                                 linetype = "solid", colour = "black", width=int.width) 
        else
          pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                     ymax=.data[[!!yupp]]),   
                                                 linetype = "solid", colour = cbPalette[5], width=int.width) 
        pred.plot <- pred.plot + 
          geom_text(data = annot, label = "Error bars are", 
                    hjust=1, vjust=-1.3, size = 2) +
          geom_text(data = annot, label = labend, 
                    hjust=1, vjust=-0.3, size = 2)
        # if (interval.annotate & low.parts[2] == "halfLeastSignificant" && !is.na(meanLSD)) 
        # { annot <- data.frame(-Inf, -Inf, 
        #                       factor(non.x.lev[1], levels = non.x.lev))
        #   names(annot) <- c(x.var, y, non.x.terms)
        #   pred.plot <- pred.plot + 
        #                 geom_text(data = annot, hjust=-0.01, vjust=-0.3, size = 2, 
        #                          label = paste("mean LSD = ",signif(meanLSD, digits=3)))
        # }
      }
    }
      if (n.non.x == 2)
      { 
        pred.plot <- pred.plot + 
          aes(x = .data[[!!x.var]], y = .data[[!!y]]) +
          scale_colour_manual(values=cbPalette) + 
          scale_shape_manual(values=symb) + 
          geom_point(shape = symb[7]) +
          facet_grid(as.formula(paste(vars[2]," ~ ",vars[1],sep="")))
        if (scheme.opt == "black")
          pred.plot <- pred.plot + geom_line(colour = "black")
        else
          pred.plot <- pred.plot + geom_line(colour = cbPalette[1])
        if (nos.lev[1] > 6 || nos.lev[2] > 6)
          pred.plot <- pred.plot +
            theme(strip.text.x = element_text(size = 8),
                  strip.text.y = element_text(size = 8, angle = -90))
        if (intervals)
        { 
          non.x.lev1 <- levels(data[[vars[1]]])
          non.x.lev2 <- levels(data[[vars[2]]])
          annot <- data.frame(Inf, -Inf, 
                              factor(non.x.lev1[length(non.x.lev1)], levels = non.x.lev1),
                              factor(non.x.lev2[length(non.x.lev2)], levels = non.x.lev2))
          names(annot) <- c(x.var, y, non.x.terms)
          if (scheme.opt == "black")
            pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                       ymax=.data[[!!yupp]]),   
                                                   linetype = "solid", colour = "black", width=int.width) 
          else
            pred.plot <- pred.plot + geom_errorbar(aes(ymin=.data[[!!ylow]], 
                                                       ymax=.data[[!!yupp]]),   
                                                   linetype = "solid", colour = cbPalette[5], width=int.width) 
          
          if (interval.annotate)
          {
            pred.plot <- pred.plot + 
              geom_text(data = annot, label = "Error bars are", 
                        hjust=1, vjust=-1.3, size = 2) +
              geom_text(data = annot, label = labend, 
                        hjust=1, vjust=-0.3, size = 2)
            # if (low.parts[2] == "halfLeastSignificant" && !is.na(meanLSD)) 
            # { annot <- data.frame(-Inf, -Inf, 
            #                       factor(non.x.lev1[1], levels = non.x.lev1),
            #                       factor(non.x.lev2[length(non.x.lev2)], levels = non.x.lev2))
            # names(annot) <- c(x.var, y, non.x.terms)
            # pred.plot <- pred.plot + 
            #   geom_text(data = annot, hjust=-0.01, vjust=-0.3, size = 2, 
            #             label = paste("mean LSD = ",signif(meanLSD, digits=3)))
            # }
          }
        }
      }
    }
  }
  if (!is.null(ggplotFuncs))
    for (f in ggplotFuncs)
      pred.plot <- pred.plot + f
  print(pred.plot)
  #Automate saving of plots in files
  if (!is.null(filestem))
  { 
    filename <- paste(filestem,".",gsub(":",".",classify, fixed= TRUE),sep="")
    if (barplot)
      filename <- paste(filename,"Bar", sep=".")
    else
      filename <- paste(filename,"Line", sep=".")
    if (intervals)
      filename <- paste(filename,abbrev,sep=".")
    filename <- paste(filename,".png",sep="")
    savePlot(filename = filename, type = "png")
    cat("\n#### Plot saved in ", filename,"\n")
  }
  invisible(pred.plot)
}

"predictPresent.asreml" <- function(asreml.obj, terms, inestimable.rm = TRUE, 
                                    linear.transformation = NULL, EGLS.linTransform = TRUE, 
                                    error.intervals = "Confidence", alpha = 0.05, 
                                    wald.tab = NULL, dDF.fault = "residual", dDF.values = NULL, 
                                    pairwise = TRUE, Vmatrix = FALSE,
                                    avsed.tolerance = 0.25, accuracy.threshold = NA, 
                                    LSDtype = "overall", LSDsupplied = NULL, LSDby = NULL, 
                                    LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
                                    x.num = NULL, x.fac = NULL, nonx.fac.order = NULL,  
                                    x.pred.values = NULL, x.plot.values = NULL, 
                                    plots = "predictions", panels = "multiple", 
                                    graphics.device = NULL, interval.annotate = TRUE,
                                    titles = NULL, colour.scheme = "colour", save.plots = FALSE, 
                                    transform.power = 1, offset = 0, scale = 1, 
                                    transform.function = "identity", 
                                    tables = "all", level.length = NA, 
                                    sortFactor = NULL, sortParallelToCombo = NULL, 
                                    sortNestingFactor = NULL, sortOrder = NULL, 
                                    decreasing = FALSE, trace = FALSE, 
                                    ggplotFuncs = NULL, ...)
#This function forms the predictions for each significant term.
#It then presents either a table or a graph based on the predicted values 
# - the decision is based on whether x.fac or x.num is in the term. 
#Might need to supply additional arguments to predict through "..."
#e.g. present
#Must have levels of x.fac in the order in which they are to be plotted
# - with dates, they should be in the form yyyymmdd
#Probably need to supply x.plot.values if x.fac is to be plotted
{ 
  tempcall <- list(...)
  if (length(tempcall)) 
  {
    #Check for deprecated argument meanLSD.type and warn
    if ("meanLSD.type" %in% names(tempcall))
      stop("meanLSD.type has been deprecated - use LSDtype")
    #Check for deprecated dDF.na argument
    if ("dDF.na" %in% names(tempcall))
      stop("The argument dDF.na has been deprecated; use dDF.fault instead.")
  }

  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  #Check that have a valid object of class asreml
  validasr <- validAsreml(asreml.obj)  
  if (is.character(validasr))
    stop(validasr)
  
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction", "supplied")
  avLSD <- AvLSD.options[check.arg.values(LSDtype, AvLSD.options)]
  if (!is.null(LSDby) &&  !is.character(LSDby))
    stop("LSDby must be a character")
  
  #check the transform.function
  link.opts <- c("identity", "log", "inverse", "sqrt", "logit", "probit", "cloglog")
  transfunc <- link.opts[check.arg.values(transform.function, link.opts)]
  
  tempcall <- list(...)
  if ("levels.length" %in% names(tempcall))
    stop("levels.length has been deprecated - use level.length")
  if ("vcov" %in% names(tempcall))
    stop("Use Vmatrix to request that the vcov matrix be saved")
  
  validwald <- validWaldTab(wald.tab)  
  if (is.character(validwald))
    stop(validwald)
  if (!is.null(x.pred.values) && !is.null(x.plot.values))
    if (length(x.pred.values) != length(x.plot.values))
    {
      if (asr4)
        stop("In analysing ",asreml.obj$formulae$fixed[[2]],
             ", length of x.pred.values and x.plot.values should be equal")
      else
        stop("In analysing ",asreml.obj$fixed.formula[[2]],
             ", length of x.pred.values and x.plot.values should be equal")
    }
  #Check options
  scheme.options <- c("black", "colour")
  scheme.opt <- scheme.options[check.arg.values(colour.scheme, scheme.options)]
  panel.options <- c("single", "multiple")
  panel.opt <- panel.options[check.arg.values(panels, panel.options)]
  plot.options <- c("none", "predictions", "backtransforms", "both")
  plot.opt <- plot.options[check.arg.values(plots, plot.options)]
  if (!is.na(avsed.tolerance) & (avsed.tolerance <0 | avsed.tolerance > 1))
    stop("avsed.tolerance should be between 0 and 1")
  int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
  int.opt <- int.options[check.arg.values(error.intervals, int.options)]
  
  #Determine if there are any model terms with x.num and any with x.fac
  if (asr4)
    model.terms <- c(attr(terms(asreml.obj$formulae$fixed), which="term.labels"),
                     attr(terms(asreml.obj$formulae$random), which="term.labels"))
  else
    model.terms <- c(attr(terms(asreml.obj$fixed.formula), which="term.labels"),
                     attr(terms(asreml.obj$random.formula), which="term.labels"))
  term.types <- trend.terms.types(terms = model.terms, 
                                  trend.num = x.num , devn.fac = x.fac)
  diff.list <- vector("list", length = length(terms))
  names(diff.list) <- gsub(":", ".", terms, fixed= TRUE)
  for (term in terms) 
    #Get the predictions
  { 
    #Make sure no functions in classify
    factors <- fac.getinTerm(term, rmfunction = TRUE)
    classify.term <- fac.formTerm(factors)
    #Check if current term involves x.num or x.fac
    if ((!is.null(x.num) && x.num %in% factors) || (!is.null(x.fac) && x.fac %in% factors))
      #Check if the set of model terms is mixed for trend and deviations
    { 
      if (term.types$has.trend.num && term.types$has.devn.fac)
        #Mixed
      { 
        if (!(x.num %in% factors))
          classify.term <- paste(classify.term, x.num, sep=":")
        else
          if (!(x.fac %in% factors))
            classify.term <- paste(classify.term, x.fac, sep=":")
      }
      diffs <- predictPlus.asreml(asreml.obj = asreml.obj, 
                                  classify = classify.term, term = term, 
                                  linear.transformation = linear.transformation, 
                                  EGLS.linTransform = EGLS.linTransform, 
                                  titles = titles, 
                                  x.num = x.num, x.fac = x.fac,  
                                  x.pred.values = x.pred.values, 
                                  x.plot.values = x.plot.values, 
                                  error.intervals = error.intervals, 
                                  avsed.tolerance = avsed.tolerance, 
                                  accuracy.threshold = accuracy.threshold, 
                                  LSDtype = LSDtype, LSDsupplied = LSDsupplied, 
                                  LSDby = LSDby, LSDstatistic = LSDstatistic, 
                                  LSDaccuracy = LSDaccuracy, 
                                  pairwise = pairwise, Vmatrix = Vmatrix, 
                                  tables = tables, 
                                  level.length = level.length, 
                                  transform.power = transform.power, 
                                  offset = offset, scale = scale, 
                                  transform.function = transfunc, 
                                  inestimable.rm = inestimable.rm, 
                                  sortFactor = sortFactor, 
                                  sortParallelToCombo = sortParallelToCombo, 
                                  sortNestingFactor = sortNestingFactor, 
                                  sortOrder = sortOrder, 
                                  decreasing = decreasing, 
                                  wald.tab = wald.tab, dDF.fault = dDF.fault, 
                                  dDF.values = dDF.values, 
                                  alpha = alpha, trace = trace, ...)
    }
    #Not mixed
    else
      #No need for x.pred.values or to modify term
      diffs <- predictPlus.asreml(asreml.obj = asreml.obj, 
                                  classify = classify.term, 
                                  linear.transformation = linear.transformation, 
                                  EGLS.linTransform = EGLS.linTransform, 
                                  titles = titles,   
                                  x.num = x.num, x.fac = x.fac,  
                                  x.plot.values = x.plot.values, 
                                  error.intervals = error.intervals, 
                                  avsed.tolerance = avsed.tolerance, 
                                  accuracy.threshold = accuracy.threshold, 
                                  LSDtype = LSDtype, LSDsupplied = LSDsupplied, 
                                  LSDby = LSDby, LSDstatistic = LSDstatistic, 
                                  LSDaccuracy = LSDaccuracy, 
                                  pairwise = pairwise, Vmatrix = Vmatrix, 
                                  tables = tables, 
                                  level.length = level.length, 
                                  transform.power = transform.power, 
                                  offset = offset, scale = scale, 
                                  transform.function = transfunc, 
                                  inestimable.rm = inestimable.rm, 
                                  sortFactor = sortFactor, 
                                  sortParallelToCombo = sortParallelToCombo, 
                                  sortNestingFactor = sortNestingFactor, 
                                  sortOrder = sortOrder, 
                                  decreasing = decreasing, 
                                  wald.tab = wald.tab, dDF.fault = dDF.fault, 
                                  dDF.values = dDF.values, 
                                  alpha = alpha, trace = trace, ...)
    new.name <- gsub(":", ".", term)
    diff.list[[new.name]] <- diffs
    #check have some predictions
    if (nrow(diffs$predictions) == 0)
      warning("No estimable predicted values")
    else
    {
      #If plot required, set up y.title
      if (plot.opt != "none")
      { 
        #Set y.title
        y.title <- as.character(attr(diffs, which = "response.title"))
        if (is.null(y.title))
          y.title <- as.character(attr(diffs, which = "response"))
        if (save.plots)
          filestem <- as.character(attr(diffs, which = "response"))
        else
          filestem <- NULL
      }
      #Plot predictions if a plot is required and data are not transformed or if 
      #plot of predictions requested
      transformed <- transform.power != 1 | offset != 0 | scale != 1 | transform.function != "identity" 
      if ((!transformed & plot.opt != "none") | plot.opt == "predictions" 
          | plot.opt == "both")
        #Plot predictions
      { 
        plotPredictions(data = diffs$predictions, classify = classify.term, 
                        y = "predicted.value", 
                        x.num = x.num, x.fac = x.fac, 
                        nonx.fac.order = nonx.fac.order, 
                        colour.scheme = colour.scheme, 
                        panels = panel.opt, graphics.device = graphics.device,
                        error.intervals = error.intervals,  
                        interval.annotate = interval.annotate, 
                        titles = titles, y.title = y.title, 
                        filestem = filestem, ggplotFuncs = ggplotFuncs, ...)
      }
      if (transformed & (plot.opt == "backtransforms"  | plot.opt == "both"))
        #Plot backtransforms
      { 
        bty.title <- paste("Backtransform of ",y.title,sep="")
        if (save.plots)
          filestem <- paste(filestem,".back",sep="")
        else
          filestem <- NULL
        plotPredictions(data = diffs$backtransforms, 
                        classify = classify.term, 
                        y = "backtransformed.predictions", 
                        x.num = x.num, x.fac = x.fac,   
                        nonx.fac.order = nonx.fac.order, 
                        colour.scheme = colour.scheme,  
                        panels = panel.opt, graphics.device = graphics.device,
                        error.intervals = error.intervals,
                        interval.annotate = interval.annotate, 
                        titles = titles, y.title = bty.title, 
                        filestem = filestem, ggplotFuncs = ggplotFuncs, ...)
      }
    }
  }
  invisible(diff.list)
}



#Calls for prediction routines:
# - predictPresent.asreml calls predictPlus.asreml & plotPredictions.data.frame
# - predictPlus.asreml calls allDifferences.data.frame