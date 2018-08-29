#alldiffs functions

#Form an alldiffs object from supplied component objects
"as.alldiffs" <- function(predictions, vcov = NULL, 
                          differences = NULL, p.differences = NULL, 
                          sed = NULL, LSD = NULL, backtransforms = NULL, 
                          response = NULL, response.title = NULL, 
                          term = NULL, classify = NULL, 
                          tdf = NULL, sortFactor = NULL, sortOrder = NULL)
{ 
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  
  #Check arguments
  #Change asreml4 names to asreml3 names
  if ("std.error" %in% colnames(predictions))
    names(predictions)[match("std.error", names(predictions))] <- "standard.error"
  if ("status" %in% colnames(predictions))
    names(predictions)[match("status", names(predictions))] <- "est.status"
  if (!is.null(sed))
    sed <- as.matrix(sed)
  if (!is.null(vcov))
    vcov <- as.matrix(vcov)
  #Check have appropriate columns
  if (!("predicted.value" %in% colnames(predictions)) || 
      !("standard.error" %in% colnames(predictions)) || !("est.status" %in% colnames(predictions))) 
    warning("Predictions argument does not include the expected column names (e.g. predicted.value)")
  npred <- nrow(predictions)
  if ((!is.null(differences) && !("matrix" %in% class(differences))) ||
      (!is.null(p.differences) && !("matrix" %in% class(p.differences))) || 
      (!is.null(sed) && !("matrix" %in% class(sed))) || 
      (!is.null(vcov) && !("matrix" %in% class(vcov))))
    warning("At least one of differences, p.differences, sed and vcov is not of type matrix")
  if (!is.null(differences) && !is.null(p.differences) && !is.null(sed))
  { 
    dimens <- c(nrow(differences), nrow(p.differences), nrow(sed), 
                ncol(differences), ncol(p.differences), ncol(sed))
    if (any(npred != dimens))
      stop("At least one of differences, p.differences or sed is not conformable with predictions")
  }
  if (!is.null(vcov))
  {
    #check that vcov conforms to predictions
    if (any(npred != c(nrow(vcov), ncol(vcov))))
      stop("vcov is not conformable with predictions")
  }
  if (!is.null(backtransforms))
  { 
    if (!("backtransformed.predictions" %in% colnames(backtransforms))) 
      warning("Backtransforms argument does not include a column named backtransformed.predictions")
    if (npred != nrow(backtransforms))
      stop("Backtransforms do not contain the same number of rows as the predictions")
  }
  #ensure diag of sed is NA
  if (!is.null(sed))
    diag(sed) <- NA
  meanLSD <- NULL
  if (!is.null(LSD))
    attr(predictions, which = "meanLSD") <- LSD$meanLSD
  p <- list(predictions = predictions, vcov = vcov, 
            differences = differences, p.differences = p.differences, sed = sed, 
            LSD = LSD, backtransforms = backtransforms)
  attr(p, which = "response") <- response
  attr(p, which = "response.title") <- response.title
  attr(p, which = "term") <- term
  attr(p, which = "classify") <- classify
  attr(p, which = "tdf") <- tdf
  attr(p, which = "sortFactor") <- sortFactor
  attr(p, which = "sortOrder") <- sortOrder
  class(p) <- "alldiffs"
  return(p)
}

"print.alldiffs" <- function(x, which = "all", ...)
{ 
  options <- c("predictions", "backtransforms", "vcov", 
               "differences", "p.differences", "sed", "LSD", "all")
  opt <- options[unlist(lapply(which, check.arg.values, options=options))]
  title <- attr(x, which = "response.title")
  if (is.null(title))
    title <- as.character(attr(x, which = "response"))
  term <- attr(x, which = "term")
  if (!is.null(term))
    title <- paste(title, " from ", term, sep="")
  #Print predictions and/or LSDs
  if ("all" %in% opt || "predictions" %in% opt || "LSD" %in% opt)
  { 
    if ("all" %in% opt || "predictions" %in% opt)
    { 
      if (!is.null(title))
      {
        tt <- paste("\n\n#### Predictions for ", title, "\n\n",sep="")
        cat(tt)
      }
      print(x$predictions)
    } 
    if (!is.null(x$LSD))
    { 
      sed.range <- abs(x$LSD$minLSD - x$LSD$maxLSD) / x$LSD$meanLSD
      if (length(x$LSD$meanLSD > 1))
      {
        cat("\n\nLSD values \n\n")
        cat("minimum LSD = ",x$LSD$minLSD,"\n")
        cat("\nmean LSD = ",x$LSD$meanLSD,"\n")
        cat("\nmaximum LSD = ",x$LSD$maxLSD,"\n")
        cat("\n(sed range / mean sed = ",signif(sed.range, digits=3),")\n\n")
      } else
      {
        cat("\n\nLSD values \n\n")
        cat("minimum LSD = ",x$LSD$minLSD,  "  mean LSD = ",x$LSD$meanLSD,
            "  maximum LSD = ",x$LSD$maxLSD,
            "\n(sed range / mean sed = ",signif(sed.range, digits=3),")\n\n")
      }
    } else
    {
      if (opt %in% c("all", "LSD"))
      {
        cat("\n\nLSD values \n\n")
        print(x$LSD)      
      }
    }
  }
  if ("all" %in% opt || "vcov" %in% opt)
  { 
    cat("\n\nVariance matrix of the predicted values \n\n")
    if (!is.null(x$vcov))
      print(zapsmall(x$vcov, 4))
    else
      print(x$vcov)
  }
  
  if ("all" %in% opt || "differences" %in% opt)
  { 
    cat("\n\nAll pairwise differences between predicted values \n\n")
    print(x$differences, digits=4)
  }
  if ("all" %in% opt || "p.differences" %in% opt)
  { 
    cat("\n\np values for all pairwise differences between predicted values \n\n")
    print(formatC(x$p.differences, digits=3, format="f"), quote=FALSE)
  }
  if ("all" %in% opt || "sed" %in% opt)
  { 
    cat("\n\nStandard errors of differences between predicted values \n\n")
    print(zapsmall(x$sed, 4))
  }
  if (("all" %in% opt & !is.null(x$backtransforms)) || "backtransforms" %in% opt)
  { 
    if (!is.null(title))
      cat("\n\n#### Backtransforms of predictions for ", title, "\n\n")
    print(x$backtransforms)
  }
  invisible()
}

makePredictionLabels <- function(predictions, classify, response = NULL, 
                                 x.num = NULL, x.fac = NULL, 
                                 level.length = NA)
{
  #determine factors for row and column names
  #Make sure no functions in classify
  factors <- fac.getinTerm(classify, rmfunction = TRUE)
  classify <- fac.formTerm(factors)
  nfac <- length(factors)
  #Check all factors in classify are in predictions
  if (length(setdiff (factors, names(predictions))) != 0)
  { 
    if (!is.null(response))
      stop("For ",response,
           ", there are factors in the classify argument that do not have columns in alldiffs.obj$predictions")
    else
      stop("There are factors in the classify argument that do not have columns in alldiffs.obj$predictions")
  }
  #Make sure only one of the numeric and factor that are parallel
  if ((!is.null(x.num) && x.num %in% factors) && (!is.null(x.fac) && x.fac %in% factors))
  { 
    k <- match(x.num, names(predictions))
    predictions <- predictions[, -k]
    nfac <- nfac - 1
  }
  
  #Generate row and column names as the combinations of the levels of factors
  if (nfac > 1)
  { 
    pred.faclist <- vector("list", length=nfac)
    nclassify <- ncol(predictions) - 3
    pred.names <- names(predictions)
    kk <- 0
    for (k in 1:nclassify)
    { 
      if (pred.names[k] %in% factors)
      { 
        kk <- kk + 1
        pred.faclist[[kk]] <- predictions[[k]]
        if (is.numeric(pred.faclist[[kk]]))
          pred.faclist[[kk]] <- factor(pred.faclist[[kk]])
        names(pred.faclist)[kk] <- pred.names[k]
      }
    }
    pred.faclist <- lapply(pred.faclist, 
                           FUN=function(ff)
                           {
                             if (is.character(levels(ff)) && !is.na(level.length))
                               ff <- factor(ff, labels=substr(levels(ff), start=1, 
                                                              stop=level.length))
                             invisible(ff)
                           })
    pred.lev <- as.character(fac.combine(pred.faclist, combine.levels=TRUE))
  }
  else
  { 
    k <- match(factors[[1]], names(predictions))
    pred.fac <- predictions[[k]]
    if (is.numeric(pred.fac))
      pred.fac <- factor(pred.fac)
    pred.lev <- levels(pred.fac)
    if (is.character(pred.lev) && !is.na(level.length))
      pred.lev <- substr(pred.lev, start=1, stop=level.length)
  }
  return(list(predictions = predictions, pred.lev = pred.lev))
}

facCombine.alldiffs <- function(object, factors, order="standard", combine.levels=TRUE, 
                                sep="_", level.length = NA,  ...)
{
  if (any(!(factors %in% names(object$predictions))))
    stop("Some factors are not in the predictions component of object")
  if (length(factors) <= 1)
    stop("Need at least two factors to combine")
  comb.fac <- fac.combine(object$predictions[factors], 
                          order = order, combine.levels = combine.levels, sep = sep)
  newfac <- paste(factors, collapse = sep)
  fstfac <- match(factors[1], names(object$predictions))
  object$predictions[fstfac] <- comb.fac
  names(object$predictions)[fstfac] <- newfac
  object$predictions <- object$predictions[-c(match(factors[-1], names(object$predictions)))]
  
  if (!is.null(object$backtransforms))
  {
    fstfac <- match(factors[1], names(object$backtransforms))
    object$backtransforms[fstfac] <- comb.fac
    names(object$backtransforms)[fstfac] <- newfac
    object$backtransforms <- object$backtransforms[-c(match(factors[-1], 
                                                            names(object$backtransforms)))]
  }
  
  classify <- attr(object, which = "classify")
  if (is.null(classify))
    stop("The alldiffs object does not have the classify attrtibute set")
  class.facs <- fac.getinTerm(classify, rmfunction = TRUE)
  class.facs[fstfac] <- newfac
  class.facs <- class.facs[-c(match(factors[-1], class.facs))]
  classify <- fac.formTerm(class.facs)
  attr(object, which = "classify") <- classify
  response <- attr(object, which = "response")
  pred.labs <- makePredictionLabels(object$predictions, classify, response)
  pred.lev <- pred.labs$pred.lev
  
  
  if (!is.null(object$vcov))
  {
    colnames(object$vcov) <- rownames(object$vcov) <- pred.lev
  }
  if (!is.null(object$differences))
  {
    colnames(object$differences) <- rownames(object$differences) <- pred.lev
  }
  if (!is.null(object$p.differences))
  {
    colnames(object$p.differences) <- rownames(object$p.differences) <- pred.lev
  }
  if (!is.null(object$sed))
  {
    colnames(object$sed) <- rownames(object$sed) <- pred.lev
  }
  return(object)
}

subset.alldiffs <- function(x, subset, rmClassifyVars = NULL, ...)
{
  #Deal with unsupported parameters
  tempcall <- list(...)
  if (length(tempcall)) 
  {
    if ("select" %in% names(tempcall))
      stop("select is not supported in subset.alldiffs")
    if ("drop" %in% names(tempcall))
      stop("drop is not supported in subset.alldiffs")
  }
  
  if (missing(subset) & is.null(rmClassifyVars)) 
    warning("neither subset or rmClassifyVars have been set")
  if (!missing(subset)) 
  {
    expr <- substitute(subset)
    cond <- eval(expr, x$predictions, parent.frame())
    if (!is.logical(cond)) 
      stop("'subset' must be logical")
    x$predictions <- x$predictions[cond & !is.na(cond),]
    x$predictions <- as.data.frame(lapply(x$predictions, 
                                          function(x)
                                          {
                                            if (class(x) == "factor")
                                              x <- factor(x)
                                            return(x)
                                          }), stringsAsFactors = FALSE)
    if (!is.null(x$vcov))
    {
      x$vcov <- x$vcov[cond, cond]
    }
    if (!is.null(x$backtransforms))
    {
      x$backtransforms <- x$backtransforms[cond & !is.na(cond),]
    }
    if (!is.null(x$differences))
    {
      x$differences <- x$differences[cond,cond]
    }
    if (!is.null(x$p.differences))
    {
      x$p.differences <- x$p.differences[cond,cond]
    }
    if (!is.null(x$sed))
    {
      x$sed <- x$sed[cond,cond]
      x <- recalcLSD(x, ...)
    }
  }
  if (!is.null(rmClassifyVars))
  {
    classify <- attr(x, which = "classify")
    if (is.null(classify))
      stop("The alldiffs object does not have the classify attrtibute set")
    rmfac <- fac.combine(as.list(x$predictions[rmClassifyVars]))
    if (length(unique(rmfac)) > 1)
      stop("The classify variables to be removed have more than one combination in the predictions; \nthe predictions cannot be unambiguosly identified.")
    class.facs <- fac.getinTerm(classify, rmfunction = TRUE)
    class.facs <- class.facs[!(class.facs %in% rmClassifyVars)]
    classify <- fac.formTerm(class.facs)
    x$predictions <- x$predictions[, -c(match(rmClassifyVars, names(x$predictions)))]
    if (!is.null(x$backtransforms))
    {
      x$backtransforms <- x$backtransforms[, -c(match(rmClassifyVars, 
                                                      names(x$backtransforms)))]
    }
    attr(x, which = "classify") <- classify
    response <- attr(x, which = "response")
    if (is.null(response))
      stop("The alldiffs object does not have the response attrtibute set")
    pred.labs <- makePredictionLabels(x$predictions, classify, response)
    pred.lev <- pred.labs$pred.lev
    
    if (!is.null(x$vcov))
    {
      colnames(x$vcov) <- rownames(x$vcov) <- pred.lev
    }
    if (!is.null(x$differences))
    {
      colnames(x$differences) <- rownames(x$differences) <- pred.lev
    }
    if (!is.null(x$p.differences))
    {
      colnames(x$p.differences) <- rownames(x$p.differences) <- pred.lev
    }
    if (!is.null(x$sed))
    {
      colnames(x$sed) <- rownames(x$sed) <- pred.lev
    }
  }
  return(x)
}

#Function to sort the components of an alldiffs object according to the predicted values 
#indexed by a factor/numeric.
#If there are more than a single factor/numeric indexing the predictions, one can
#(1) sort within each of the variables other than the sortFactor, or
#(2) sort in parallel in which case either the first value in the predictions 
#component or single nominated values, for all but the sortFactor, can be used 
#to identify a subset to be used to determine the order for which all other subsets 
#are sorted.
sort.alldiffs <- function(x, decreasing = FALSE, classify = NULL, 
                          sortFactor = NULL, sortWithinVals = NULL, 
                          sortOrder = NULL, ...)
{
  if (!class(x) == "alldiffs")
    stop("x must be of class alldiffs")
  
  if (is.null(classify))
    classify <- attr(x, which = "classify")
  else
  {
    cl <- attr(x, which = "classify")
    if (!is.null(cl))
      if (classify != cl)
        warning(paste("Supplied classify is not the same as the classify attribute of",
                      "the alldiffs object"))
  }
  class.names <-  fac.getinTerm(classify)
  if (!all(class.names %in% names(x$predictions)))
    stop(paste("The predictions data.frame does not have a column for each variable", 
               "in the classify stored with alldiffs", sep = " "))
  nclassify <- length(class.names)
  if (nclassify < 1)
    stop("Cannot find the classify variables")
  if (!is.null(sortFactor))
  {
    if (!is.factor(x$predictions[[sortFactor]]))
      stop("sortFactor must be a factor")
    if (length(sortFactor) != 1)
      stop("Can only supply one sortFactor name")
  }
  
  #Get the order in which to sort the predictions
  if (nclassify == 1)
  {
    if (is.null(sortFactor))
    {
      sortFactor <- class.names[1]
      if (!is.null(sortFactor))
        if (!is.factor(sortFactor))
          stop("Single variable in classify is not a factor")
    }
    if (is.null(sortOrder))
    {
      val.ord <- data.frame(X1 = x$predictions[[sortFactor]], 
                            Ord = order(x$predictions$predicted.value, 
                                        decreasing=decreasing))
      names(val.ord)[1] <- sortFactor
    } else #have sortOrder
    {
      old.levs <- levels(x$predictions[[sortFactor]])
      if (length(sortOrder) != length(old.levs))
        stop("The number of values in sortOrder must equal the number of levels of sortFactor")
      val.ord <- data.frame(X1 = x$predictions[[sortFactor]],
                            Ord = match(as.character(sortOrder), old.levs))
      if (any(is.na(val.ord$Ord)))
        stop("Not all the values in sortOrder are levels in SortFactor")
    }
    tmp <- val.ord
  } else #classify > 1
  {
    if (is.null(sortFactor))
      stop("The classify for the predictions has multiple variables - need to set sortFactor")
    other.vars <- class.names[-(match(sortFactor, class.names))]
    
    #Work out order for sortFactor
    if (is.null(sortOrder))
    {
      if (is.null(sortWithinVals))
      {
        sortWithinVals <- as.list(x$predictions[1, other.vars])
        sortWithinVals <- lapply(sortWithinVals, 
                                 function(var)
                                 {if (is.factor(var)) {var <- as.character(var)}; return(var)})
        names(sortWithinVals) <- other.vars
      } else
      {
        if (length(sortWithinVals) != length(other.vars))
          stop("Need a value for each classify variable except sortFactor")
        else
        {
          if (!all(other.vars %in% names(sortWithinVals)))
            stop("The names in sortWithinVals do not match the names in the classify")
          if (!all(unlist(lapply(sortWithinVals, function(el){length(el)==1}))))
            stop("Each component of sortWithinVals can have only a single value")
        }
      }
      subs <- TRUE
      for (name in names(sortWithinVals))
      {
        subs <- subs & x$predictions[name] == sortWithinVals[[name]]
      }
      val.ord <- data.frame(X1 = x$predictions[subs, sortFactor], 
                            Ord = order(x$predictions[subs, "predicted.value"], 
                                        decreasing=decreasing))
      names(val.ord)[1] <- sortFactor
    } else #have sortOrder
    {
      old.levs <- levels(x$predictions[[sortFactor]])
      if (length(sortOrder) != length(old.levs))
        stop("The number of values in sortOrder must equal the number of levels of sortFactor")
      val.ord <- data.frame(X1 = factor(old.levs, levels = old.levs),
                            Ord = match(as.character(sortOrder), old.levs))
      names(val.ord)[1] <- sortFactor
      if (any(is.na(val.ord$Ord)))
        stop("Not all the values in sortOrder are levels in SortFactor")
    }
    
    #Get the order for the full set of predictions
    # - this work for unequal replication of sortFactor
    tmp <- x$predictions[c(other.vars, sortFactor)]
    rownames(tmp) <- as.character(1:nrow(tmp))
    tmp <- tmp[do.call(order,tmp),]
    reord <- as.numeric(rownames(tmp))
    tmp <- suppressMessages(plyr::join(tmp, val.ord))
    tmp <- split(tmp, as.list(tmp[other.vars]))
    tmp <- lapply(tmp, 
                  function(d)
                  { d$Ord <- c(1:nrow(d))[order(order(d$Ord))]; return(d)})
    tmp <- do.call(rbind, tmp)
    starts <- reshape::melt(table(tmp[c(other.vars)]))
    if (length(other.vars) == 1)
      names(starts)[1] <- other.vars
    names(starts)[match("value", names(starts))] <- "start"
    starts$start <- c(0,cumsum(starts$start[1:(nrow(starts)-1)]))
    tmp <- suppressMessages(plyr::join(tmp, starts))
    tmp$Ord <- with(tmp, reord[Ord + start])
    if (any(is.na(tmp[sortFactor])))
      warning("The order for some predicted values could not be determined.\n",
              "Make sure that all values of sortFactor occur for the subset defined by sortWithinVals.")
  }
  
  #Order the components that are present
  if (is.factor(x$predictions[[sortFactor]]))
  {
    newlevs <- as.character(levels(x$predictions[[sortFactor]]))[val.ord$Ord]
    x$predictions[sortFactor] <- factor(x$predictions[[sortFactor]], 
                                        levels = newlevs)
    if (!is.null(x$backtransforms))
      x$backtransforms[sortFactor] <- factor(x$backtransforms[[sortFactor]], 
                                             levels = newlevs)
  }
  x$predictions <- x$predictions[tmp$Ord,]
  if (!is.null(x$vcov))
    x$vcov <- x$vcov[tmp$Ord, tmp$Ord]
  if (!is.null(x$backtransforms))
    x$backtransforms <- x$backtransforms[tmp$Ord,]
  if (!is.null(x$differences))
    x$differences <- x$differences[tmp$Ord, tmp$Ord]
  if (!is.null(x$p.differences))
    x$p.differences <- x$p.differences[tmp$Ord, tmp$Ord]
  if (!is.null(x$sed))
    x$sed <- x$sed[tmp$Ord, tmp$Ord]
  
  #Set attributes
  attr(x, which = "sortFactor") <- sortFactor
  attr(x, which = "sortOrder") <- newlevs
  
  return(x)
}

recalcLSD.alldiffs <- function(alldiffs.obj, meanLSD.type = "overall", LSDby = NULL, 
                               alpha = 0.05, ...)
{
  kattr <- attributes(alldiffs.obj)
  alldiffs.obj <- allDifferences(alldiffs.obj$predictions, 
                                 classify = attr(alldiffs.obj, which = "classify"), 
                                 vcov = alldiffs.obj$vcov, 
                                 differences = alldiffs.obj$differences, 
                                 p.differences = alldiffs.obj$p.differences,
                                 sed = alldiffs.obj$sed, 
                                 tdf = attr(alldiffs.obj, which = "tdf"),
                                 meanLSD.type = meanLSD.type, LSDby = LSDby, ...)
  newattr <- attributes(alldiffs.obj)
  #Find missing atTributes in new alldiffs.obj and add them back in 
  kattr <- kattr[names(kattr)[!(names(kattr) %in% names(newattr))]]
  if (length(kattr) > 0)
  {
    newattr <- c(newattr,kattr)
    attributes(alldiffs.obj) <- newattr
  }
  return(alldiffs.obj)
}

redoErrorIntervals.alldiffs <- function(alldiffs.obj, error.intervals = "Confidence", 
                                        alpha = 0.05, avsed.tolerance = 0.25, 
                                        meanLSD.type = "overall", LSDby = NULL, ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction")
  avLSD <- AvLSD.options[check.arg.values(meanLSD.type, AvLSD.options)]
  if (!is.null(LSDby) &&  !is.character(LSDby))
    stop("LSDby must be a character")
  
  if (!is.na(avsed.tolerance) & (avsed.tolerance <0 | avsed.tolerance > 1))
    stop("avsed.tolerance should be between 0 and 1")
  int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
  int.opt <- int.options[check.arg.values(error.intervals, int.options)]
  
  denom.df <- attr(alldiffs.obj, which = "tdf")
  
  #Remove any intervals
  if (any(grepl("upper", names(alldiffs.obj$predictions), fixed = TRUE)))
  {
    alldiffs.obj$predictions <- alldiffs.obj$predictions[, -pmatch(c("lower.", "upper."), 
                                                                   names(alldiffs.obj$predictions))]
  }
  
  #Add lower and upper uncertainty limits
  if (int.opt != "none")
  { 
    revert <- FALSE
    if (is.null(denom.df) && c("Confidence", "halfLeastSignificant") %in% int.opt)
    {
      warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                    "- reverting to Standard Error"))
      int.opt <- "StandardError"
    }
    if (int.opt == "StandardError")
      alldiffs.obj$predictions <- within(alldiffs.obj$predictions, 
                                         { lower.StandardError.limit <- 
                                           alldiffs.obj$predictions[["predicted.value"]] - 
                                           alldiffs.obj$predictions[["standard.error"]]
                                         upper.StandardError.limit <- 
                                           alldiffs.obj$predictions[["predicted.value"]] + 
                                           alldiffs.obj$predictions[["standard.error"]]
                                         })
    else
    {
      t.value = qt(1-alpha/2, denom.df)
    }
    if (int.opt == "halfLeastSignificant" && (nrow(alldiffs.obj$predictions) != 1))
    { 
      overall.meanLSD <- sqrt(mean(alldiffs.obj$sed*alldiffs.obj$sed, na.rm = TRUE))
      overall.sed.range <- abs((max(alldiffs.obj$sed*alldiffs.obj$sed, na.rm = TRUE) - 
                                  min(alldiffs.obj$sed*alldiffs.obj$sed, na.rm = TRUE))) /overall.meanLSD
      overall.meanLSD <- t.value * overall.meanLSD
      nLSD <- length(alldiffs.obj$LSD$meanLSD)
      sed.range <- abs(alldiffs.obj$LSD$minLSD - alldiffs.obj$LSD$maxLSD) /  alldiffs.obj$LSD$meanLSD
      if (!is.na(avsed.tolerance) & overall.sed.range <= avsed.tolerance) #always plot overall LSD
      {
        alldiffs.obj$predictions <- within(alldiffs.obj$predictions, 
                                           { 
                                             lower.halfLeastSignificant.limit <- 
                                               alldiffs.obj$predictions[["predicted.value"]] - 
                                               0.5 * overall.meanLSD
                                             upper.halfLeastSignificant.limit <- 
                                               alldiffs.obj$predictions[["predicted.value"]] + 
                                               0.5 * overall.meanLSD
                                           })
      } else #process for each meanLSD.type option
      {
        if (avLSD == "overall")
        {
          if (nLSD != 1)
            stop("There is not just one LSD for meanLSD.type overall")
          warning("The avsed.tolerance is exceeded - reverting to confidence intervals")
          revert = TRUE
        } else
        {
          if (avLSD == "factor.combinations")
          {
            if (is.list(LSDby))
            {
              if (any(unlist(lapply(LSDby, function(x) class(x)!="factor"))))
                stop("Some components of the LSDby list are not factors")
              fac.list <- LSDby
            } else
            {
              if (class(LSDby) == "factor")
                fac.list <- list(LSDby)
              else
              {
                if (is.character(LSDby))
                  fac.list <- as.list(alldiffs.obj$predictions[LSDby])
                else
                  stop("by is not one of the allowed class of inputs")
              }
            }
            #Form levels combination for mean LSDs
            levs <- levels(fac.combine(as.list(alldiffs.obj$predictions[LSDby]), combine.levels = TRUE))
            #Check have got the correct LSDs
            if (is.null(rownames(alldiffs.obj$LSD)) | nLSD != length(levs) | 
                any(levs != rownames(alldiffs.obj$LSD)))
              stop(paste("For meanLSD.type factor.combinations, the LSD component of the alldiffs.obj", 
                         "must be a named vector of the LSDs for each combination of the factors in LSDby", 
                         sep = " "))
            if (any(na.omit(sed.range) > avsed.tolerance))
            {
              warning("The avsed.tolerance is exceeded for factor combinations - reverting to confidence intervals")
              revert <- TRUE
            } else #plot factor.combination LSD
            {
              levs <- strsplit(levs, ",", fixed = TRUE)
              nfac <- length(levs[[1]])
              LSD.dat <- as.data.frame(do.call(cbind,lapply(1:nfac, 
                                                            function(k)
                                                            {
                                                              unlist(lapply(levs, 
                                                                            function(lev,k )lev[k], 
                                                                            k = k))
                                                            })))
              LSD.dat <- cbind(LSD.dat, alldiffs.obj$LSD$meanLSD)
              names(LSD.dat) <- c(LSDby, "meanLSD")
              alldiffs.obj$predictions <- merge(alldiffs.obj$predictions, LSD.dat, 
                                                all.x = TRUE, sort = FALSE)
              alldiffs.obj$predictions <- within(alldiffs.obj$predictions, 
                                                 { 
                                                   lower.halfLeastSignificant.limit <- 
                                                     alldiffs.obj$predictions[["predicted.value"]] - 
                                                     0.5 * alldiffs.obj$predictions$meanLSD
                                                   upper.halfLeastSignificant.limit <- 
                                                     alldiffs.obj$predictions[["predicted.value"]] + 
                                                     0.5 * alldiffs.obj$predictions$meanLSD
                                                 })
              alldiffs.obj$predictions <- alldiffs.obj$predictions[, -match("meanLSD",
                                                                            names(alldiffs.obj$predictions))]
            }
          } else
          {
            if (avLSD == "per.prediction")
            {
              if (nLSD != nrow(alldiffs.obj$predictions))
                stop("The numbers of LSDs and predicted values are not equal for meanLSD.type per.prediction")
              if (any(na.omit(sed.range) > avsed.tolerance))
              {
                warning("The avsed.tolerance is exceeded for one or more predictions - reverting to confidence intervals")
                revert <- TRUE
              } else #plot per predictions LSD
                alldiffs.obj$predictions <- within(alldiffs.obj$predictions, 
                                                   { 
                                                     lower.halfLeastSignificant.limit <- 
                                                       alldiffs.obj$predictions[["predicted.value"]] - 
                                                       0.5 * alldiffs.obj$LSD$meanLSD
                                                     upper.halfLeastSignificant.limit <- 
                                                       alldiffs.obj$predictions[["predicted.value"]] + 
                                                       0.5 * alldiffs.obj$LSD$meanLSD
                                                   })
            } 
          } 
        }
      }
    }
    if (int.opt == "Confidence" || revert)
      alldiffs.obj$predictions <- within(alldiffs.obj$predictions, 
                                         { lower.Confidence.limit <- alldiffs.obj$predictions[["predicted.value"]] - 
                                           qt(1-alpha/2, denom.df) * alldiffs.obj$predictions[["standard.error"]]
                                         upper.Confidence.limit <- alldiffs.obj$predictions[["predicted.value"]] + 
                                           qt(1-alpha/2, denom.df) * alldiffs.obj$predictions[["standard.error"]]
                                         })
    ks <- NA
    if ("est.status" %in% names(alldiffs.obj$predictions))
      ks <- match("est.status", names(alldiffs.obj$predictions))
    else
    {
      if ("status" %in% names(alldiffs.obj$predictions))
        ks <- match("status", names(alldiffs.obj$predictions))
    }
    if (!is.na(ks) && ks != length(names(alldiffs.obj$predictions)))
      alldiffs.obj$predictions <- alldiffs.obj$predictions[, c(1:(ks-1), (ks+1), (ks+2), ks)]
  }
  return(alldiffs.obj)
}


"allDifferences.data.frame" <- function(predictions, classify, vcov = NULL, 
                                        differences = NULL, p.differences = NULL, 
                                        sed = NULL, LSD = NULL, meanLSD.type = "overall", 
                                        LSDby = NULL, backtransforms = NULL, 
                                        response = NULL, response.title = NULL, 
                                        term = NULL, tdf = NULL, 
                                        x.num = NULL, x.fac = NULL,level.length = NA, 
                                        pairwise = TRUE, alpha = 0.05,
                                        inestimable.rm = TRUE, 
                                        sortFactor = NULL, sortWithinVals = NULL, 
                                        sortOrder = NULL, decreasing = FALSE, 
                                        ...)
#a function to do the calculations to form an alldiffs object
#takes a table of asreml predictions and forms associated statistics
#  for all pairwise differences
{ 
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = TRUE)
  
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction")
  avLSD <- AvLSD.options[check.arg.values(meanLSD.type, AvLSD.options)]
  if (!is.null(LSDby) &&  !is.character(LSDby))
    stop("LSDby must be a character")
  
  tempcall <- list(...)
  if ("levels.length" %in% names(tempcall))
    stop("levels.length has been deprecated - use level.length")
  
  alldiffs.obj <- as.alldiffs(predictions = predictions, 
                              vcov = vcov,
                              differences = differences, 
                              p.differences = p.differences, 
                              sed = sed, LSD = LSD, 
                              backtransforms = backtransforms, 
                              response = response, 
                              response.title = response.title, 
                              term = term, classify = classify, 
                              tdf = tdf)
  
  #Check alldiffs.obj
  if (pairwise && is.null(alldiffs.obj$sed))
    stop(paste("No sed supplied in alldiffs.obj \n",
               "- can obtain using sed=TRUE in predict.asreml"))
  predictions <- alldiffs.obj$predictions
  rownames(predictions) <- NULL
  #Retain only estimable predictions
  if (asr4)
    which.estim <- (predictions$est.status == "Estimable")
  else
    which.estim <- (predictions$est.status == "Estimable")
  if (inestimable.rm & sum(which.estim) != nrow(predictions))
  { 
    predictions <- predictions[which.estim, ]
    if (nrow(predictions) == 0)
      warning("There are no estimable predictions")
    #Make sure all factors have only observed levels
    predictions[1:ncol(predictions)] <- 
      lapply(1:ncol(predictions), 
             function(k, data)
             { if (is.factor(data[[k]]))
               data[[k]] <- factor(data[[k]])
             return(data[[k]])
             }, predictions)
    rownames(predictions) <- NULL
    alldiffs.obj$predictions <- predictions
    if (!is.null(alldiffs.obj$vcov))
    { 
      if (inestimable.rm)
        alldiffs.obj$vcov <- alldiffs.obj$vcov[which.estim, which.estim]
    }
    if (!is.null(alldiffs.obj$sed))
    { 
      if (inestimable.rm)
        alldiffs.obj$sed <- alldiffs.obj$sed[which.estim, which.estim]
      diag(alldiffs.obj$sed) <- NA
    }
    
    #Reset the other components to NULL
    alldiffs.obj <- as.alldiffs(predictions = alldiffs.obj$predictions, 
                                vcov = alldiffs.obj$vcov, 
                                differences = NULL, 
                                p.differences = NULL, 
                                sed = alldiffs.obj$sed, LSD = NULL, 
                                backtransforms = backtransforms, 
                                response = response, 
                                response.title = response.title, 
                                term = term, classify = classify, 
                                tdf = tdf)
    attr(alldiffs.obj, which = "meanLSD") <- NULL
  }
  response <- as.character(attr(alldiffs.obj, which = "response"))
  
  #Sort if sortFactor set
  if (!is.null(sortFactor))
    alldiffs.obj <- sort(alldiffs.obj, decreasing = decreasing, sortFactor = sortFactor, 
                         sortWithinVals = sortWithinVals, sortOrder = sortOrder)
  
  #Ensure that predictions and other components are in standard order for the classify
  class <- unlist(strsplit(classify, ":", fixed = TRUE))
  if (!all(class == names(alldiffs.obj$predictions)[1:length(class)]))
  {
    rest <- names(alldiffs.obj$predictions)[(length(class)+1):ncol(alldiffs.obj$predictions)]
    alldiffs.obj$predictions <- cbind(alldiffs.obj$predictions[class],
                                      alldiffs.obj$predictions[rest])
    rownames(alldiffs.obj$predictions) <- NULL
  }
  ord <- do.call(order, alldiffs.obj$predictions)
  alldiffs.obj$predictions <- alldiffs.obj$predictions[ord,]
  predictions <- alldiffs.obj$predictions
  pred.labs <- makePredictionLabels(alldiffs.obj$predictions, classify, response,
                                    x.num = x.num, x.fac = x.fac, 
                                    level.length = level.length)
  #alldiffs.obj$predictions <- pred.labs$predictions
  pred.lev <- pred.labs$pred.lev
  if (!is.null(alldiffs.obj$backtransforms))
    alldiffs.obj$backtransforms <- alldiffs.obj$backtransforms[ord,]
  if (!is.null(alldiffs.obj$differences))
  {
    alldiffs.obj$differences <- alldiffs.obj$differences[ord,ord]
    colnames(alldiffs.obj$differences) <- rownames(alldiffs.obj$differences) <- pred.lev
  }
  if (!is.null(alldiffs.obj$p.differences))
  {
    alldiffs.obj$p.differences <- alldiffs.obj$p.differences[ord,ord]
    colnames(alldiffs.obj$p.differences) <- rownames(alldiffs.obj$p.differences) <- pred.lev
  }
  if (!is.null(alldiffs.obj$vcov))
  {
    alldiffs.obj$vcov <- alldiffs.obj$vcov[ord,ord]
    colnames(alldiffs.obj$vcov) <- rownames(alldiffs.obj$vcov) <- pred.lev
  }
  if (!is.null(alldiffs.obj$sed) && length(pred.lev) > 1)
  {
    alldiffs.obj$sed <- alldiffs.obj$sed[ord,ord]
    colnames(alldiffs.obj$sed) <- rownames(alldiffs.obj$sed) <- pred.lev
  }
  
  #Retain variance matrix, if vcov is not NULL
  if (!is.null(alldiffs.obj$vcov))
  { 
    if (ncol(alldiffs.obj$vcov) != length(pred.lev) | nrow(alldiffs.obj$vcov) != length(pred.lev))
      stop(paste("Dimensions of variance matrix not equal to \n",
                 "the number of observed levels combinations of the factors"))
    dimnames(alldiffs.obj$vcov) <- list(pred.lev, pred.lev)
  } 
  
  #Form all pairwise differences, if not present and store
  if (is.null(alldiffs.obj$differences) & pairwise)
  { 
    pred.diff <- outer(predictions$predicted.value, predictions$predicted.value, "-")
    if (nrow(alldiffs.obj$sed) != nrow(pred.diff) | 
        ncol(alldiffs.obj$sed) != ncol(pred.diff))
      stop("The matrix of pairwise differences and sed are not conformable")
    if (ncol(alldiffs.obj$sed) != length(pred.lev) | nrow(alldiffs.obj$sed) != length(pred.lev))
      stop(paste("Dimensions of differences and sed not equal to \n",
                 "the number of observed levels combinations of the factors"))
    dimnames(pred.diff) <- list(pred.lev, pred.lev)
    dimnames(alldiffs.obj$sed) <- list(pred.lev, pred.lev)
  } else
  { 
    if (!pairwise)
    { 
      alldiffs.obj["differences"] <- list(NULL)
      alldiffs.obj["p.differences"] <- list(NULL)
    }
  }
  
  #Check if tdf available
  denom.df <- attr(alldiffs.obj, which = "tdf")
  if (is.null(denom.df))
    warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                  "- p-values and LSDs not calculated"))
  else
  { 
    #calculate p-values, if not present
    if (is.null(alldiffs.obj$p.differences) & pairwise)
    { 
      p.diff <- abs(pred.diff)/alldiffs.obj$sed
      p.diff <- 2*pt(p.diff, df = denom.df, lower.tail = FALSE)
      alldiffs.obj$differences <- pred.diff
      alldiffs.obj$p.differences <- p.diff
    }
    
    #calculate LSDs, if not present
    if (is.null(alldiffs.obj$LSD) && pairwise && (nrow(alldiffs.obj$predictions) != 1))
    { 
      t.value = qt(1-alpha/2, denom.df)
      if (avLSD == "overall")
      {
        minLSD <- t.value * min(alldiffs.obj$sed, na.rm = TRUE)
        maxLSD <- t.value * max(alldiffs.obj$sed, na.rm = TRUE)
        meanLSD <- t.value * sqrt(mean(alldiffs.obj$sed * alldiffs.obj$sed, 
                                       na.rm = TRUE))
      } else 
      {
        if (avLSD == "factor.combinations") #factor.combinations
        {
          if (is.null(LSDby))
            stop("Need to specify factors using LSDby for meanLSD.typ = factor.combinations")
          LSDs <- sliceLSDs(alldiffs.obj, by = LSDby, t.value = t.value, alpha = alpha)
          meanLSD <- LSDs$meanLSD
          names(meanLSD) <- rownames(LSDs)
          minLSD <- LSDs$minLSD
          names(minLSD) <- rownames(LSDs)
          maxLSD <- LSDs$maxLSD
          names(maxLSD) <- rownames(LSDs)
        } else #per.prediction
        {
          meanLSD <- t.value * sqrt(apply(alldiffs.obj$sed*alldiffs.obj$sed, 
                                          FUN = mean, MARGIN = 1, na.rm = TRUE))
          maxLSD <- t.value * apply(alldiffs.obj$sed, FUN = max, MARGIN = 1, na.rm = TRUE)
          minLSD <- t.value * apply(alldiffs.obj$sed, FUN = min, MARGIN = 1, na.rm = TRUE)
        }
      }
      alldiffs.obj$LSD <- data.frame(minLSD  = minLSD, 
                                     meanLSD = meanLSD, 
                                     maxLSD = maxLSD)
      attr(alldiffs.obj, which = "meanLSD") <- meanLSD
    }
  }
  
  return(alldiffs.obj)
}

"addBacktransforms.alldiffs" <- function(alldiffs.obj, 
                                        transform.power = 1, offset = 0, scale = 1)
{  
  #Add backtransforms if there has been a transformation
  backtransforms <- NULL
  if (nrow(alldiffs.obj$predictions) > 0 && (transform.power != 1 || offset != 0 || scale != 1))
  { 
    denom.df <- attr(alldiffs.obj, which = "tdf")
    if (is.null(denom.df))
      warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                    "- p-values and LSDs not calculated"))
    backtransforms <- alldiffs.obj$predictions
    kp <- match("predicted.value", names(backtransforms))
    #Check if LSD used for predictions and so need to compute CIs
    if ((strsplit(names(backtransforms)[kp+2], ".", 
                  fixed=TRUE))[[1]][2] == "halfLeastSignificant")
    { 
      names(backtransforms)[kp+2] <- "lower.Confidence.limit" 
      names(backtransforms)[kp+3] <- "upper.Confidence.limit" 
      backtransforms <- within(backtransforms, 
                               { 
                                 lower.Confidence.limit <- alldiffs.obj$predictions[["predicted.value"]] - 
                                   qt(1-alpha/2, denom.df) * alldiffs.obj$predictions[["standard.error"]]
                                 upper.Confidence.limit <- alldiffs.obj$predictions[["predicted.value"]] + 
                                   qt(1-alpha/2, denom.df) * alldiffs.obj$predictions[["standard.error"]]
                               })
    }
    names(backtransforms)[match("predicted.value", names(backtransforms))] <- 
      "backtransformed.predictions"
    #Backtransform predictions and intervals for power transformation
    if (transform.power == 0)
    { 
      backtransforms$backtransformed.predictions <- 
                                 exp(backtransforms$backtransformed.predictions)
      backtransforms[[kp+2]] <- exp(backtransforms[[kp+2]])
      backtransforms[[kp+3]] <- exp(backtransforms[[kp+3]])
    } else
      if (transform.power != 1)
      { 
        backtransforms$backtransformed.predictions <- 
          backtransforms$backtransformed.predictions^(1/transform.power)
        backtransforms[[kp+2]] <- backtransforms[[kp+2]]^(1/transform.power)
        backtransforms[[kp+3]] <- backtransforms[[kp+3]]^(1/transform.power)
      } 
    #Backtransform for offset and scale
    if (offset !=0 || scale != 1)
    { 
      backtransforms$backtransformed.predictions <- 
        (backtransforms$backtransformed.predictions - offset)/scale
      backtransforms[[kp+2]] <- (backtransforms[[kp+2]] - offset)/scale
      backtransforms[[kp+3]] <- (backtransforms[[kp+3]] - offset)/scale
    }
    #Set standard.error to missing if a power transformation has been used
    if (transform.power != 1)
    {
      ks <- match("standard.error", names(backtransforms))
      backtransforms[[ks]] <- NA
    } else
    {
      if (scale != 1)
      {
        ks <- match("standard.error", names(backtransforms))
        backtransforms[[ks]] <- backtransforms[[ks]] / scale
      }
    }
    alldiffs.obj$backtransforms <- backtransforms
  }
  return(alldiffs.obj)
}

"reorderClassify.alldiffs" <- function(alldiffs.obj, newclassify, 
                                       sortFactor = NULL, sortWithinVals = NULL, 
                                       sortOrder = NULL, decreasing = FALSE, 
                                       ...)
{
  alldiffs.obj <- allDifferences(alldiffs.obj$predictions, classify = newclassify, 
                                 vcov = alldiffs.obj$vcov,
                                 differences = alldiffs.obj$differences, 
                                 p.differences = alldiffs.obj$p.differences, 
                                 sed = alldiffs.obj$sed,
                                 LSD = alldiffs.obj$LSD, 
                                 backtransforms = alldiffs.obj$bakctransforms,
                                 response = attr(alldiffs.obj, which = "response"), 
                                 response.title = attr(alldiffs.obj, 
                                                       which = "response.title"),
                                 term = attr(alldiffs.obj, which = "term"), 
                                 tdf = attr(alldiffs.obj, which = "tdf"),
                                 sortFactor = sortFactor, sortOrder = sortOrder, 
                                 sortWithinVals = sortWithinVals,
                                 decreasing = decreasing, ...)
  return(alldiffs.obj)
}

"linTransform.alldiffs" <- function(alldiffs.obj, classify = NULL, term = NULL, 
                                    linear.transformation = NULL, Vmatrix = FALSE, 
                                    error.intervals = "Confidence", avsed.tolerance = 0.25, 
                                    meanLSD.type = "overall", LSDby = NULL, 
                                    response = NULL, response.title = NULL, 
                                    x.num = NULL, x.fac = NULL, 
                                    tables = "all", level.length = NA, 
                                    pairwise = TRUE, alpha = 0.05,
                                    transform.power = 1, offset = 0, scale = 1, 
                                    inestimable.rm = TRUE, 
                                    ...)
{
  #Check if want a linear transformation
  if (is.null(linear.transformation))
    warning("A linear transformation has not been specified")
  else
  {
    #Get error.intervals
    int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
    int.opt <- int.options[check.arg.values(error.intervals, int.options)]
    #Check have vcov
    if (int.opt != "none" && is.null(alldiffs.obj$vcov))
      stop("Need to have stored the variance matrix of the predictions in alldiffs.obj")
    
    #Get table option and  check if must form pairwise differences
    tab.options <- c("none", "predictions", "vcov", "backtransforms", 
                     "differences", "p.differences", "sed", "LSD", "all")
    table.opt <- tab.options[unlist(lapply(tables, check.arg.values, options=tab.options))]
    if ("all" %in% table.opt || "differences" %in% table.opt)
      pairwise <- TRUE
    if (inherits(linear.transformation, what = "matrix"))
      lintrans.type <- "matrix"
    else 
    {
      if (inherits(linear.transformation, what = "formula"))
        lintrans.type <- "submodel"
      else
        stop("linear.transformation should be either a matrix or a model")
    }
    
    #get attributes from alldiffs object
    if (is.null(response))
      response <- attr(alldiffs.obj, which = "response")
    if (is.null(response.title))
    {
      response.title <- attr(alldiffs.obj, which = "response.title")
      response.title <- paste(response.title, "transform(s)", sep = " ")
    }
    if (is.null(term))
      term <- attr(alldiffs.obj, which = "term")
    if (is.null(classify))
      classify <- attr(alldiffs.obj, which = "classify")
    denom.df <- attr(alldiffs.obj, which = "tdf")
    if (is.null(denom.df))
      warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                    "- p-values and LSDs not calculated"))
    
    #Project predictions on submodel, if required
    if (lintrans.type == "submodel")
    {
      #Form projector on predictions for submodel
      suppressWarnings(Q <- pstructure(linear.transformation, grandMean = TRUE, 
                                       orthogonalize = "eigen", aliasing.print = FALSE, 
                                       data = alldiffs.obj$predictions)$Q)
      Q.submod <- Q[[1]]
      if (length(Q) > 1)
        for (k in 2:length(Q))
          Q.submod <- Q.submod + Q[[k]]
      Q.submod <- projector(Q.submod)
      
      #Check that submodel is a subspace of the classify space
      Q <- pstructure(as.formula(paste("~", classify, sep = " ")), 
                      grandMean = TRUE, data = alldiffs.obj$predictions)$Q
      if (any(abs(Q.submod %*% projector(Q[[1]] + Q[[2]]) - Q.submod) > 1e-08))
        stop("Model space for", linear.transformation, 
             " is not a subspace of the space for the classify ", classify)
      
      #Form predictions projected onto submodel
      lintrans <- alldiffs.obj$predictions
      lintrans$predicted.value <- as.vector(Q.submod %*% lintrans$predicted.value)
      
      # Calculate standard errors andthe variance matrix for differences between predictions
      if (!is.null(alldiffs.obj$vcov))
      {
        lintrans.vcov <- Q.submod %*% alldiffs.obj$vcov %*% Q.submod
        lintrans$standard.error <- as.vector(sqrt(diag(lintrans.vcov)))
        n <- nrow(lintrans.vcov)
        lintrans.sed <- matrix(rep(diag(lintrans.vcov), each = n), nrow = n) + 
          matrix(rep(diag(lintrans.vcov), times = n), nrow = n) - 
          2 * lintrans.vcov
        lintrans.sed <- sqrt(lintrans.sed)  
      } else
      {
        lintrans$standard.error <- NA
        lintrans.sed <- NULL  
      }
      
      #Form alldiffs object for linear transformation
      if (!Vmatrix)
        lintrans.vcov <- NULL
      diffs <- allDifferences(predictions = lintrans, vcov = lintrans.vcov, 
                              sed = lintrans.sed, meanLSD.type = meanLSD.type, LSDby = LSDby, 
                              response = response, response.title =  response.title, 
                              term = term, classify = classify, 
                              tdf = denom.df, 
                              x.num = x.num, x.fac = x.fac,
                              level.length = level.length, 
                              pairwise = pairwise, 
                              inestimable.rm = inestimable.rm, 
                              alpha = alpha)
    } else  #Form estimated linear combinations from predictions using a matrix
    {
      #check that the number of predictions conform
      if ((ncol(linear.transformation) != nrow(alldiffs.obj$predictions)))
        stop("The number of columns in linear.transformation is not equal to ", 
             "the number of estimable predictions")
      if (!is.null(rownames(linear.transformation)))
        lintrans <- data.frame(X = rownames(linear.transformation))
      else
        lintrans <- data.frame(X = 1:nrow(linear.transformation))
      names(lintrans) <- "Combination"
      lintrans$Combination <- factor(lintrans$Combination, levels = lintrans$Combination)
      lintrans$predicted.value <- as.vector(linear.transformation %*% 
                                              alldiffs.obj$predictions$predicted.value)
      lintrans$est.status <- "Estimable"
      lintrans$est.status[is.na(lintrans$predicted.value)] <- "Aliased"
      
      # Calculate standard errors andthe variance matrix for differences between predictions
      if (!is.null(alldiffs.obj$vcov))
      {
        lintrans.vcov <- linear.transformation %*% alldiffs.obj$vcov %*% t(linear.transformation)
        lintrans$standard.error <- as.vector(sqrt(diag(lintrans.vcov)))
        n <- nrow(lintrans.vcov)
        lintrans.sed <- matrix(rep(diag(lintrans.vcov), each = n), nrow = n) + 
          matrix(rep(diag(lintrans.vcov), times = n), nrow = n) - 2 * lintrans.vcov
        lintrans.sed <- sqrt(lintrans.sed)  
      } else
      {
        lintrans$standard.error <- NA
        lintrans.sed <- NULL  
      }
      
      #Form alldiffs object for linear transformation
      if (!Vmatrix)
        lintrans.vcov <- NULL
      diffs <- allDifferences(predictions = lintrans, vcov = lintrans.vcov, 
                              sed = lintrans.sed, meanLSD.type = meanLSD.type, LSDby = LSDby, 
                              response = response, response.title =  response.title, 
                              term = classify, classify = "Combination", 
                              tdf = denom.df, 
                              x.num = x.num, x.fac = x.fac,
                              level.length = level.length, 
                              pairwise = pairwise, 
                              inestimable.rm = inestimable.rm, 
                              alpha = alpha)
    }
    
    #Add lower and upper uncertainty limits
    diffs <- redoErrorIntervals.alldiffs(diffs, error.intervals = error.intervals,
                                         alpha = alpha, avsed.tolerance = avsed.tolerance,
                                         meanLSD.type = meanLSD.type, LSDby = LSDby)
    
    #Add backtransforms if there has been a transformation
    diffs <- addBacktransforms.alldiffs(alldiffs.obj = diffs, 
                                        transform.power = transform.power, 
                                        offset = offset, scale = scale)

    #Outut tables according to table.opt
    if (!("none" %in% table.opt))
      print(diffs, which = table.opt)
  }
  return(diffs)
}
