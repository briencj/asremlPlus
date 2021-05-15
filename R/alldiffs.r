#alldiffs functions

"is.predictions.frame" <- function(object)
{
  inherits(object, "predictions.frame") && inherits(object, "data.frame")
}


"validPredictionsFrame" <- function(object)
{
  ispredframe <- TRUE 
  #Check that is a data.frame
  if (!is.data.frame(object))
  {
    ispredframe[1] <- FALSE
    ispredframe <- c(ispredframe, 
                     "\n  Predictions.frame is not a data.frame")
  }
  #Check have appropriate columns
  if (!any(c("predicted.value", "backtransformed.predictions") %in% colnames(object)) || 
      !all(c("standard.error", "est.status") %in% colnames(object)))
  {
    ispredframe[1] <- FALSE
    ispredframe <- c(ispredframe, 
                     "\n  Predictions.frame does not include the expected column names",
                     paste("\n  (must be predicted.value or backtransformed.predictions, standard.error and est.status\n",
                           "- not std.error and status)"))
  }    
  if (length(ispredframe) > 1)
    ispredframe[1] <- "Error in validPredictionsFrame : "
  return(ispredframe)
}

"as.predictions.frame" <- function(data, predictions = NULL, se = NULL, 
                                   est.status = NULL, interval.type = NULL, 
                                   interval.names = NULL)
{
  #Check interval.typ argument
  int.type <-NULL
  if (!is.null(interval.type))
  {
    options <- c("CI", "SE", "halfLSD")
    int.type <- options[check.arg.values(interval.type, options)]
  }
  
  ## Modify data to be compatible with a predictions.frame
  if (!is.null(predictions))
  {
    if (!any(c("predicted.value", "backtransformed.predictions") %in% names(data)))
      names(data)[match(predictions, names(data))] <- c("predicted.value")
  }
  if (!is.null(se))
  {
    if (!("standard.error" %in% names(data)))
      names(data)[match(se, names(data))] <- c("standard.error")
  }
  if (!is.null(est.status))
  {
    if (!("est.status" %in% names(data)))
      names(data)[match(est.status, names(data))] <- c("est.status")
  }
  if (!is.null(int.type))
  {
    if (length(interval.names) != 2)
      stop("The number of interval names does not equal 2")
    if (int.type == "SE")
    {
      if (!all( c("lower.StandardError.limit", 
                  "upper.StandardError.limit") %in% names(data)))
        int.names <- c("lower.StandardError.limit", "upper.StandardError.limit")
    } else if (int.type == "CI")
    {
      if (!all( c("lower.Confidence.limit", "upper.Confidence.limit") %in% names(data)))
        int.names <- c("lower.Confidence.limit", "upper.Confidence.limit")
    } else if (int.type == "halfLSD")
    {
      if (!all( c("lower.halfLeastSignificant.limit", 
                  "upper.halfLeastSignificant.limit") %in% names(data)))
        int.names <- c("lower.halfLeastSignificant.limit", "upper.halfLeastSignificant.limit")
    } 
    names(data)[match(interval.names, names(data))] <- int.names
  }
  if (!("est.status" %in% names(data)))
  {
    data$est.status <- "Estimable"
    if ("predicted.value" %in% names(data))
      data$est.status[is.na(data$predicted.value)] <- "Aliased"
    else if ("backtransformed.predictions" %in% names(data))
    {
      data$est.status[is.na(data$backtransformed.predictions)] <- "Aliased"
    }
  }
  class(data) <- c("predictions.frame", "asreml.predict", "data.frame")
  # if ("asreml.predict" %in% class(data))
  #   class(data) <- c("predictions.frame", "asreml.predict", "data.frame")
  # else
  #   class(data) <- c("predictions.frame", "data.frame")
  
  #Check that have valid predictions.frame
  validpframe <- validPredictionsFrame(data)  
  if (is.character(validpframe))
    stop(validpframe)
  
  return(data)
}

setOldClass("predictions.frame")

#Form an alldiffs object from supplied component objects
#A function that constructs an alldiffs object without the validity check
"makeAlldiffs" <- function(predictions, vcov = NULL, 
                           differences = NULL, p.differences = NULL, 
                           sed = NULL, LSD = NULL, backtransforms = NULL, 
                           response = NULL, response.title = NULL, 
                           term = NULL, classify = NULL, 
                           tdf = NULL, sortFactor = NULL, sortOrder = NULL)
{
  #Check arguments
  if (!is.null(sed))
    sed <- as.matrix(sed)
  if (!is.null(vcov))
    vcov <- as.matrix(vcov)
  npred <- nrow(predictions)
  if ((!is.null(differences) && !("matrix" %in% class(differences))) ||
      (!is.null(p.differences) && !("matrix" %in% class(p.differences))) || 
      (!is.null(sed) && !("matrix" %in% class(sed))) || 
      (!is.null(vcov) && !("matrix" %in% class(vcov))))
    warning("At least one of differences, p.differences, sed and vcov is not of type matrix")
  #Check dimensions
  if (!all(unlist(lapply(list(differences, p.differences, sed, vcov), 
                         function(comp, npred) 
                         {
                           dimsOK <- TRUE
                           if (is.null(comp))
                           {
                             if (any(dim(comp) != npred))
                               dimsOK <- FALSE
                           }
                           return(dimsOK)
                         }, npred = npred))))
    stop("At least one of differences, p.differences, sed or vcov is not conformable with predictions")
  #Check backtransforms
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
  if (!is.null(LSD))
    attr(predictions, which = "meanLSD") <- LSD$meanLSD
  else
    attr(predictions, which = "meanLSD") <- NA
  
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


"as.alldiffs" <- function(predictions, vcov = NULL, 
                          differences = NULL, p.differences = NULL, 
                          sed = NULL, LSD = NULL, backtransforms = NULL, 
                          response = NULL, response.title = NULL, 
                          term = NULL, classify = NULL, 
                          tdf = NULL, sortFactor = NULL, sortOrder = NULL)
{ 
  #Change asreml4 names to asreml3 names
  predictions <- as.predictions.frame(predictions, se = "std.error", est.status = "status")

  p <- makeAlldiffs(predictions = predictions, vcov = vcov, 
                    differences = differences, p.differences = p.differences, 
                    sed = sed, LSD = LSD, backtransforms = backtransforms, 
                    response = response, response.title = response.title, 
                    term = term, classify = classify, 
                    tdf = tdf, sortFactor = sortFactor, sortOrder = sortOrder)
  
  return(p)
}

"is.alldiffs" <- function(object)
{
  inherits(object, "alldiffs")
}

"validAlldiffs" <- function(object)
{
  isalldiff <- TRUE 
  #Check have only legal attributes
  if (!all(names(attributes(object)) %in% c("names", "class", "response", "response.title",
                                            "classify", "term", "tdf",
                                            "sortFactor", "sortOrder", 
                                            "meanLSD", "meanLSD.type", "LSDby")))
  {
    isalldiff[1] <- FALSE
    isalldiff <- c(isalldiff, 
                   paste("\n  An unexpected attribute is present in"), 
                   deparse(substitute(object)))
  }
  #Check class
  if (!is.alldiffs(object))
  {
    isalldiff[1] <- FALSE 
    isalldiff <- c(isalldiff, paste("\n ", deparse(substitute(object)),
                                    "is not of class alldiffs",
                                    "(use class<- to assign to alldiffs class)"))
  }
  #Check have all components
  if (!all(c("predictions", "vcov", "differences", "p.differences", "sed", "LSD",
             "backtransforms") %in% names(object)))
  {
    warning("Not all of predictions, vcov, differences, p.differences, sed, LSD and
             backtransforms are present in ", deparse(substitute(object)))
   }
  #Check predictions frame
  valpred <- validPredictionsFrame(object$predictions)
  if (is.character(valpred))
  {
    isalldiff[1] <- FALSE
    isalldiff <- c(isalldiff, 
                   valpred[2:length(valpred)],
                   paste("\n  Predictions frame is the predictions component of", 
                         deparse(substitute(object))))
  }    
  #Check that classify variables are the first variables in classify order 
  #in the predictions.frame
  classify <- attr(object, which = "classify")
  if (is.null(classify))
  {
    isalldiff[1] <- FALSE
    isalldiff <- c(isalldiff, 
                   paste("\n ",deparse(substitute(object)),
                         "does not have a classify attribute",
                         "\n  (if appropriate, specify classify argument when creating",
                         deparse(substitute(object)),")"))
  } else
  {
    class <- unlist(strsplit(classify, ":", fixed = TRUE))
    if (!all(class == names(object$predictions)[1:length(class)]))
    {
      isalldiff[1] <- FALSE
      isalldiff <- c(isalldiff, 
                     paste("\n  Initial columns of predictions component in", 
                           deparse(substitute(object)), "are not the classify",
                           "variables in the same order as in the classify"))
    }
  }
  #Check components that should be matrices
  npred <- nrow(object$predictions)
  if ((!is.null(object$differences) && !("matrix" %in% class(object$differences))) ||
      (!is.null(object$p.differences) && !("matrix" %in% class(object$p.differences))) || 
      (!is.null(object$sed) && !("matrix" %in% class(object$sed))) || 
      (!is.null(object$vcov) && !("matrix" %in% class(object$vcov))))
  {
    isalldiff[1] <- FALSE 
    isalldiff <- c(isalldiff, 
                   paste("\n  At least one of differences, p.differences, sed and vcov in", 
                         deparse(substitute(object)), "is not of type matrix"))
  }
  #Check dimensions
  if (!all(unlist(lapply(c("differences", "p.differences", "sed", "vcov"), 
                         function(comp, object, npred) 
                         {
                           dimsOK <- TRUE
                           if (is.null(object[[comp]]))
                           {
                             if (any(dim(comp) != npred))
                               dimsOK <- FALSE
                           }
                           return(dimsOK)
                         }, object = object, npred = npred))))
  {
    isalldiff[1] <- FALSE 
    isalldiff <- c(isalldiff, 
                   paste("\n  At least one of differences, p.differences, sed and vcov in", 
                         deparse(substitute(object)), 
                         "is not conformable with predictions"))
  }
  #Check backtransforms, if present
  if (!is.null(object$backtransforms))
  { 
    if (!("backtransformed.predictions" %in% colnames(object$backtransforms))) 
    {
      isalldiff[1] <- FALSE
      isalldiff <- c(isalldiff, paste("\n  Backtransforms argument does not include a column",
                                      "named backtransformed.predictions"))
    }
    if (npred != nrow(object$backtransforms))
    {
      isalldiff[1] <- FALSE
      isalldiff <- c(isalldiff, paste("\n  Backtransforms do not contain the same number of rows",
                                      "as the predictions"))
    }
  }
  #ensure diag of sed is NA
  if (!is.null(object$sed))
    if (!all(is.na(diag(object$sed))))
    {
      isalldiff[1] <- FALSE 
      isalldiff <- c(isalldiff, 
                     paste("\n  Not all diagonal elements of the sed component of", 
                           deparse(substitute(object)), "are NA"))
    }
  if (length(isalldiff) > 1)
    isalldiff[1] <- paste("Error(s) in validAlldiffs(", deparse(substitute(object)), ") : ")
  return(isalldiff)
}

setOldClass("alldiffs")

"print.predictions.frame" <- function(x, title = NULL,  
                                      which.predictions = c("title", "heading", "table"), 
                                      colourise = FALSE, ...)
{
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = FALSE)
  
  options <- c("title", "heading", "table", "all")
  opt <- options[unlist(lapply(which.predictions, check.arg.values, options=options))]
  
  if (any(c("title", "all") %in% opt))
  {
    if (!is.null(title))
      cat(title)
    else
      cat("\n\n#### Predictions\n\n")
  }
  if (!any(c("table", "all") %in% opt))
  {
    if ("heading" %in% opt)
    {
      hd <- attr(x, which = "heading")
      for (i in 1:length(hd))
        cat(hd[i],"\n")
    }
  } else #print out the table, possibly with the heading
  {
    if (any(c("heading", "all") %in% opt) && !is.null(asr4) && asr4)
    {
      asr.col <- asreml::asreml.options()$colourise
      # if (xor(colourise,asr.col))
      #   asreml::asreml.options(colourise = colourise)
      if (length(asr.col) == 0)
      {
        if (colourise) 
          asreml::asreml.options(colourise = colourise)
      } else 
        if (xor(colourise,asr.col))
          asreml::asreml.options(colourise = colourise)
      class(x) <- c("asreml.predict", "data.frame")
      print(x, ...)
      asreml::asreml.options(colourise = asr.col)
    } else
    {
      if (!any(c("heading", "all") %in% opt))
      {
        class(x) <- class(x)[-match("asreml.predict", class(x))]
      } else
      {
        hd <- attr(x, which = "heading")
        for (i in 1:length(hd))
          cat(hd[i],"\n")
      }
      if (any(c("table", "all") %in% opt))
        print.data.frame(x, ...)
    }
  }
  invisible()
}

"print.alldiffs" <- function(x, which = "all", colourise = FALSE, ...)
{ 
  asr4 <- isASRemlVersionLoaded(4, notloaded.fault = FALSE)

  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(x)  
  if (is.character(validalldifs))
   stop(validalldifs)
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
      } else
        tt <- "\n\n#### Predictions\n\n"
      print.predictions.frame(x$predictions, title = tt, colourise = colourise, ...)
    } 
    if (!is.null(x$LSD))
    { 
      sed.range <- abs(x$LSD$minLSD - x$LSD$maxLSD) / x$LSD$meanLSD
      sed.range[is.nan(sed.range)] <- 0
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
  if (classify == "(Intercept)")
  {
    factors <- classify
  } else
  {
    factors <- fac.getinTerm(classify, rmfunction = TRUE)
    classify <- fac.formTerm(factors)
  }
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
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  if (is.character(validalldifs))
    stop(validalldifs)
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
    stop("The alldiffs object does not have the classify attribute set")
  class.facs <- fac.getinTerm(classify, rmfunction = TRUE)
  class.facs[fstfac] <- newfac
  class.facs <- class.facs[-c(match(factors[-1], class.facs))]
  classify <- fac.formTerm(class.facs)
  attr(object, which = "classify") <- classify
  response <- attr(object, which = "response")
  pred.labs <- makePredictionLabels(object$predictions, classify, response)
  pred.lev <- pred.labs$pred.lev
  #Set meanLSD attribute of predictions component
  predictions <- object$predictions
  if (is.null(object$LSD))
    attr(predictions, which = "meanLSD") <- NA
  else
    attr(predictions, which = "meanLSD") <- object$LSD$meanLSD
  object$predictions <- predictions
  
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

facRecast.alldiffs <- function(object, factor, newlevels = NULL, newlabels, ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  if (is.character(validalldifs))
    stop(validalldifs)
  fac <- factor
  if (any(!(fac %in% names(object$predictions))))
    stop("Some factors are not in the predictions component of object")
  if (length(fac) != 1)
    stop("Only one factor at a time")
  if (length(unique(newlevels))< length(newlevels))
    stop("The set of newlevels must be unique")
  
  #Recast the factor
  if (is.null(newlevels))
    newlevels <- levels(object$predictions[[fac]])
  object$predictions[fac] <- factor(object$predictions[[fac]], levels = newlevels, 
                                      labels = newlabels, ...)
  
  #revise the alldiffs component
  if (!is.null(object$backtransforms))
  {
    object$backtransforms[fac] <- object$predictions[fac]
  }
  
  classify <- attr(object, which = "classify")
  if (is.null(classify))
    stop("The alldiffs object does not have the classify attribute set")
  response <- attr(object, which = "response")
  pred.labs <- makePredictionLabels(object$predictions, classify, response)
  pred.lev <- pred.labs$pred.lev
  #Set meanLSD attribute of predictions component
  predictions <- object$predictions
  if (is.null(object$LSD))
    attr(predictions, which = "meanLSD") <- NA
  else
    attr(predictions, which = "meanLSD") <- object$LSD$meanLSD
  object$predictions <- predictions
  
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
  if (!is.null(newlevels))
    object <- renewClassify(object, newclassify = attr(object, which = "classify"))
  return(object)
}

facRecode.alldiffs <- function(object, factor, newlevels,  ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  if (is.character(validalldifs))
    stop(validalldifs)
  fac <- factor
  if (any(!(fac %in% names(object$predictions))))
    stop("Some factors are not in the predictions component of object")
  if (length(fac) != 1)
    stop("Only one factor at a time")
  object$predictions[fac] <- fac.recode(object$predictions[[fac]], newlevels = newlevels, ...)
  
  if (!is.null(object$backtransforms))
  {
    object$backtransforms[fac] <- object$predictions[fac]
  }
  
  classify <- attr(object, which = "classify")
  if (is.null(classify))
    stop("The alldiffs object does not have the classify attribute set")
  response <- attr(object, which = "response")
  pred.labs <- makePredictionLabels(object$predictions, classify, response)
  pred.lev <- pred.labs$pred.lev
  #Set meanLSD attribute of predictions component
  predictions <- object$predictions
  if (is.null(object$LSD))
    attr(predictions, which = "meanLSD") <- NA
  else
    attr(predictions, which = "meanLSD") <- object$LSD$meanLSD
  object$predictions <- predictions
  
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

facRename.alldiffs <- function(object, factor.names, newnames,  ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  if (is.character(validalldifs))
    stop(validalldifs)
  if (any(!(factor.names %in% names(object$predictions))))
    stop("Some factors are not in the predictions component of object")
  if (length(factor.names) != length(newnames))
    stop("The number of factor.names and newnames must be the same")
  names(object$predictions)[match(factor.names, 
                                  names(object$predictions))] <- newnames
  
  if (!is.null(object$backtransforms))
  {
    names(backtransforms$predictions)[match(newnames, 
                                            names(object$backtransforms))] <- newnames
  }
  
  classify <- attr(object, which = "classify")
  if (is.null(classify))
    stop("The alldiffs object does not have the classify attribute set")
  class.facs <- fac.getinTerm(classify, rmfunction = TRUE)
  class.facs[match(factor.names, class.facs)] <- newnames
  classify <- fac.formTerm(class.facs)
  attr(object, which = "classify") <- classify
  response <- attr(object, which = "response")
  #Set meanLSD attribute of predictions component
  predictions <- object$predictions
  if (is.null(object$LSD))
    attr(predictions, which = "meanLSD") <- NA
  else
    attr(predictions, which = "meanLSD") <- object$LSD$meanLSD
  object$predictions <- predictions
  
  return(object)
}

subset.alldiffs <- function(x, subset = rep(TRUE, nrow(x$predictions)), 
                            rmClassifyVars = NULL, ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(x)  
  if (is.character(validalldifs))
    stop(validalldifs)
  #Save attributes
  x.attr <- attributes(x)
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
    #Set meanLSD attribute of predictions component
    predictions <- x$predictions
    if (is.null(x$LSD))
      attr(predictions, which = "meanLSD") <- NA
    else
      attr(predictions, which = "meanLSD") <- x$LSD$meanLSD
    x$predictions <- predictions
    
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
    classify <- x.attr[["classify"]]
    if (is.null(classify))
      stop("The alldiffs object does not have the classify attrtibute set")
    class.vars <- fac.getinTerm(classify, rmfunction = TRUE)
    if (!(all(rmClassifyVars %in% class.vars)))
      stop("Not all the rmClassifyVars are in the classify for the alldiffs.object")
    newclass <- setdiff(class.vars, rmClassifyVars)
    if (any(table(x$predictions[newclass]) > 1))
      stop("The classify variables remaining after rmClassifyVars excluded do not uniquely index the predictions")
    # rmfac <- fac.combine(as.list(x$predictions[rmClassifyVars]))
    # if (length(unique(rmfac)) > 1)
    #   stop("The classify variables to be removed have more than one combination in the predictions; \nthe predictions cannot be unambiguosly identified.")
    # class.facs <- fac.getinTerm(classify, rmfunction = TRUE)
    # class.facs <- class.facs[!(class.facs %in% rmClassifyVars)]
    # classify <- fac.formTerm(class.facs)
    classify <- fac.formTerm(newclass)
    x$predictions <- x$predictions[, -c(match(rmClassifyVars, names(x$predictions)))]
    if (!is.null(x$backtransforms))
    {
      x$backtransforms <- x$backtransforms[, -c(match(rmClassifyVars, 
                                                      names(x$backtransforms)))]
    }
    x.attr["classify"] <- classify
    response <- x.attr[["response"]]
    if (is.null(response))
      stop("The alldiffs object does not have the response attribute set")
    pred.labs <- makePredictionLabels(x$predictions, classify, response)
    pred.lev <- pred.labs$pred.lev
    attributes(x) <- x.attr
    
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
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(x)  
  if (is.character(validalldifs))
    stop(validalldifs)
  
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
      if (nrow(val.ord) < length(levels(val.ord$X1)))
      {
        fac.levs <- levels(val.ord$X1)
        val.ord <- rbind(val.ord,
                         data.frame(X1 = fac.levs[!(fac.levs %in% val.ord$X1)],
                                    Ord = (nrow(val.ord)+1):length(fac.levs)))
      }
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
  }  
  #Get the order for the full set of predictions
  # - this work for unequal replication of sortFactor
  newlevs <-  val.ord[val.ord$Ord, sortFactor]
  x <- facRecast(x, sortFactor, newlevels = newlevs)

  #   tmp <- x$predictions[c(other.vars, sortFactor)]
  #   rownames(tmp) <- as.character(1:nrow(tmp))
  #   tmp <- tmp[do.call(order,tmp),]
  #   reord <- as.numeric(rownames(tmp))
  #   tmp <- suppressMessages(plyr::join(tmp, val.ord))
  #   tmp <- split(tmp, as.list(tmp[other.vars]))
  #   tmp <- lapply(tmp, 
  #                 function(d)
  #                 { d$Ord <- c(1:nrow(d))[order(order(d$Ord))]; return(d)})
  #   tmp <- do.call(rbind, tmp)
  #   starts <- reshape::melt(table(tmp[c(other.vars)]))
  #   if (length(other.vars) == 1)
  #     names(starts)[1] <- other.vars
  #   names(starts)[match("value", names(starts))] <- "start"
  #   starts$start <- c(0,cumsum(starts$start[1:(nrow(starts)-1)]))
  #   tmp <- suppressMessages(plyr::join(tmp, starts))
  #   tmp$Ord <- with(tmp, reord[Ord + start])
  #   if (any(is.na(tmp[sortFactor])))
  #     warning("The order for some predicted values could not be determined.\n",
  #             "Make sure that all values of sortFactor occur for the subset defined by sortWithinVals.")
  # }
  # 
  # #Order the components that are present
  # if (is.factor(x$predictions[[sortFactor]]))
  # {
  #   newlevs <- as.character(levels(x$predictions[[sortFactor]]))[val.ord$Ord]
  #   x$predictions[sortFactor] <- factor(x$predictions[[sortFactor]], 
  #                                       levels = newlevs)
  #   if (!is.null(x$backtransforms))
  #     x$backtransforms[sortFactor] <- factor(x$backtransforms[[sortFactor]], 
  #                                            levels = newlevs)
  # }
  # x$predictions <- x$predictions[tmp$Ord,]
  # #Set meanLSD attribute of predictions component
  # predictions <- x$predictions
  # if (is.null(x$LSD))
  #   attr(predictions, which = "meanLSD") <- NA
  # else
  #   attr(predictions, which = "meanLSD") <- x$LSD$meanLSD
  # x$predictions <- predictions
  # 
  # if (!is.null(x$vcov))
  #   x$vcov <- x$vcov[tmp$Ord, tmp$Ord]
  # if (!is.null(x$backtransforms))
  #   x$backtransforms <- x$backtransforms[tmp$Ord,]
  # if (!is.null(x$differences))
  #   x$differences <- x$differences[tmp$Ord, tmp$Ord]
  # if (!is.null(x$p.differences))
  #   x$p.differences <- x$p.differences[tmp$Ord, tmp$Ord]
  # if (!is.null(x$sed))
  #   x$sed <- x$sed[tmp$Ord, tmp$Ord]
  
  #Set attributes
#  if (is.null(sortFactor)) sortFactor <- NA
  attr(x, which = "sortFactor") <- sortFactor
  attr(x, which = "sortOrder") <- newlevs
  
  return(x)
}

recalcLSD.alldiffs <- function(alldiffs.obj, meanLSD.type = "overall", LSDby = NULL, 
                               alpha = 0.05, ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
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
  #Find missing attributes in new alldiffs.obj and add them back in 
  kattr <- kattr[names(kattr)[!(names(kattr) %in% names(newattr))]]
  if (length(kattr) > 0)
  {
    newattr <- c(newattr,kattr)
    attributes(alldiffs.obj) <- newattr
  }
  return(alldiffs.obj)
}

"LSDstats" <- function(sed, t.value, zero.tolerance = 1e-10)
{
  zero.tolerance <- zero.tolerance
  ksed <- as.vector(sed)
  ksed <- na.omit(ksed *ksed)
  max.ksed <- max(ksed)
  #Retain only nonzero variances
  if (max.ksed > zero.tolerance && sum(ksed/max.ksed > zero.tolerance) > 0)
    ksed <- ksed[ksed/max.ksed > zero.tolerance]
  else if (max.ksed < zero.tolerance)
    ksed <- 0
  minLSD <- t.value * sqrt(min(ksed))
  maxLSD <- t.value * sqrt(max(ksed))
  meanLSD <- t.value * sqrt(mean(ksed))
#  stats <- cbind(minLSD, meanLSD, maxLSD)
  stats <- data.frame(minLSD = minLSD, meanLSD = meanLSD, maxLSD = maxLSD)
  return(stats)
}

#Function to calculate the LSDs for combinations of the levels of the by factor(s)
sliceLSDs <- function(alldiffs.obj, by, t.value, alpha = 0.05, zero.tolerance = 1E-04)
{
  classify <- attr(alldiffs.obj, which = "classify")
  if (!all(unlist(lapply(by, grepl, x = classify, fixed = TRUE))))
    stop("One of the elements of LSDby is not in the classify")
  
  sed <- alldiffs.obj$sed
  denom.df <- attr(alldiffs.obj, which = "tdf")
  if (is.null(denom.df))
  {
    warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                  "- p-values and LSDs not calculated"))
    LSDs <- NULL
  } else
  {
    t.value = qt(1-alpha/2, denom.df)
    #Process the by argument
    if (is.list(by))
    {
      fac.list <- by
    } else
    {
      if (class(by) == "factor")
        fac.list <- list(by)
      else
      {
        if (is.character(by))
          fac.list <- as.list(alldiffs.obj$predictions[by])
        else
          stop("by is not one of the allowed class of inputs")
      }
    }
    #Convert any non-factors and form levels combination for which a mean LSD is required
    fac.list <- lapply(fac.list, 
                       function(x) 
                       {
                         if (class(x)!="factor")
                           x <- factor(x)
                         return(x)
                       })
    fac.comb <- fac.combine(fac.list, combine.levels = TRUE)
    if (length(fac.comb) != nrow(sed))
      stop("Variable(s) in LSDby argument are not the same length as the order of the sed matrix")
    levs <- levels(fac.comb)
    combs <- strsplit(levs, ",", fixed = TRUE)
    #Get the LSDs
    LSDs <- lapply(levs, 
                   function(lev, sed, t.value)
                   {
                     krows <- lev == fac.comb
                     if (length(fac.comb[krows]) == 1)
                     {
                       warning(paste("LSD calculated for a single prediction",
                                     "- applies to two independent predictions with the same standard error"))
                       stats <- rep(t.value * sqrt(2) * 
                                       alldiffs.obj$predictions$standard.error[krows] /2,
                                     3)
                       names(stats) <- c("minLSD", "meanLSD", "maxLSD")
                     } else
                     {
                       ksed <- sed[krows, krows]
                       stats <- LSDstats(ksed, t.value)
                     }
                     return(stats)
                   }, sed = sed, t.value = t.value)
    if (!is.null(LSDs))
    {
      LSDs <- cbind(levs,
                    as.data.frame(do.call(rbind, LSDs)))
      rownames(LSDs) <- levs
    }
  }  
  return(LSDs)
}

#transform info is only passed through redoErrorIntervals; if backtransforms are required 
#then redoErrorIntervals must be called with appropriate transform info
redoErrorIntervals.alldiffs <- function(alldiffs.obj, error.intervals = "Confidence", 
                                        alpha = 0.05, avsed.tolerance = 0.25, 
                                        meanLSD.type = NULL, LSDby = NULL, ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction")
  avLSD <- AvLSD.options[check.arg.values(meanLSD.type, AvLSD.options)]
  if (length(avLSD) != 1)
    avLSD <- NULL
  if (!is.null(LSDby) &&  !is.character(LSDby))
    stop("LSDby must be a character")
  
  if (!is.na(avsed.tolerance) & (avsed.tolerance <0 | avsed.tolerance > 1))
    stop("avsed.tolerance should be between 0 and 1")
  int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
  int.opt <- int.options[check.arg.values(error.intervals, int.options)]
  
  denom.df <- attr(alldiffs.obj, which = "tdf")
  preds.hd <- attr(alldiffs.obj$predictions, which = "heading")
  
  #Remove any intervals
  if (any(grepl("lower", names(alldiffs.obj$predictions), fixed = TRUE)
          | grepl("upper", names(alldiffs.obj$predictions), fixed = TRUE)))
  {
    cols <- pmatch(c("lower.", "upper."), names(alldiffs.obj$predictions))
    cols <- cols[!is.na(cols)]
    alldiffs.obj$predictions <- alldiffs.obj$predictions[, -cols]
  }
  
  #Add lower and upper uncertainty limits
  if (int.opt != "none")
  { 
    revert <- FALSE
    if (is.na(denom.df) && c("Confidence", "halfLeastSignificant") %in% int.opt)
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
      #Make sure that the correct type of LSDs are available
      if (is.null(avLSD))
      {
        avLSD <- attr(alldiffs.obj, which = "meanLSD.type")
        if (is.null(avLSD))
        {
          avLSD <- "overall"
        }
        if (avLSD == "factor.combinations" && is.null(LSDby))
        {
          LSDby <- attr(alldiffs.obj, which = "LSDby")
          if (is.null(LSDby))
            stop("meanLSD.type is factor.combinations, but LSDby is not set")
        }
      }
      if (avLSD != "factor.combinations")
        LSDby <- NULL

      #Determine if no LSD component or the avLSD and LSDby do not match the attributes of alldiff.obj
      avLSD.diff <- attr(alldiffs.obj, which = "meanLSD.type")
      LSDby.diff <- attr(alldiffs.obj, which = "LSDby")
      avLSD.same <- TRUE
      LSDby.same <- TRUE
      if (!(is.null(avLSD.diff) & is.null(avLSD)))
      {
        if (!is.null(avLSD.diff) & !is.null(avLSD))
        {
          avLSD.same <- avLSD.diff == avLSD
          if (avLSD.same & avLSD == "factor.combinations")
          {
            if (is.null(LSDby.diff) & is.null(LSDby))
              stop("LSDby not set for meanLSD.type set to factor.combinations")
            else 
              if (!is.null(LSDby.diff) & !is.null(LSDby))
              {
                LSDby.same <- all(LSDby.diff == LSDby)
              } else
                LSDby.same <- FALSE
          }
        }
      }
      #If no LSD component or not match recalcLSDs
      if (!is.null(alldiffs.obj$LSD) | !avLSD.same | !LSDby.same)
        alldiffs.obj <- recalcLSD(alldiffs.obj, meanLSD.type = avLSD, LSDby = LSDby, 
                                  alpha = alpha, ...)
      #Calculate overall and individual sed ranges and overall mean LSD
      overall.LSDs <- LSDstats(alldiffs.obj$sed, t.value = t.value)
      rownames(overall.LSDs) <- "overall"
      overall.sed.range <- unlist(abs(overall.LSDs["maxLSD"] - overall.LSDs["minLSD"]) / 
                                    overall.LSDs["meanLSD"])
      if (is.nan(overall.sed.range))
        overall.sed.range <- 0
      overall.meanLSD <- unlist(overall.LSDs["meanLSD"])
      nLSD <- length(alldiffs.obj$LSD$meanLSD)
      sed.range <- abs(alldiffs.obj$LSD$minLSD - alldiffs.obj$LSD$maxLSD) /  alldiffs.obj$LSD$meanLSD
      sed.range[is.nan(sed.range)] <- 0
      # if (!is.na(avsed.tolerance) & overall.sed.range <= avsed.tolerance) #always plot overall LSD
      # {
      #   alldiffs.obj$predictions <- within(alldiffs.obj$predictions, 
      #                                      { 
      #                                        lower.halfLeastSignificant.limit <- 
      #                                          alldiffs.obj$predictions[["predicted.value"]] - 
      #                                          0.5 * overall.meanLSD
      #                                        upper.halfLeastSignificant.limit <- 
      #                                          alldiffs.obj$predictions[["predicted.value"]] + 
      #                                          0.5 * overall.meanLSD
      #                                      })
      # } else #process for each meanLSD.type option
      # {
      #   if (avLSD == "overall")
      #   {
      #     if (nLSD != 1)
      #       stop("There is not just one LSD for meanLSD.type overall")
      #     warning("The avsed.tolerance is exceeded - reverting to confidence intervals")
      #     revert = TRUE
      #   } else
      #process for each meanLSD.type option
      {
        if (avLSD == "overall")
        {
          if (nLSD != 1)
            stop("There is not just one LSD for meanLSD.type overall")
          rownames(alldiffs.obj$LSD) <- "overall"
          if (!is.na(avsed.tolerance) & overall.sed.range <= avsed.tolerance)
          {
            alldiffs.obj$predictions <- within(alldiffs.obj$predictions,
                                               {
                                                 if (overall.meanLSD == 0)
                                                 {
                                                   lower.halfLeastSignificant.limit <- NA
                                                   upper.halfLeastSignificant.limit <- NA
                                                 } else
                                                 {
                                                   lower.halfLeastSignificant.limit <-
                                                     alldiffs.obj$predictions[["predicted.value"]] -
                                                     0.5 * overall.meanLSD
                                                   upper.halfLeastSignificant.limit <-
                                                     alldiffs.obj$predictions[["predicted.value"]] +
                                                     0.5 * overall.meanLSD
                                                 }
                                               })
          } else
          {              
            warning("The avsed.tolerance is exceeded - reverting to confidence intervals")
            revert = TRUE
          } 
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
                  stop("LSDby is not one of the allowed class of inputs")
              }
            }
            #Form levels combination for mean LSDs
            fac.list <- lapply(fac.list, 
                               function(x) 
                               {
                                 if (class(x)!="factor")
                                   x <- factor(x)
                                 return(x)
                               })
            levs <- levels(fac.combine(fac.list, combine.levels = TRUE))
            #Check have got the correct LSDs
            if (is.null(rownames(alldiffs.obj$LSD)) | nLSD != length(levs) | 
                any(levs != rownames(alldiffs.obj$LSD)))
              stop(paste("For meanLSD.type factor.combinations, the LSD component of the alldiffs.obj", 
                         "must be a named vector of the LSDs for each combination of the factors in LSDby", 
                         sep = " "))
            if (any(na.omit(sed.range) > avsed.tolerance))
            {
              warning("The avsed.tolerance is exceeded for the factor combinations - reverting to confidence intervals")
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
              #save currrent order of predictions and use to restore after the merge
              preds <- alldiffs.obj$predictions
              preds.attr <- attributes(alldiffs.obj$predictions)
              preds$rows <- 1:nrow(preds)
              preds.nam <- names(preds)
              preds <- merge(preds, LSD.dat, all.x = TRUE, sort = FALSE)
              preds <- preds[c(preds.nam, setdiff(names(preds), preds.nam))]
              if (any(diff(preds$rows != 1)))
                preds[preds$rows, ] <- preds
              preds <- preds[, -match("rows", names(preds))]
              preds <- within(preds,                                                  
                              { 
                                lower.halfLeastSignificant.limit <- 
                                  preds[["predicted.value"]] - 0.5 * preds$meanLSD
                                upper.halfLeastSignificant.limit <- 
                                  preds[["predicted.value"]] + 0.5 * preds$meanLSD
                                if (any(preds$meanLSD == 0))
                                {
                                  lower.halfLeastSignificant.limit[preds$meanLSD == 0] <- NA
                                  upper.halfLeastSignificant.limit[preds$meanLSD == 0] <- NA
                                } 
                                  
                              })
              alldiffs.obj$predictions <- preds[, -match("meanLSD", names(preds))]
              attr(alldiffs.obj$predictions, which = "heading") <- preds.attr$heading
              class(alldiffs.obj$predictions) <- preds.attr$class
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
    klen <- length(names(alldiffs.obj$predictions))
    if (!is.na(ks) && ks != klen)
      alldiffs.obj$predictions <- alldiffs.obj$predictions[, c(1:(ks-1), (ks+1):klen, ks)]
  }
  #Set meanLSD attribute of predictions component
  predictions <- alldiffs.obj$predictions
  attributes(predictions) <- attributes(alldiffs.obj$predictions)
  if (is.null(alldiffs.obj$LSD))
    attr(predictions, which = "meanLSD") <- NA
  else
    attr(predictions, which = "meanLSD") <- alldiffs.obj$LSD$meanLSD
  alldiffs.obj$predictions <- predictions
  attributes(alldiffs.obj$predictions) <- attributes(predictions)
  attr(alldiffs.obj$predictions, which = "heading")  <- preds.hd
  
  #Add backtransforms if there has been a transformation
  if (is.null(alldiffs.obj$backtransforms))
  {
    transform.power = 1; offset <- 0; scale <- 1
    tempcall <- list(...)
    if (!("transform.power"  %in% names(tempcall)))
      transform.power <- 1
    else
      transform.power <- tempcall$transform.power
    if (!("offset"  %in% names(tempcall)))
      offset = 0
    else
      offset <- tempcall$offset
    if (!("scale"  %in% names(tempcall)))
      scale = 1
    else
      scale <- tempcall$scale
    
  } else
  {
    transform.power = attr(alldiffs.obj$backtransforms, which = "transform.power")
    offset = attr(alldiffs.obj$backtransforms, which = "offset")
    scale = attr(alldiffs.obj$backtransforms, which = "scale")
  }
  alldiffs.obj <- addBacktransforms.alldiffs(alldiffs.obj = alldiffs.obj, 
                                             transform.power = transform.power, 
                                             offset = offset, scale = scale)

  return(alldiffs.obj)
}

#allDifferences does not change Error.Intervals, 
#but adds backtransforms depending on transform info
"allDifferences.data.frame" <- function(predictions, classify, vcov = NULL, 
                                        differences = NULL, p.differences = NULL, 
                                        sed = NULL, LSD = NULL, meanLSD.type = "overall", 
                                        LSDby = NULL, backtransforms = NULL, 
                                        response = NULL, response.title = NULL, 
                                        term = NULL, tdf = NULL, 
                                        x.num = NULL, x.fac = NULL, level.length = NA, 
                                        pairwise = TRUE, alpha = 0.05,
                                        transform.power = 1, offset = 0, scale = 1, 
                                        inestimable.rm = TRUE, 
                                        sortFactor = NULL, sortWithinVals = NULL, 
                                        sortOrder = NULL, decreasing = FALSE, 
                                        ...)
#a function to do the calculations to form an alldiffs object
#takes a table of asreml predictions and forms associated statistics
#  for all pairwise differences
{ 
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction")
  avLSD <- AvLSD.options[check.arg.values(meanLSD.type, AvLSD.options)]
  if (!is.null(LSDby) &&  !is.character(LSDby))
    stop("LSDby must be a character")
  
  tempcall <- list(...)
  if ("levels.length" %in% names(tempcall))
    stop("levels.length has been deprecated - use level.length")
  
  #Change asreml4 names to asreml3 names
  predictions <- as.predictions.frame(predictions, se = "std.error", est.status = "status")
  
  alldiffs.obj <- makeAlldiffs(predictions = predictions, 
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
  if (pairwise && is.null(alldiffs.obj$sed) && is.null(alldiffs.obj$vcov))
    stop(paste("No sed or vcov supplied in alldiffs.obj \n",
               "- can obtain using sed=TRUE or vcov=TRUE in predict.asreml"))
  predictions <- alldiffs.obj$predictions
  preds.attr <- attributes(alldiffs.obj$predictions)
  rownames(predictions) <- NULL
  #Retain only estimable predictions
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
    attr(alldiffs.obj$predictions, which = "heading") <- preds.attr$heading
    class(alldiffs.obj$predictions) <- preds.attr$class
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
    alldiffs.obj <- makeAlldiffs(predictions = alldiffs.obj$predictions, 
                                 vcov = alldiffs.obj$vcov, 
                                 differences = NULL, 
                                 p.differences = NULL, 
                                 sed = alldiffs.obj$sed, LSD = NULL, 
                                 backtransforms = backtransforms, 
                                 response = response, 
                                 response.title = response.title, 
                                 term = term, classify = classify, 
                                 tdf = tdf)
    predictions <- alldiffs.obj$predictions
    preds.attr <- attributes(alldiffs.obj$predictions)
    attr(predictions, which = "meanLSD") <- NA
    alldiffs.obj$predictions <- predictions
    attributes(alldiffs.obj$predictions) <- preds.attr
  }
  response <- as.character(attr(alldiffs.obj, which = "response"))
  
  #Deal with case when have vcov, but not sed
  if (pairwise && !is.null(alldiffs.obj$vcov) && is.null(alldiffs.obj$sed))
  {
    alldiffs.obj$sed <- alldiffs.obj$vcov
    n <- nrow(alldiffs.obj$sed)
    dvcov <- diag(alldiffs.obj$sed)
    alldiffs.obj$sed <- matrix(rep(dvcov, each = n), nrow = n) + 
      matrix(rep(dvcov, times = n), nrow = n) - 2 * alldiffs.obj$sed
    alldiffs.obj$sed <- sqrt(alldiffs.obj$sed)
    diag(alldiffs.obj$sed) <- NA_real_
  }
  
  #Ensure that the columns of predictions are in the same order as the classify 
  class <- unlist(strsplit(classify, ":", fixed = TRUE))
  if (!all(class == names(alldiffs.obj$predictions)[1:length(class)]))
  {
    rest <- names(alldiffs.obj$predictions)[(length(class)+1):ncol(alldiffs.obj$predictions)]
    preds.attr <- attributes(alldiffs.obj$predictions)
    alldiffs.obj$predictions <- cbind(alldiffs.obj$predictions[class],
                                      alldiffs.obj$predictions[rest])
    rownames(alldiffs.obj$predictions) <- NULL
    attr(alldiffs.obj$predictions, which = "heading") <- preds.attr$heading
    class(alldiffs.obj$predictions) <- preds.attr$class
  }

  #Sort if sortFactor set
  if (!is.null(sortFactor))
    alldiffs.obj <- sort(alldiffs.obj, decreasing = decreasing, sortFactor = sortFactor, 
                         sortWithinVals = sortWithinVals, sortOrder = sortOrder)
  
  #Make sure that the predictions and other components are in standard order for the classify
  ord <- do.call(order, alldiffs.obj$predictions)
  alldiffs.obj$predictions <- alldiffs.obj$predictions[ord,]
  predictions <- alldiffs.obj$predictions
  pred.labs <- makePredictionLabels(alldiffs.obj$predictions, classify, response,
                                    x.num = x.num, x.fac = x.fac, 
                                    level.length = level.length)
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
    
    #Set meanLSD attribute of predictions
    zero.tolerance <- 1e-12
    if (pairwise && (nrow(alldiffs.obj$predictions) != 1))
    { 
      #calculate LSDs, if not present
      if (is.null(alldiffs.obj$LSD))
      {
        t.value = qt(1-alpha/2, denom.df)
        if (avLSD == "overall")
        {
          LSDs<- LSDstats(alldiffs.obj$sed, t.value)
          rownames(LSDs) <- "overall"
          minLSD <- LSDs["minLSD"]
          maxLSD <- LSDs["maxLSD"]
          meanLSD <- LSDs["meanLSD"]
        } else 
        {
          if (avLSD == "factor.combinations") #factor.combinations
          {
            if (is.null(LSDby))
              stop("Need to specify factors using LSDby for meanLSD.typ = factor.combinations")
            LSDs <- sliceLSDs(alldiffs.obj, by = LSDby, t.value = t.value, alpha = alpha, 
                              zero.tolerance = zero.tolerance)
            meanLSD <- LSDs$meanLSD
            names(meanLSD) <- rownames(LSDs)
            minLSD <- LSDs$minLSD
            names(minLSD) <- rownames(LSDs)
            maxLSD <- LSDs$maxLSD
            names(maxLSD) <- rownames(LSDs)
          } else #per.prediction
          {
            zero.tolerance = 1E-04
            max.var <- max(alldiffs.obj$sed*alldiffs.obj$sed, na.rm = TRUE)
            if (max.var > zero.tolerance || 
                sum(alldiffs.obj$sed*alldiffs.obj$sed/max.var < zero.tolerance) > 0)
              alldiffs.obj$sed[alldiffs.obj$sed*alldiffs.obj$sed/max.var < zero.tolerance] <- NA
            meanLSD <- t.value * sqrt(apply(alldiffs.obj$sed*alldiffs.obj$sed, 
                                            FUN = mean, MARGIN = 1, na.rm = TRUE))
            maxLSD <- t.value * apply(alldiffs.obj$sed, FUN = max, MARGIN = 1, na.rm = TRUE)
            minLSD <- t.value * apply(alldiffs.obj$sed, FUN = min, MARGIN = 1, na.rm = TRUE)
          }
        }
        alldiffs.obj$LSD <- data.frame(minLSD  = minLSD, 
                                       meanLSD = meanLSD, 
                                       maxLSD = maxLSD)
        attr(alldiffs.obj, which = "meanLSD.type") <- avLSD
        attr(alldiffs.obj, which = "LSDby") <- LSDby
      } 
    } 
  }
  #Set meanLSD attribute of predictions component
  predictions <- alldiffs.obj$predictions
  attributes(predictions) <- attributes(alldiffs.obj$predictions)
  if (is.null(alldiffs.obj$LSD))
    attr(predictions, which = "meanLSD") <- NA
  else
    attr(predictions, which = "meanLSD") <- alldiffs.obj$LSD$meanLSD
  alldiffs.obj$predictions <- predictions
  attributes(alldiffs.obj$predictions) <- attributes(predictions)
  
  #Add backtransforms if there has been a transformation
  alldiffs.obj <- addBacktransforms.alldiffs(alldiffs.obj, 
                                             transform.power = transform.power, 
                                             offset = offset, scale = scale)
  
  #Check that have a valid alldiffs object
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  
  return(alldiffs.obj)
}

"addBacktransforms.alldiffs" <- function(alldiffs.obj, transform.power = 1, 
                                         offset = 0, scale = 1, ...)
{  
  #Add backtransforms if there has been a transformation
  if (nrow(alldiffs.obj$predictions) > 0 && (transform.power != 1 || offset != 0 || scale != 1))
  { 
    denom.df <- attr(alldiffs.obj, which = "tdf")
    if (is.null(denom.df))
      warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                    "- p-values and LSDs not calculated"))
    backtransforms <- alldiffs.obj$predictions
    kp <- match("predicted.value", names(backtransforms))
    kpl <- pmatch("lower.", names(backtransforms))
    kpu <- pmatch("upper.", names(backtransforms))
    ## As of 3/4/2019 I am allowing backtransformed halfLSD intervals
    #Check if LSD used for predictions and so need to compute CIs
    # if ((strsplit(names(backtransforms)[kp+2], ".", 
    #               fixed=TRUE))[[1]][2] == "halfLeastSignificant")
    # { 
    #   names(backtransforms)[kp+2] <- "lower.Confidence.limit" 
    #   names(backtransforms)[kp+3] <- "upper.Confidence.limit" 
    #   backtransforms <- within(backtransforms, 
    #                            { 
    #                              lower.Confidence.limit <- alldiffs.obj$predictions[["predicted.value"]] - 
    #                                qt(1-alpha/2, denom.df) * alldiffs.obj$predictions[["standard.error"]]
    #                              upper.Confidence.limit <- alldiffs.obj$predictions[["predicted.value"]] + 
    #                                qt(1-alpha/2, denom.df) * alldiffs.obj$predictions[["standard.error"]]
    #                            })
    # }
    names(backtransforms)[match("predicted.value", names(backtransforms))] <- 
      "backtransformed.predictions"
    kpl <- pmatch("lower.", names(backtransforms))
    kpu <- pmatch("upper.", names(backtransforms))
    err.int <- TRUE
    if (is.na(kpl) || is.na(kpu))
      err.int <- FALSE
    #Backtransform predictions and intervals for power transformation
    if (transform.power == 0)
    { 
      backtransforms$backtransformed.predictions <- 
                                 exp(backtransforms$backtransformed.predictions)
      if (err.int)
      {
        backtransforms[[kpl]] <- exp(backtransforms[[kpl]])
        backtransforms[[kpu]] <- exp(backtransforms[[kpu]])
      }
    } else
      if (transform.power != 1)
      { 
        backtransforms$backtransformed.predictions <- 
          backtransforms$backtransformed.predictions^(1/transform.power)
        if (err.int)
        {
          backtransforms[[kpl]] <- backtransforms[[kpl]]^(1/transform.power)
          backtransforms[[kpu]] <- backtransforms[[kpu]]^(1/transform.power)
        }  
      } 
    #Backtransform for offset and scale
    if (offset !=0 || scale != 1)
    { 
      backtransforms$backtransformed.predictions <- 
        (backtransforms$backtransformed.predictions - offset)/scale
      if (err.int)
      {
        backtransforms[[kpl]] <- (backtransforms[[kpl]] - offset)/scale
        backtransforms[[kpu]] <- (backtransforms[[kpu]] - offset)/scale
      }
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
    #Set meanLSD attribute of predictions component
    attr(backtransforms, which = "meanLSD") <- NA
    attr(backtransforms, which = "transform.power") <- transform.power
    attr(backtransforms, which = "offset") <- offset
    attr(backtransforms, which = "scale") <- scale
    alldiffs.obj$backtransforms <- backtransforms
  }
  return(alldiffs.obj)
}

"renewClassify.alldiffs" <- function(alldiffs.obj, newclassify, 
                                     sortFactor = NULL, sortWithinVals = NULL, 
                                     sortOrder = NULL, decreasing = FALSE, 
                                     ...)
{
  kattr <- attributes(alldiffs.obj)
  alldiffs.obj <- allDifferences(alldiffs.obj$predictions, classify = newclassify, 
                                 vcov = alldiffs.obj$vcov,
                                 differences = alldiffs.obj$differences, 
                                 p.differences = alldiffs.obj$p.differences, 
                                 sed = alldiffs.obj$sed,
                                 LSD = alldiffs.obj$LSD, 
                                 backtransforms = alldiffs.obj$backtransforms,
                                 response = attr(alldiffs.obj, which = "response"), 
                                 response.title = attr(alldiffs.obj, 
                                                       which = "response.title"),
                                 term = attr(alldiffs.obj, which = "term"), 
                                 tdf = attr(alldiffs.obj, which = "tdf"),
                                 sortFactor = sortFactor, sortOrder = sortOrder, 
                                 sortWithinVals = sortWithinVals,
                                 decreasing = decreasing, ...)
  newattr <- attributes(alldiffs.obj)
  #Find missing attributes in new alldiffs.obj and add them back in 
  kattr <- kattr[names(kattr)[!(names(kattr) %in% names(newattr))]]
  if (length(kattr) > 0)
  {
    newattr <- c(newattr,kattr)
    attributes(alldiffs.obj) <- newattr
  }
  #Check that the newclassify uniquely indexes the predictions
  newclass.vars <- fac.getinTerm(newattr$classify, rmfunction = TRUE)
  if (any(table(alldiffs.obj$predictions[newclass.vars]) > 1))
    stop("The newclassify variables do not uniquely index the predictions")
 
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
    
    #get attributes from predictions
    preds.attr <- attributes(alldiffs.obj$predictions)
    
    #get attributes from backtransforms
    transform.power = 1; offset <- 0; scale <- 1
    if (!is.null(alldiffs.obj$backtransforms))
    {
      transform.power = attr(alldiffs.obj$backtransforms, which = "transform.power")
      offset = attr(alldiffs.obj$backtransforms, which = "offset")
      scale = attr(alldiffs.obj$backtransforms, which = "scale")
    } 
    
    #Project predictions on submodel, if required
    if (lintrans.type == "submodel")
    {
      #Check that factors in LSDby are in the formula
      term.obj <- as.terms.object(linear.transformation, alldiffs.obj)
      lintrans.fac <- rownames(attr(term.obj, which = "factor"))
      if (!all(LSDby %in% lintrans.fac))
        warning("Some factors in the LSDby are not in the linear.transformation submodel")
      
      #Form projector on predictions for submodel
      suppressWarnings(Q <- pstructure(linear.transformation, grandMean = TRUE, 
                                       orthogonalize = "eigen", #aliasing.print = FALSE, 
                                       data = alldiffs.obj$predictions)$Q)
      Q.submod <- Q[[1]]
      if (length(Q) > 1)
        for (k in 2:length(Q))
          Q.submod <- Q.submod + Q[[k]]
      Q.submod <- projector(Q.submod)
      
      #Process the classify to ensure there is a separate term for covariates
      vars <- fac.getinTerm(classify, rmfunction = TRUE)
      facs <- covs <- list()
      for (var in vars)
      {
        if (is.numeric(alldiffs.obj$predictions[[var]]))
          covs <- c(covs, list(var))
        else
          facs <- c(facs, list(var))
      }
      if (length(facs) > 0)
      {
        full.mod <- fac.formTerm(facs)
        if (length(covs) > 0)
        {
          covs <- paste(unlist(covs), collapse = " + ")
          full.mod <- paste0(full.mod,"/(",covs,")")
        }
      } else #no facs
      {
        if (length(covs) == 0)
          stop("Did not find any factors or covariates in the classify")
        full.mod <- paste(unlist(covs), collapse = " + ")
      }
      full.mod <- as.formula(paste0("~ ", full.mod))

      #Check that submodel is a subspace of the classify space
      Q <- pstructure(full.mod, grandMean = TRUE, data = alldiffs.obj$predictions)$Q
      Q.class <- Q[[1]]
      if (length(Q) > 1)
        for (k in 2:length(Q))
          Q.class <- Q.class + Q[[k]]
      Q.class <- projector(Q.class)
      
      if (any(abs(Q.submod %*% Q.class - Q.submod) > 1e-08))
        stop("Model space for ", linear.transformation, ", with ", degfree(Q.submod), 
             " DF, is not a subspace of the space for the classify ", classify, 
             ", with ", degfree(Q.class), " DF.")
      
      #Form predictions projected onto submodel
      lintrans <- alldiffs.obj$predictions
      lintrans$predicted.value <- as.vector(Q.submod %*% lintrans$predicted.value)
      
      # Calculate standard errors and the variance matrix for differences between predictions
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
      preds.attr$heading <- c(paste("The original predictions, obtained as described below, have",
                                    "\nbeen linearly transformed to form estimated marginal means.", 
                                    sep = ""),
                              preds.attr$heading)
      attr(lintrans, which = "heading") <- preds.attr$heading
      class(lintrans) <- preds.attr$class
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
      
      # Calculate standard errors and the variance matrix for differences between predictions
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
      preds.attr$heading <- c(paste("The original predictions, obtained as described below, have",
                                    "\nbeen linearly transformed.", sep = ""),
                              preds.attr$heading)
      attr(lintrans, which = "heading") <- preds.attr$heading
      class(lintrans) <- preds.attr$class
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
                                         meanLSD.type = meanLSD.type, LSDby = LSDby,
                                         transform.power = transform.power, 
                                         offset = offset, scale = scale)
    
    #Outut tables according to table.opt
    if (!("none" %in% table.opt))
      print(diffs, which = table.opt)
  }
  return(diffs)
}
