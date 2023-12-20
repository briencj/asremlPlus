#alldiffs functions
#The prime alldiffs functions are allDifferences.data.frame, redoErrorIntervals, renewClassify and addBackTransforms;
#They are basic to building an alldiffs object; most other functions are utility functions.
#
#     redoErrorIntervals Passes LSD arguments to recalcLSD; calculates intervals; sets LSD attributes for predictions; 
#             |          calls addBacktransforms that sets backtransforms attributes
#             |          (should not be called without transform arguments; calls)
#             |
#         recalcLSD      Passes LSD parameters to allDifferences; sets object attributes
#             |
#             |
#      allDifferences    Calculates LSDs using LSDstats and sliceLSDs; sets object LSD attributes; 
#                        sets predictions attribute, but not predictions LSD attr; 
#                        calls addBacktransforms that sets backtransforms attributes

"is.LSD.frame" <- function(object)
{
  inherits(object, "LSD.frame") && inherits(object, "data.frame")
}


"validLSDFrame" <- function(object)
{
  isLSDframe <- TRUE 
  #Check that is a data.frame
  if (!is.data.frame(object))
  {
    isLSDframe[1] <- FALSE
    isLSDframe <- c(isLSDframe, 
                     "\n  LSD.frame is not a data.frame")
  }
  #Check have appropriate columns
  if (!all(c("minLSD", "meanLSD", "maxLSD", "assignedLSD", "accuracyLSD") %in% colnames(object)))
  {
    isLSDframe[1] <- FALSE
    isLSDframe <- c(isLSDframe, 
                     "\n  LSD.frame does not include the expected column names",
                     paste("\n  (must be minLSD, meanLSD, maxLSD, assignedLSD and accuracyLSD"))
  }    
  if (length(isLSDframe) > 1)
    isLSDframe[1] <- "Error in validLSDFrame : "
  return(isLSDframe)
}

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

"as.predictions.frame" <- function(data, classify = NULL, predictions = NULL, se = NULL, 
                                   est.status = NULL, interval.type = NULL, 
                                   interval.names = NULL)
{
  #Check interval.type argument
  int.type <-NULL
  if (!is.null(interval.type))
  {
    options <- c("CI", "SE", "halfLSD")
    int.type <- options[check.arg.values(interval.type, options)]
  }
  
  ## Modify data to be compatible with a predictions.frame
  # Convert any characters to factors
  if (!is.null(classify))
  {
    vars <- fac.getinTerm(classify, rmfunction=TRUE)
    data[vars] <- cbind(lapply(vars, 
                               function(x, dat)
                               {
                                 var <- dat[[x]]
                                 if (is.character(var))
                                   var <- factor(var, levels = unique(var))
                                 return(var)
                               }, dat = data))
  }
  #If necessary, set the name of the column containing the predictions
  if (!is.null(predictions))
  {
    if (!any(c("predicted.value", "backtransformed.predictions") %in% names(data)))
      names(data)[match(predictions, names(data))] <- c("predicted.value")
  }
  #Check name of se column
  if (!is.null(se))
  {
    if (!("standard.error" %in% names(data)))
      names(data)[match(se, names(data))] <- c("standard.error")
  }
  #Check for est.status column
  if (!is.null(est.status))
  {
    if (!("est.status" %in% names(data)))
      names(data)[match(est.status, names(data))] <- c("est.status")
  }
  #Deal with columns containing limits for error intervals
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

  #Check that have valid predictions.frame
  validpframe <- validPredictionsFrame(data)  
  if (is.character(validpframe))
    stop(validpframe)
  
  return(data)
}

setOldClass("predictions.frame")

"print.LSDdata" <- function(x,  which.print = c("statistics", "false.pos", "false.neg"), ...)
{
  options <- c("frequencies", "distinct.vals", "statistics", "accuracy", "false.pos", "false.neg", 
               "per.pred.accuracy", "LSDmatrix", "summary", "all")
  opt <- options[unlist(lapply(which.print, check.arg.values, options=options))]
  if (all(c("summary", "all") %in% opt))
    stop("Can only specify one of summary and all for which argument")
  
  #make change to control printing
  class(x) <- c("LSDdata", "data.frame")
  
  if (any(c("frequencies", "all") %in% opt))
  {
    cat("\n\n####  Frequency distribution of LSDs \n\n")
    fr <- as.data.frame(x$frequencies)
    fr <- cbind(rownames(fr),fr)
    rownames(fr) <- NULL
    names(fr) <- c("midpoint", "frequency")
    print(fr, ...)
  }
  
  if (any(c("distinct.vals", "summary", "all") %in% opt))
  {
    cat("\n\n####  Distinct LSD values \n\n")
    print(x$distinct.vals, ...)
  }
  
  if (any(c("statistics", "summary", "all") %in% opt))
  {
    cat("\n\n####  Statistics calculated from LSD values \n\n")
    print(x$statistics, ...)
  }
  
  if (any(c("accuracy", "all") %in% opt))
  {
    cat(paste0("\n\n####  Accuracy (", attr(x, which = "LSDaccuracy"), 
               ") of statistics calculated from LSD values \n\n"))
    print(x$accuracy, ...)
  }
  
  if (any(c("false.pos", "summary", "all") %in% opt))
  {
    cat(paste0("\n\n####  False positives resulting from the use of various LSD statistics\n\n"))
    print(x$false.pos, ...)
  }
  
  if (any(c("false.neg", "summary", "all") %in% opt))
  {
    cat(paste0("\n\n####  False negatives resulting from the use of various LSD statistics\n\n"))
    print(x$false.neg, ...)
  }
  
  if (any(c("per.pred.accuracy", "all") %in% opt))
  {
    cat(paste0("\n\n####  Accuracy (", attr(x, which = "LSDaccuracy"), 
               ") for each prediction if LSD statistics are used \n\n"))
    print(x$per.pred.accuracy, ...)
  }
  
  if (any(c("LSDmatrix", "all") %in% opt))
  {
    cat("\n\n####  Matrix of all LSD values \n\n")
    print(x$LSD, ...)
  }
  
  invisible()
}

#Form an alldiffs object from supplied component objects
#A function that constructs an alldiffs object without the validity check
"makeAlldiffs" <- function(predictions, vcov = NULL, 
                           differences = NULL, p.differences = NULL, 
                           sed = NULL, LSD = NULL, backtransforms = NULL, 
                           response = NULL, response.title = NULL, 
                           term = NULL, classify = NULL, 
                           tdf = NULL, alpha = 0.05,
                           sortFactor = NULL, sortOrder = NULL)
{
  #Check arguments
  if (!is.null(sed))
    sed <- as.matrix(sed)
  if (!is.null(vcov))
    vcov <- as.matrix(vcov)
  npred <- nrow(predictions)
  if ((!is.null(differences) && !(inherits(differences, "matrix"))) ||
      (!is.null(p.differences) && !(inherits(p.differences, "matrix"))) || 
      (!is.null(sed) && !(inherits(sed, "matrix"))) || 
      (!is.null(vcov) && !(inherits(vcov, "matrix"))))
    warning("At least one of differences, p.differences, sed and vcov is not of type matrix")
  #Check dimensions
  if (!all(unlist(lapply(list(differences, p.differences, sed, vcov), 
                         function(comp, npred) 
                         {
                           dimsOK <- TRUE
                           if (!is.null(comp))
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

  p <- list(predictions = predictions, vcov = vcov, 
            differences = differences, p.differences = p.differences, sed = sed, 
            LSD = LSD, backtransforms = backtransforms)
  attr(p, which = "response") <- response
  attr(p, which = "response.title") <- response.title
  attr(p, which = "term") <- term
  attr(p, which = "classify") <- classify
  attr(p, which = "tdf") <- tdf
  attr(p, which = "alpha") <- alpha
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
                          tdf = NULL, alpha = 0.05, sortFactor = NULL, sortOrder = NULL)
{ 
  #Change asreml4 names to asreml3 names
  predictions <- as.predictions.frame(predictions, classify = classify, 
                                      se = "std.error", est.status = "status")

  p <- makeAlldiffs(predictions = predictions, vcov = vcov, 
                    differences = differences, p.differences = p.differences, 
                    sed = sed, LSD = LSD, backtransforms = backtransforms, 
                    response = response, response.title = response.title, 
                    term = term, classify = classify, 
                    tdf = tdf, alpha = alpha, sortFactor = sortFactor, sortOrder = sortOrder)
  return(p)
}

"is.alldiffs" <- function(object)
{
  inherits(object, "alldiffs")
}

#Function to rename meanLSD.type attribute to LSDtype, if it is present
renameAttr <- function(object, change.attribs)
{
  oldattribs <- names(change.attribs)
  for (oldattr in oldattribs)
  {
    if (!is.null(attr(object, which = oldattr)))
    {
      attr(object, which = change.attribs[[oldattr]]) <- attr(object, which = oldattr)
      attr(object, which = oldattr) <- NULL
    }
    
  }
  return(object)
}
renameDiffsAttr <- function(object)
{
  if (is.null(attr(object, which = "LSDtype")) && !is.null(attr(object, which = "meanLSD")))
    attr(object, which = "meanLSD") <- NULL
  object <- renameAttr(object, 
                       change.attribs = list(meanLSD.type = "LSDtype", meanLSD = "LSDvalues"))
  
  predictions <- object$predictions
  if (is.null(attr(predictions, which = "LSDtype")) && !is.null(attr(predictions, which = "meanLSD")))
    attr(predictions, which = "meanLSD") <- NULL
  predictions <- renameAttr(predictions, 
                            change.attribs = list(meanLSD.type = "LSDtype", meanLSD = "LSDvalues"))
  object$predictions <- predictions
  
  if (!is.null(attr(object, which = "meanLSD")))
    attr(object, which = "meanLSD") <- NULL
  if (!is.null(object$backtransforms))
  {
    backtransforms <- object$backtransforms
    backtransforms <- renameAttr(backtransforms, 
                                 change.attribs = list(meanLSD.type = "LSDtype"))
    object$backtransforms <- backtransforms
  }
  return(object)
}

#Function collect attributes from alldiffs object and its predictions and backtransforms components
getAllAttr.alldiffs <- function(alldiffs.obj)
{  
  kattr <- list(object = attributes(alldiffs.obj), preds = attributes(alldiffs.obj$predictions), 
                back = attributes(alldiffs.obj$backtransforms))
  return(kattr)
}

#Function compare attributes from alldiffs object and its predictions and backtransforms components with
# a previously saved set of attributes and add in those that are missing
addMissingAttr.alldiffs <- function(alldiffs.obj, kattr)
{
  newattr <- getAllAttr.alldiffs(alldiffs.obj)
  #Find missing attributes in new alldiffs.obj and add them back in 
  newattr <- mapply(function(newatt, katt)
  {
    katt <- katt[names(katt)[!(names(katt) %in% names(newatt))]]
    if (length(katt) > 0)
      newatt <- c(newatt,katt)
    return(newatt)
  }, newatt = newattr, katt = kattr, SIMPLIFY = FALSE)
  attributes(alldiffs.obj) <- newattr$object
  attributes(alldiffs.obj$predictions) <- newattr$preds
  if (!is.null(alldiffs.obj$backtransforms))
    attributes(alldiffs.obj$backtransforms) <- newattr$back
  return(alldiffs.obj)
}

#Function to check that the classify variables are the initial variables in a predictions.frame 
"getClassifyVars" <- function(classify)
{ 
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = FALSE)
  if (!is.null(asr4.2) && asr4.2)
    class.nam <- fac.getinTerm(classify) #removes parentheses from "(Intercept)"
  else
    class.nam <- unlist(strsplit(classify, "\\:"))
  return(class.nam)
}

checkClassifyVars.predictions.frame <- function(predictions, classify.names)
{
  if (!all(classify.names %in% names(predictions)))
    stop("The predictions data.frame does not have a column for each of the following variables", 
               "in the classify \n", paste0(setdiff(classify.names, predictions), collapse = ", "))
  invisible()
}

"validAlldiffs" <- function(object)
{
  isalldiff <- TRUE 
  #Check have only legal attributes
  if (!all(names(attributes(object)) %in% c("names", "class", "response", "response.title",
                                            "classify", "term", "tdf", "alpha", 
                                            "sortFactor", "sortOrder", 
                                            "meanLSD", "meanLSD.type", 
                                            "LSDtype", "LSDby", "LSDstatistic", "LSDaccuracy")))
  {
    isalldiff[1] <- FALSE
    isalldiff <- c(isalldiff, 
                   paste("\n  An unexpected attribute is present in "), 
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
    class.names <- getClassifyVars(classify)
    if (!all(class.names == names(object$predictions)[1:length(class.names)]))
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
    if ((!is.null(object$differences) && !(inherits(object$differences, "matrix"))) ||
      (!is.null(object$p.differences) && !(inherits(object$p.differences, "matrix"))) || 
      (!is.null(object$sed) && !(inherits(object$sed, "matrix"))) || 
      (!is.null(object$vcov) && !(inherits(object$vcov, "matrix"))))
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
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = FALSE)
  
  #determine factors for row and column names
  #Make sure no functions in classify
  if (classify == "(Intercept)")
  {
    if (asr4.2)
      factors <- "Intercept"
    else
      factors <- classify
  } else
  {
    factors <- fac.getinTerm(classify, rmfunction = TRUE)
    classify <- fac.formTerm(factors)
  }
  nfac <- length(factors)
  
  #Check all factors in classify are in predictions
  if (length(setdiff(factors, names(predictions))) != 0)
  { 
    if (!is.null(response))
      stop("For ",response,
           ", the following  factors in the classify argument do not have columns in alldiffs.obj$predictions:\n",
           paste0(setdiff(factors, names(predictions)), collapse = ", "))
    else
      stop("The following  factors in the classify argument do not have columns in alldiffs.obj$predictions:\n",
           paste0(setdiff(factors, names(predictions)), collapse = ", "))
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
        if (is.numeric(pred.faclist[[kk]]) || is.character(pred.faclist[[kk]]))
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
  #Deal with unsupported parameters
  tempcall <- list(...)  
  if (length(tempcall)) 
    if (any(c("LSDtype", "LSDby", "LSDstatistic", "LSDaccuracy", 
              "avsed.tolerance", "accuracy.threshold", "alpha") %in% names(tempcall)))
      stop("To change any of LSDtype, LSDby, LSDstatistic, LSDaccuracy, ",
           "avsed.tolerance, accuracy.threshold and alpha use redoErrorIntervals.alldiffs")
  
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  if (is.character(validalldifs))
    stop(validalldifs)
  object <- renameDiffsAttr(object)
  if (any(!(factors %in% names(object$predictions))))
    stop("Some factors are not in the predictions component of object")
  if (length(factors) <= 1)
    stop("Need at least two factors to combine")

  LSDtype <- attr(object, which = "LSDtype")
  LSDby <- attr(object, which = "LSDby")
  LSDstatistic <- attr(object, which = "LSDstatistic")
  LSDaccuracy <- attr(object, which = "LSDaccuracy")
  avsed.tolerance <- attr(object$predictions, which = "avsed.tolerance")
  if (is.null(avsed.tolerance))
    avsed.tolerance = 0.25
  accuracy.threshold <- attr(object$predictions, which = "accuracy.threshold")
  if (is.null(accuracy.threshold))
    accuracy.threshold <- NA
  alpha <- attr(object, which = "alpha")
  
  if (!is.null(attr(object, which = "LSDby")))
  {
    if (any(factors %in% LSDby))
    {
      if (all(factors %in% LSDby) && length(factors) == length(LSDby))
      {
        LSDby <- paste(factors, collapse = sep)
        attr(object, which = "LSDby") <- LSDby
        if (!is.null(attr(object$predictions, which = "LSDby")))
          attr(object$predictions, which = "LSDby") <- LSDby
        if (!is.null(attr(object$backtransforms, which = "LSDby")))
          attr(object$backtransforms, which = "LSDby") <- LSDby
      } else
      {
        LSDby <- setdiff(LSDby, factors) 
        if (length(LSDby) == 0)
        {
          LSDtype <- "overall"
          LSDby <- NULL
          warning("Factor(s) in the set of factors to be combined have been removed from the LSDby, reducing LSDby to NULL")
        } else
          warning("Factor(s) in the set of factors to be combined have been removed from the LSDby, reducing LSDby to ",
                  paste(LSDby, collapse = ","))
        if (any(grepl("StandardError.limit", names(object$predictions), fixed = TRUE)))
          error.intervals <- "StandardError"
        if (any(grepl("Confidence.limit", names(object$predictions), fixed = TRUE)))
          error.intervals <- "Confidence"
        if (any(grepl("halfLeastSignificant.limit", names(object$predictions), fixed = TRUE)))
          error.intervals <- "halfLeastSignificant"
        if (!is.null(LSDtype) && LSDtype != "supplied")
          object <- redoErrorIntervals(object, error.intervals = error.intervals,
                                       avsed.tolerance = avsed.tolerance, accuracy.threshold = accuracy.threshold, 
                                       LSDtype = LSDtype, LSDby = LSDby, 
                                       LSDstatistic = LSDstatistic, LSDaccuracy = LSDaccuracy, 
                                       alpha = alpha, ...)
      } 
    }
  }
  
  #Combine the factors
  pred.attr <- attributes(object$predictions)
  comb.fac <- fac.combine(object$predictions[factors], 
                          order = order, combine.levels = combine.levels, sep = sep)
  newfac <- paste(factors, collapse = sep)
  fstfac <- match(factors[1], names(object$predictions))
  object$predictions[fstfac] <- comb.fac
  names(object$predictions)[fstfac] <- newfac
  predictions <- object$predictions[-c(match(factors[-1], names(object$predictions)))]
  pred.attr$names <- attr(predictions, which = "names")
  attributes(predictions) <- pred.attr
  object$predictions <- predictions
  
  if (!is.null(object$backtransforms))
  {
    back.attr <- attributes(object$backtransforms)
    fstfac <- match(factors[1], names(object$backtransforms))
    object$backtransforms[fstfac] <- comb.fac
    names(object$backtransforms)[fstfac] <- newfac
    backtransforms <- object$backtransforms[-c(match(factors[-1], 
                                                            names(object$backtransforms)))]
    back.attr$names <- attr(backtransforms, which = "names")
    attributes(backtransforms) <- back.attr
    object$backtransforms <- backtransforms
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
  
  if (any(grepl("StandardError.limit", names(object$predictions), fixed = TRUE)))
    error.intervals <- "StandardError"
  if (any(grepl("Confidence.limit", names(object$predictions), fixed = TRUE)))
    error.intervals <- "Confidence"
  if (any(grepl("halfLeastSignificant.limit", names(object$predictions), fixed = TRUE)))
    error.intervals <- "halfLeastSignificant"
  if (!is.null(LSDtype) && LSDtype != "supplied")
    object <- redoErrorIntervals(object, error.intervals = error.intervals,
                                 avsed.tolerance = avsed.tolerance, accuracy.threshold = accuracy.threshold, 
                                 LSDtype = LSDtype, LSDby = LSDby, 
                                 LSDstatistic = LSDstatistic, LSDaccuracy = LSDaccuracy,
                                 alpha = alpha, ...)
  return(object)
}

facRecast.alldiffs <- function(object, factor, levels.order = NULL, newlabels = NULL, ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  if (is.character(validalldifs))
    stop(validalldifs)
  object <- renameDiffsAttr(object)

  fac <- factor
  if (any(!(fac %in% names(object$predictions))))
    stop("Some factors are not in the predictions component of object")
  if (length(fac) != 1)
    stop("Only one factor at a time")
  oldlevs <- levels(object$predictions[[fac]])
  if (!is.null(levels.order))
    if (!all(oldlevs %in% levels.order) | length(unique(levels.order)) < length(levels.order))
      stop("The set of levels.order must be unique and contain all the levels in the factor being recast")
  if (!is.null(newlabels))
    if (length(newlabels) != length(oldlevs) | length(unique(newlabels)) < length(newlabels))
      stop("The newlabels must be a set of unique values equal in length to the number of levels in the factor being recast")
  
  #Prepare for recasting
  classify <- attr(object, which = "classify")
  if (is.null(classify))
    stop("The alldiffs object does not have the classify attribute set")
  response <- attr(object, which = "response")
  
  #Recast levels order
  if (!is.null(levels.order))
  {
    object$predictions[fac] <- factor(object$predictions[[fac]], levels = levels.order, ...)
    object$predictions <- object$predictions[do.call(order, object$predictions),]
    
    #revise the alldiffs component
    if (!is.null(object$backtransforms))
    {
      object$backtransforms[fac] <- factor(object$backtransforms[[fac]], levels = levels.order, ...)
      object$backtransforms <- object$backtransforms[do.call(order, object$backtransforms),]
    }
    
    pred.labs <- makePredictionLabels(object$predictions, classify, response)
    pred.lev <- pred.labs$pred.lev

    #revise the alldiffs component
    if (!is.null(object$vcov))
      object$vcov <- object$vcov[pred.lev, pred.lev]

    if (!is.null(object$differences))
      object$differences <- object$differences[pred.lev, pred.lev]

    if (!is.null(object$p.differences))
      object$p.differences <- object$p.differences[pred.lev, pred.lev]

    if (!is.null(object$sed))
      object$sed <- object$sed[pred.lev, pred.lev]
    
  }
  
  #Recast the labels
  if (!is.null(newlabels))
  {
    object$predictions[fac] <- factor(object$predictions[[fac]], labels = newlabels, ...)
    
    #revise the alldiffs component
    if (!is.null(object$backtransforms))
      object$backtransforms[fac] <- object$predictions[fac]

    pred.labs <- makePredictionLabels(object$predictions, classify, response)
    pred.lev <- pred.labs$pred.lev

    #revise the alldiffs component
    if (!is.null(object$vcov))
      rownames(object$vcov) <- colnames(object$vcov) <- pred.lev

    if (!is.null(object$differences))
      rownames(object$differences) <- colnames(object$differences) <- pred.lev

    if (!is.null(object$p.differences))
      rownames(object$p.differences) <- colnames(object$p.differences) <- pred.lev
    
    if (!is.null(object$sed))
      rownames(object$sed) <- colnames(object$sed) <- pred.lev
  }
  
  #Modify the classify
  object <- renewClassify(object, newclassify = attr(object, which = "classify"))

  #deal with LSD component and LSDvalues
  if (!is.null(object$LSD))
  {
    LSDtype <- attr(object, which = "LSDtype")
    if (is.null(LSDtype)) LSDtype <- "overall"
    LSDstatistic <- attr(object, which = "LSDstatistic")
    if (is.null(LSDstatistic)) LSDstatistic <- "mean"
    alpha <- attr(object, which = "alpha")
    if (is.null(alpha)) alpha <- 0.05
    if (!is.null(LSDtype) && LSDtype != "supplied")
      object <- recalcLSD(object, 
                          LSDtype = LSDtype, 
                          LSDby = attr(object, which = "LSDby"),
                          LSDstatistic = LSDstatistic,
                          alpha = alpha)
    if (!is.null(attr(object$predictions, which = "LSDvalues")))
    {
      LSDstat <- attr(object$predictions, which = "LSDstatistic")
      LSDname <- LSDstat2name(LSDstat)
      attr(object$predictions, which = "LSDvalues") <- object$LSD[[LSDname]]
      names(attr(object$predictions, which = "LSDvalues")) <- rownames(object$LSD)
    }
  }
  return(object)
}

facRename.alldiffs <- function(object, factor.names, newnames,  ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  object <- renameDiffsAttr(object)
  if (is.character(validalldifs))
    stop(validalldifs)
  if (any(!(factor.names %in% names(object$predictions))))
    stop("Some factors are not in the predictions component of object")
  if (length(factor.names) != length(newnames))
    stop("The number of factor.names and newnames must be the same")
  pred.attr <- attributes(object$predictions)
  predictions <- object$predictions
  names(predictions)[match(factor.names, names(predictions))] <- newnames
  pred.attr$names <- attr(predictions, which = "names")

  if (!is.null(object$LSD) & !is.null(attr(object, which = "LSDby")))
  {
    LSDby <- attr(object, which = "LSDby")
    if (any(factor.names %in% LSDby))
    {
      which.facs <- LSDby %in% factor.names
      LSDby[match(factor.names[which.facs], LSDby)] <- newnames[which.facs]
      attr(object, which = "LSDby") <- LSDby
      pred.attr$LSDby <- LSDby
    }
  } else
    LSDby <- NULL
  attributes(predictions) <- pred.attr
  object$predictions <- predictions
  
  if (!is.null(object$backtransforms))
  {
    back.attr <- attributes(object$backtransforms)
    backtransforms <- object$backtransforms
    names(backtransforms)[match(factor.names, names(backtransforms))] <- newnames
    back.attr$names <- attr(backtransforms, which = "names")
    back.attr$LSDby <- LSDby
    attributes(backtransforms) <- back.attr
    object$backtransforms <- backtransforms
  }
  
  classify <- attr(object, which = "classify")
  if (is.null(classify))
    stop("The alldiffs object does not have the classify attribute set")
  class.facs <- fac.getinTerm(classify, rmfunction = TRUE)
  class.facs[match(factor.names, class.facs)] <- newnames
  classify <- fac.formTerm(class.facs)
  attr(object, which = "classify") <- classify
  response <- attr(object, which = "response")
  return(object)
}

subset.alldiffs <- function(x, subset = rep(TRUE, nrow(x$predictions)), 
                            rmClassifyVars = NULL, ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(x)  
  if (is.character(validalldifs))
    stop(validalldifs)
  x <- renameDiffsAttr(x)

  #Save attributes
  x.attr <- attributes(x)
  old.attrs <- getAllAttr.alldiffs(x)

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
                                            if (inherits(x, "factor"))
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
      LSDtype <- attr(x, which = "LSDtype")
      if (is.null(LSDtype)) LSDtype <- "overall"
      LSDstatistic <- attr(x, which = "LSDstatistic")
      if (is.null(LSDstatistic)) LSDstatistic <- "mean"
      alpha <- attr(x, which = "alpha")
      if (is.null(alpha)) alpha <- 0.05
      if (!is.null(LSDtype) && LSDtype != "supplied")
        x <- recalcLSD(x, 
                       LSDtype = LSDtype, 
                       LSDby = attr(x, which = "LSDby"),
                       LSDstatistic = LSDstatistic,
                       alpha = alpha)
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
  x <- addMissingAttr.alldiffs(x, old.attrs)
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
sort.predictions.frame <- function(x, decreasing = FALSE, classify, sortFactor = NULL, 
                                   sortParallelToCombo = NULL, sortNestingFactor = NULL, 
                                   sortOrder = NULL, ...)
{
  #Deal with deprecated sortWithinVals argument
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("sortWithinVals" %in% names(tempcall))
      stop("sortWithinVals has been deprecated in sort.predictions.frame - use sortParallelToCombo")
  
  #Change asreml4 names to asreml3 names, if necessary
  if (inherits(x, what = "asreml.predict") && !inherits(x, what = "predictions.frame"))
      x <- as.predictions.frame(x, classify = classify, 
                                se = "std.error", est.status = "status")
  
  #Check that a valid predictions 
  validPredictionsFrame <- validPredictionsFrame(x)  
  if (is.character(validPredictionsFrame))
    stop(validPredictionsFrame)
  
  class.names <-  getClassifyVars(classify)
  if (!all(class.names %in% names(x)))
    stop("The predictions data.frame does not have a column for the following variables ", 
         "in the classify stored with alldiffs\n",  
         paste0(setdiff(class.names, names(x)), collapse = ", "))
  nclassify <- length(class.names)
  if (nclassify < 1)
    stop("Cannot find the classify variables")
  if (!is.null(sortFactor))
  {
    if (length(sortFactor) != 1)
      stop("Can only supply one sortFactor name")
    if (!is.factor(x[[sortFactor]]))
      stop("sortFactor must be a factor")
  }
  if (!is.null(sortNestingFactor))
  {
    if (length(sortNestingFactor) != 1)
      stop("Can only supply one sortNestingFactor name")
    if (!is.factor(x[[sortNestingFactor]]))
      stop("sortNestingFactor must be a factor")
  }
  
  #Save the current attributs
  x.attr <- attributes(x)
  #Get the order in which to sort the predictions
  if (nclassify == 1 || all(class.names %in% c(sortFactor, sortNestingFactor)))
  {
    if (is.null(sortFactor))
    {
      sortFactor <- class.names[1]
      if (!is.null(sortFactor))
        if (!is.factor(sortFactor))
          stop("Single variable in classify is not a factor")
    }
    if (is.null(sortOrder))#sorting order for the sortFactor
    {
      if (is.null(sortNestingFactor))
      {
        val.ord <- data.frame(X1 = x[[sortFactor]], 
                              Ord = order(x$predicted.value, 
                                          decreasing=decreasing))
        names(val.ord)[1] <- sortFactor
      } else #nested sorting order for the sortFactor
      {
        if (sortFactor == sortNestingFactor)
          stop("sortFactor and sortNestingFactor must be different factors")
        tmp <- x
        repln <- table(tmp[[sortNestingFactor]])
        tmp <- split(tmp, tmp[[sortNestingFactor]])
        val.ord <- unlist(lapply(tmp, 
                                 function(pred) ord <-order(pred$predicted.value, 
                                                            decreasing = decreasing)))
        val.ord <- val.ord + rep(cumsum(c(0, repln[1:(length(repln)-1)])), repln)
        val.ord <- data.frame(X1 = x[[sortFactor]],
                              Ord = val.ord)
        names(val.ord)[1] <- sortFactor
      }
    } else #have sortOrder
    {
      old.levs <- levels(x[[sortFactor]])
      if (length(sortOrder) != length(old.levs))
        stop("The number of values in sortOrder must equal the number of levels of sortFactor")
      val.ord <- data.frame(X1 = x[[sortFactor]],
                            Ord = match(as.character(sortOrder), old.levs))
      names(val.ord)[1] <- sortFactor
      if (any(is.na(val.ord$Ord)))
        stop("Not all the values in sortOrder are levels in SortFactor")
    }
    tmp <- val.ord
  } else #classify > 1
  {
    if (is.null(sortFactor))
      stop("The classify for the predictions has multiple variables - need to set sortFactor")
    
    other.vars <- setdiff(class.names, c(sortFactor, sortNestingFactor))
    #Work out order for sortFactor
    if (is.null(sortOrder))
    {
      if (is.null(sortParallelToCombo))
      {
        sortParallelToCombo <- as.list(x[1, other.vars])
        sortParallelToCombo <- lapply(sortParallelToCombo, 
                                      function(var)
                                      {if (is.factor(var)) {var <- as.character(var)}; return(var)})
        names(sortParallelToCombo) <- other.vars
      } else
      {
        if (length(sortParallelToCombo) != length(other.vars))
          stop("Need a value for each classify variable except sortFactor (and sortNestingFactor)")
        else
        {
          if (!all(other.vars %in% names(sortParallelToCombo)))
            stop("The names in sortParallelToCombo do not match the names in the classify")
          if (!all(unlist(lapply(sortParallelToCombo, function(el){length(el)==1}))))
            stop("Each component of sortParallelToCombo can have only a single value")
        }
      }
      subs <- TRUE
      for (name in names(sortParallelToCombo))
      {
        subs <- subs & x[name] == sortParallelToCombo[[name]]
      }
      if (is.null(sortNestingFactor))
      {
        val.ord <- data.frame(X1 = x[subs, sortFactor], 
                              Ord = order(x[subs, "predicted.value"], 
                                          decreasing=decreasing))
        if (nrow(val.ord) < length(levels(val.ord$X1)))
        {
          fac.levs <- levels(val.ord$X1)
          val.ord <- rbind(val.ord,
                           data.frame(X1 = fac.levs[!(fac.levs %in% val.ord$X1)],
                                      Ord = (nrow(val.ord)+1):length(fac.levs)))
        }
        names(val.ord)[1] <- sortFactor
      } else #do a nested sort of the sortFactor
      {
        if (sortFactor == sortNestingFactor)
          stop("sortFactor and sortNestingFactor must be different factors")
        tmp <- x[subs, ]
        repln <- table(tmp[[sortNestingFactor]])
        tmp <- split(tmp, tmp[[sortNestingFactor]])
        val.ord <- unlist(lapply(tmp, 
                                   function(pred) ord <-order(pred$predicted.value, 
                                                              decreasing = decreasing)))
        val.ord <- val.ord + rep(cumsum(c(0, repln[1:(length(repln)-1)])), repln)
        val.ord <- data.frame(X1 = x[subs, sortFactor],
                                Ord = val.ord)
        names(val.ord)[1] <- sortFactor
      }
    } else #have sortOrder
    {
      old.levs <- levels(x[[sortFactor]])
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
  # - this works for unequal replication of sortFactor
  newlevs <-  as.character(val.ord[val.ord$Ord, sortFactor])
  if (!all(levels(x[[sortFactor]] %in% newlevs)))
    stop("Not all levels of ",sortFactor, "occur in the sorted set of predicted values")
  x[sortFactor] <- dae::fac.recast(x[[sortFactor]], levels.order = newlevs)
  #Make sure that the predictions and other components are in standard order for the classify
  x <- x[do.call(order, x),]

  #Set attributes
  attributes(x) <- x.attr
  attr(x, which = "sortFactor") <- sortFactor
  attr(x, which = "sortOrder") <- newlevs
  return(x)
}

sort.alldiffs <- function(x, decreasing = FALSE, classify = NULL, sortFactor = NULL, 
                          sortParallelToCombo = NULL, sortNestingFactor = NULL, 
                          sortOrder = NULL, ...)
{
  #Deal with deprecated sortWithinVals argument
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("sortWithinVals" %in% names(tempcall))
      stop("sortWithinVals has been deprecated in sort.alldiffs - use sortParallelToCombo")
  x <- renameDiffsAttr(x)
  
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
  
  #Sort the predictions frame
  pred <- sort.predictions.frame(x = x$predictions, decreasing = decreasing, 
                                 classify = classify, sortFactor = sortFactor, 
                                 sortParallelToCombo = sortParallelToCombo, 
                                 sortNestingFactor = sortNestingFactor, 
                                 sortOrder = sortOrder, ...)
  x$predictions <- pred

  #Get the order for the full set of predictions
  # - this works for unequal replication of sortFactor
  sortOrder <-  attr(pred, which = "sortOrder")
  x <- facRecast(x, factor = sortFactor, levels.order = sortOrder)
  
  #Set attributes
  attr(x, which = "classify") <- classify
  attr(x, which = "sortFactor") <- sortFactor
  attr(x, which = "sortOrder") <- sortOrder
  return(x)
}

#recalcLSD uses allDIfferences to recalculate the LSD component so cannot use ... to pass parameters 
#  that are in the allDifferences call
recalcLSD.alldiffs <- function(alldiffs.obj, 
                               LSDtype = "overall", LSDsupplied = NULL, 
                               LSDby = NULL, LSDstatistic = "mean", 
                               LSDaccuracy = "maxAbsDeviation", 
                               alpha = 0.05, ...)
{
  #Check for deprecated argument meanLSD.type and warn
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("meanLSD.type" %in% names(tempcall))
      stop("meanLSD.type has been deprecated - use LSDtype")
  if (any(c("transform.power", "offset", "scale", "transform.function")  %in% names(tempcall)))
    stop(cat("Including transform.power, offset, scale or transform.function in the call is invalid \n",
             "- they are obtained from the backtransform component\n"))
  
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction", "supplied")
  avLSD <- AvLSD.options[check.arg.values(LSDtype, AvLSD.options)]
  if (length(avLSD) != 1)
    avLSD <- NULL
  
  LSDstat <- getLSDstatOpt(LSDstatistic = LSDstatistic, avLSD = avLSD, LSDby = LSDby)
  
  LSDacc.options <- c("maxAbsDeviation", "maxDeviation", "q90Deviation", "RootMeanSqDeviation")
  LSDacc <- LSDacc.options[check.arg.values(LSDaccuracy, LSDacc.options)]
  if (length(LSDacc) == 0)
    LSDacc <- "maxAbsDeviation"
  
  #determine transform arguments
  if (is.null(alldiffs.obj$backtransforms))
  {
    transform.power <- 1; offset <- 0; scale <- 1; transform.function <- "identity"
  } else
  {
    transform.power = attr(alldiffs.obj$backtransforms, which = "transform.power")
    offset = attr(alldiffs.obj$backtransforms, which = "offset")
    scale = attr(alldiffs.obj$backtransforms, which = "scale")
    transform.function = attr(alldiffs.obj$backtransforms, which = "transform.function")
    if (is.null(transform.function))
      transform.function <- identity
  }

  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  alldiffs.obj <- renameDiffsAttr(alldiffs.obj)
  kattr <- getAllAttr.alldiffs(alldiffs.obj)
  alldiffs.obj <- allDifferences(alldiffs.obj$predictions, 
                                 classify = attr(alldiffs.obj, which = "classify"), 
                                 vcov = alldiffs.obj$vcov, 
                                 differences = alldiffs.obj$differences, 
                                 p.differences = alldiffs.obj$p.differences,
                                 sed = alldiffs.obj$sed, 
                                 backtransforms = alldiffs.obj$backtransforms,
                                 transform.power = transform.power, 
                                 offset = offset, 
                                 scale = scale, 
                                 transform.function = transform.function, 
                                 tdf = attr(alldiffs.obj, which = "tdf"), alpha = alpha,
                                 LSDtype = avLSD, LSDsupplied = LSDsupplied, 
                                 LSDby = LSDby, LSDstatistic = LSDstat, 
                                 LSDaccuracy = LSDacc, ...)
  alldiffs.obj <- addMissingAttr.alldiffs(alldiffs.obj, kattr)
  return(alldiffs.obj)
}

#Function to explore the LSD values
exploreLSDs.alldiffs <- function(alldiffs.obj,  LSDtype = "overall", LSDby = NULL, 
                                 LSDaccuracy = "maxAbsDeviation", alpha = 0.05, digits = 3, 
                                 retain.zeroLSDs = FALSE, zero.tolerance = .Machine$double.eps ^ 0.5,
                                 plotHistogram = FALSE, ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  alldiffs.obj <- renameDiffsAttr(alldiffs.obj)
  
  AvLSD.options <- c("overall", "factor.combinations")
  avLSD <- AvLSD.options[check.arg.values(LSDtype, AvLSD.options)]
  if (length(avLSD) != 1)
    avLSD <- NULL
  
  LSDacc.options <- c("maxAbsDeviation", "maxDeviation", "q90Deviation", "RootMeanSqDeviation")
  LSDacc <- LSDacc.options[check.arg.values(LSDaccuracy, LSDacc.options)]
  if (length(LSDacc) == 0)
    LSDacc <- "maxAbsDeviation"
  
  LSDstat.hdr <- c("min", "quant10", "quant25", "mean", "median", "quant75", "quant90", "max")

  #Deal with case when have vcov, but not sed
  if (is.null(alldiffs.obj$sed))
  {
    if (!is.null(alldiffs.obj$vcov))
    {
      alldiffs.obj$sed <- alldiffs.obj$vcov
      n <- nrow(alldiffs.obj$sed)
      dvcov <- diag(alldiffs.obj$sed)
      alldiffs.obj$sed <- matrix(rep(dvcov, each = n), nrow = n) + 
        matrix(rep(dvcov, times = n), nrow = n) - 2 * alldiffs.obj$sed
      alldiffs.obj$sed <- sqrt(alldiffs.obj$sed)
      diag(alldiffs.obj$sed) <- NA_real_
    } else
      stop("Neither the vcov or sed components are present in the alldiffs.obj")
  }
  
  
  if (!all(LSDby %in% names(alldiffs.obj$predictions)))
    stop("At least one element of LSDby is not in the predictions component of the alldiffs object\n")
  # classify <- attr(alldiffs.obj, which = "classify")
  # if (!all(unlist(lapply(LSDby, grepl, x = classify, fixed = TRUE))))
  #   stop("One of the elements of LSDby is not in the classify")
  
  denom.df <- attr(alldiffs.obj, which = "tdf")
  if (is.null(denom.df))
    stop(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
               "- LSDs cannot be calculated"))
  t.value = qt(1-alpha/2, denom.df) 
  LSDs <- t.value * alldiffs.obj$sed
  
  #Prepare for frequencies
  LSD.dat <- as.data.frame(getUpperTri(LSDs))
  names(LSD.dat) <- "LSD"
  freq <- hist(LSD.dat$LSD, plot = FALSE, include.lowest = TRUE)
  breaks <- freq$breaks
  
  if (avLSD == "overall")
  {
    #Get distinct  values
    distinct <- sort(unique(signif(na.omit(getUpperTri(LSDs)), digits = digits)))

    #Remove NAs and zero values
    rm.list <- rm.nazero(LSD.dat$LSD, getUpperTri(alldiffs.obj$differences), 
                         retain.zeroLSDs = retain.zeroLSDs, 
                         zero.tolerance = zero.tolerance)
    kLSDs.vec <- rm.list$ksed
    kdifs.vec <- rm.list$kdif
    
    #Get statistics
    allstats <- LSDallstats(kLSDs.vec, kdifs.vec, t.value = 1, LSDaccuracy = LSDacc, 
                            retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance)
    
    #Get per.pred.accuracy
    predacc <- do.call(cbind, lapply(LSDstat.hdr, 
                                     function(LSDstatistic, LSDs, allstats, LSDaccuracy, t.value, 
                                              retain.zeroLSDs, zero.tolerance)
                                     { 
                                       acc <- LSDpred.acc(LSDs, 
                                                          assignedLSD = allstats$statistics[[LSDstatistic]], 
                                                          LSDaccuracy = LSDaccuracy, t.value = t.value, 
                                                          retain.zeroLSDs = retain.zeroLSDs, 
                                                          zero.tolerance = zero.tolerance)
                                       acc <- as.data.frame(acc)
                                       names(acc) <- LSDstatistic
                                       return(acc)
                                     }, 
                                     LSDs = LSDs, allstats = allstats, LSDaccuracy = LSDacc, t.value = 1, 
                                     retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance))
    rownames(predacc) <- rownames(LSDs)
    counts <- freq$counts
    names(counts) <- as.character(freq$mids)
    LSD.list <- list(frequencies = counts, distinct.vals = distinct, 
                     statistics = allstats$statistics, accuracy = allstats$accuracy, 
                     false.pos = allstats$false.pos, false.neg = allstats$false.neg, 
                     per.pred.accuracy = predacc, LSD = LSDs)
    if (plotHistogram)
      print(ggplot(LSD.dat, aes(x = .data[["LSD"]])) + geom_histogram(breaks = breaks) + theme_bw())
  } else #factor.combinations
  {
    LSD.list <- sliceAll(alldiffs.obj, by = LSDby, t.value = t.value, LSDaccuracy = LSDacc, 
                         breaks = breaks, digits = digits, plotHistogram = plotHistogram, 
                         retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance)
    LSD.list <- c(LSD.list, list(LSD = LSDs))
  }
  
  #Set attributes on the lsd.list
  class(LSD.list) <- c("LSDdata", "list")
  attr(LSD.list, which = "LSDtype") <- avLSD
  attr(LSD.list, which = "LSDby") <- LSDby
  attr(LSD.list, which = "LSDaccuracy") <- LSDacc
  attr(LSD.list, which = "alpha") <- alpha
  attr(LSD.list, which = "retain.zeroLSDs") <- retain.zeroLSDs
  return(LSD.list)
}

pickLSDstatistics.alldiffs <- function(alldiffs.obj, 
                                       LSDtype = "overall", LSDby = NULL, 
                                       alpha = 0.05, digits = 3, 
                                       false.pos.wt = NULL, retain.zeroLSDs = FALSE, 
                                       zero.tolerance = .Machine$double.eps ^ 0.5, 
                                       ...)
{
  lsd.errors <- exploreLSDs(alldiffs.obj, LSDtype = LSDtype, LSDby = LSDby, 
                       retain.zeroLSDs = retain.zeroLSDs, 
                       zero.tolerance = zero.tolerance, 
                       ...)
  lsd.errors <- c(lsd.errors["false.pos"], lsd.errors["false.neg"])
  lsd.errors$false.pos <- lsd.errors$false.pos[,-1]
  lsd.errors$false.neg <- lsd.errors$false.neg[,-1]
  nfalserows <- nrow(lsd.errors$false.neg) 
  lsdstats <- sapply(1:nfalserows,
                    function(krow, lsd)
                    {
                      #Determine whether there are no pos or neg errors 
                      no.errors <- lsd$false.pos[krow,] == 0 & lsd$false.neg[krow,] == 0
                      if (any(no.errors)) #no error
                        klsd <- colnames(no.errors)[min(which(no.errors))]
                      else 
                      {  
                        if (is.null(false.pos.wt)) # get the LSD with the min false.pos and, amongst these, the min false.neg
                        {
                          no.pos <- which(lsd$false.pos[krow, ]  == min(lsd$false.pos[krow, ]))
                          min.neg <- lsd$false.neg[krow, ][no.pos]
                          klsd <- names(min.neg)[min(which(min.neg == min(min.neg)))]
                        } else # get the LSD with the min weight sum of false.pos and false.neg
                        {
                          false.no <- lsd$false.pos[krow, ] * false.pos.wt + lsd$false.neg[krow, ]
                          klsd <- names(false.no)[min(which(false.no == min(false.no)))]
                        }
                      }
                      klsd <- gsub("quant", "q", klsd)
                      return(klsd)
                    }, lsd = lsd.errors)
  return(lsdstats)
}


#redoErrorIntervals calls recalcLSD to compute LSD component, which in turn calls allDifferences, but ... does not pass parameters to recalcLSD; 
#redoErrorIntervals is responsible for setting LSD attributes for alldiffs.obj and predictions.frame;
#at the end redoErrorIntervals also calls addBackTransforms, but ... does not pass parameters to addBackTransforms
redoErrorIntervals.alldiffs <- function(alldiffs.obj, error.intervals = "Confidence", alpha = 0.05, 
                                        avsed.tolerance = 0.25, accuracy.threshold = NA, 
                                        LSDtype = NULL, LSDsupplied = NULL, 
                                        LSDby = NULL, LSDstatistic = "mean", 
                                        LSDaccuracy = "maxAbsDeviation", 
                                        retain.zeroLSDs = FALSE, 
                                        zero.tolerance = .Machine$double.eps ^ 0.5, ...)
{
  tempcall <- list(...)
  if ("meanLSD.type" %in% names(tempcall))
    stop("meanLSD.type has been deprecated - use LSDtype")
  if (any(c("transform.power", "offset", "scale", "transform.function")  %in% names(tempcall)))
    stop(cat("Including transform.power, offset, scale or transform.function in the call is invalid \n",
             "- they are obtained from the backtransform component\n"))
  if (any(c("LSDtype", "LSDby", "LSDstatistic", "LSDaccuracy", 
            "avsed.tolerance", "accuracy.threshold", "alpha") %in% names(tempcall)))
    stop("Do not try to change LSDtype, LSDby, LSDstatistic, LSDaccuracy, ",
         "avsed.tolerance, accuracy.threshold or alpha using ...")
  
  #determine transform arguments
  if (is.null(alldiffs.obj$backtransforms))
  {
    transform.power <- 1; offset <- 0; scale <- 1; transform.function <- "identity"
  } else
  {
    transform.power <- attr(alldiffs.obj$backtransforms, which = "transform.power")
    offset <- attr(alldiffs.obj$backtransforms, which = "offset")
    scale <- attr(alldiffs.obj$backtransforms, which = "scale")
    transform.function <- attr(alldiffs.obj$backtransforms, which = "transform.function")
    if (is.null(transform.function))
      transform.function <- "identity"
  }
  
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  alldiffs.obj <- renameDiffsAttr(alldiffs.obj)
  
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction", "supplied")
  avLSD <- AvLSD.options[check.arg.values(LSDtype, AvLSD.options)]
  if (length(avLSD) != 1)
    avLSD <- NULL
 
  LSDstat <- getLSDstatOpt(LSDstatistic = LSDstatistic, avLSD = avLSD, LSDby = LSDby)
  
  LSDname <- paste0(gsub("imum", "", LSDstat, fixed = TRUE), "LSD")
  if (!is.null(LSDby) &&  !is.character(LSDby))
    stop("LSDby must be a character")

  LSDacc.options <- c("maxAbsDeviation", "maxDeviation", "q90Deviation", "RootMeanSqDeviation")
  LSDacc <- LSDacc.options[check.arg.values(LSDaccuracy, LSDacc.options)]
  if (length(LSDacc) == 0)
    LSDacc <- "maxAbsDeviation"
  
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
  
  if (int.opt == "halfLeastSignificant" && is.null(alldiffs.obj$sed))
  {
    if (is.null(alldiffs.obj$vcov))
      stop("cannot compute LSDs because there is no sed or vcov component")
    else #compute sed component
      alldiffs.obj <- makeSED(alldiffs.obj)
  }
  
  #Deal with LSD component
  #Make sure that the correct type of LSDs are available
  if (is.null(avLSD))
  {
    avLSD <- attr(alldiffs.obj, which = "LSDtype")
    if (is.null(avLSD))
    {
      avLSD <- "overall"
    }
    if (avLSD == "factor.combinations" && is.null(LSDby))
    {
      LSDby <- attr(alldiffs.obj, which = "LSDby")
      if (is.null(LSDby))
        stop("LSDtype is factor.combinations, but LSDby is not set")
    }
  }
  if (!(avLSD %in% c("factor.combinations", "supplied")))
    LSDby <- NULL
  
  
  #Determine if no LSD component or the avLSD and LSDby do not match the attributes of alldiffs.obj
  avLSD.diff <- attr(alldiffs.obj, which = "LSDtype")
  LSDby.diff <- attr(alldiffs.obj, which = "LSDby")
  avLSD.same <- TRUE
  LSDby.same <- TRUE
  if (!(is.null(avLSD.diff) && is.null(avLSD)))
  {
    if (!is.null(avLSD.diff) && !is.null(avLSD))
    {
      avLSD.same <- avLSD.diff == avLSD
      if (avLSD.same & avLSD == "factor.combinations")
      {
        if (is.null(LSDby.diff) && is.null(LSDby))
          stop("LSDby not set for LSDtype set to factor.combinations")
        else 
          if (!is.null(LSDby.diff) && !is.null(LSDby))
          {
            LSDby.same <- all(LSDby.diff == LSDby)
          } else
            LSDby.same <- FALSE
      }
    }
  }
  
  if (!is.null(alldiffs.obj$sed))
  {
    #If no LSD component or not match, call recalcLSDs
    if (!is.null(alldiffs.obj$LSD) || !avLSD.same || !LSDby.same)
      alldiffs.obj <- recalcLSD(alldiffs.obj, 
                                LSDtype = avLSD, LSDsupplied = LSDsupplied, 
                                LSDby = LSDby, LSDstatistic = LSDstat, LSDaccuracy = LSDacc,
                                alpha = alpha, ...)
    else #check for assigned LSD
      alldiffs.obj <- checkLSD(alldiffs.obj)
    
    #Calculate overall sed ranges
    t.value = qt(1-alpha/2, denom.df)
    ksed <- getUpperTri(alldiffs.obj$sed)
    kdif <- getUpperTri(alldiffs.obj$differences)
    #remove NA and zero values
    rm.list <- rm.nazero(ksed, kdif, retain.zeroLSDs = retain.zeroLSDs, 
                         zero.tolerance = zero.tolerance)
    ksed <- rm.list$ksed
    kdif <- rm.list$kdif
    overall.LSDs <- LSDstats(ksed, kdif, t.value = t.value)
    rownames(overall.LSDs) <- "overall"
    overall.sed.range <- unlist(abs(overall.LSDs["maxLSD"] - overall.LSDs["minLSD"]) / 
                                  overall.LSDs["meanLSD"])
    if (is.nan(overall.sed.range))
      overall.sed.range <- 0
    
    #Calculate sed ranges from the LSD component
    nLSD <- length(alldiffs.obj$LSD$meanLSD)
    sed.range <- abs(alldiffs.obj$LSD$minLSD - alldiffs.obj$LSD$maxLSD) /  alldiffs.obj$LSD$meanLSD
    sed.range[is.nan(sed.range)] <- 0
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
    
    #Standard errors
    if (int.opt == "StandardError")
      alldiffs.obj$predictions <- within(alldiffs.obj$predictions, 
                                         { lower.StandardError.limit <- 
                                           alldiffs.obj$predictions[["predicted.value"]] - 
                                           alldiffs.obj$predictions[["standard.error"]]
                                         upper.StandardError.limit <- 
                                           alldiffs.obj$predictions[["predicted.value"]] + 
                                           alldiffs.obj$predictions[["standard.error"]]
                                         })

    #half-LSIs
    if (int.opt == "halfLeastSignificant" && (nrow(alldiffs.obj$predictions) != 1))
    { 
       #process for each LSDtype option
      {
        LSDvalues <- NULL
        if (avLSD == "overall" || (avLSD == "supplied" && all(rownames(alldiffs.obj$LSD) == "overall")))
        {
          if (nLSD != 1)
            stop("There is not just one LSD for LSDtype overall")
          rownames(alldiffs.obj$LSD) <- "overall"
          if (!is.na(avsed.tolerance) & overall.sed.range > avsed.tolerance)
          {
            warning("The avsed.tolerance is exceeded - reverting to confidence intervals")
            revert <- TRUE
          } else #plot factor.combination LSD
          {
            alldiffs.obj$predictions <- within(alldiffs.obj$predictions,
                                               {
                                                 if (alldiffs.obj$LSD$assignedLSD < zero.tolerance)
                                                 {
                                                   lower.halfLeastSignificant.limit <- NA
                                                   upper.halfLeastSignificant.limit <- NA
                                                 } else
                                                 {
                                                   lower.halfLeastSignificant.limit <-
                                                     alldiffs.obj$predictions[["predicted.value"]] -
                                                     0.5 * alldiffs.obj$LSD$assignedLSD #[[LSDname]]
                                                   upper.halfLeastSignificant.limit <-
                                                     alldiffs.obj$predictions[["predicted.value"]] +
                                                     0.5 * alldiffs.obj$LSD$assignedLSD #[[LSDname]]
                                                 }
                                               })
            if (!is.na(accuracy.threshold))
            {
              #calculate the accuracy of individual observations
              alldiffs.obj$predictions$LSDwarning <- LSDpred.acc(alldiffs.obj$sed, 
                                                                 assignedLSD = alldiffs.obj$LSD$assignedLSD, 
                                                                 LSDaccuracy = LSDacc, t.value = t.value,
                                                                 retain.zeroLSDs = retain.zeroLSDs, 
                                                                 zero.tolerance = zero.tolerance) > accuracy.threshold
            }
            LSDvalues <- alldiffs.obj$LSD$assignedLSD #overallLSDstat
          }
        } else
        {
          if (avLSD == "factor.combinations" || (avLSD == "supplied" && !all(rownames(alldiffs.obj$LSD) == "overall")))
          {
            #Form levels combination for  LSDs
            levs <- levels(fac.LSDcombs.alldiffs(alldiffs.obj, LSDby))

            #Check have got the correct LSDs in the LSD component
            if (is.null(rownames(alldiffs.obj$LSD)) | nLSD != length(levs) | 
                any(levs != rownames(alldiffs.obj$LSD)))
              stop(paste("For LSDtype factor.combinations, the LSD component of the alldiffs.obj", 
                         "must be a named vector of the LSDs for each combination of the factors in LSDby", 
                         sep = " "))
            if (!is.na(avsed.tolerance) & any(na.omit(sed.range) > avsed.tolerance))
            {
              warning("The avsed.tolerance is exceeded for the factor combinations - reverting to confidence intervals")
              revert <- TRUE
            } else # factor.combination LSD
            {
              LSD.dat <- addByFactorsToLSD.alldiffs(alldiffs.obj, LSDby = LSDby)$LSD
              LSD.dat <- LSD.dat[c(LSDby, "assignedLSD")]
              #save currrent order of predictions and use to restore after the merge
              preds <- alldiffs.obj$predictions
              preds.knam <- names(preds)
              preds.attr <- attributes(alldiffs.obj$predictions)
              preds$rows <- 1:nrow(preds)
              preds.nam <- names(preds)
              preds <- dplyr::left_join(preds, LSD.dat)
              preds <- preds[c(preds.nam, setdiff(names(preds), preds.nam))]
              if (any(diff(preds$rows != 1)))
                preds[preds$rows, ] <- preds
              preds <- preds[, -match("rows", names(preds))]
              preds <- within(preds,                                                  
                              { 
                                lower.halfLeastSignificant.limit <- 
                                  preds[["predicted.value"]] - 0.5 * preds$assignedLSD
                                upper.halfLeastSignificant.limit <- 
                                  preds[["predicted.value"]] + 0.5 * preds$assignedLSD
                                if (any(preds$assignedLSD < zero.tolerance))
                                {
                                  lower.halfLeastSignificant.limit[preds$assignedLSD < zero.tolerance] <- NA
                                  upper.halfLeastSignificant.limit[preds$assignedLSD < zero.tolerance] <- NA
                                } 
                              })
              if (!is.na(accuracy.threshold))
                preds <- within(preds, 
                                {
                                  LSDwarning <- sliceAccs(alldiffs.obj, by = LSDby, 
                                                          LSDstatistic = LSDstat, LSDaccuracy = LSDacc, 
                                                          alpha = alpha)
                                  LSDwarning <- LSDwarning > accuracy.threshold
                                })
              LSDvalues <-  alldiffs.obj$LSD$assignedLSD
              names(LSDvalues) <- rownames(alldiffs.obj$LSD)
              alldiffs.obj$predictions <- preds[, -match("assignedLSD", names(preds))]
              attr(alldiffs.obj$predictions, which = "heading") <- preds.attr$heading
              class(alldiffs.obj$predictions) <- preds.attr$class
            }
          } else
          {
            #per-prediction
            if (avLSD == "per.prediction")
            {
              if (nLSD != nrow(alldiffs.obj$predictions))
                stop("The numbers of LSDs and predicted values are not equal for LSDtype per.prediction")
              if (!is.na(avsed.tolerance) & any(na.omit(sed.range) > avsed.tolerance))
              {
                warning("The avsed.tolerance is exceeded for one or more predictions - reverting to confidence intervals")
                revert <- TRUE
              } else #plot per predictions LSD
              {
                alldiffs.obj$predictions <- within(alldiffs.obj$predictions, 
                                                   { 
                                                     lower.halfLeastSignificant.limit <- 
                                                       alldiffs.obj$predictions[["predicted.value"]] - 
                                                       0.5 * alldiffs.obj$LSD$assignedLSD
                                                     upper.halfLeastSignificant.limit <- 
                                                       alldiffs.obj$predictions[["predicted.value"]] + 
                                                       0.5 * alldiffs.obj$LSD$assignedLSD
                                                   })
                if (!is.na(accuracy.threshold))
                  alldiffs.obj$predictions$LSDwarning <- alldiffs.obj$LSD$accuracyLSD > accuracy.threshold
                LSDvalues <- alldiffs.obj$LSD$assignedLSD
              }
            } 
          } 
        }
        if (!revert)
        {
          #Add LSD attributes, except LSDvalues, to predictions frame
          preds <- alldiffs.obj$predictions
          attributes(preds) <- attributes(alldiffs.obj$predictions)
          attr(preds, which = "LSDtype") <- avLSD
          attr(preds, which = "LSDby") <- LSDby
          attr(preds, which = "LSDstatistic") <- LSDstat
          attr(preds, which = "LSDaccuracy") <- LSDacc
          attr(preds, which = "avsed.tolerance") <- avsed.tolerance
          attr(preds, which = "accuracy.threshold") <- accuracy.threshold
          attr(alldiffs.obj, which = "alpha") <- alpha
          attr(preds, which = "LSDvalues") <- LSDvalues
          alldiffs.obj$predictions <- preds
          attributes(alldiffs.obj$predictions) <- attributes(preds)
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
    {
      pr.attr <- attributes(alldiffs.obj$predictions)
      alldiffs.obj$predictions <- alldiffs.obj$predictions[, c(1:(ks-1), (ks+1):klen, ks)]
      pr.attr$names <- attr(alldiffs.obj$predictions, which = "names")
      attributes(alldiffs.obj$predictions) <- pr.attr
    }
  }
  #Add LSD attributes to the alldiffs object
  attr(alldiffs.obj, which = "LSDtype") <- avLSD
  attr(alldiffs.obj, which = "LSDby") <- LSDby
  attr(alldiffs.obj, which = "LSDstatistic") <- LSDstat
  attr(alldiffs.obj, which = "LSDaccuracy") <- LSDacc
  attr(alldiffs.obj, which = "alpha") <- alpha
  attr(alldiffs.obj$predictions, which = "heading")  <- preds.hd
  
  #Add backtransforms if there has been a transformation
  alldiffs.obj <- addBacktransforms.alldiffs(alldiffs.obj = alldiffs.obj, 
                                             transform.power = transform.power, 
                                             offset = offset, scale = scale, 
                                             transform.function = transform.function)
  return(alldiffs.obj)
}

makeSED <- function(alldiffs.obj)
{ 
  alldiffs.obj$sed <- alldiffs.obj$vcov
  n <- nrow(alldiffs.obj$sed)
  dvcov <- diag(alldiffs.obj$sed)
  alldiffs.obj$sed <- matrix(rep(dvcov, each = n), nrow = n) + 
    matrix(rep(dvcov, times = n), nrow = n) - 2 * alldiffs.obj$sed
  alldiffs.obj$sed <- sqrt(alldiffs.obj$sed)
  diag(alldiffs.obj$sed) <- NA_real_
  return(alldiffs.obj)
}

#allDifferences does not change Error.Intervals, 
#but adds backtransforms depending on transform info
#allDifferences is responsible calculating the LSD component; it used sliceLSDs and LSDstats
#allDifferences calls addBackTransforms to modify backtransforms object
"allDifferences.data.frame" <- function(predictions, classify, vcov = NULL, 
                                        differences = NULL, p.differences = NULL, 
                                        sed = NULL, LSD = NULL,
                                        LSDtype = "overall", LSDsupplied = NULL, 
                                        LSDby = NULL, LSDstatistic = "mean", 
                                        LSDaccuracy = "maxAbsDeviation", 
                                        retain.zeroLSDs = FALSE,
                                        zero.tolerance = .Machine$double.eps ^ 0.5, 
                                        backtransforms = NULL, 
                                        response = NULL, response.title = NULL, 
                                        term = NULL, tdf = NULL, 
                                        x.num = NULL, x.fac = NULL, level.length = NA, 
                                        pairwise = TRUE, alpha = 0.05,
                                        transform.power = 1, offset = 0, scale = 1, 
                                        transform.function = "identity", 
                                        inestimable.rm = TRUE, 
                                        sortFactor = NULL, sortParallelToCombo = NULL, 
                                        sortNestingFactor = NULL, sortOrder = NULL, 
                                        decreasing = FALSE, ...)
#a function to do the calculations to form an alldiffs object
#takes a table of asreml predictions and forms associated statistics
#  for all pairwise differences
{ 
  asr4.2 <- isASReml4_2Loaded(4.2, notloaded.fault = FALSE)
  
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction", "supplied")
  avLSD <- AvLSD.options[check.arg.values(LSDtype, AvLSD.options)]
  if (!is.null(LSDby) &&  !is.character(LSDby))
    stop("LSDby must be a character")
  if (!is.null(LSDsupplied) && avLSD != "supplied")
    warning("LSDsupplied is not NULL and LSDtype in not set to supplied - LSDsupplied will be ignored.")

  LSDstat <- getLSDstatOpt(LSDstatistic = LSDstatistic, avLSD = avLSD, LSDby = LSDby)
  
  LSDname <- paste0(gsub("imum", "", LSDstat, fixed = TRUE), "LSD")
 
  LSDacc.options <- c("maxAbsDeviation", "maxDeviation", "q90Deviation", "RootMeanSqDeviation")
  LSDacc <- LSDacc.options[check.arg.values(LSDaccuracy, LSDacc.options)]
  if (length(LSDacc) == 0)
    LSDacc <- "maxAbsDeviation"
  
  
  tempcall <- list(...)
  if ("levels.length" %in% names(tempcall))
    stop("levels.length has been deprecated - use level.length")
  if ("meanLSD.type" %in% names(tempcall))
    stop("meanLSD.type has been deprecated - use LSDtype")
  
  #determine transform arguments from backtransforms, if set
  if (!is.null(backtransforms))
  {
    transform.power = attr(backtransforms, which = "transform.power")
    offset = attr(backtransforms, which = "offset")
    scale = attr(backtransforms, which = "scale")
    transform.function = attr(backtransforms, which = "transform.function")
    if (is.null(transform.function))
      transform.function <- "identity"
  }
  
  #Change asreml4 names to asreml3 names
  preds.attr <- attributes(predictions)
  predictions <- as.predictions.frame(predictions, classify = classify, 
                                      se = "std.error", est.status = "status")
  preds.attr$names <- attr(predictions, which = "names")
  attributes(predictions) <- preds.attr
  predictions <- renameAttr(predictions, 
                            change.attribs = list(meanLSD.type = "LSDtype", meanLSD = "LSDvalues"))

  alldiffs.obj <- makeAlldiffs(predictions = predictions, 
                               vcov = vcov,
                               differences = differences, 
                               p.differences = p.differences, 
                               sed = sed, LSD = LSD, 
                               backtransforms = backtransforms, 
                               response = response, 
                               response.title = response.title, 
                               term = term, classify = classify, 
                               tdf = tdf, alpha = alpha)
  
  #Check alldiffs.obj
  if (pairwise && is.null(alldiffs.obj$sed) && is.null(alldiffs.obj$vcov))
    stop(paste("No sed or vcov supplied in alldiffs.obj \n",
               "- can obtain using sed=TRUE or vcov=TRUE in predict.asreml"))
  alldiffs.obj <- renameDiffsAttr(alldiffs.obj)
  predictions <- alldiffs.obj$predictions
  rownames(predictions) <- NULL
  preds.attr <- attributes(predictions)

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
                                 backtransforms = alldiffs.obj$backtransforms, 
                                 response = response, 
                                 response.title = response.title, 
                                 term = term, classify = classify, 
                                 tdf = tdf, alpha = alpha)
  }
  response <- as.character(attr(alldiffs.obj, which = "response"))
  
  #Deal with case when have vcov, but not sed
  if (pairwise && !is.null(alldiffs.obj$vcov) && is.null(alldiffs.obj$sed))
    alldiffs.obj <- makeSED(alldiffs.obj)

  #Check that differences are consistent with predictions
  pred.diff <- outer(predictions$predicted.value, predictions$predicted.value, "-")
  if (any(na.omit(abs(pred.diff-alldiffs.obj$differences)) > .Machine$double.eps ^ 0.5))
    stop("The differnces in the differences component of the alldiffs object ",
         "are inconsistent with the predictions in the predictions component")
  
  #Ensure that the columns of predictions are in the same order as the classify 
  class.names <- getClassifyVars(classify)
  if (!all(class.names == names(alldiffs.obj$predictions)[1:length(class.names)]))
  {
    rest <- names(alldiffs.obj$predictions)[(length(class.names)+1):ncol(alldiffs.obj$predictions)]
    preds.attr <- attributes(alldiffs.obj$predictions)
    alldiffs.obj$predictions <- cbind(alldiffs.obj$predictions[class.names],
                                      alldiffs.obj$predictions[rest])
    rownames(alldiffs.obj$predictions) <- NULL
    attr(alldiffs.obj$predictions, which = "heading") <- preds.attr$heading
    class(alldiffs.obj$predictions) <- preds.attr$class
  }
  
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
  
  #Sort if sortFactor set
  # if (!is.null(sortFactor))
  #   alldiffs.obj <- sort(alldiffs.obj, decreasing = decreasing, sortFactor = sortFactor,
  #                        sortParallelToCombo = sortParallelToCombo,
  #                        sortNestingFactor = sortNestingFactor, sortOrder = sortOrder)
  
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
  if (is.null(denom.df) && avLSD != "supplied")
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
    
    #Set LSD component
    {
      if (pairwise && (nrow(alldiffs.obj$predictions) != 1))
      { 
        #calculate LSDs, if not present - it seems that they are never present as set to NULL in makeAlldiffs call
        if (is.null(alldiffs.obj$LSD))
        {
          minLSD <- maxLSD <- meanLSD <- NULL
          t.value = qt(1-alpha/2, denom.df)
          if (avLSD == "overall" || (avLSD == "supplied" && is.null(LSDby)))
          {
            ksed <- getUpperTri(alldiffs.obj$sed)
            kdif <- getUpperTri(alldiffs.obj$differences)
            #remove NA and zero values
            rm.list <- rm.nazero(ksed, kdif, retain.zeroLSDs = retain.zeroLSDs, 
                                 zero.tolerance = zero.tolerance)
            ksed <- rm.list$ksed
            kdif <- rm.list$kdif
            LSDs<- LSDstats(ksed, kdif, t.value, LSDstatistic = LSDstat, LSDaccuracy = LSDacc)
            rownames(LSDs) <- "overall"
          } 
          if (avLSD == "factor.combinations" || (avLSD == "supplied" && !is.null(LSDby))) #factor.combinations
          {
            if (is.null(LSDby))
              stop("Need to specify factors using LSDby for LSDtype = factor.combinations")
            LSDs <- sliceLSDs(alldiffs.obj, by = LSDby, LSDstatistic = LSDstat, LSDaccuracy = LSDacc, 
                              t.value = t.value, alpha = alpha,
                              retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance)
          } 
          if (avLSD == "per.prediction")  #per.prediction
          {
            #set up LSD data.frame
            LSDs <- LSDpred.stats(alldiffs.obj$sed, alldiffs.obj$differences, t.value = t.value, 
                                  LSDstatistic = LSDstat, 
                                  LSDaccuracy = LSDacc, retain.zeroLSDs = retain.zeroLSDs, 
                                  zero.tolerance = zero.tolerance)
          }
          alldiffs.obj$LSD <- LSDs
          if (avLSD == "supplied")
          {
            alldiffs.obj <- addLSDsupplied(alldiffs.obj, LSDsupplied = LSDsupplied, LSDby = LSDby, 
                                           denom.df = denom.df, alpha = alpha)
            #Calculate the accuracy of the supplied LSDs
            if (is.null(LSDby))
            {
              ksed <- getUpperTri(alldiffs.obj$sed)
              kdif <- getUpperTri(alldiffs.obj$differences)
              #remove NA and zero values
              rm.list <- rm.nazero(ksed, kdif, retain.zeroLSDs = retain.zeroLSDs, 
                                   zero.tolerance = zero.tolerance)
              ksed <- rm.list$ksed
              kdif <- rm.list$kdif
              alldiffs.obj$LSD$accuracyLSD <- LSDaccmeas(ksed = ksed, 
                                                         assignedLSD = alldiffs.obj$LSD$assignedLSD, 
                                                         t.value = t.value, LSDaccuracy = LSDacc)
              #Recalculate the false significances
              falsesig <- falseSignif(ksed = ksed, kdif = kdif, assignedLSD = alldiffs.obj$LSD$assignedLSD, 
                                      t.value = t.value)
              alldiffs.obj$LSD$falsePos <- falsesig["false.pos"]
              alldiffs.obj$LSD$falseNeg <- falsesig["false.neg"]
            }
            else
            {
              slLSD <- sliceLSDs(alldiffs.obj, by = LSDby, t.value = t.value, 
                                                        LSDstatistic = LSDstat, LSDaccuracy = LSDacc, 
                                                        alpha = alpha, which.stats = "evalLSD", 
                                                        retain.zeroLSDs = retain.zeroLSDs, 
                                                        zero.tolerance = zero.tolerance)
              alldiffs.obj$LSD[c("accuracyLSD", "falsePos", "falseNeg")] <- 
                slLSD[c("accuracyLSD", "false.pos", "false.neg")]
            }
          }

          attr(alldiffs.obj, which = "LSDtype") <- avLSD
          attr(alldiffs.obj, which = "LSDby") <- LSDby
          attr(alldiffs.obj, which = "LSDstatistic") <- LSDstat
          attr(alldiffs.obj, which = "LSDaccuracy")
        } 
      } 
    }
  }
  #Add backtransforms if there has been a transformation
  alldiffs.obj <- addBacktransforms.alldiffs(alldiffs.obj, 
                                             transform.power = transform.power, 
                                             offset = offset, scale = scale, 
                                             transform.function = transform.function)
  
  #Sort if sortFactor set
  if (!is.null(sortFactor))
    alldiffs.obj <- sort(alldiffs.obj, decreasing = decreasing, sortFactor = sortFactor,
                         sortParallelToCombo = sortParallelToCombo,
                         sortNestingFactor = sortNestingFactor, sortOrder = sortOrder)

  #Check that have a valid alldiffs object
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  
  return(alldiffs.obj)
}

#addBackTransforms adds or recalculates the backtransforms component of an alldiffs.object
#It is invalid to use ... to pass transform.power, offset and scale to it
"addBacktransforms.alldiffs" <- function(alldiffs.obj, 
                                         transform.power = 1, offset = 0, scale = 1, 
                                         transform.function = "identity", ...)
{  
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  alldiffs.obj <- renameDiffsAttr(alldiffs.obj)
  
  #check the transform.function
  link.opts <- c("identity", "log", "inverse", "sqrt", "logit", "probit", "cloglog")
  transfunc <- link.opts[check.arg.values(transform.function, link.opts)]
  
  trans <- transform.power != 1 | offset != 0 | scale != 1
  if (transfunc != "identity" && trans)
    stop("Cannot have both transform.function and one or more of transform.power, ",
         "offset and scale not equal to their defaults")
  
  #'### Save approx.se of backtransformed predictions, if these are available or have been saved
  approx.se <- NULL
  if ("approx.se" %in% names(alldiffs.obj$predictions))
    approx.se <- alldiffs.obj$predictions$approx.se
  else 
  {
    backtransfunc <- attr(alldiffs.obj$backtransforms, which = "transform.function")
    if (is.null(backtransfunc)) backtransfunc <- "identity"
    if (!is.null(alldiffs.obj$backtransforms) && transfunc == backtransfunc)
      approx.se <- alldiffs.obj$backtransforms$standard.error
  }
  
  #Add backtransforms if there has been a transformation 
  if (nrow(alldiffs.obj$predictions) > 0 && (trans || transfunc != "identity"))
  { 
    denom.df <- attr(alldiffs.obj, which = "tdf")
    if (is.null(denom.df))
      warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                    "- p-values and LSDs not calculated"))
    backtransforms <- alldiffs.obj$predictions
    attr(backtransforms, which = "LSDvalues") <- NULL
    kp <- match("predicted.value", names(backtransforms))
    kpl <- pmatch("lower.", names(backtransforms))
    kpu <- pmatch("upper.", names(backtransforms))
    ## As of 3/4/2019 I am allowing backtransformed halfLSD intervals
    #Check if LSD used for predictions and so need to compute CIs
    names(backtransforms)[match("predicted.value", names(backtransforms))] <- 
      "backtransformed.predictions"
    kpl <- pmatch("lower.", names(backtransforms))
    kpu <- pmatch("upper.", names(backtransforms))
    err.int <- TRUE
    if (is.na(kpl) || is.na(kpu))
      err.int <- FALSE
    #Backtransform predictions and intervals for power transformation or transformation function
    if (transform.power == 0 || transfunc == "log")
    { 
      backtransforms$backtransformed.predictions <- 
        exp(backtransforms$backtransformed.predictions)
      if (err.int)
      {
        backtransforms[[kpl]] <- exp(backtransforms[[kpl]])
        backtransforms[[kpu]] <- exp(backtransforms[[kpu]])
      }
    } else
    { 
      if (transform.power != 1 || transfunc %in% c("inverse", "sqrt"))
      { 
        tpower <- transform.power
        if (transfunc == "inverse")
          tpower <- -1
        else
        {
          if (transfunc == "sqrt")
            tpower <- 0.5
        }
        backtransforms$backtransformed.predictions <- 
          backtransforms$backtransformed.predictions^(1/transform.power)
        if (err.int)
        {
          backtransforms[[kpl]] <- backtransforms[[kpl]]^(1/transform.power)
          backtransforms[[kpu]] <- backtransforms[[kpu]]^(1/transform.power)
        }  
      } else
      {
        if (transfunc == "logit")
        { 
          backtransforms <- within(backtransforms, 
                                   {
                                     backtransformed.predictions <-exp(backtransformed.predictions)
                                     backtransformed.predictions <- backtransformed.predictions / 
                                       (1+backtransformed.predictions)
                                   })
          if (err.int)
          {
            backtransforms[[kpl]] <- exp(backtransforms[[kpl]])
            backtransforms[[kpl]] <- backtransforms[[kpl]]/(1 + backtransforms[[kpl]])
            backtransforms[[kpu]] <- exp(backtransforms[[kpu]])/(1 - backtransforms[[kpu]])
            backtransforms[[kpu]] <- backtransforms[[kpu]]/(1 + backtransforms[[kpu]])
          }  
        } else
        {
          if (transfunc == "cloglog")
          { 
            backtransforms <- within(backtransforms, 
                                     {
                                       backtransformed.predictions <-
                                         1 - exp(-exp(backtransformed.predictions))
                                     })
            if (err.int)
            {
              backtransforms[[kpl]] <- 1 - exp(-exp(backtransforms[[kpl]]))
              backtransforms[[kpu]] <- 1 - exp(-exp(backtransforms[[kpu]]))
            }  
          } else
          {  
            if (transfunc == "probit")
            { 
              backtransforms <- within(backtransforms, 
                                       {
                                         backtransformed.predictions <-
                                           pnorm(backtransformed.predictions)
                                       })
              if (err.int)
              {
                backtransforms[[kpl]] <- pnorm(backtransforms[[kpl]])
                backtransforms[[kpu]] <- pnorm(backtransforms[[kpu]])
              } 
            } 
          }
        }
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
    #Set standard.error to missing if a power transformation or transform.function other than identity has been used
    if (transfunc != "identity")
    {
      #deal with backtransformed columns inserted by asreml
      if ("transformed.value" %in% names(alldiffs.obj$predictions))  
      {
        if (!(all((backtransforms$backtransformed.predictions - backtransforms$transformed.value) 
                  < .Machine$double.eps ^ 0.5)))
        {
          warning("The column transformed value inserted by asreml is not the same as the column ",
                  "backtransformed prediction obtained using the inverse transform.function")
          ks <- match("standard.error", names(backtransforms))
          backtransforms[[ks]] <- NA
        } else #have duplicate columns and approx.se
        {
          #remove columns from predictions and back transforms
          #         if ("approx.se" %in% names(backtransforms))
          if (!is.null(approx.se))
            backtransforms$standard.error <- backtransforms$approx.se
          komit <-na.omit(-match(c("transformed.value", "approx.se"), 
                                 names(backtransforms)))
          if (length(komit) > 0)
            backtransforms <- backtransforms[, komit]
          komit <-na.omit(-match(c("transformed.value", "approx.se"), 
                                 names(alldiffs.obj$predictions)))
          if (length(komit) > 0)
            alldiffs.obj$predictions <- alldiffs.obj$predictions[, komit]
        }          
      } else
      {  
        if ("standard.error" %in% names(backtransforms))
        {
          if (!is.null(approx.se))
            backtransforms$standard.error <- approx.se
          else
            backtransforms$standard.error <- NA
        }
      }
    } else
    {
      if (trans) #there has been a transformation
      { 
        if (transform.power != 1)
        {
          ks <- match("standard.error", names(backtransforms))
          backtransforms[[ks]] <- NA
        }
        if (scale != 1)
        {
          ks <- match("standard.error", names(backtransforms))
          backtransforms[[ks]] <- backtransforms[[ks]] / scale
        }
      }
    }
    #Set attributes of backtransform component
    attr(backtransforms, which = "LSDtype") <- attr(alldiffs.obj$predictions, which = "LSDtype")
    attr(backtransforms, which = "LSDby") <- attr(alldiffs.obj$predictions, which = "LSDby")
    attr(backtransforms, which = "LSDstatistic") <- attr(alldiffs.obj$predictions, which = "LSDstatistic")
    attr(backtransforms, which = "LSDaccuracy") <- attr(alldiffs.obj$predictions, which = "LSDaccuracy")
    attr(backtransforms, which = "avsed.tolerance") <- attr(alldiffs.obj$predictions, which = "avsed.tolerance")
    attr(backtransforms, which = "accuracy.threshold") <- attr(alldiffs.obj$predictions, which = "accuracy.threshold")
    attr(backtransforms, which = "alpha") <- attr(alldiffs.obj$predictions, which = "alpha")
    attr(backtransforms, which = "transform.power") <- transform.power
    attr(backtransforms, which = "offset") <- offset
    attr(backtransforms, which = "scale") <- scale
    attr(backtransforms, which = "transform.function") <- transfunc
    alldiffs.obj$backtransforms <- backtransforms
  }
  return(alldiffs.obj)
}

"renewClassify.alldiffs" <- function(alldiffs.obj, newclassify, 
                                     sortFactor = NULL, sortParallelToCombo = NULL, 
                                     sortNestingFactor = NULL, sortOrder = NULL, 
                                     decreasing = FALSE, ...)
{
  tempcall <- list(...)
  if (any(c("transform.power", "offset", "scale", "transform.function")  %in% names(tempcall)))
    stop(cat("Including transform.power, offset, scale or transform.function in the call is invalid \n",
             "- they are obtained from the backtransform component\n"))
  
  #Check for meanLSD.type and, if found, rename to LSDtype
  alldiffs.obj <- renameDiffsAttr(alldiffs.obj)
  
  #determine transform arguments
  if (is.null(alldiffs.obj$backtransforms))
  {
    transform.power <- 1; offset <- 0; scale <- 1; transform.function <- "identity"
  } else
  {
    transform.power <- attr(alldiffs.obj$backtransforms, which = "transform.power")
    offset <- attr(alldiffs.obj$backtransforms, which = "offset")
    scale <- attr(alldiffs.obj$backtransforms, which = "scale")
    transform.function <- attr(alldiffs.obj$backtransforms, which = "transform.function")
    if (is.null(transform.function))
      transform.function <- "identity"
  }
  
  kattr <- getAllAttr.alldiffs(alldiffs.obj)
  alldiffs.obj <- allDifferences(alldiffs.obj$predictions, classify = newclassify, 
                                 vcov = alldiffs.obj$vcov,
                                 differences = alldiffs.obj$differences, 
                                 p.differences = alldiffs.obj$p.differences, 
                                 sed = alldiffs.obj$sed,
                                 LSD = alldiffs.obj$LSD, 
                                 backtransforms = alldiffs.obj$backtransforms,
                                 transform.power = transform.power, 
                                 offset = offset, 
                                 scale = scale,
                                 transform.function = transform.function, 
                                 response = attr(alldiffs.obj, which = "response"), 
                                 response.title = attr(alldiffs.obj, 
                                                       which = "response.title"),
                                 term = attr(alldiffs.obj, which = "term"), 
                                 tdf = attr(alldiffs.obj, which = "tdf"),
                                 alpha = attr(alldiffs.obj, which = "alpha"),
                                 sortFactor = sortFactor, sortOrder = sortOrder, 
                                 sortParallelToCombo = sortParallelToCombo,
                                 sortNestingFactor = sortNestingFactor, 
                                 decreasing = decreasing, ...)
  alldiffs.obj <- addMissingAttr.alldiffs(alldiffs.obj, kattr)
  #Check that the newclassify uniquely indexes the predictions
  newclass.vars <- fac.getinTerm(attr(alldiffs.obj, which = "classify"), rmfunction = TRUE)
  if (any(table(alldiffs.obj$predictions[newclass.vars]) > 1))
    stop("The newclassify variables do not uniquely index the predictions")
 
  return(alldiffs.obj)
}

#calls allDiferences.data.frame, but cannot use ... to pass arguments to allDifferences
"linTransform.alldiffs" <- function(alldiffs.obj, classify = NULL, term = NULL, 
                                    linear.transformation = NULL, Vmatrix = FALSE, 
                                    error.intervals = "Confidence", 
                                    avsed.tolerance = 0.25, accuracy.threshold = NA, 
                                    LSDtype = "overall", LSDsupplied = NULL, 
                                    LSDby = NULL, LSDstatistic = "mean", 
                                    LSDaccuracy = "maxAbsDeviation",
                                    zero.tolerance = .Machine$double.eps ^ 0.5, 
                                    response = NULL, response.title = NULL, 
                                    x.num = NULL, x.fac = NULL, 
                                    tables = "all", level.length = NA, 
                                    pairwise = TRUE, alpha = 0.05, 
                                    inestimable.rm = TRUE, 
                                    ...)
{
  #Check for deprecated argument meanLSD.type and warn
  tempcall <- list(...)
  if (length(tempcall)) 
    if ("meanLSD.type" %in% names(tempcall))
      stop("meanLSD.type has been deprecated - use LSDtype")

  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  alldiffs.obj <- renameDiffsAttr(alldiffs.obj)
  
  AvLSD.options <- c("overall", "factor.combinations", "per.prediction", "supplied")
  avLSD <- AvLSD.options[check.arg.values(LSDtype, AvLSD.options)]
  if (length(avLSD) != 1)
    avLSD <- NULL
  
  LSDstat <- getLSDstatOpt(LSDstatistic = LSDstatistic, avLSD = avLSD, LSDby = LSDby)
  
  LSDacc.options <- c("maxAbsDeviation", "maxDeviation", "q90Deviation", "RootMeanSqDeviation")
  LSDacc <- LSDacc.options[check.arg.values(LSDaccuracy, LSDacc.options)]
  if (length(LSDacc) == 0)
    LSDacc <- "maxAbsDeviation"

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
      
    if (is.null(classify))
      classify <- attr(alldiffs.obj, which = "classify")
    if (is.null(classify))
      stop("The classify has been neither set nor is an attribute of the alldiffs.obj")
    denom.df <- attr(alldiffs.obj, which = "tdf")
    if (is.null(denom.df))
      warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                    "- p-values and LSDs not calculated"))
    
    #get and set term
    if (is.null(term))
    {
      term <- attr(alldiffs.obj, which = "term")
      if (is.null(term))
      {
        term <- attr(alldiffs.obj, which = "classify")
        attr(alldiffs.obj, which = "term") <- term
      }
    }
    #get attributes from predictions
    preds.attr <- attributes(alldiffs.obj$predictions)
    
    #get transform attributes from backtransforms
    transform.power <- 1; offset <- 0; scale <- 1; transform.function <- "identity"
    if (!is.null(alldiffs.obj$backtransforms))
    {
      transform.power = attr(alldiffs.obj$backtransforms, which = "transform.power")
      offset = attr(alldiffs.obj$backtransforms, which = "offset")
      scale = attr(alldiffs.obj$backtransforms, which = "scale")
      transform.function = attr(alldiffs.obj$backtransforms, which = "transform.function")
      if (is.null(transform.function))
        transform.function <- identity
    } 
    
    #Project predictions on submodel, if required
    if (lintrans.type == "submodel")
    {
      colnam <- names(alldiffs.obj$predictions)
      
      #Check that factors in LSDby are in the formula
      term.obj <- as.terms.object(linear.transformation, alldiffs.obj)
      lintrans.fac <- rownames(attr(term.obj, which = "factor"))
      if (!is.null(lintrans.fac) && !all(lintrans.fac %in% colnam))
        stop("Some factors in the linear.transformation are not in the predictions component of the alldiffs object\n")
      if (!is.null(lintrans.fac) && !all(LSDby %in% lintrans.fac))
        warning("Some factors in the LSDby are not in the linear.transformation submodel")
      
      #Form projector on predictions for submodel
      suppressWarnings(Q <- dae::pstructure(linear.transformation, grandMean = TRUE, 
                                            orthogonalize = "eigen", 
                                            data = alldiffs.obj$predictions)$Q)
      Q.submod <- Q[[1]]
      if (length(Q) > 1)
        for (k in 2:length(Q))
          Q.submod <- Q.submod + Q[[k]]
      Q.submod <- dae::projector(Q.submod)
      
      #Check classify variables
      vars <- fac.getinTerm(classify, rmfunction = TRUE)
      if (!all(vars %in% colnam))
        stop("Not all of the variables in the classify are in the predictions component of the alldiffs object\n")
      #Process the classify to ensure there is a separate term for covariates
      tmp <- alldiffs.obj$predictions
      facs <- covs <- list()
      for (var in vars)
      {
        if (is.numeric(tmp[[var]]))
          covs <- c(covs, list(var))
        else
          facs <- c(facs, list(var))
      }
      if (length(facs) > 0)
      {
        tmp$fac.comb <- fac.combine(as.list(tmp[unlist(facs)]))
        full.mod <- "fac.comb"
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
      Q <- dae::pstructure(full.mod, grandMean = TRUE, data = tmp)$Q
      Q.class <- Q[[1]]
      if (length(Q) > 1)
        for (k in 2:length(Q))
          Q.class <- Q.class + Q[[k]]
      Q.class <- dae::projector(Q.class)

      if (any(abs(Q.submod %*% Q.class - Q.submod) > 1e-08))
        stop("Model space for ", linear.transformation, ", with ", degfree(Q.submod), 
             " DF, is not a subspace of the space for the classify ", classify, 
             ", with ", degfree(Q.class), " DF.")
      
      #Form predictions projected onto submodel
      lintrans <- alldiffs.obj$predictions
      lintrans$predicted.value <- as.vector(Q.submod %*% lintrans$predicted.value)
      zeroes <- abs(lintrans$predicted.value) < zero.tolerance
      if (any(zeroes))
        lintrans$predicted.value[zeroes] <- 0
      
      # Calculate standard errors and the variance matrix for differences between predictions
      if (!is.null(alldiffs.obj$vcov))
      {
        lintrans.vcov <- Q.submod %*% alldiffs.obj$vcov %*% Q.submod
        lintrans.vcov <- setToZero(lintrans.vcov, zero.tolerance = zero.tolerance)
        lintrans$standard.error <- as.vector(sqrt(diag(lintrans.vcov)))
        n <- nrow(lintrans.vcov)
        lintrans.sed <- matrix(rep(diag(lintrans.vcov), each = n), nrow = n) + 
          matrix(rep(diag(lintrans.vcov), times = n), nrow = n) - 
          2 * lintrans.vcov
        lintrans.sed <- setToZero(lintrans.sed, zero.tolerance = zero.tolerance)
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
                              sed = lintrans.sed, 
                              LSDtype = LSDtype, LSDsupplied = LSDsupplied, 
                              LSDby = LSDby, LSDstatistic = LSDstat, 
                              LSDaccuracy = LSDacc, 
                              response = response, response.title =  response.title, 
                              term = term, classify = classify, 
                              tdf = denom.df, 
                              transform.power = transform.power, 
                              offset = offset, scale = scale, 
                              transform.function = "identity", 
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
      lintrans$predicted.value <- setToZero(lintrans$predicted.value, 
                                            zero.tolerance = zero.tolerance)
      lintrans$est.status <- "Estimable"
      lintrans$est.status[is.na(lintrans$predicted.value)] <- "Aliased"
      
      # Calculate standard errors and the variance matrix for differences between predictions
      if (!is.null(alldiffs.obj$vcov))
      {
        lintrans.vcov <- linear.transformation %*% alldiffs.obj$vcov %*% t(linear.transformation)
        lintrans.vcov <- setToZero(lintrans.vcov, zero.tolerance = zero.tolerance)
        lintrans$standard.error <- as.vector(sqrt(diag(lintrans.vcov)))
        n <- nrow(lintrans.vcov)
        lintrans.sed <- matrix(rep(diag(lintrans.vcov), each = n), nrow = n) + 
          matrix(rep(diag(lintrans.vcov), times = n), nrow = n) - 2 * lintrans.vcov
        lintrans.sed <- setToZero(lintrans.sed, zero.tolerance = zero.tolerance)
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
                              sed = lintrans.sed, 
                              response = response, response.title =  response.title, 
                              term = classify, classify = "Combination", 
                              tdf = denom.df, 
                              transform.power = transform.power, 
                              offset = offset, scale = scale, 
                              transform.function = transform.function, 
                              x.num = x.num, x.fac = x.fac,
                              level.length = level.length, 
                              pairwise = pairwise, 
                              inestimable.rm = inestimable.rm, 
                              alpha = alpha)
    }
    
    #Add lower and upper uncertainty limits
    diffs <- redoErrorIntervals.alldiffs(diffs, error.intervals = error.intervals, alpha = alpha, 
                                         avsed.tolerance = avsed.tolerance, 
                                         accuracy.threshold = accuracy.threshold,
                                         LSDtype = LSDtype, LSDsupplied = LSDsupplied, 
                                         LSDby = LSDby, LSDstatistic = LSDstat,
                                         LSDaccuracy = LSDacc, zero.tolerance = zero.tolerance, ...)
 
    #Outut tables according to table.opt
    if (!("none" %in% table.opt))
      print(diffs, which = table.opt)
  }
  return(diffs)
}

preserveAttributes.alldiffs <- function(alldiffs.obj)
{
  if (!is.null(alldiffs.obj$predictions))
  {
    alldiffs.obj$predictions <- sticky::sticky(alldiffs.obj$predictions)
    cl.vals <- class(alldiffs.obj$predictions)
    class(alldiffs.obj$predictions) <- cl.vals[c(2:length(cl.vals),1)]
  }
  if (!is.null(alldiffs.obj$LSD))
  {
    alldiffs.obj$LSD <- sticky::sticky(alldiffs.obj$LSD)
    cl.vals <- class(alldiffs.obj$LSD)
    class(alldiffs.obj$LSD) <- cl.vals[c(2:length(cl.vals),1)]
  }
  if (!is.null(alldiffs.obj$backtransforms))
  {
    alldiffs.obj$backtransforms <- sticky::sticky(alldiffs.obj$backtransforms)
    cl.vals <- alldiffs.obj$backtransforms
    class(alldiffs.obj$backtransforms) <- cl.vals[c(2:length(cl.vals),1)]
  }
  return(alldiffs.obj)
}

