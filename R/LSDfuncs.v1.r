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
    
    #Is there only one value for the sed and this is zero and all preds are equal?
    if (length(unique(kLSDs.vec)) == 1 && all(kLSDs.vec == 0) && 
        length(unique(kdifs.vec)) == 1  && 
        diff(range(alldiffs.obj$predictions$standard.error)) < zero.tolerance)
    {
      lsd1 <- t.value * sqrt(2) * alldiffs.obj$predictions$standard.error[1]
      LSDs <- putUpperTri(LSDs, lsd1)
      LSDs <- putLowerTri(LSDs, lsd1)
      
      #Prepare for frequencies
      LSD.dat <- as.data.frame(getUpperTri(LSDs))
      names(LSD.dat) <- "LSD"
      freq <- hist(LSD.dat$LSD, plot = FALSE, include.lowest = TRUE)
      breaks <- freq$breaks
      
      #Get distinct  values
      distinct <- sort(unique(signif(na.omit(LSD.dat$LSD), digits = digits)))
      
      #Remove NAs and zero values
      rm.list <- rm.nazero(LSD.dat$LSD, getUpperTri(alldiffs.obj$differences), 
                           retain.zeroLSDs = retain.zeroLSDs, 
                           zero.tolerance = zero.tolerance)
      kLSDs.vec <- rm.list$ksed
      kdifs.vec <- rm.list$kdif
      allstats <- LSDallstats(kLSDs.vec, kdifs.vec, t.value = 1, LSDaccuracy = LSDacc, 
                              retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance)
      allstats$statistics$c <- 0
      allstats <- c(allstats[1], 
                    lapply(allstats[-1], 
                           function(stat) 
                           {
                             stat <- c(0, rep(NA_real_, length(stat)-1))
                             names(stat) <- c("c", LSDstat.hdr)
                             return(stat)
                           }))
      
      #Get per.pred.accuracy
      predacc <- do.call(cbind, lapply(LSDstat.hdr, 
                                       function(LSDstatistic, LSDs, allstats, LSDaccuracy, 
                                                t.value, retain.zeroLSDs, zero.tolerance)
                                       { 
                                         acc <- LSDpred.acc(LSDs, 
                                                            assignedLSD = allstats$statistics[[LSDstatistic]], 
                                                            LSDaccuracy = LSDaccuracy, 
                                                            t.value = t.value, 
                                                            retain.zeroLSDs = retain.zeroLSDs, 
                                                            zero.tolerance = zero.tolerance)
                                         acc <- as.data.frame(acc)
                                         names(acc) <- LSDstatistic
                                         return(acc)
                                       }, 
                                       LSDs = LSDs, allstats = allstats, LSDaccuracy = LSDacc, t.value = 1, 
                                       retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance))
      rownames(predacc) <- rownames(LSDs)
      predacc <- as.data.frame(lapply(predacc, function(x) x <- rep(NA_real_, length(x))), 
                               row.names = rownames(predacc))
    } else #multiple values
    { 
      #Get statistics
      allstats <- LSDallstats(kLSDs.vec, kdifs.vec, t.value = 1, LSDaccuracy = LSDacc, 
                              retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance)
      
      #Get per.pred.accuracy
      predacc <- do.call(cbind, lapply(LSDstat.hdr, 
                                       function(LSDstatistic, LSDs, allstats, LSDaccuracy, 
                                                t.value, retain.zeroLSDs, zero.tolerance)
                                       { 
                                         acc <- LSDpred.acc(LSDs, 
                                                            assignedLSD = allstats$statistics[[LSDstatistic]], 
                                                            LSDaccuracy = LSDaccuracy, 
                                                            t.value = t.value, 
                                                            retain.zeroLSDs = retain.zeroLSDs, 
                                                            zero.tolerance = zero.tolerance)
                                         acc <- as.data.frame(acc)
                                         names(acc) <- LSDstatistic
                                         return(acc)
                                       }, 
                                       LSDs = LSDs, allstats = allstats, LSDaccuracy = LSDacc, t.value = 1, 
                                       retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance))
      rownames(predacc) <- rownames(LSDs)
    }
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
  lsd.errors$false.pos <- lsd.errors$false.pos[,-1] #remove c
  lsd.errors$false.neg <- lsd.errors$false.neg[,-1] #remove c
  nfalserows <- nrow(lsd.errors$false.neg) 
  lsdstats <- sapply(1:nfalserows,
                     function(krow, lsd)
                     {
                       #Determine whether there are no pos or neg errors 
                       no.errors <- lsd$false.pos[krow,] == 0 & lsd$false.neg[krow,] == 0
                       if (all(is.na(no.errors)))
                         klsd <- "min"
                       else
                       {                      
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
                             klsd <- names(false.no)[which(false.no == min(false.no))]
                             if (length(klsd) > 1) #if several, select one with min false negatives
                               klsd <- klsd[min(which(lsd$false.pos[krow, klsd]  == 
                                                        min(lsd$false.pos[krow, klsd])))]
                           }
                         }
                       }
                       klsd <- gsub("quant", "q", klsd)
                       return(klsd)
                     }, lsd = lsd.errors)
  return(lsdstats)
}

#Function to find LSD values with the minimum errors
findLSDminerrors.alldiffs <- function(alldiffs.obj, 
                                      LSDtype = "overall", LSDby = NULL, 
                                      alpha = 0.05, 
                                      false.pos.wt = 10, nvalues = 100,
                                      retain.zeroLSDs = FALSE, 
                                      zero.tolerance = .Machine$double.eps ^ 0.5, 
                                      trace = FALSE, ...)
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
  
  denom.df <- attr(alldiffs.obj, which = "tdf")
  if (is.null(denom.df))
    stop(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
               "- LSDs cannot becalculated"))
  
  t.value = qt(1-alpha/2, denom.df)
  # Calculate the LSDs and difs   
  LSD.list <- sliceLSDmat(alldiffs.obj, type = avLSD, by = LSDby, 
                          t.value = t.value,  alpha = alpha, 
                          retain.zeroLSDs = retain.zeroLSDs, 
                          zero.tolerance = zero.tolerance)
  if (trace) print(LSD.list)
  if (length(false.pos.wt) == 1)
    false.pos.wt <- rep(false.pos.wt, length(LSD.list))
  else
  {
    if (length(false.pos.wt) != length(LSD.list))
      stop(paste0("false.pos.wt must have 1 or the number of combinations (", length(LSD.list), 
                  ") of the LSDby factors."))
  }
  
  #Search for minimum LSD
  stepsize <- 1/nvalues
  optLSDs <- mapply(function(kLSDs, kpos.wt)
  {
    #Search for optimal LSD
    LSDrange <- range(kLSDs$lsd)
    if (trace) {cat("\n\n#### New set\n"); print(LSDrange)}
    if (diff(LSDrange) > zero.tolerance)
    { 
      #Do a grid search for the minLSD
      testLSDs <- lapply(seq(0, 1, stepsize),
                         function(step, LSDrange, kLSDs, kpos.wt)
                         {
                           testlsd <- LSDrange[1] + step * diff(LSDrange)
                           false.vals <- c(testlsd, 
                                           falseErrorNums(testlsd, kLSDs, 
                                                          kpos.wt))
                           names(false.vals)[1] <- "LSD"
                           return(false.vals)
                         }, LSDrange, kLSDs, kpos.wt)
      if (trace) {cat("\n#### Initial grid search\n") ; print(testLSDs)}
      testLSDs <- as.data.frame(do.call(rbind, testLSDs))
      rownames(testLSDs) <- NULL
      min.criterion <- min(testLSDs["false.criterion"], na.rm = TRUE)
      which.min <- testLSDs["false.criterion"] == min.criterion
      testLSDs <- testLSDs[which.min,]
      #select those with the minimum false positives
      minpos <- min(testLSDs["false.pos"], na.rm = TRUE)
      testLSDs <- testLSDs[testLSDs["false.pos"] == minpos, ]
      optLSD <- testLSDs[1,]
      #Check in the neighbourhood below the minlsd for a smaller lsd
      if (optLSD$LSD > LSDrange[1])
      { 
        testlsd <- optLSD$LSD
        substepsize <- 0.1
        subminlsd <- testlsd - stepsize
        stepval <- substepsize*stepsize
        if (trace) cat("\n#### Finer grid search\n")
        for (step in seq(1-substepsize, 0, -substepsize))
        {
          testlsd <-testlsd - stepval 
          false.nums <- falseErrorNums(testlsd, kLSDs, kpos.wt)
          if (trace) print(c(testlsd, false.nums))
          if (optLSD["false.criterion"] != false.nums["false.criterion"] || 
              optLSD["false.pos"] != false.nums["false.pos"]) break
        }
        optlsd <- testlsd + stepval
        optLSD <- c(optlsd, falseErrorNums(optlsd, kLSDs, kpos.wt))
      }
    } else
    {
      falsesig <- rep(NA, 3)
      names(falsesig) <- c("false.pos", "false.neg", "false.criterion")
      optLSD <- c(LSDrange[1], falsesig)
    }
    names(optLSD)[1] <- "LSD"
    return(optLSD)
  }, LSD.list, false.pos.wt, SIMPLIFY = FALSE)
  optLSDs <- as.data.frame(do.call(rbind, optLSDs))
  
  #Set attributes on the LSD.list
  attr(optLSDs, which = "LSDtype") <- avLSD
  attr(optLSDs, which = "LSDby") <- LSDby
  attr(optLSDs, which = "alpha") <- alpha
  attr(optLSDs, which = "retain.zeroLSDs") <- retain.zeroLSDs
  return(optLSDs)
}

falseErrorNums <- function(testLSD.val, kLSDs, false.pos.wt)
{
  sig.actual <- kLSDs$sig.actual
  sig.approx <- abs(kLSDs$dif) >= testLSD.val 
  falsesig <- c(sum(!sig.actual & sig.approx, na.rm = TRUE), 
                sum(sig.actual & !sig.approx, na.rm = TRUE))
  falsesig <- c(falsesig, 
                falsesig[1] * false.pos.wt + falsesig[2])
  names(falsesig) <- c("false.pos", "false.neg", "false.criterion")
  return(falsesig)
}

falseSignifCriterion <- function(testLSD.val, kLSDs, false.pos.wt)
{
  falsesig <- falseErrorNums(testLSD.val, kLSDs, false.pos.wt)
  return(falsesig["false.criterion"])
}

#Function to produce the LSDs for combinations of the levels of the by factor(s)
sliceLSDmat <- function(alldiffs.obj, type, by, 
                        t.value,  alpha = 0.05, 
                        retain.zeroLSDs = FALSE, 
                        zero.tolerance = .Machine$double.eps ^ 0.5)
{
  if (!all(by %in% names(alldiffs.obj$predictions)))
    stop("At least one element of LSDby is not in the predictions component of the alldiffs object\n")
  #  classify <- attr(alldiffs.obj, which = "classify")
  #  if (!all(unlist(lapply(by, grepl, x = classify, fixed = TRUE))))
  #    stop("One of the elements of LSDby is not in the classify")
  
  lsd <- t.value * alldiffs.obj$sed
  dif <- alldiffs.obj$differences
  diag(dif) <- NA_real_
  
  if (type == "overall")
  {
    #Remove NAs and zero values
    rm.list <- rm.nazero(getUpperTri(lsd), getUpperTri(dif), 
                         retain.zeroLSDs = retain.zeroLSDs, 
                         zero.tolerance = zero.tolerance)
    names(rm.list) <- c("lsd", "dif")
    if (length(rm.list$lsd) == 1 && length(rm.list$dif) == 1  && rm.list$lsd == 0 && 
        diff(range(alldiffs.obj$predictions$standard.error)) < zero.tolerance)
      rm.list$lsd <- t.value * sqrt(2) * alldiffs.obj$predictions$standard.error[1]
    rm.list$sig.actual <- abs(rm.list$dif) >= rm.list$lsd
    LSDs <- list(overall = rm.list)
  } else
  {    
    #Get the LSDs
    fac.comb <- fac.LSDcombs.alldiffs(alldiffs.obj, by)
    levs <- levels(fac.comb)
    #loop over LSDby combinations
    LSDs <- lapply(levs, 
                   function(lev, lsd, dif, t.value, alldiffs.obj)
                   {
                     krows <- lev == fac.comb
                     if (!any(krows) || length(fac.comb[krows]) == 1) 
                     {
                       #have a single prediction
                       rm.list <- list(lsd = t.value * sqrt(2) * 
                                         alldiffs.obj$predictions$standard.error[krows],
                                       dif = 0,
                                       sig.actual = NA)
                     } else  #have several predictions
                     {
                       klsd <- getUpperTri(lsd[krows, krows])
                       kdif <- getUpperTri(dif[krows, krows])
                       rm.list <- rm.nazero(klsd, kdif, retain.zeroLSDs = retain.zeroLSDs,
                                            zero.tolerance = zero.tolerance)
                       #Is there only one value for the sed and this is zero?
                       if (all(abs(klsd) < zero.tolerance) && 
                           diff(range(alldiffs.obj$predictions$standard.error[krows])) < zero.tolerance)
                       {
                         rm.list <- list(lsd = t.value * sqrt(2) *
                                           alldiffs.obj$predictions$standard.error[krows][1],
                                         dif = 0)
                       } else
                       {
                         #remove NA and zero values
                         rm.list <- rm.nazero(klsd, kdif, retain.zeroLSDs = retain.zeroLSDs,
                                              zero.tolerance = zero.tolerance)
                         names(rm.list) <- c("lsd", "dif")
                       }
                       rm.list$sig.actual <- abs(rm.list$dif) >= rm.list$lsd
                     }
                     return(rm.list)
                   }, lsd = lsd, dif = dif, t.value = t.value, alldiffs.obj = alldiffs.obj)
    names(LSDs) <- levs
  }
  return(LSDs)
}

##Function to add letters indicating pairwise significances between predictions based on t-tests
##Adapted from Clayton Forknall's LSDsubscripts script by Chris Brien
addPairwiseLetters.alldiffs <- function(alldiffs.obj, within = NULL, 
                                        alpha = 0.05, ...)
{
  #Remove a pairwiseSignific column, if there is one in the predictions
  if ("pairwiseSignif" %in% names(alldiffs.obj$predictions))
  { 
    preds <- alldiffs.obj$predictions
    preds <- preds[, -match("pairwiseSignif", names(preds))]
    alldiffs.obj$predictions <- preds
  }
  
  ##First grab preds and append the classify labels
  preds <- alldiffs.obj$predictions
  nfetch <- match("predicted.value", names(preds))
  preds <- preds[1:nfetch]
  preds$labs <- row.names(alldiffs.obj$p.differences)
  
  ##order from largest to smallest prediction
  preds <- preds[order(preds$predicted.value, decreasing = TRUE),]
  
  ##pull out the p.differences matrix
  pdiff.mat <- alldiffs.obj$p.differences
  
  ##order the pdiff.mat to be consistent with preds
  pdiff.mat <- pdiff.mat[preds$labs, preds$labs]
  preds <- subset(preds, select = -labs)
  
  ##turn the pdiff.mat into a binary matrix of 0s and 1s depending on whether the p-value is <= 0.05
  pdiff.bin <- pdiff.mat
  pdiff.bin[pdiff.mat > alpha] <- 0
  pdiff.bin[pdiff.mat <= alpha] <- 1
  
  #ensure that all diagonal elements of pdiff.bin are zero
  diag(pdiff.bin) <- 0

  #also ensure that NaNs are 0s too - this is for cases where we are using 
  #                                   linTransform output and chosen models
  pdiff.bin[is.nan(pdiff.bin)] <- 0

  #Call calcLetters to compute the significance letters.  
  if (is.null(within))
  { 
    subs <- calcLetters(pdiff.bin)
    preds$pairwiseSignif <- subs
  }
  else
  { 
    if (!all(within %in% names(alldiffs.obj$predictions)))
      stop("At least one element of within is not in the predictions component of the alldiffs object\n")
    subs <- slicepMat(pdiff.bin, within = within, preds = preds)
    preds <- dplyr::left_join(preds, subs)
  }
  
  #Add letters to predictions.frame
  alldiffs.obj$predictions <- dplyr::left_join(alldiffs.obj$predictions, preds)
  alldiffs.obj <- renewClassify(alldiffs.obj, 
                                   newclassify = attr(alldiffs.obj, which = "classify"))
  
  #Returning the alldiffs.obj
  return(alldiffs.obj)
}

calcLetters <- function(pdiff.bin)
{
  ###################################################################
  #Determine subscripts#
  #Adpated from LSDcalc script by Clayton Forknall
  ###################################################################
  
  #putting together index of letters and symbols
  letter.array <- c(letters,toupper(letters),"0","1","2","3","4","5","6","7","8","9","!","@","#","$","%","^","&","*","+","=")
  
  #defining variables to be used in this functiions
  nrows <- nrow(pdiff.bin)
  nminus1 <- nrows-1
  rowind <- 1:nrows
  
  #start the subscripting procedure
  subs <- NA
  subs[1] <- letter.array[1] #making "a" the first subscript
  subs[2:nrows] <- 1 #initalising the vector of subscripts
  
  #initalise some useful variables
  test <- 1
  testlast <- 1
  top <- 1 #the position of the column of pdiff.bin being considered
  topm1 <- 1
  imaxsub <- 1 #the an index of the maximum subscript that has been used
  total <- NULL
  total[1:nrows] <- 1 #used to test significance between one column of pdiff.bin to the next
  testdiff <- 1
  nextmean <- 1 #this is a trigger variable used throughout the loop below
  nsrow <- NULL
  nsrow1 <- NULL
  nsrow2 <- NULL
  bi <- NULL
  
  #begin determining the subscripts
  for (I in 2:nrows){
    top <- I-1
    nextmean <- 2
    tmp <- pdiff.bin[1:I,I]
    total[I] <- sum(tmp)
    
    if (total[I]==0){ #mean is not significant from previous -> give same subscript
      subs[I] <- paste(letter.array[1:imaxsub],collapse="")
    }
    else if(total[I]==top){ #mean is significant from all previous means -> give new subscript
      imaxsub <- imaxsub+1
      subs[I] <- letter.array[imaxsub]
    }
    else { #mean is significant from at least one and not significant from at least one mean
      for (J in top:1){ #test for exact equality in terms of sig differences
        tmp1 <- pdiff.bin[1:I,I]
        tmp2 <- pdiff.bin[1:I,J]
        testdiff <- sum(abs(tmp1-tmp2))
        
        if (testdiff==0){ #if exact equality in significant differences -> use same subscript and move onto next mean
          subs[I] <- subs[J]
          nextmean <- 1
          break
        }
      }
      #check at this point if nextmean ==1, if so nothing more to do
      #if nextmean==2, then get new subscript + assign new subscript/s to all nonsignificant
      if (nextmean!=1){
        k <- 0
        for (J in 1:top){ #record the index of the rows we are interested in
          sigdiff <- pdiff.bin[J,I]
          if (sigdiff==0){ #if not significant (ie =0)
            k <- k+1
            nsrow[k] <- J
          }
        }
        
        mnewsub <- k
        for (jns in 1:mnewsub){ #loop at most mnewsub times
          imaxsub <- imaxsub+1
          if (jns==1){
            subs[I] <- letter.array[imaxsub]
          }
          else {
            subs[I] <- paste(subs[I],letter.array[imaxsub], sep="")
          }
          k1 <- 0
          k2 <- 0
          for (j in 1:k){
            J <- nsrow[j]
            if (j == 1){ #know we can add new subscript without conflict
              k1 <- k1+1
              nsrow1[k1] <- J
              subs[J] <- paste(subs[J],letter.array[imaxsub], sep="")
            }
            else {
              thistime <- 1
              for (j1 in 1:k1){
                J1 <- nsrow1[j1]
                test <- pdiff.bin[J1,J]
                if (test == 1){
                  thistime <- 0
                  k2 <- k2+1
                  nsrow2[k2] <- J
                }
                if (thistime == 0){
                  break
                }
              }
              if (thistime ==1){ #ok to add subscript
                k1 <- k1+1
                nsrow1[k1] <- J
                subs[J] <- paste(subs[J],letter.array[imaxsub], sep="")
              }
            }
          }
          if (k1 == k){
            break
          }
          k <- k2
          for (k3 in 1:k){
            nsrow[k3] <- nsrow2[k3]
          }
        } # end of allocation of new letter(s) for this row
      }
    }
  }
  
  #remove any redundant letters
  a <- matrix(NA,nrow=nrows,ncol=imaxsub) #number of letters in use form cols, while means form rows
  b <- matrix(NA,nrow=nrows,ncol=imaxsub)
  
  for (j in 1:imaxsub){
    a[,j] <- 0
    b[,j] <- 0
  }
  nb <- 0
  
  for (i in 1:nrows){ #extract pattern from subs
    for (j in 1:imaxsub){
      if(length(grep(letter.array[j],subs[i]))==1){
        a[i,j] <- 1
      }
    }
  }
  
  redund <- array(0,dim=imaxsub)
  redundi <- 0
  
  for (j2 in 1:imaxsub){ #identify any redundant columns
    j <- imaxsub-j2+1
    for (j3 in 1:j2){
      j1 <- imaxsub-j3+1
      if((redund[j1]==0) && (j1!=j)){
        testv <- a[,j1] & a[,j]
        #testv[TRUE] <- 1
        #testv[FALSE] <- 0
        test1 <- sum(abs(a[,j1]-testv))
        test <- sum(abs(a[,j]-testv))
        if(test1 == 0){
          redund[j1] <- 1
          redundi <- 1
        }
        else if(test == 0){
          a[,j] <- a[,j1]
          redund[j1] <- 1
          redundi <- 1
        }
      }
    }
  }
  
  if (redundi==1){ #redundant columns were found, so need to redo subscripts
    j1 <- 0
    for (j in 1:imaxsub){ #pack up subscript indicators into b[]
      if (redund[j]==0){
        j1 <- j1+1
        b[,j1] <- a[,j]
      }
    }
    imaxsub <- j1
    for (j in 1:imaxsub){
      a[,j] <- b[,j]
      redund[j] <- 0
    }
  }
  
  j1 <- 0
  for (i in 1:nrows){ #sort columns so letters will start at the top
    for (j in 1:imaxsub){
      test <- a[i,j]
      if (redund[j]==0 & test==1){
        j1 <- j1+1
        b[,j1] <- a[,j]
        redund[j] <- 1
        if (j1 == imaxsub){
          break
        }
      }
    }
    if (j1==imaxsub){
      break
    }
  }
  
  #checking final result against input significant, just to be sure
  algfhead <- 0
  brow <- matrix(NA,ncol=nrows,nrow=imaxsub)
  for (i in 1:nrows){
    for (j in 1:imaxsub){
      brow[j,i] <- b[i,j]
    }
  }
  for (i in 1:nrows){
    for (j in 1:i){
      if (j != i){
        bsig <- sum(brow[,i] & brow[,j])==0
        if (bsig != pdiff.bin[i,j]){
          if(algfhead==0){
            print("Warning: subscript lettering approximate")
            algfhead <- 1
          }
          if (bsig==1){
            print.temp <- paste0("Difference between means ",j," and ",i," is not significant")
          }
          else {
            print.temp <- paste0("Difference between means ",j," and ",i," is significant")
          }
          #print(print.temp)
        }
      }
    }
  }
  
  #regenerate subscripts
  imaxsub1 <- imaxsub+1
  letter.array[imaxsub1] <- "_"
  for (i in 1:nrows){
    j1 <- 0
    for (j in 1:imaxsub){
      test <- b[i,j]
      if (test==1){
        j1 <- j1+1
        bi[j1] <- j
      }
    }
    ttmp <- paste(letter.array[bi[1:j1]], collapse="")
    subs[i] <- ttmp
  }
  return(subs)
}

#Function to produce the LSDs for combinations of the levels of the by factor(s)
slicepMat <- function(pdiff.bin, within, preds)
{
  #Get the letters
  fac.comb <- fac.Predscombs.predictions.frame(preds, by = within)
  levs <- levels(fac.comb)
  #loop over within combinations
  subs <- lapply(levs, 
                 function(lev, pdiff.bin, preds)
                 {
                   krows <- lev == fac.comb
                   subs.lev <- preds[krows, ]
                   if (!any(krows) || length(fac.comb[krows]) == 1) 
                   {
                     #have a single prediction
                     subs.lev <- NULL
                   } else  #have several predictions
                   {
                     pdiff.lev <- pdiff.bin[krows, krows]
                     subs.lev$pairwiseSignif <- calcLetters(pdiff.lev)
                   }
                   return(subs.lev)
                 }, pdiff.bin = pdiff.bin, preds = preds)
  subs <- do.call(rbind, subs)

  return(subs)
}

