#### LSD functions
#Function to convert a LSDstatistic option into a column name in the LSD component
LSDstat2name <- function(LSDstat)
  LSDname <- paste0(gsub("imum", "", LSDstat, fixed = TRUE), "LSD")

getUpperTri <- function(x)
  x <- x[upper.tri(x)]

getLSDstatOpt <- function(LSDstatistic, avLSD, LSDby)
{
  LSDstat.options <- c("minimum", "q10", "q25", "mean", "median", "q75", "q90", "maximum")
  if (length(LSDstatistic) == 0)
    LSDstat <- "mean"
  else 
  {
    if (length(LSDstatistic) == 1)
      LSDstat <- LSDstat.options[check.arg.values(LSDstatistic, LSDstat.options)]
    else
    {
      if (avLSD == "factor.combinations")
      {
        if (is.null(LSDby))
          stop("Multiple values for LSDstatistic have been specified so that LSDby must not be NULL")
        LSDstat <- LSDstat.options[unlist(lapply(LSDstatistic, check.arg.values, options=LSDstat.options))]
      }
      else
        stop(paste("LSDstatistic must contain only one value for LSDtype", avLSD))
    }
  }  
  return(LSDstat)
}

#Function to make a combined factor of the LSDby factor
#Convert any non-factors and form levels combination for which a LSDs are required
fac.LSDcombs.alldiffs <- function(alldiffs.obj, by)
{
  if (!is.character(by))
    stop("by is not a character")
  fac.list <- as.list(alldiffs.obj$predictions[by])
  fac.list <- lapply(fac.list, 
                     function(x) 
                     {
                       if (!inherits(x, "factor"))
                         x <- factor(x)
                       return(x)
                     })
  fac.comb <- fac.combine(fac.list, combine.levels = TRUE)
  if (length(fac.comb) != nrow(alldiffs.obj$sed))
    stop("Variable(s) in (LSD)by argument are not the same length as the order of the sed matrix")
  return(fac.comb)
}

#Check LSDsupplied and make sure in format for inclusion as an LSD component
addLSDsupplied <- function(alldiffs.obj, LSDsupplied, LSDby, denom.df, alpha)
{
  if (is.null(LSDsupplied))
    stop("Must set LSDsupplied for LSDtype set to supplied")
  
  if (inherits(LSDsupplied, what = "data.frame"))
  {
    if (ncol(LSDsupplied) > 1)
    {
      if (ncol(LSDsupplied) != length(LSDby)+1 || any(names(LSDsupplied)[1:length(LSDby)] != LSDby))
        stop("If LSD supplied as a data.frame with more than one column ", 
             "it should have a column for each element of LSDby followed by a column of LSD values")
      fac.comb <- fac.combine(as.list(alldiffs.obj$predictions[LSDby]), combine.levels = TRUE)
      LSDnam <- names(LSDsupplied)[ncol(LSDsupplied)]
      LSDsupplied <- LSDsupplied[[LSDnam]]
      names(LSDsupplied) <- levels(fac.comb)
    } else #one-column LSDsupplied
    {
      if (is.null(rownames(LSDsupplied)))
        stop("Must supply LSDby combination as row names to LSDsupplied")
      #Convert any non-factors and form levels combination for which a LSDs are required
      if (is.null(LSDby))
      {
        if (rownames(LSDsupplied) != "overall")
          stop("A single row should have row name equal to overall")
        LSDsupplied <- LSDsupplied["overall",1]
      } else
      {
        fac.comb <- fac.LSDcombs.alldiffs(alldiffs.obj, by = LSDby)
        levs <- levels(fac.comb)
        if (any(levs != rownames(LSDsupplied)))
          stop("The combinations of the LSDby variables for LSDsupplied and those in the predictions do not match")
        LSDsupplied <- LSDsupplied[,1]
        names(LSDsupplied) <- levels(fac.comb)
      }
    }
  } else #a numeric
  {
    if (is.null(names(LSDsupplied)))
      stop("Must supply LSDby combinations or overall as names for LSDsupplied")
    #Convert any non-factors and form levels combination for which a LSDs are required
    if (is.null(LSDby))
    {
      if (names(LSDsupplied) != "overall")
        stop("A single value should have name overall")
      # LSDsupplied <- data.frame(meanLSD = LSDsupplied, row.names = "overall")
      names(LSDsupplied) <- "overall"
    } else
    {
      levs <- levels(fac.LSDcombs.alldiffs(alldiffs.obj, by = LSDby))
      if (any(levs != rownames(LSDsupplied)))
        stop("The combinations of the LSDby variables for LSDsupplied and those in the predictions do not match")
    }
  }
  #Put the supplied LSD into an LSD data.frame
  alldiffs.obj$LSD$assignedLSD <- LSDsupplied
  return(alldiffs.obj)
}

addByFactorsToLSD.alldiffs <- function(alldiffs.obj, LSDby = NULL)
{
  if (is.null(LSDby))
    LSDby <- attr(alldiffs.obj, which = "LSDby")
  
  LSD.facsvars <- alldiffs.obj$predictions[LSDby]
  names(LSD.facsvars) <- LSDby
  
  #Determine if any numerics
  any.num <- lapply(LSD.facsvars, 
                    function(fac) {is.numeric(fac)})
  
  #Convert numerics to factors
  if (any(unlist(any.num)))
  {
    nams.num <- names(any.num)
    LSD.facsvars <- lapply(LSD.facsvars, 
                           function(fac)
                           {
                             if (is.numeric(fac))
                               fac <- factor(fac)
                             return(fac)
                           }) 
    
  } 
  #Generate factors and combine
  combs <- dae::fac.genfactors(LSD.facsvars)
  combs <- cbind(fac.comb = fac.combine(combs, combine.levels = TRUE), combs)

  #Convert numerics back from factors
  if (any(unlist(any.num)))
    combs[nams.num] <- lapply(combs[nams.num], 
                              function(fac) {fac <- as.numfac(fac); return(fac)}) 
  
  #Add LSDby factors to LSD object
  combs <- dplyr::left_join(data.frame(fac.comb = rownames(alldiffs.obj$LSD)),
                            combs)
  combs <- combs[-match("fac.comb", names(combs))]
  alldiffs.obj$LSD <- cbind(combs, alldiffs.obj$LSD)
  return(alldiffs.obj)
}

checkLSD <- function(alldiffs.obj)
{
  if (!"assignedLSD" %in% names(alldiffs.obj$LSD))
  { 
    alldiffs.obj$LSD$assignedLSD <- alldiffs.obj$LSD$meanLSD
    alldiffs.obj$LSD$accuracyLSD <- NA_real_
  }
  return(alldiffs.obj)
}

#Function to find nonzero values in a vector using zero.tolerance to determine if zero
findNonzero <- function(x, zero.tolerance = .Machine$double.eps ^ 0.5)
{
  max.x <- max(x)
  which.nonzero <- NULL
  #Retain only nonzero variances
  if (max.x > zero.tolerance && sum(x/max.x > zero.tolerance) > 0)
    which.nonzero <- x/max.x > zero.tolerance
  return(which.nonzero)
}

#Function to remove nas & zero values from a vector using zero.tolerance to determine if zero
rm.nazero <- function(ksed, kdif = NULL, retain.zeroLSDs = FALSE, zero.tolerance = .Machine$double.eps ^ 0.5)
{
  #remove NA and zero values
  which.na <- is.na(ksed)
  ksed <- ksed[!which.na]
  if (!is.null(kdif))
    kdif <- kdif[!which.na]
  if (!retain.zeroLSDs)
  {
    which.nonzero <- findNonzero(ksed, zero.tolerance = zero.tolerance)
    if (is.null(which.nonzero))
    {
      ksed <- 0
      if (!is.null(kdif))
        kdif <- 0
    } else
    {
      ksed <- ksed[which.nonzero]
      if (!is.null(kdif))
        kdif <- kdif[which.nonzero]
    }
  }
  return(list(ksed = ksed, kdif = kdif)) 
}

#sed should be a vector that has had NAs and zeroes removed
"LSDaccmeas" <- function(ksed, assignedLSD, t.value, LSDaccuracy = "maxAbsDeviation")
{
  #Determine the accuracy of the assigned LSD 
  if (is.na(assignedLSD) || all(is.na(ksed)))
    accuracyLSD <- NA_real_
  else
  {
    if (LSDaccuracy == "maxAbsDeviation")
      accuracyLSD <- max(abs(t.value*ksed - assignedLSD))/assignedLSD
    else
    {
      if (LSDaccuracy == "maxDeviation")
        accuracyLSD <- max(t.value*ksed - assignedLSD)/assignedLSD
      else
      {
        if (LSDaccuracy == "q90Deviation")
          accuracyLSD <- quantile(t.value*ksed - assignedLSD, 0.90)/assignedLSD
        else
        {
          if (LSDaccuracy == "RootMeanSqDeviation")
            accuracyLSD <- sqrt(mean((t.value*ksed - assignedLSD)*
                                       (t.value*ksed - assignedLSD))) / assignedLSD 
          else
            accuracyLSD <- NULL
        }
      }
    }
  }
  return(accuracyLSD)
}

#ksed and kdif should be vectors that have had NAs and zeroes removed
falseSignif <- function(ksed, kdif, assignedLSD, t.value)
{
  sig.actual <- abs(kdif) >= t.value * ksed
  sig.approx <- abs(kdif) >= assignedLSD
  falsesig <- c(sum(!sig.actual & sig.approx, na.rm = TRUE), 
                sum(sig.actual & !sig.approx, na.rm = TRUE))
  names(falsesig) <- c("false.pos", "false.neg")
  return(falsesig)
}

#sed should be a vector that has had NAs and zeroes removed
"LSDstats" <- function(ksed, kdif, t.value, LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation") 
{
  #calculate LSD statistics
  stats <- data.frame(c = length(ksed),
                      minLSD = t.value * min(ksed),
                      meanLSD = t.value * sqrt(mean(ksed*ksed)),
                      maxLSD = t.value * max(ksed),
                      assignedLSD = NA_real_,
                      accuracyLSD = NA_real_,
                      falsePos = NA_real_,
                      falseNeg = NA_real_)
  #Set assigned LSD for use with halfLSIs
  if (LSDstatistic == "median")
    stats$assignedLSD <- t.value * median(ksed)
  else
  {
    if (LSDstatistic == "q10")
      stats$assignedLSD <- t.value * quantile(ksed, probs = 0.10)
    else
    {
      if (LSDstatistic == "q25")
        stats$assignedLSD <- t.value * quantile(ksed, probs = 0.25)
      else
      {  
        if (LSDstatistic == "q75")
          stats$assignedLSD <- t.value * quantile(ksed, probs = 0.75)
        else
          {
            if (LSDstatistic == "q90")
              stats$assignedLSD <- t.value * quantile(ksed, probs = 0.90)
            else
              stats$assignedLSD <- stats[[LSDstat2name(LSDstatistic)]]
          }
        }
    }
  }
   
  #Determine the accuracy of the assigned LSD 
  stats$accuracyLSD <- LSDaccmeas(ksed = ksed, assignedLSD = stats$assignedLSD, 
                                  t.value = t.value, LSDaccuracy = LSDaccuracy)
  
  #Calculate the number of false positives and negatives
  falsesig <- falseSignif(ksed = ksed, kdif = kdif, assignedLSD = stats$assignedLSD, 
                          t.value = t.value)
  stats$falsePos <- falsesig["false.pos"]
  stats$falseNeg <- falsesig["false.neg"]
  return(stats)
}

#Function to calculate the LSDs for combinations of the levels of the by factor(s)
sliceLSDs <- function(alldiffs.obj, by, t.value, LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
                      alpha = 0.05, which.stats = "all", 
                      retain.zeroLSDs = FALSE, zero.tolerance = .Machine$double.eps ^ 0.5)
{
  if (!all(by %in% names(alldiffs.obj$predictions)))
    stop("At least one element of LSDby is not in the predictions component of the alldiffs object\n")
#  classify <- attr(alldiffs.obj, which = "classify")
#  if (!all(unlist(lapply(by, grepl, x = classify, fixed = TRUE))))
#    stop("One of the elements of LSDby is not in the classify")
  
  denom.df <- attr(alldiffs.obj, which = "tdf")
  if (is.null(denom.df))
  {
    warning(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
                  "- p-values and LSDs not calculated"))
    LSDs <- NULL
  } else
  {
    t.value = qt(1-alpha/2, denom.df)
    sed <- alldiffs.obj$sed
    dif <- alldiffs.obj$differences
    diag(dif) <- NA_real_
    
    #Get the LSDs
    fac.comb <- fac.LSDcombs.alldiffs(alldiffs.obj, by)
    levs <- levels(fac.comb)
    if (length(LSDstatistic) ==  1)
      LSDstatistic <- rep(LSDstatistic, length(levs))
    if (length(levs) != length(LSDstatistic))
      stop("The length of LSDstatistic should be the same as the number of combinations of the LSDby variables")
    names(LSDstatistic) <- levs
    #loop over LSDby combinations
    LSDs <- lapply(levs, 
                   function(lev, sed, dif, t.value)
                   {
                     krows <- lev == fac.comb
                     if (length(fac.comb[krows]) == 1) #have a single prediction
                     {
                       warning(paste("LSD calculated for a single prediction",
                                     "- applies to two independent predictions with the same standard error"))
                       if (which.stats == "all")
                       {
                         stats <- c(0, rep(t.value * sqrt(2) * 
                                          alldiffs.obj$predictions$standard.error[krows], times = 4), NA_real_, NA_real_, NA_real_)
                         names(stats) <- c("c", "minLSD", "meanLSD", "maxLSD", "assignedLSD", "accuracyLSD", 
                                           "false.pos", "false.neg")
                       } else
                         if (which.stats == "evalLSD")
                           stats <- NA_real_
                         else
                           stop("Unknown which.stats option in SliceLSDs")
                     } else  #have several predictions
                     {
                       ksed <- getUpperTri(sed[krows, krows])
                       kdif <- getUpperTri(dif[krows, krows])
                       #remove NA and zero values
                       rm.list <- rm.nazero(ksed, kdif, retain.zeroLSDs = retain.zeroLSDs, 
                                            zero.tolerance = zero.tolerance)
                       ksed <- rm.list$ksed
                       kdif <- rm.list$kdif
                       if (which.stats == "all")
                         stats <- LSDstats(ksed = ksed, kdif = kdif, t.value, 
                                           LSDstatistic = LSDstatistic[lev], LSDaccuracy = LSDaccuracy)
                       else
                       {  
                         if (which.stats == "evalLSD")
                         {
                           stats <- LSDaccmeas(ksed, 
                                               assignedLSD = 
                                                 alldiffs.obj$LSD$assignedLSD[
                                                   rownames(alldiffs.obj$LSD) == lev], 
                                               t.value = t.value, LSDaccuracy = LSDaccuracy)
                           names(stats) <- "accuracyLSD"
                           #Calculate the number of false positives and negatives
                           falsesig <- falseSignif(ksed = ksed, kdif = kdif, 
                                                   assignedLSD = 
                                                     alldiffs.obj$LSD$assignedLSD[rownames(alldiffs.obj$LSD) 
                                                                                == lev], 
                                                   t.value = t.value)
                           stats <- c(stats,falsesig)
                         }
                         else
                           stop("Unknown which.stats option in SliceLSDs")
                       }
                     }
                     return(stats)
                   }, sed = sed, dif = dif, t.value = t.value)
    if (!is.null(LSDs))
    {
      LSDs <- as.data.frame(do.call(rbind, LSDs))
      rownames(LSDs) <- levs
    }
  }  
  return(LSDs)
}

#Function to calculate the Accuracy for each prediction
sliceAccs <- function(alldiffs.obj, by, LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
                      alpha = 0.05, retain.zeroLSDs = FALSE, zero.tolerance = .Machine$double.eps ^ 0.5)
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
    Acc <- NULL
  } else
  {
    t.value = qt(1-alpha/2, denom.df)
    
    #Get the Acc
    fac.comb <- fac.LSDcombs.alldiffs(alldiffs.obj, by)
    levs <- levels(fac.comb)
    Acc <- lapply(levs, 
                  function(lev, sed, t.value)
                  {
                    krows <- lev == fac.comb
                    assignedLSD <- alldiffs.obj$LSD$assignedLSD[rownames(alldiffs.obj$LSD)  == lev]
                    if (length(fac.comb[krows]) == 1)
                    {
                      warning(paste("LSD calculated for a single prediction",
                                    "- applies to two independent predictions with the same standard error"))
                      acc <- NA_real_
                      names(acc) <- lev
                    } else
                    {
                      ksed <- sed[krows, krows]
                      acc <-apply(ksed, MARGIN = 1, FUN =  
                                    function(sedrow, assignedLSD, t.value, LSDaccuracy, 
                                             retain.zeroLSDs, zero.tolerance)
                                    {
                                      sedrow <- as.vector(sedrow)
                                      sedrow <- rm.nazero(sedrow, retain.zeroLSDs = retain.zeroLSDs,
                                                           zero.tolerance = zero.tolerance)$ksed
                                      acc <- LSDaccmeas(sedrow, t.value = t.value, assignedLSD, 
                                                        LSDaccuracy = LSDaccuracy)
                                      return(acc)
                                    }, 
                                  assignedLSD = assignedLSD, t.value = t.value, LSDaccuracy = LSDaccuracy, 
                                  retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance)
                      acc <- unlist(acc)
                    }
                    return(acc)
                  }, sed = sed, t.value = t.value)
    if (!is.null(Acc))
      Acc <- unlist(Acc)
  }  
  return(Acc)
}
#Function to produce the per.prediction LSD stats 
LSDpred.stats <- function(LSD.mat, dif.mat, LSDstatistic, LSDaccuracy, t.value = 1, 
                          retain.zeroLSDs, zero.tolerance)
{
  perpred <- lapply(1:nrow(LSD.mat), FUN =  
                      function(k, LSD.mat, dif.mat, LSDstatistic, LSDaccuracy, t.value = 1, 
                               retain.zeroLSDs, zero.tolerance)
                      {
                        LSDrow <- as.vector(LSD.mat[k,])
                        difrow <- as.vector(dif.mat[k,])
                        #remove NA and zero values
                        rm.list <- rm.nazero(LSDrow, difrow, retain.zeroLSDs = retain.zeroLSDs, 
                                             zero.tolerance = zero.tolerance)
                        LSDrow <- rm.list$ksed
                        difrow <- rm.list$kdif
                        x <- LSDstats(LSDrow, difrow, t.value = t.value, LSDstatistic = LSDstatistic, 
                                      LSDaccuracy = LSDaccuracy)
                        return(x)
                      }, LSD.mat = LSD.mat, dif.mat = dif.mat, 
                    LSDstatistic = LSDstatistic, LSDaccuracy = LSDaccuracy, t.value = t.value, 
                    retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance)
  perpred <- do.call(rbind, perpred)
  return(perpred)
}

#Function to produce just the per.prediction accuracy for an assigned LSD value 
LSDpred.acc <- function(LSD.mat, assignedLSD, LSDaccuracy, t.value = 1, 
                        retain.zeroLSDs, zero.tolerance)
{
  perpred <- apply(LSD.mat, MARGIN = 1, FUN =  
                     function(LSDrow, assignedLSD, LSDaccuracy, t.value = 1, 
                              retain.zeroLSDs, zero.tolerance)
                     {
                       LSDrow <- as.vector(LSDrow)
                       LSDrow <- rm.nazero(LSDrow, retain.zeroLSDs = retain.zeroLSDs, 
                                           zero.tolerance = zero.tolerance)$ksed
                       x <- LSDaccmeas(LSDrow, t.value = t.value, assignedLSD = assignedLSD, 
                                       LSDaccuracy = LSDaccuracy)
                       return(x)
                     }, 
                   assignedLSD = assignedLSD, LSDaccuracy = LSDaccuracy, t.value = t.value, 
                   retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance)
  perpred <- unlist(perpred)
  return(perpred)
}

#Function to get all stats and accuracies for all kLSDs and kdifs supplied in vectors
"LSDallstats" <- function(kLSDs, kdifs, LSDaccuracy = "maxAbsDeviation", t.value = 1, 
                          retain.zeroLSDs = FALSE, 
                          zero.tolerance = .Machine$double.eps ^ 0.5) 
{
  LSDstat.labs <- c("min", "quant10", "quant25", "mean", "median", "quant75", "quant90", "max")
  
  #calculate LSD statistics
  c <- length(kLSDs)
  quants <- quantile(kLSDs, c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 1))
  stats <- c(c, quants[1:3], sqrt(mean(kLSDs*kLSDs)), quants[4:7])
  names(stats) <- c("c", LSDstat.labs)
  stats <- as.data.frame(as.list(stats))

  #Calculate the number of false positives and negatives
  falsesig <- lapply(LSDstat.labs, 
                     function(LSDstat, kLSDs, kdifs, stats, t.value)
                     {
                       fsig <- falseSignif(ksed = kLSDs, kdif = kdifs, assignedLSD = stats[[LSDstat]], 
                                   t.value = t.value)
                     }, kLSDs = kLSDs, kdifs = kdifs, stats = stats, t.value = t.value)
  if (!is.null(falsesig))
  {
    false.pos <- c(list(c = c), lapply(falsesig, function(comp) comp["false.pos"]))
    false.pos <- as.data.frame(do.call(cbind, false.pos))
    names(false.pos) <- c("c", LSDstat.labs)
    false.neg <- c(list(c = c), lapply(falsesig, function(comp) comp["false.neg"]))
    false.neg <- as.data.frame(do.call(cbind, false.neg))
    names(false.neg) <- c("c", LSDstat.labs)
  }
  
  #Determine the accuracy of the assigned LSD 
  Acc <- lapply(LSDstat.labs, 
                function(LSDstat, kLSDs, stats, LSDaccuracy)
                {
                  acc <- LSDaccmeas(ksed = kLSDs, assignedLSD = stats[[LSDstat]], 
                                    t.value = 1, LSDaccuracy = LSDaccuracy)
                }, kLSDs = kLSDs, stats = stats, LSDaccuracy = LSDaccuracy)
  Acc <- as.data.frame(c(list(length(kLSDs)), Acc))
  names(Acc) <- c("c", LSDstat.labs)
  return(list(statistics = stats, accuracy = Acc, false.pos = false.pos, false.neg = false.neg))
}

#Function to get stats, accuracy, per.pred.acc for sets of LSDs specified by the by argument
sliceAll <- function(alldiffs.obj, by, t.value, LSDaccuracy = "maxAbsDeviation", breaks, 
                     digits = 3, plotHistogram = FALSE, 
                     retain.zeroLSDs = FALSE, zero.tolerance = .Machine$double.eps ^ 0.5)
{
  LSDstat.labs <- c("min", "quant10", "quant25", "mean", "median", "quant75", "quant90", "max")
  
  if (!all(by %in% names(alldiffs.obj$predictions)))
    stop("At least one element of LSDby is not in the predictions component of the alldiffs object\n")
#  classify <- attr(alldiffs.obj, which = "classify")
#  if (!all(unlist(lapply(by, grepl, x = classify, fixed = TRUE))))
#    stop("One of the elements of LSDby is not in the classify")
  
  LSDs <-t.value * alldiffs.obj$sed
  difs <- alldiffs.obj$differences
  
  #Get the LSDs
  fac.comb <- fac.LSDcombs.alldiffs(alldiffs.obj, by)
  levs <- levels(fac.comb)
  #loop over LSDby combinations
  LSDsall <- lapply(levs, 
                    function(lev, LSDs, t.value)
                    {
                      krows <- lev == fac.comb
                      if (length(fac.comb[krows]) == 1) #have a single prediction
                      {
                        warning(paste("LSD calculated for a single prediction",
                                      "- applies to two independent predictions with the same standard error"))
                        stats <- c(0, rep(t.value * sqrt(2) * alldiffs.obj$predictions$standard.error[krows], 
                                     times = length(LSDstat.labs)))
                        names(stats) <- c("c", LSDstat.labs)
                        stats <- as.data.frame(as.list(stats))
                        acc <- c(0, rep(NA_real_, times = length(LSDstat.labs)))
                        names(acc) <- c("c", LSDstat.labs)
                        acc <- as.data.frame(as.list(acc))
                        fneg <- fpos <- acc
                        per.pred <- acc[,-1]
                        names(per.pred) <- LSDstat.labs
                        stats <- list(statistics = stats, accuracy = acc, false.pos = fpos, false.neg = fneg, 
                                      per.predictions = per.pred)
                      } else  #have several predictions
                      {
                        kLSDs <- LSDs[krows, krows]
                        kdifs <- difs[krows, krows]
                        #remove NA and zero values
                        rm.list <- rm.nazero(getUpperTri(kLSDs), getUpperTri(kdifs), 
                                             retain.zeroLSDs = retain.zeroLSDs, 
                                             zero.tolerance = zero.tolerance)
                        kLSDs.vec <- rm.list$ksed
                        kdifs.vec <- rm.list$kdif
                        distinct <- sort(unique(signif(kLSDs.vec, digits = digits)))
                        stats <- LSDallstats(kLSDs.vec, kdifs.vec, t.value = 1, LSDaccuracy = LSDaccuracy)
                        
                        #get per.prediction accuracies
                        per.pred <- lapply(LSDstat.labs, 
                                           function(LSDstatistic, kLSDs, stats, LSDaccuracy, t.value, 
                                                    retain.zeroLSDs, zero.tolerance) 
                                           {
                                             predacc <- LSDpred.acc(kLSDs, 
                                                                    assignedLSD = stats$statistics[[LSDstatistic]], 
                                                                    LSDaccuracy = LSDaccuracy, t.value = t.value,
                                                                    retain.zeroLSDs = retain.zeroLSDs, 
                                                                    zero.tolerance = zero.tolerance)
                                           }, kLSDs = kLSDs, stats = stats, LSDaccuracy = LSDaccuracy, t.value = 1,
                                           retain.zeroLSDs = retain.zeroLSDs, zero.tolerance = zero.tolerance)
                        per.pred <- as.data.frame(do.call(cbind, per.pred))
                        names(per.pred) <- LSDstat.labs

                        #Get the frequencies
                        kLSD.dat <- as.data.frame(kLSDs[upper.tri(kLSDs)])
                        names(kLSD.dat) <- "LSD"
                        freq <- hist(kLSD.dat$LSD, plot = FALSE, include.lowest = TRUE, breaks = breaks)
                        stats <- c(stats, list(distinct = distinct, frequencies = freq, per.prediction = per.pred))
                      }
                      return(stats)
                    }, LSDs = LSDs, t.value = t.value)
  names(LSDsall) <- levs
  if (!is.null(LSDsall))
  {
    stats <- lapply(LSDsall, function(comp) comp$statistics)
    stats <- as.data.frame(do.call(rbind, stats))
    acc <- lapply(LSDsall, function(comp) comp$accuracy)
    acc <- as.data.frame(do.call(rbind, acc))
    false.pos <- lapply(LSDsall, function(comp) comp$false.pos)
    false.pos <- as.data.frame(do.call(rbind, false.pos))
    false.neg <- lapply(LSDsall, function(comp) comp$false.neg)
    false.neg <- as.data.frame(do.call(rbind, false.neg))
    freq.labs <- as.character(LSDsall[[1]]$frequencies$mids)
    freq <- lapply(LSDsall, function(comp) comp$frequencies$counts)
    freq <- as.data.frame(do.call(rbind, freq))
    names(freq) <- freq.labs
    distinct <- lapply(LSDsall, function(comp) comp$distinct)
    per.pred <- lapply(LSDsall, function(comp) comp$per.prediction)
    per.pred <- as.data.frame(do.call(rbind, per.pred))
    rownames(per.pred) <- rownames(LSDs)
  }
  
  if (plotHistogram)
  {
    LSD.dat <- lapply(levs, 
                      function(lev, LSDs)
                      {
                        krows <- lev == fac.comb
                        kLSDs <- LSDs[krows, krows]
                        kLSDs <- as.data.frame(kLSDs[upper.tri(kLSDs)])
                        kLSDs <- cbind(lev, kLSDs)
                        names(kLSDs) <- c("Combination", "LSD")
                        return(kLSDs)
                      }, LSDs = LSDs)
    LSD.dat <- do.call(rbind, LSD.dat)
    plt <- ggplot(LSD.dat, aes(x = .data[["LSD"]])) + 
      geom_histogram(breaks = breaks) + 
      theme_bw() + 
      facet_grid(rows = eval(parse(text="vars(Combination)")))
    print(plt)    
  }
  
  return(list(frequencies = freq, distinct.vals = distinct, statistics = stats, accuracy = acc, 
              false.pos = false.pos, false.neg = false.neg, per.pred.accuracy = per.pred))
}
