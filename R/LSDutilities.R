#### LSD functions
#Functionto convert a LSDstatistic option into a column name in the LSD component
LSDstat2name <- function(LSDstat)
  LSDname <- paste0(gsub("imum", "", LSDstat, fixed = TRUE), "LSD")

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
                       if (class(x)!="factor")
                         x <- factor(x)
                       return(x)
                     })
  fac.comb <- fac.combine(fac.list, combine.levels = TRUE)
  if (length(fac.comb) != nrow(alldiffs.obj$sed))
    stop("Variable(s) in (LSD)by argument are not the same length as the order of the sed matrix")
  return(fac.comb)
}

#Check LSDsupplied and make sure in format for inclusion as an LSD component
# makeLSDsupplied <- function(alldiffs.obj, LSDsupplied, LSDby, denom.df, alpha, zero.tolerance)
# {
#   if (is.null(LSDsupplied))
#     stop("Must set LSDsupplied for LSDtype set to supplied")
#   
#   if (inherits(LSDsupplied, what = "data.frame"))
#   {
#     if (ncol(LSDsupplied) > 1)
#     {
#       if (ncol(LSDsupplied) != length(LSDby)+1 || names(LSDsupplied)[1:length(LSDby)] != LSDby)
#         stop("If LSD supplied as a data.frame with more than one column ", 
#              "it should have a column for each element of LSDby followed by a column of LSD values")
#       fac.comb <- fac.combine(as.list(alldiffs.obj$predictions[LSDby]), combine.levels = TRUE)
#       LSDnam <- names(LSDsupplied)[ncol(LSDsupplied)]
#       LSDsupplied <- LSDsupplied[[LSDnam]]
#       names(LSDsupplied) <- levels(fac.comb)
#     } else #one-column LSDsupplied
#     {
#       if (is.null(rownames(LSDsupplied)))
#         stop("Must supply LSDby combination as row names to LSDsupplied")
#       #Convert any non-factors and form levels combination for which a LSDs are required
#       if (is.null(LSDby))
#       {
#         if (rownames(LSDsupplied) != "overall")
#           stop("A single row should have row name equal to overall")
#       } else
#       {
#         fac.comb <- fac.LSDcombs.alldiffs(alldiffs.obj, by = LSDby)
#         levs <- levels(fac.comb)
#         if (any(levs != rownames(LSDsupplied)))
#           stop("The combinations of the LSDby variables for LSDsupplied and those in the predictions do not match")
#         LSDsupplied <- LSDsupplied[,1]
#         names(LSDsupplied) <- levels(fac.comb)
#       }
#     }
#   } else #a numeric
#   {
#     if (is.null(names(LSDsupplied)))
#       stop("Must supply LSDby combinations or overall as names for LSDsupplied")
#     #Convert any non-factors and form levels combination for which a LSDs are required
#     if (is.null(LSDby))
#     {
#       if (names(LSDsupplied) != "overall")
#         stop("A single value should have name overall")
#      # LSDsupplied <- data.frame(meanLSD = LSDsupplied, row.names = "overall")
#       names(LSDsupplied) <- "overall"
#     } else
#     {
#       levs <- levels(fac.LSDcombs.alldiffs(alldiffs.obj, by = LSDby))
#       if (any(levs != rownames(LSDsupplied)))
#         stop("The combinations of the LSDby variables for LSDsupplied and those in the predictions do not match")
#     }
#   }
#   
#   #Put the supplied LSD into an LSD data.frame
#   #Check if tdf available
#   if (is.null(denom.df))
#   {
#     warning("Only supplied LSD values stored in LSD component")
#     minLSD <- rep(NA, nrow(LSDsupplied))
#     maxLSD <- rep(NA, nrow(LSDsupplied))
#     meanLSD <- LSDsupplied
#   } else
#   {
#     t.value = qt(1-alpha/2, denom.df)
#     if (is.null(LSDby))
#     {
#       LSDs<- LSDstats(alldiffs.obj$sed, t.value)
#       rownames(LSDs) <- "overall"
#     }
#     else
#       LSDs <- sliceLSDs(alldiffs.obj, by = LSDby, t.value = t.value, alpha = alpha,
#                         zero.tolerance = zero.tolerance)
#     minLSD <- LSDs["minLSD"]
#     maxLSD <- LSDs["maxLSD"]
#     meanLSD <- LSDsupplied
#   }
#   alldiffs.obj$LSD <- data.frame(minLSD  = minLSD,
#                                  meanLSD = LSDsupplied,
#                                  maxLSD = maxLSD)
#   rownames(alldiffs.obj$LSD) <- rownames(LSDs)
#   attr(alldiffs.obj, which = "LSDtype") <- "supplied"
#   attr(alldiffs.obj, which = "LSDby") <- LSDby
#   attr(alldiffs.obj, which = "LSDstatistic") <- "mean"
#   return(alldiffs.obj)
# }

#Check LSDsupplied an make sure in format for inclusion as an LSD component
addLSDsupplied <- function(alldiffs.obj, LSDsupplied, LSDby, denom.df, alpha, zero.tolerance)
{
  if (is.null(LSDsupplied))
    stop("Must set LSDsupplied for LSDtype set to supplied")
  
  if (inherits(LSDsupplied, what = "data.frame"))
  {
    if (ncol(LSDsupplied) > 1)
    {
      if (ncol(LSDsupplied) != length(LSDby)+1 || names(LSDsupplied)[1:length(LSDby)] != LSDby)
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
  combs <- dae::fac.genfactors(alldiffs.obj$predictions[LSDby])
  combs <- cbind(fac.comb = fac.combine(combs, combine.levels = TRUE), combs)
  combs <- merge(data.frame(fac.comb = rownames(alldiffs.obj$LSD)),
                 combs, all.x = TRUE, sort = FALSE)
  combs <- combs[-match("fac.comb", names(combs))]
  alldiffs.obj$LSD <- cbind(combs, alldiffs.obj$LSD)
  return(alldiffs.obj)
}

checkLSD <- function(alldiffs.obj)
{
  if (!"assignedLSD" %in% names(alldiffs.obj$LSD))
  { 
    alldiffs.obj$LSD$assignedLSD <- alldiffs.obj$LSD$meanLSD
    alldiffs.obj$LSD$accuracyLSD <- NA
  }
  return(alldiffs.obj)
}

"LSDaccmeas" <- function(ksed, assignedLSD, t.value, LSDaccuracy = "maxAbsDeviation")
{
  ksed <- na.omit(as.vector(ksed))
  #Determine the accuracy of the assigned LSD 
  if (is.na(assignedLSD) || all(is.na(ksed)))
    accuracyLSD <- NA
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
        if (LSDaccuracy == "90Deviation")
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

"LSDstats" <- function(sed, t.value, LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
                       zero.tolerance = 1e-10) 
{
  ksed <- na.omit(as.vector(sed))
#  ksed <- na.omit(ksed *ksed)
  max.ksed <- max(ksed)
  #Retain only nonzero variances
  if (max.ksed > zero.tolerance && sum(ksed/max.ksed > zero.tolerance) > 0)
    ksed <- ksed[ksed/max.ksed > zero.tolerance]
  else if (max.ksed < zero.tolerance)
    ksed <- 0
  
  #calculate LSD statistics
  stats <- data.frame(minLSD = t.value * min(ksed),
                      meanLSD = t.value * sqrt(mean(ksed*ksed)),
                      maxLSD = t.value * max(ksed),
                      assignedLSD = NA,
                      accuracyLSD = NA)
  medianLSD <- t.value * median(ksed)
  LSDname <- LSDstat2name(LSDstatistic)
  #Set assigned LSD for use with halfLSIs
  if (LSDstatistic != "median")
    stats$assignedLSD <- stats[[LSDname]]
  else
  {
    stats$assignedLSD <- medianLSD
  }
  
  #Determine the accuracy of the assigned LSD 
  stats$accuracyLSD <- LSDaccmeas(ksed = ksed, assignedLSD = stats$assignedLSD, 
                                  t.value = t.value, LSDaccuracy = LSDaccuracy)
  return(stats)
}

#Function to calculate the LSDs for combinations of the levels of the by factor(s)
sliceLSDs <- function(alldiffs.obj, by, t.value, LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
                      alpha = 0.05, which.stats = "all", zero.tolerance = 1E-04)
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
    
    #Get the LSDs
    fac.comb <- fac.LSDcombs.alldiffs(alldiffs.obj, by)
    levs <- levels(fac.comb)
    #loop over LSDbt combinations
    LSDs <- lapply(levs, 
                   function(lev, sed, t.value)
                   {
                     krows <- lev == fac.comb
                     if (length(fac.comb[krows]) == 1) #have a single prediction
                     {
                       warning(paste("LSD calculated for a single prediction",
                                     "- applies to two independent predictions with the same standard error"))
                       if (which.stats == "all")
                       {
                         stats <- c(rep(t.value * sqrt(2) * 
                                          alldiffs.obj$predictions$standard.error[krows] /2, times = 4), NA)
                         names(stats) <- c("minLSD", "meanLSD", "maxLSD", "assignedLSD", "accuracyLSD")
                       } else
                         if (which.stats == "accuracyLSD")
                           stats <- NA
                         else
                           stop("Unknown which.stats option in SliceLSDs")
                     } else  #have several predictions
                     {
                       ksed <- sed[krows, krows]
                       ksed <- na.omit(as.vector(ksed))
                       if (which.stats == "all")
                         stats <- LSDstats(ksed, t.value, LSDstatistic = LSDstatistic, LSDaccuracy = LSDaccuracy)
                       else
                       {  
                         if (which.stats == "accuracyLSD")
                           stats <- LSDaccmeas(ksed, 
                                               assignedLSD = alldiffs.obj$LSD$assignedLSD[rownames(alldiffs.obj$LSD) 
                                                                                          == lev], 
                                               t.value = t.value, LSDaccuracy = LSDaccuracy)
                         else
                           stop("Unknown which.stats option in SliceLSDs")
                       }
                     }
                     return(stats)
                   }, sed = sed, t.value = t.value)
    if (!is.null(LSDs))
    {
      # LSDs <- cbind(levs,
      #               as.data.frame(do.call(rbind, LSDs)))
      LSDs <- as.data.frame(do.call(rbind, LSDs))
      rownames(LSDs) <- levs
    }
  }  
  return(LSDs)
}

#Function to calculate the LSDs for combinations of the levels of the by factor(s)
sliceAccs <- function(alldiffs.obj, by, LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
                      alpha = 0.05, zero.tolerance = 1E-04)
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
                      acc <- NA
                      names(acc) <- lev
                    } else
                    {
                      ksed <- sed[krows, krows]
                      acc <-apply(ksed, MARGIN = 1, FUN =  
                                    function(sedrow, assignedLSD, t.value, LSDstatistic, LSDaccuracy)
                                    {
                                      acc <- LSDaccmeas(sedrow, assignedLSD, t.value, LSDaccuracy = LSDaccuracy)
                                      return(acc)
                                    }, 
                                  assignedLSD = assignedLSD, t.value = t.value, LSDaccuracy = LSDaccuracy)
                      
                      acc <- unlist(acc)
                    }
                    return(acc)
                  }, sed = sed, t.value = t.value)
    if (!is.null(Acc))
      Acc <- unlist(Acc)
    #Acc <- do.call(rbind, Acc)$accuracyLSD
  }  
  return(Acc)
}

