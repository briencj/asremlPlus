#Function to explore the LSD values
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
                      minlsd <- min(kLSDs$lsd)
                      maxlsd <- max(kLSDs$lsd)
                      if (trace) {cat("\n\n#### New set\n"); print(c(minlsd,maxlsd))}
                      range <- maxlsd - minlsd
                      if (range > 0)
                      { 
                        #Do a grid search for the minLSD
                        testLSDs <- lapply(seq(0, 1, stepsize),
                                           function(step, range, kLSDs, kpos.wt)
                                           {
                                             testlsd <- minlsd + step * range
                                             false.vals <- c(testlsd, 
                                                             falseErrorNums(testlsd, kLSDs, 
                                                                            kpos.wt))
                                             names(false.vals)[1] <- "LSD"
                                             return(false.vals)
                                           }, range, kLSDs, kpos.wt)
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
                        if (optLSD$LSD > minlsd)
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
                        optLSD <- c(minlsd, falseErrorNums(minlsd, kLSDs, kpos.wt))
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
                     if (length(fac.comb[krows]) == 1) #have a single prediction
                     {
                       rm.list <- list (klsd = t.value * sqrt(2) * 
                                          alldiffs.obj$predictions$standard.error[krows],
                                        kdif = 0)
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
                         # rm.list <- list (lsd = t.value * sqrt(2) * 
                         #                    alldiffs.obj$predictions$standard.error[krows][1],
                         #                  dif = 0)
                         rm.list <- list (lsd = 0, dif = 0)
                       } else
                       {
                         #remove NA and zero values
                         rm.list <- rm.nazero(klsd, kdif, retain.zeroLSDs = retain.zeroLSDs,
                                              zero.tolerance = zero.tolerance)
                         names(rm.list) <- c("lsd", "dif")
                       }
                     }
                     rm.list$sig.actual <- abs(rm.list$dif) >= rm.list$lsd
                     return(rm.list)
                   }, lsd = lsd, dif = dif, t.value = t.value, alldiffs.obj = alldiffs.obj)
    names(LSDs) <- levs
  }
  return(LSDs)
}
