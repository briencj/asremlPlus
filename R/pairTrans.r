#'## Compute CIs for prediction ratios
#'### Functions for computing Ratio CIs 
FiellerRatioCI <- function(a, b, V, t)
{
  theta <- a/b
  g <- (t^2)*V[2,2]/b^2
  se <- sqrt(V[1,1] - 2*theta*V[1,2] + theta^2 * V[2,2] - g*(V[1,1]-V[1,2]^2/V[2,2]))
  LCL <- (1/(1-g))*(theta- g*V[1,2]/V[2,2] - t/b * se)
  UCL <- (1/(1-g))*(theta- g*V[1,2]/V[2,2] + t/b * se)
  if (!is.na(UCL) & !is.na(LCL) & UCL < LCL) {
    UCL <- NA; LCL <- NA
  }
  return(list(pred.ratio=theta, upper.Confidence.limit = UCL, lower.Confidence.limit = LCL))
}

makeContrastMatrix <- function(preds, indx, pairs.factor, first.level, second.level)
{
  X <- model.matrix(~ Combinations-1, data = preds)
  ## Form the contrast matrix
  #Determine the columns of X for each contrast
  colsets <- split(colnames(X), preds[indx], lex.order = TRUE)
  prediv <- split(preds, preds[indx], lex.order = TRUE)
  locns <- lapply(colsets, 
                  function(x, cols = colnames(X)) 
                  { pair <-  na.omit(match(x, cols)); names(pair) <- x; return(pair) })
  contr <- mapply(function(pred, cols, pairs.factor, first.level, second.level) 
  {
    names(cols) <- gsub("Combinations", "", names(cols), fixed = TRUE)
    #    k1 <- (1:length(cols))[grepl(first.level, names(cols), fixed = TRUE)]
    #    k2 <- (1:length(cols))[grepl(second.level, names(cols), fixed = TRUE)]
    k1 <- (1:length(cols))[pred[pairs.factor] == first.level]
    k2 <- (1:length(cols))[pred[pairs.factor] == second.level]

    if (any(c(length(k1),length(k2)) == 0)) #no contrast
      cols <- NA
    else
    {
      if (length(k1) != 1)
        stop("There is more than one observation with the level ", first.level, " in combination ", names(cols))
      else
      {
        if (length(k2) != 1)
          stop("There is more than one observation with the level ", second.level, " in combination ", names(cols))
        else
        {
          cols[1:length(cols)] <- 0
          cols[k1] <- 1
          cols[k2] <- -1
        }
      }
    }
    return(cols)
  }, prediv, locns, MoreArgs = list(pairs.factor = pairs.factor, 
                                   first.level = first.level, second.level = second.level), 
  SIMPLIFY = FALSE)
  
  
  #Substitute the contrast in one column and NA in the other for each contrast
  for (k in 1:length(locns))
  {
    for (j in 1:length(locns[[k]]))
    {
      X[locns[[k]][j], locns[[k]][1]] <- contr[[k]][j]
      if (j != 1)
        X[locns[[k]][j], locns[[k]][j]] <- 0
    }
  }
  
  #Form contrast matrix
  X <- X[, unlist(locns)[order(unlist(locns))]]
  which.cols <- unlist(lapply(1:ncol(X), 
                              function(k) 
                              {
                                if   (any(is.na(as.vector(X[,k])))) retain <- FALSE 
                                else retain <- !all(X[,k] == 0)
                              }))
  L <- t(X[, which.cols])
#  lev <- paste0("Combinations",second.level)
  lev <- paste0("Combinations",levels(preds[[pairs.factor]])[1])
  rownames(L) <- gsub(paste0(",",lev,","), "", rownames(L), fixed = TRUE)
  rownames(L) <- gsub(paste0(",",lev), "", rownames(L), fixed = TRUE)
  rownames(L) <- gsub(paste0(lev,","), "", rownames(L), fixed = TRUE)
  return(L)
}


fac.split <- function(factor, new.factors, sep = ",")
{
  fac.sep <- strsplit(factor, split = sep)
  #name the values in each elements of fac.sep 
  fac.sep <- lapply(fac.sep, 
                    function(obs, new.factors)
                    {
                      names(obs) <- names(new.factors); return(obs)
                    }, new.factors = new.factors)
  #form separate columns for the values in each element of fac.sep
  new.facs <- lapply(names(new.factors), 
                     function(elem)
                     {
                       unlist(lapply(fac.sep, function(obs, elem) obs[elem], 
                                     elem = elem))
                     })
  names(new.facs) <- names(new.factors)
  df <- as.data.frame(new.facs, stringsAsFactors = FALSE)
  for ( fac in names(new.factors))
  {
    df[fac] <- factor(df[[fac]], levels = new.factors[[fac]])
    df[fac] <- factor(df[[fac]]) #make sure no unused levels
  }
  return(df)
}


ratioTransform.alldiffs <- function(alldiffs.obj, ratio.factor, 
                                    numerator.levels, denominator.levels, 
                                    method = "Fieller", alpha = 0.05, 
                                    response = NULL, response.title = NULL, 
                                    tables = "predictions", 
                                    ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  alldiffs.obj <- renameDiffsAttr(alldiffs.obj)
  
  
  meth.options <- c("Fieller")
  method.opt <- meth.options[check.arg.values(method, meth.options)]
  options <- c("none", "predictions", "backtransforms", "vcov", 
               "differences", "p.differences", "sed", "LSD", "all")
  opt <- options[unlist(lapply(tables, check.arg.values, options=options))]

  #Get a t-value (use z if no tdf attribute)
  tdf <- attr(alldiffs.obj, which = "tdf")
  if (!is.null(tdf))
  {
    t <- qt(1 - alpha/2, df = tdf)
  } else
    t <- qnorm(1 - alpha/2)
  
  if (is.null(response))
    response <- attr(alldiffs.obj, which = "response")
  if (is.null(response.title))
    response.title <- attr(alldiffs.obj, which = "response.title")
  classify <- attr(alldiffs.obj, which = "classify")

  facs <- fac.getinTerm(classify, rmfunction = TRUE)
  alldiffs.obj <- renewClassify(alldiffs.obj, newclassify = classify) #make sure in standard order
  tmp <- alldiffs.obj$predictions
  tmp$Combinations <- dae::fac.combine(as.list(tmp[facs]), combine.levels = TRUE, sep = ",")
  if (length(ratio.factor) != 1 | !(ratio.factor %in% facs))
    stop("ratio factor must specify a single factor that is in the classify attribute of the alldiffs object")
  if (!all(c(numerator.levels, denominator.levels) %in% levels(tmp[[ratio.factor]])))
    stop("Not all numerator.levels and denominator.levels are levels in ", ratio.factor, " in the alldiffs object")
  sortFactor <- attr(alldiffs.obj, which = "sortFactor")
  
  indx <- setdiff(facs, ratio.factor)

  # if (length(indx) == 1)
  # {
  #   ratios.dat <- as.data.frame(ratios.dat, row.names = NULL)
  #   names(ratios.dat)[1] <- indx
  # }
  tmp <- split(tmp, tmp[indx], sep = ",", lex.order = TRUE)
  Ratios <- lapply(denominator.levels, function(denom.lev,tmp) 
  {
    lapply(numerator.levels, function(num.lev, denom.lev, tmp) 
    { 
      
      if (!(num.lev %in% alldiffs.obj$predictions[[ratio.factor]]) || 
          !(denom.lev %in% alldiffs.obj$predictions[[ratio.factor]]) ||
          num.lev == denom.lev)
        ratioCIs <-NULL
      else
      {
        ratioCIs <- do.call(rbind, lapply(tmp, function(preds, num.lev, denom.lev, V, tval)
        {
          indx.vals <- (lapply(preds[indx], function(x) as.character(x)[1])) #current levels of indx factors
          ratio.pair <- c(indx.vals, rep(NA,3))
          names(ratio.pair) <- c(indx, "pred.ratio", "upper.Confidence.limit", "lower.Confidence.limit")
          if (all(c(num.lev, denom.lev) %in% preds[[ratio.factor]]))
          {
            a <- preds[preds[ratio.factor] == num.lev, "predicted.value"]
            b <- preds[preds[ratio.factor] == denom.lev, "predicted.value"]
            labs <- preds$Combinations
            posns <- unlist(lapply(c(num.lev,denom.lev), grep, x =labs, fixed = TRUE))
            lab <- labs[posns]
            Vpair <- V[lab,lab]
            ratio.pair <- as.data.frame(c(indx.vals, (FiellerRatioCI(a, b, Vpair, t = t))))
          }
          return(ratio.pair)
        }, V = alldiffs.obj$vcov, num.lev = num.lev, denom.lev = denom.lev))
        
        rownames(ratioCIs) <- NULL
        names(ratioCIs)[match("pred.ratio", names(ratioCIs))] <- "predicted.value"
        ratioCIs$standard.error <- NA
        ratioCIs$est.status <- "Estimable"
        ratio.names <- names(ratioCIs)
        n <- match("standard.error", ratio.names)
        k <- match("upper.Confidence.limit", ratio.names)
        if (n-k+1 != 3)
          stop("Standard.error and confidence limits not together in ratios data frame")
        ratio.names[k:n] <- c("standard.error", "upper.Confidence.limit", "lower.Confidence.limit")
        ratioCIs <- ratioCIs[ratio.names]
        ratioCIs[indx] <-  lapply(indx, 
                                  function(fac) new.fac <- factor(ratioCIs[[fac]], 
                                                                  levels = levels(alldiffs.obj$predictions[[fac]])))
        ratioCIs <- as.predictions.frame(ratioCIs)
        class(ratioCIs) <- c("predictions.frame", "data.frame")
        attr(ratioCIs, which = "response") <- response
        attr(ratioCIs, which = "response.title") <- response.title
        
        #sort the predictions.frame, if it was previously sorted
        if (!is.null(sortFactor))
        {
          if (sortFactor != ratio.factor)
            ratioCIs <- sort(ratioCIs, sortFactor = sortFactor, classify = paste(indx, collapse = ":"),
                             sortOrder = attr(alldiffs.obj, which = "sortOrder"))
        }
      }
      return(ratioCIs)
    }, tmp = tmp, denom.lev = denom.lev)
  }, tmp = tmp)
  Ratios <- unlist(Ratios, recursive = FALSE)
  nams <- outer(numerator.levels, denominator.levels, paste, sep = ",")
  names(Ratios) <- nams
  
  
  if ("predictions" == opt)
  {
    cat("\n\n#### Table(s) of ratio predictions and confidence limits for ", response, "\n\n")
    print(Ratios)
  }
  invisible(Ratios)                      
}
  
pairdiffsTransform.alldiffs <- function(alldiffs.obj, pairs.factor, first.levels, second.levels, 
                                        Vmatrix = FALSE, 
                                        error.intervals = "Confidence", 
                                        avsed.tolerance = 0.25, accuracy.threshold = NA, 
                                        LSDtype = "overall", LSDsupplied = NULL, LSDby = NULL, 
                                        LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
                                        response = NULL, response.title = NULL, tables = "all", 
                                        pairwise = TRUE, alpha = 0.05, ...)
{
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(alldiffs.obj)  
  if (is.character(validalldifs))
    stop(validalldifs)
  alldiffs.obj <- renameDiffsAttr(alldiffs.obj)
  attr(alldiffs.obj, which = "alpha") <- alpha

  #get transform attributes from backtransforms
  transform.power = 1; offset <- 0; scale <- 1
  if (!is.null(alldiffs.obj$backtransforms))
  {
    transform.power = attr(alldiffs.obj$backtransforms, which = "transform.power")
    offset = attr(alldiffs.obj$backtransforms, which = "offset")
    scale = attr(alldiffs.obj$backtransforms, which = "scale")
  } 
  
  sortFactor <- attr(alldiffs.obj, which = "sortFactor")
  
  int.options <- c("none", "Confidence", "StandardError", "halfLeastSignificant")
  int.opt <- int.options[check.arg.values(error.intervals, int.options)]
  options <- c("none", "predictions", "backtransforms", "vcov", 
               "differences", "p.differences", "sed", "LSD", "all")
  opt <- options[unlist(lapply(tables, check.arg.values, options=options))]
  #Get a t-value (use z if no tdf attribute)
  tdf <- attr(alldiffs.obj, which = "tdf")
  if (!is.null(tdf))
  {
    t <- qt(1 - alpha/2, df = tdf)
  } else
    t <- qnorm(1 - alpha/2)
  
  if (is.null(response))
    response <- attr(alldiffs.obj, which = "response")
  if (is.null(response.title))
    response.title <- attr(alldiffs.obj, which = "response.title")
  classify <- attr(alldiffs.obj, which = "classify")
  facs <- fac.getinTerm(classify, rmfunction = TRUE)
  alldiffs.obj <- renewClassify(alldiffs.obj, newclassify = classify) #make sure in standard order
  tmp <- alldiffs.obj$predictions
  tmp$Combinations <- dae::fac.combine(as.list(tmp[facs]), combine.levels = TRUE, sep = ",")
  if (length(pairs.factor) != 1 | !(pairs.factor %in% facs))
    stop("pairs factor must specify a single factor that is in the classify attribute of the alldiffs object")
  if (!all(c(first.levels, second.levels) %in% levels(tmp[[pairs.factor]])))
    stop("Not all first.levels and second.levels are levels in ", pairs.factor, " in the all diffs object")
  
  indx <- setdiff(facs, pairs.factor)

  Diffs <- lapply(second.levels, function(second.lev,tmp) 
  {
    lapply(first.levels, function(first.lev, second.lev, tmp) 
    { 
      
      if (!(first.lev %in% alldiffs.obj$predictions[[pairs.factor]]) || 
          !(second.lev %in% alldiffs.obj$predictions[[pairs.factor]]))
        diffs <-NULL
      else
      {
        L <- makeContrastMatrix(tmp, indx = indx, pairs.factor = pairs.factor, first.lev, second.lev)
  
        #form the factors indexing the pair-difference predictions
        new.facs <- lapply(indx, function(fac, data)
        {
          levels(tmp[[fac]])
        }, data = tmp)
        names(new.facs) <- indx
        pairs.dat <- fac.split(rownames(L), new.factors = new.facs, sep = ",")
  
        #Calculate the differences for the current pair
        diffs <- linTransform(alldiffs.obj, classify = classify, 
                              linear.transformation = L, Vmatrix = Vmatrix, 
                              error.intervals = "Confidence",
                              avsed.tolerance = avsed.tolerance, accuracy.threshold = accuracy.threshold, 
                              response = response, response.title = response.title, 
                              pairwise = pairwise, alpha = alpha, 
                              tables = "none", ...)
        diffs$predictions <- cbind(pairs.dat, diffs$predictions)
        diffs$predictions <- diffs$predictions[, -match("Combination", names(diffs$predictions))]
        diffs <- renewClassify(diffs, newclassify = paste(indx, collapse = ":"))
        diffs <- redoErrorIntervals(diffs,  error.intervals = int.opt, alpha = alpha, 
                                    avsed.tolerance = avsed.tolerance, accuracy.threshold = accuracy.threshold, 
                                    LSDtype = LSDtype, LSDby = LSDby, LSDsupplied = LSDsupplied, 
                                    LSDstatistic = LSDstatistic, LSDaccuracy = LSDaccuracy)
        # if (!is.null(diffs$backtransforms))
        # {
        #   diffs$backtransforms <- cbind(pairs.dat, diffs$backtransforms)
        #   diffs$backtransforms <- diffs$backtransforms[, -match("Combination", names(diffs$backtransforms))]
        #   #Find missing attributes in new alldiffs.obj and add them back in 
        #   newattr <- attributes(diffs$backtransforms)
        #   back.attr <- attributes(alldiffs.obj$backtransforms)
        #   back.attr <- back.attr[names(back.attr)[!(names(back.attr) %in% names(newattr))]]
        #   if (length(back.attr) > 0)
        #   {
        #     newattr <- c(newattr,back.attr)
        #     attributes(diffs$backtransforms) <- newattr
        #   }
        # }
        
        #sort the alldiffs if it was previously sorted
        if (!is.null(sortFactor))
        {
          if (sortFactor != pairs.factor)
            diffs <- sort(diffs, sortFactor = sortFactor,
                             sortOrder = attr(diffs, which = "sortOrder"))
        }
        
        if (opt != "none")
          print(diffs, which = opt)
        return(diffs)
      }
    }, tmp = tmp, second.lev = second.lev)
  }, tmp = tmp)
  Diffs <- unlist(Diffs, recursive = FALSE)
  nams <- outer(first.levels, second.levels, paste, sep = ",")
  names(Diffs) <- nams
  
  invisible(Diffs)                      
}
