#'### Functions to plot LSD-values for a set of pairwise differences
#+
"plotLSDs.data.frame" <- function(object, LSD = "LSDs",  x, y, alpha = 0.05,
                                  triangles = "both", gridspacing = 0,  
                                  title = NULL, axis.labels = NULL, axis.text.size = 12, 
                                  colours = RColorBrewer::brewer.pal(3, "Set2"), 
                                  ggplotFuncs = NULL, printPlot = TRUE, ...)
  #Plots a data.frame of LSD-values that has columns x, y, LSD
{ 
  show.sig = FALSE #leave here in case wnat to indicate significant differences 
                   # - would need to add p.differences to data.frame.
  if (!all(c(x,y,LSD) %in% names(object)))
  {
    stop(paste(c(x,y,LSD)[!(c(x,y,LSD) %in% names(object))], collapse = ", "), " is/are not in object")
  }

    options <- c("both", "upper", "lower")
  tri.opt <- options[check.arg.values(triangles, options)]
  # if (!(alpha %in% c(0.05, 0.10)))
  #   stop("alpha must be either 0.05 or 0.10")
  plt <- NULL
  
  if (all(is.na(object[[LSD]])))
    warning("All LSD-values are NA for this plot")
  else
  {
    object <- na.omit(object)
    object[x] <- factor(object[[x]])
    object[y] <- factor(object[[y]])
    labs <- sort(levels(object[[x]]))
    #  if (any(labs != sort(levels(object[[y]]))))
    #    stop("The row and column labels of differences are not the same")
    plt <- ggplot(object, aes(x = .data[[!!x]], y = .data[[!!y]], fill=.data[[!!LSD]])) +
      geom_tile() +
      scale_fill_gradientn(colours=colours) + #, 
 #                          values = c(0, 0.001, 0.01, 0.05, 0.10, 1), 
#                           limits=c(0,1)) +
      labs(x=axis.labels, y=axis.labels, title=title) + 
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=axis.text.size),
            axis.text.y=element_text(size=axis.text.size),
            plot.title=element_text(face="bold"),
            panel.grid = element_blank(),
            legend.position = "right", 
            legend.key.height=unit(2,"lines"), legend.key.width=unit(1.5,"lines"), 
            aspect.ratio = 1) +
      guides(fill=guide_colourbar(title = "LSD", nbin=50))
#     if (show.sig)
#     { 
#       object$sig <- ifelse(object[p] > alpha, "",
# #                           ifelse(object[p] > 0.05, ".",
#                                   ifelse(object[p] > 0.01, "*",
#                                          ifelse(object[p] > 0.001, "**",
#                                                 "***")))#)
#       if (alpha == 0.1)
#         if (any(object[p] <= alpha & object[p] > 0.05))
#           object$sig[(object[p] <= alpha & object[p] > 0.05)] <- "."
#       plt <- plt + geom_text(data=object, aes(label=.data[["sig"]]), size=3)
#     }
    
    if (gridspacing[1] > 0)
    {
      if (length(gridspacing) > 1)
      {
        grids <- cumsum(gridspacing)+0.5
      } else
      {
        nlabs <- length(labs)
        grids <- seq(gridspacing + 0.5, nlabs, gridspacing)
      }
      if (tri.opt == "lower")
        plt <- plt + geom_hline(yintercept = grids) + geom_vline(xintercept = grids - 1)
      else
      {
        if (tri.opt == "upper")
          plt <- plt + geom_hline(yintercept = grids - 1) + geom_vline(xintercept = grids)
        else
          plt <- plt + geom_hline(yintercept = grids) + geom_vline(xintercept = grids)
      }
    }
    
    if (!is.null(plt) && !is.null(ggplotFuncs))
      for (f in ggplotFuncs)
        plt <- plt + f
    
    if (!is.null(plt) && printPlot)
      print(plt)
  }
  invisible(plt)
}

"plotLSDs.alldiffs" <- function(object, alpha = 0.05, sections = NULL, 
                                gridspacing = 0, factors.per.grid = 0, 
                                triangles = "both", 
                                title = NULL, axis.labels = TRUE, axis.text.size = 12, 
                                sep=",", colours = RColorBrewer::brewer.pal(3, "Set2"), 
                                ggplotFuncs = NULL, printPlot = TRUE, 
                                sortFactor = NULL, sortParallelToCombo = NULL, 
                                sortNestingFactor = NULL, sortOrder = NULL, decreasing = FALSE, 
                                ...)
  #Plots a matrix of LSD-values or, when  predictions are for combinations of two 
  #factors, produces a plot for each levels combination of  nominated section factors 
  #object is an all.diffs object with sed from which to calculate the LSDs
  #show.sig is a logical indicating whether to put stars onto plot
  #title is a character string giving the plot main title
{ 
  show.sig = FALSE #leave here in case want to indicate significant differences 
                   # - would need to add p.differences to data.frame.
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  if (is.character(validalldifs))
    stop(validalldifs)
  object <- renameDiffsAttr(object)
  
  if (is.null(object$sed))
    stop("The sed component of object cannot be NULL")
  if (all(is.na(object$sed)))
    stop("All sed values are NA")
  options <- c("both", "upper", "lower")
  tri.opt <- options[check.arg.values(triangles, options)]
  
  
  #Sort alldiffs components, if required
  if (!is.null(sortFactor))
    object <- sort(object, decreasing = decreasing, sortFactor = sortFactor, 
         sortParallelToCombo = sortParallelToCombo, sortNestingFactor = sortNestingFactor, 
         sortOrder = sortOrder)
  
  classify <- attr(object, which = "classify")
  
  #Get LSDs
  denom.df <- attr(object, which = "tdf")
  if (is.null(denom.df))
    stop(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
               "- LSDs cannot be calculated"))
  t.value = qt(1-alpha/2, denom.df) 
  LSDs <- t.value * object$sed

  #Convert LSDs to a data.frame
  if (tri.opt == "upper") 
    object[lower.tri(LSDs)] <- NA
  else 
  {
    if (tri.opt == "lower")
      LSDs[upper.tri(LSDs)] <- NA
  }
  LSD.dat <- within(reshape2::melt(LSDs), 
                    { 
                      Var1 <- factor(Var1, labels=dimnames(object$sed)[[1]])
                      Var2 <- factor(Var2, labels=levels(Var1))
                    })
  names(LSD.dat) <- c("Rows","Columns","LSDs")
  
  #prepare for plotting
  n <- nrow(LSD.dat)
  facs <- fac.getinTerm(classify)
  if (any(is.na(match(facs, names(object$predictions)))))
    stop("Some factors in the classify for object are not in the predictions")
  else
    facs[match(facs, names(object$predictions))] <- facs
  nfacs <- length(facs)
  #Function to calculate gridspacing
  autogridspace <- function(object, plotfacs, factors.per.grid = 0)
  {
    gridspace = 0
    if (length(plotfacs) > 1) #only compute gridspace for more than one factor per section
    {
      gridfacs <- plotfacs[order(1:(length(plotfacs)-factors.per.grid), decreasing = TRUE)]
      gridspace <- as.vector(table(object[gridfacs]))
      gridspace <- gridspace[-length(gridspace)]
    }
    return(gridspace)
  }
  
  #Do plots
  plts <- list()
  if (is.null(sections))
  { 
    pairname <- NULL
    if (axis.labels)
    {
      pairname <-  fac.getinTerm(classify)
      pairname <- paste(pairname, collapse = ", ")
      
    }
    #Do single plot
    if (factors.per.grid > 0)
      gridspacing <- autogridspace(object = object$predictions, plotfacs = facs, 
                                   factors.per.grid = factors.per.grid)
    if (is.null(title))
      plts[[1]]  <- plotLSDs.data.frame(object = LSD.dat, x = "Rows", y = "Columns", alpha = alpha, 
                                        gridspacing = gridspacing, #show.sig = show.sig, 
                                        triangles = triangles, 
                                        axis.labels = pairname, colours = colours, 
                                        printPlot = printPlot, 
                                        axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
    else
      plts[[1]] <- plotLSDs.data.frame(object = LSD.dat, x = "Rows", y = "Columns", alpha = alpha, 
                                       gridspacing = gridspacing, #show.sig = show.sig, 
                                       triangles = triangles, 
                                       title = title, axis.labels = pairname, 
                                       colours = colours, printPlot = printPlot, 
                                       axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
  } else #have sections
  { 
    #Prepare for sectioning
    sec.pos <- match(sections, facs)
    pairdiffs <- facs[-match(sections, facs)]
    pd.pos <- match(pairdiffs, facs)
    facOrg <- function(LSD.fac, LSD.dat, object, facs, sections, pairdiffs, sep)
    {
      facs.levs <- strsplit(as.character(LSD.dat[[LSD.fac]]), split=sep, fixed=TRUE)
      facs.levs <- as.data.frame(do.call(rbind, facs.levs), stringsAsFactors = FALSE)
      names(facs.levs) <- facs
      facs.levs <- as.data.frame(lapply(facs, 
                                        function(fac, facs.levs, predictions)
                                        {
                                          levs <- levels(factor(predictions[[fac]]))
                                          new.fac <- factor(facs.levs[[fac]], levels = levs)
                                          return(new.fac)
                                        }, 
                                        facs.levs = facs.levs, 
                                        predictions = object$predictions
      ), stringsAsFactors = FALSE)
      names(facs.levs) <- facs
      LSD.dat$sections <- fac.combine(facs.levs[sections], combine.levels = TRUE)
      LSD.dat[LSD.fac] <- fac.combine(facs.levs[pairdiffs], combine.levels = TRUE)
      return(LSD.dat)
    }
    
    LSD.dat <- facOrg("Rows", LSD.dat, object, facs, sections = sections, pairdiffs = pairdiffs, sep = sep)
    names(LSD.dat)[match(c("sections"), names(LSD.dat))] <- "sections1"
    LSD.dat <- facOrg("Columns", LSD.dat, object, facs, sections = sections, pairdiffs = pairdiffs, sep = sep)
    names(LSD.dat)[match(c("sections"), names(LSD.dat))] <- "sections2"
    
    sect.lev <- levels(LSD.dat$sections1)
    if (!all(sect.lev == levels(LSD.dat$sections2)))
      stop("Sectioning levels in rows and columns of sed do not match")
    
    #Do plot for each section
    sec.name <- pairname <- NULL
    if (axis.labels)
    {
      secname <- paste(sections, collapse = ', ')
      pairname <- paste(pairdiffs, collapse = ', ')
    }
    if (factors.per.grid > 0)
      object$predictions$sections <- fac.combine(object$predictions[sections], 
                                                 combine.levels = TRUE)
    for (j in sect.lev)
    { 
      psect <- LSD.dat[LSD.dat$sections1==j & LSD.dat$sections2==j, ]
      if (factors.per.grid > 0)
      {
        objsect <- object$predictions[object$predictions$sections == j,]
        gridspacing <- autogridspace(object = objsect, plotfacs = pairdiffs, 
                                     factors.per.grid = factors.per.grid)
      }
      if (is.null(title))
        plts[[j]] <- plotLSDs.data.frame(object = psect, x = "Rows", y = "Columns", alpha = alpha, 
                                         gridspacing = gridspacing, #show.sig = show.sig, 
                                         triangles = triangles, 
                                         title = paste("Plot of LSD.dat for ",
                                                       secname," = ",j, sep = ""),
                                         axis.labels = pairname, colours = colours, 
                                         printPlot = printPlot, 
                                         axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
      else
        plts[[j]] <- plotLSDs.data.frame(object = psect, x = "Rows", y = "Columns", alpha = alpha, 
                                         gridspacing = gridspacing, #show.sig = show.sig, 
                                         triangles = triangles, 
                                         title = paste(title," - ", secname," = ", j, sep=""),
                                         axis.labels = pairname, colours = colours, 
                                         printPlot = printPlot, 
                                         axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
    }  
  }
  invisible(list(LSDs = LSD.dat, plots = plts))
}

#'### Functions to plot the reliability of LSD results for a set of pairwise differences
#+
"plotLSDerrors.data.frame" <- function(object, LSDresults = "LSDresults",  x, y, alpha = 0.05,
                                            triangles = "both", gridspacing = 0,  
                                            title = NULL, axis.labels = NULL, axis.text.size = 12, 
                                            colours = c("white","blue","red","grey"), 
                                            ggplotFuncs = NULL, printPlot = TRUE, ...)
  #Plots a data.frame of LSD-results that has columns x, y, 
{ 
  if (!all(c(LSDresults,x,y) %in% names(object)))
  {
    print(names(object))
    stop(paste(c(x,y,LSDresults)[!(c(x,y,LSDresults) %in% names(object))], collapse = ", "), " is/are not in object")
  }
  options <- c("both", "upper", "lower")
  tri.opt <- options[check.arg.values(triangles, options)]
  plt <- NULL
  
  if (all(is.na(object[[LSDresults]])))
    warning("All LSD-results are NA for this plot")
  else
  {
    object <- na.omit(object)
    object[x] <- factor(object[[x]])
    object[y] <- factor(object[[y]])
    cols <- colours[c("Ok","FN","FP", "na") %in% unique(object$LSDresults)]
    object$LSDresults <- factor(object$LSDresults, levels = c("Ok","FN","FP", "na"), 
                                labels = c("Ok", "False\nNegative", "False\nPositive", "n.a."))
    labs <- sort(levels(object[[x]]))
    #  if (any(labs != sort(levels(object[[y]]))))
    #    stop("The row and column labels of differences are not the same")
    plt <- ggplot(object, aes(x = .data[[!!x]], y = .data[[!!y]], fill=.data[[!!LSDresults]])) +
      geom_tile() +
      scale_fill_manual(values=cols) + 
      labs(x=axis.labels, y=axis.labels, title=title) + 
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=axis.text.size),
            axis.text.y=element_text(size=axis.text.size),
            plot.title=element_text(face="bold"),
            panel.grid = element_blank(),
            legend.position = "right", 
            legend.key = element_rect(colour = "black", linewidth = 1),
            legend.key.height=unit(2,"lines"), legend.key.width=unit(1.5,"lines"), 
            aspect.ratio = 1) +
      guides(fill=guide_legend(title = "LSD results"))
    
    if (gridspacing[1] > 0)
    {
      if (length(gridspacing) > 1)
      {
        grids <- cumsum(gridspacing)+0.5
      } else
      {
        nlabs <- length(labs)
        grids <- seq(gridspacing + 0.5, nlabs, gridspacing)
      }
      if (tri.opt == "lower")
        plt <- plt + geom_hline(yintercept = grids) + geom_vline(xintercept = grids - 1)
      else
      {
        if (tri.opt == "upper")
          plt <- plt + geom_hline(yintercept = grids - 1) + geom_vline(xintercept = grids)
        else
          plt <- plt + geom_hline(yintercept = grids) + geom_vline(xintercept = grids)
      }
    }
    
    if (!is.null(plt) && !is.null(ggplotFuncs))
      for (f in ggplotFuncs)
        plt <- plt + f
    
    if (!is.null(plt) && printPlot)
      print(plt)
  }
  invisible(plt)
}

"plotLSDerrors.alldiffs" <- function(object, alpha = 0.05, useIntervals = FALSE, 
                                     sections = NULL, 
                                     gridspacing = 0, factors.per.grid = 0, 
                                     triangles = "both", 
                                     title = NULL, axis.labels = TRUE, axis.text.size = 12, 
                                     sep=",", colours = c("white","blue","red","grey"), 
                                     ggplotFuncs = NULL, printPlot = TRUE, 
                                     sortFactor = NULL, sortParallelToCombo = NULL, 
                                     sortNestingFactor = NULL, sortOrder = NULL, decreasing = FALSE, 
                                     ...)
  #Plots a matrix of LSD-values or, when  predictions are for combinations of two 
  #factors, produces a plot for each levels combination of  nominated section factors 
  #
  #object is an all.diffs object with sed from which to calculate the LSDs
  #title is a character string giving the plot main title
{ 
  #The seperator between the levels of the factors that form the rownames of the LSD component
  #is assumed to be a "," - see fac.LSDcombs.alldiffs.
  #Set sep.levs in case it is decided to allow the user to specify the seperator
  sep.levs = ","
  
  show.sig = FALSE #leave here in case want to indicate significant differences 
  # - would need to add p.differences to data.frame.
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  if (is.character(validalldifs))
    stop(validalldifs)
  object <- renameDiffsAttr(object)
  
  if (is.null(object$sed))
    stop("The sed component of object cannot be NULL")
  if (all(is.na(object$sed)))
    stop("All sed values are NA")
  options <- c("both", "upper", "lower")
  tri.opt <- options[check.arg.values(triangles, options)]
  
  
  #Sort alldiffs components, if required
  if (!is.null(sortFactor))
    object <- sort(object, decreasing = decreasing, sortFactor = sortFactor, 
                   sortParallelToCombo = sortParallelToCombo, sortNestingFactor = sortNestingFactor, 
                   sortOrder = sortOrder)
  
  classify <- attr(object, which = "classify")
  
  #Form LSD results
  denom.df <- attr(object, which = "tdf")
  if (is.null(denom.df))
    stop(paste("The degrees of freedom of the t-distribtion are not available in alldiffs.obj\n",
               "- LSDs cannot be calculated"))
  sig.actual <- object$p.differences <= alpha
  
  if (!useIntervals)
  {  
    LSDapprox <- object$LSD["assignedLSD"]
    #Expand LSDapprox into an LSDresults matrix that is the same size as sed 
    if (rownames(LSDapprox)[1] == "overall")
    {
      LSDresults <- matrix(LSDapprox[1,1], nrow = nrow(object$sed), ncol = ncol(object$sed))
      rownames(LSDresults) <- colnames(LSDresults) <- rownames(object$sed)
    } else
    {
      nr <- nrow(object$sed)
      #Retrieve the LSDby attribute
      LSDby <- attr(object, which = "LSDby")
      if (is.null(LSDby))
        stop("The LSDby attribute of the supplied alldiffs.object is NULL")

      fac.combs.sed <- fac.combine(as.list(object$predictions[LSDby]), combine.levels = TRUE, 
                                   sep =sep)
      LSDresults <- matrix(nrow = nr, ncol = ncol(object$sed))
      
      for (lev in rownames(LSDapprox))
      {
        rownames(LSDresults) <- colnames(LSDresults) <- rownames(object$sed)
        kcells <- fac.combs.sed == lev
        if (!any(kcells))
          stop("No elements of the sed component of object correspond to the LSD for ",lev)
        LSDresults[kcells, kcells] <- as.vector(LSDapprox[rownames(LSDapprox) == lev,1])
      }
    }
    diag(LSDresults) <- NA
    sig.approx <- abs(object$differences) >= LSDresults
  } else
  {
    #Find the intervals
    preds <- object$predictions
    lims <- names(preds)[grep(".limit", names(preds))]
    if (length(lims) > 2)
      stop("There are more than two columns whose names end with .limit")
    lims <- c(lower = lims[grep("lower.", lims)], upper = lims[grep("upper", lims)])
    #Calculate significances using the intervals
    sig.approx <- apply(preds, MARGIN = 1, 
                        FUN = function(irow, preds, lims) 
                        {
                          apply(preds, MARGIN = 1, 
                                FUN = function(jrow, irow, lims)
                                {
                                  reslt <- (irow[lims["lower"]] >= jrow[lims["upper"]] | 
                                              irow[lims["upper"]] <= jrow[lims["lower"]])
                                  return(reslt)
                                }, irow = irow, lims = lims)
                        }, 
                        preds = preds, lims = lims, simplify = FALSE)
    sig.approx <- do.call(rbind, sig.approx)
    diag(sig.approx) <- NA
    rownames(sig.approx) <- colnames(sig.approx) <- rownames(object$p.differences)
    LSDresults <- sig.approx
  }
  
  #Determine the veracity of using the assigned LSDs
  falsepos <- !sig.actual & sig.approx
  falseneg <- sig.actual & !sig.approx
  LSDresults[!is.na(LSDresults)] <- "Ok"
  LSDresults[is.na(LSDresults)] <- "na"
  LSDresults[falseneg] <- "FN"
  LSDresults[falsepos] <-"FP"
  
  #Convert LSDresults to a data.frame
  if (tri.opt == "upper") 
    object[lower.tri(LSDresults)] <- NA
  else 
  {
    if (tri.opt == "lower")
      LSDresults[upper.tri(LSDresults)] <- NA
  }
  LSDres.dat <- within(reshape2::melt(LSDresults), 
                       { 
                         Var1 <- factor(Var1, labels=dimnames(object$sed)[[1]])
                         Var2 <- factor(Var2, labels=levels(Var1))
                       })
  names(LSDres.dat) <- c("Rows","Columns","LSDresults")
  
  #prepare for plotting
  n <- nrow(LSDres.dat)
  facs <- fac.getinTerm(classify)
  if (any(is.na(match(facs, names(object$predictions)))))
    stop("Some factors in the classify for object are not in the predictions")
  else
    facs[match(facs, names(object$predictions))] <- facs
  nfacs <- length(facs)
  #Function to calculate gridspacing
  autogridspace <- function(object, plotfacs, factors.per.grid = 0)
  {
    gridspace = 0
    if (length(plotfacs) > 1) #only compute gridspace for more than one factor per section
    {
      gridfacs <- plotfacs[order(1:(length(plotfacs)-factors.per.grid), decreasing = TRUE)]
      gridspace <- as.vector(table(object[gridfacs]))
      gridspace <- gridspace[-length(gridspace)]
    }
    return(gridspace)
  }
  
  #Do plots
  plts <- list()
  if (is.null(sections))
  { 
    pairname <- NULL
    if (axis.labels)
    {
      pairname <-  fac.getinTerm(classify)
      pairname <- paste(pairname, collapse = ", ")
      
    }
    #Do single plot
    if (factors.per.grid > 0)
      gridspacing <- autogridspace(object = object$predictions, plotfacs = facs, 
                                   factors.per.grid = factors.per.grid)
    if (is.null(title))
      plts[[1]] <- plotLSDerrors.data.frame(object = LSDres.dat, x = "Rows", y = "Columns", alpha = alpha, 
                                            gridspacing = gridspacing, 
                                            triangles = triangles, 
                                            axis.labels = pairname, colours = colours, 
                                            printPlot = printPlot, 
                                            axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
    else
      plts[[1]] <-   plotLSDerrors.data.frame(object = LSDres.dat, x = "Rows", y = "Columns", alpha = alpha, 
                                              gridspacing = gridspacing, 
                                              triangles = triangles, 
                                              title = title, axis.labels = pairname, 
                                              colours = colours, printPlot = printPlot, 
                                              axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
  } else #have sections
  { 
    #Prepare for sectioning
    sec.pos <- match(sections, facs)
    pairdiffs <- facs[-match(sections, facs)]
    pd.pos <- match(pairdiffs, facs)
    facOrg <- function(LSD.fac, LSDres.dat, object, facs, sections, pairdiffs, sep)
    {
      facs.levs <- strsplit(as.character(LSDres.dat[[LSD.fac]]), split=sep, fixed=TRUE)
      facs.levs <- as.data.frame(do.call(rbind, facs.levs), stringsAsFactors = FALSE)
      names(facs.levs) <- facs
      facs.levs <- as.data.frame(lapply(facs, 
                                        function(fac, facs.levs, predictions)
                                        {
                                          levs <- levels(factor(predictions[[fac]]))
                                          new.fac <- factor(facs.levs[[fac]], levels = levs)
                                          return(new.fac)
                                        }, 
                                        facs.levs = facs.levs, 
                                        predictions = object$predictions
      ), stringsAsFactors = FALSE)
      names(facs.levs) <- facs
      LSDres.dat$sections <- fac.combine(facs.levs[sections], combine.levels = TRUE)
      LSDres.dat[LSD.fac] <- fac.combine(facs.levs[pairdiffs], combine.levels = TRUE)
      return(LSDres.dat)
    }
    
    LSDres.dat <- facOrg("Rows", LSDres.dat, object, facs, sections = sections, pairdiffs = pairdiffs, sep = sep)
    names(LSDres.dat)[match(c("sections"), names(LSDres.dat))] <- "sections1"
    LSDres.dat <- facOrg("Columns", LSDres.dat, object, facs, sections = sections, pairdiffs = pairdiffs, sep = sep)
    names(LSDres.dat)[match(c("sections"), names(LSDres.dat))] <- "sections2"
    
    sect.lev <- levels(LSDres.dat$sections1)
    if (!all(sect.lev == levels(LSDres.dat$sections2)))
      stop("Sectioning levels in rows and columns of sed do not match")
    
    #Do plot for each section
    sec.name <- pairname <- NULL
    if (axis.labels)
    {
      secname <- paste(sections, collapse = ', ')
      pairname <- paste(pairdiffs, collapse = ', ')
    }
    if (factors.per.grid > 0)
      object$predictions$sections <- fac.combine(object$predictions[sections], 
                                                 combine.levels = TRUE)
    for (j in sect.lev)
    { 
      psect <- LSDres.dat[LSDres.dat$sections1==j & LSDres.dat$sections2==j, ]
      if (factors.per.grid > 0)
      {
        objsect <- object$predictions[object$predictions[sections] == j,]
        gridspacing <- autogridspace(object = objsect, plotfacs = pairdiffs, 
                                     factors.per.grid = factors.per.grid)
      }
      if (is.null(title))
        plts[[j]] <- plotLSDerrors.data.frame(object = psect, x = "Rows", y = "Columns", alpha = alpha, 
                                              gridspacing = gridspacing, 
                                              triangles = triangles, 
                                              title = paste("Plot of LSD errors for ",
                                                            secname," = ",j, sep = ""),
                                              axis.labels = pairname, colours = colours, 
                                              printPlot = printPlot, 
                                              axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
      else
        plts[[j]] <- plotLSDerrors.data.frame(object = psect, x = "Rows", y = "Columns", alpha = alpha, 
                                              gridspacing = gridspacing, 
                                              triangles = triangles, 
                                              title = paste(title," - ", secname," = ", j, sep=""),
                                              axis.labels = pairname, colours = colours, 
                                              printPlot = printPlot, 
                                              axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
      
    }  
  }
  invisible(list(LSDresults = LSDres.dat, plots = plts))
}

