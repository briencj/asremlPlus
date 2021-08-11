#'### Functions to plot LSD-values for a set of pairwise differences
#+
"plotLSDs.data.frame" <- function(object, LSD = "LSD",  x, y, alpha = 0.05,
                                  triangles = "both", gridspacing = 0,  
                                  title = NULL, axis.labels = NULL, axis.text.size = 12, 
                                  colours = RColorBrewer::brewer.pal(3, "Set2"), 
                                  ggplotFuncs = NULL, printPlot = TRUE, ...)
  #Plots a data.frame of LSD-values that has columns x, y, LSD
{ 
  show.sig = FALSE #leave here in case wnat to indicate significant differences 
                   # - would need to add p.differences to data.frame.
  if (!all(c(LSD,x,y) %in% names(object)))
    stop("One or more of the columns LSD, x and y are not in object")
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
    plt <- ggplot(object, aes_string(x = x, y = y, fill=LSD)) +
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
#       plt <- plt + geom_text(data=object, aes_string(label="sig"), size=3)
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
  show.sig = FALSE #leave here in case wnat to indicate significant differences 
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
  rownames(LSDs) <- colnames(LSDs) <- NULL #needed because reshape::melt throws a warning re type.convert
  LSD.dat <- within(reshape::melt(LSDs), 
                    { 
                      X1 <- factor(X1, labels=dimnames(object$sed)[[1]])
                      X2 <- factor(X2, labels=levels(X1))
                    })
  names(LSD.dat)[match("value", names(LSD.dat))] <- "LSD"
  
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
      plotLSDs.data.frame(object = LSD.dat, x = "X1", y = "X2", alpha = alpha, 
                          gridspacing = gridspacing, #show.sig = show.sig, 
                          triangles = triangles, 
                          axis.labels = pairname, colours = colours, 
                          printPlot = printPlot, 
                          axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
    else
      plotLSDs.data.frame(object = LSD.dat, x = "X1", y = "X2", alpha = alpha, 
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
    
    LSD.dat <- facOrg("X1", LSD.dat, object, facs, sections = sections, pairdiffs = pairdiffs, sep = sep)
    names(LSD.dat)[match(c("sections"), names(LSD.dat))] <- "sections1"
    LSD.dat <- facOrg("X2", LSD.dat, object, facs, sections = sections, pairdiffs = pairdiffs, sep = sep)
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
        plotLSDs.data.frame(object = psect, x = "X1", y = "X2", alpha = alpha, 
                            gridspacing = gridspacing, #show.sig = show.sig, 
                            triangles = triangles, 
                            title = paste("Plot of LSD.dat for ",
                                          secname," = ",j, sep = ""),
                            axis.labels = pairname, colours = colours, 
                            printPlot = printPlot, 
                            axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
      else
        plotLSDs.data.frame(object = psect, x = "X1", y = "X2", alpha = alpha, 
                            gridspacing = gridspacing, #show.sig = show.sig, 
                            triangles = triangles, 
                            title = paste(title," - ", secname," = ", j, sep=""),
                            axis.labels = pairname, colours = colours, 
                            printPlot = printPlot, 
                            axis.text.size = axis.text.size, ggplotFuncs = ggplotFuncs)
      
    }  
  }
  invisible(LSD.dat)
}
