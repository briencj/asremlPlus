#'### Functions to plot p-values for a set of pairwise differences
#+
"plotPvalues.data.frame" <- function(object, p = "p",  x, y, 
                                     gridspacing = 0, show.sig = FALSE, 
                                     triangles = "both", 
                                     title = NULL, axis.labels = NULL, 
                                     colours = RColorBrewer::brewer.pal(3, "Set2"), 
                                     ggplotFuncs = NULL, printPlot = TRUE, ...)
  #Plots a data.frame of p-values that has columns x, y, p
{ 
  if (!all(c(p,x,y) %in% names(object)))
    stop("One or more of the columns p, x and y are not in object")
  options <- c("both", "upper", "lower")
  tri.opt <- options[check.arg.values(triangles, options)]
  plt <- NULL
  
  if (all(is.na(object[[p]])))
    warning("All p-values are NA for this plot")
  else
  {
    object <- na.omit(object)
    object[x] <- factor(object[[x]])
    object[y] <- factor(object[[y]])
    labs <- sort(levels(object[[x]]))
    #  if (any(labs != sort(levels(object[[y]]))))
    #    stop("The row and column labels of differences are not the same")
    plt <- ggplot(object, aes_string(x = x, y = y, fill=p)) +
      geom_tile() +
      scale_fill_gradientn(colours=colours, 
                           values = c(0, 0.001, 0.01, 0.05, 0.10, 1), 
                           limits=c(0,1)) +
      labs(x=axis.labels, y=axis.labels, title=title) + 
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=12),
            axis.text.y=element_text(size=12),
            plot.title=element_text(face="bold"),
            panel.grid = element_blank(),
            legend.position = "right", 
            legend.key.height=unit(2,"lines"), legend.key.width=unit(1.5,"lines"), 
            aspect.ratio = 1) +
      guides(fill=guide_colourbar(title = "p", nbin=50))
    if (show.sig)
    { 
      object$sig <- ifelse(object[p] > 0.1, "",
                           ifelse(object[p] > 0.05, ".",
                                  ifelse(object[p] > 0.01, "*",
                                         ifelse(object[p] > 0.001, "**",
                                                "***"))))
      plt <- plt + geom_text(data=object, aes_string(label="sig"), size=3)
    }
    
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

"plotPvalues.alldiffs" <- function(object, sections = NULL, 
                                   gridspacing = 0, factors.per.grid = 0, 
                                   show.sig = FALSE, 
                                   triangles = "both", 
                                   title = NULL, axis.labels = TRUE, sep=",", 
                                   colours = RColorBrewer::brewer.pal(3, "Set2"), 
                                   ggplotFuncs = NULL, printPlot = TRUE, 
                                   sortFactor = NULL, sortWithinVals = NULL, 
                                   sortOrder = NULL, decreasing = FALSE, 
                                   ...)
  #Plots a matrix of p-values or, when  predictions are for combinations of two 
  #factors, produces a plot for each levels combination of  nominated section factors 
  #object is an all.diffs object with p.differences to be plotted
  #show.sig is a logical indicating whether to put stars onto plot
  #title is a character string giving the plot main title
{ 
  #Check that a valid object of class alldiffs
  validalldifs <- validAlldiffs(object)  
  if (is.character(validalldifs))
    stop(validalldifs)
  if (is.null(object$p.differences))
    stop("The p.differences component of object cannot be NULL")
  if (all(is.na(object$p.differences)))
    stop("All p.differences are NA")
  options <- c("both", "upper", "lower")
  tri.opt <- options[check.arg.values(triangles, options)]
  
  
  #Sort alldiffs components, if required
  if (!is.null(sortFactor))
    object <- sort(object, decreasing = decreasing, sortFactor = sortFactor, 
         sortWithinVals = sortWithinVals, sortOrder = sortOrder)
  
  classify <- attr(object, which = "classify")
  #Get differences and convert to a data.frame
  if (tri.opt == "upper") 
    object$p.differences[lower.tri(object$p.differences)] <- NA
  else 
  {
    if (tri.opt == "lower")
      object$p.differences[upper.tri(object$p.differences)] <- NA
  }
  p <- within(reshape::melt(object$p.differences), 
              { 
                X1 <- factor(X1, levels=dimnames(object$p.differences)[[1]])
                X2 <- factor(X2, levels=levels(X1))
              })
  names(p)[match("value", names(p))] <- "p"
#  names(p)[c(1,2)] <- paste(classify, c(".1",".2"), sep = "")
  
  #prepare for plotting
  n <- nrow(p)
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
      plotPvalues.data.frame(object = p, x = "X1", y = "X2", 
                             gridspacing = gridspacing, show.sig = show.sig, 
                             triangles = triangles, 
                             axis.labels = pairname, colours = colours, 
                             printPlot = printPlot, ggplotFuncs = ggplotFuncs)
    else
      plotPvalues.data.frame(object = p, x = "X1", y = "X2",  
                             gridspacing = gridspacing, show.sig = show.sig, 
                             triangles = triangles, 
                             title = title, axis.labels = pairname, 
                             colours = colours, printPlot = printPlot, 
                             ggplotFuncs = ggplotFuncs)
    
  } else #have sections
  { 
    #Prepare for sectioning
    sec.pos <- match(sections, facs)
    pairdiffs <- facs[-match(sections, facs)]
    pd.pos <- match(pairdiffs, facs)
    facOrg <- function(p.fac, p, object, facs, sections, pairdiffs, sep)
    {
      facs.levs <- strsplit(as.character(p[[p.fac]]), split=sep, fixed=TRUE)
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
                                        ))
      names(facs.levs) <- facs
      p$sections <- fac.combine(facs.levs[sections], combine.levels = TRUE)
      p[p.fac] <- fac.combine(facs.levs[pairdiffs], combine.levels = TRUE)
      return(p)
    }
    
    p <- facOrg("X1", p, object, facs, sections = sections, pairdiffs = pairdiffs, sep = sep)
    names(p)[match(c("sections"), names(p))] <- "sections1"
    p <- facOrg("X2", p, object, facs, sections = sections, pairdiffs = pairdiffs, sep = sep)
    names(p)[match(c("sections"), names(p))] <- "sections2"

    sect.lev <- levels(p$sections1)
    if (!all(sect.lev == levels(p$sections2)))
      stop("Sectioning levels in rows and columns of p.differences do not match")
    
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
      psect <- p[p$sections1==j & p$sections2==j, ]
      if (factors.per.grid > 0)
      {
        objsect <- object$predictions[object$predictions$sections == j,]
        gridspacing <- autogridspace(object = objsect, plotfacs = pairdiffs, 
                                     factors.per.grid = factors.per.grid)
      }
      if (is.null(title))
        plotPvalues.data.frame(object = psect, x = "X1", y = "X2", 
                               gridspacing = gridspacing, show.sig = show.sig, 
                               triangles = triangles, 
                               title = paste("Plot of p-values for ",
                                             secname," = ",j, sep = ""),
                               axis.labels = pairname, colours = colours, 
                               printPlot = printPlot, ggplotFuncs = ggplotFuncs)
      else
        plotPvalues.data.frame(object = psect, x = "X1", y = "X2", 
                               gridspacing = gridspacing, show.sig = show.sig, 
                               triangles = triangles, 
                               title = paste(title," - ", secname," = ", j, sep=""),
                               axis.labels = pairname, colours = colours, 
                               printPlot = printPlot, ggplotFuncs = ggplotFuncs)
      
    }  
  }
  invisible(p)
}
