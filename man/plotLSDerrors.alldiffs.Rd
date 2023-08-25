\name{plotLSDerrors.alldiffs}
\alias{plotLSDerrors.alldiffs}
\alias{plotLSDerrors}
\title{Plots a map of the errors that occur in using the computed LSD values for pairwise differences between predictions.}
\description{Produces a plot of the errors that occur in using the computed  LSD values for 
             pairwise differences predictions by comparing the result obtained from using the 
             LSDs stored in the \code{assignedLSD} column of the \code{LSD} component of 
             the \code{\link{alldiffs.object}} with those computed from the \code{sed} 
             component using the \code{t}-value for the \code{df} stored in the  \code{tdf} 
             attribute of the \code{\link{alldiffs.object}}. 
             The \code{sed} component is generally a matrix whose rows and columns 
             are labelled by the levels of one or more factors, the set of labels being 
             the same for rows and columns. The \code{sections} argument allows multiple 
             plots to be produced, one for each combination of the levels of the factors 
             listed in \code{sections}. Otherwise, a single plot is produced for all 
             observed combinations of the levels of the factors in the \code{classify} 
             attribute for the \code{\link{alldiffs.object}}. The plots are produced using 
             \code{\link{plotLSDerrors.data.frame}}. The order of plotting the levels of 
             one of the factors indexing the predictions can be modified using 
             \code{\link{sort.alldiffs}}.
}
\usage{plotLSDerrors(object, ...)
\method{plotLSDerrors}{alldiffs}(object, alpha = 0.05, useIntervals = FALSE, 
              sections = NULL, gridspacing = 0, factors.per.grid = 0, 
              triangles = "both", title = NULL, 
              axis.labels = TRUE, axis.text.size = 12, 
              sep=",", colours = c("white","blue","red","grey"), 
              ggplotFuncs = NULL, printPlot = TRUE, 
              sortFactor = NULL, sortParallelToCombo = NULL, 
              sortNestingFactor = NULL, sortOrder = NULL, 
              decreasing = FALSE, ...)}
\arguments{
 \item{object}{An \code{\link{alldiffs.object}} with both \code{LSD} and  
              \code{sed} components that are not \code{NULL}.}
 \item{alpha}{A \code{\link{numeric}} giving the significance level for the LSD.}
 \item{useIntervals}{A \code{\link{logical}} indicating whether to use the interval limits 
              stored in the \code{predictions} component of \code{object}, instead of
              the LSDs stored in the \code{LSD} component, for determining whether 
              pairs of predictions are significantly different. It allows a check of how 
              the \code{error.intervals} in the \code{predictions} component will perform 
              if they are used for all pairwise predictions comparisons, whereas the 
              comparisons to which the LSDs apply 
              may be restricted by the setting of the \code{LSDby} attribute of 
              \code{object}. There is no restriction on the \code{error.intervals} that 
              can be used. However, the limits for them must be in columns 
              in the \code{predictions} component of \code{object} and their names 
              must end with \code{.limits} and begin with \code{lower.} and \code{upper.}}
 \item{sections}{A \code{character} listing the names of the factors that are to be used 
                 to break the plot into sections. A separate plot will be produced for 
                 each observed combination of the levels of these factors.}
 \item{gridspacing}{A \code{numeric} specifying the number(s) of rows and columns 
             that form groups in the grid of differences. An alternative is to specify 
             the \code{factors.per.grid} argument to have the grid spacings automatically 
             calculated. Grids are most useful when  
             two or more factors index the rows and columns. If a single, nonzero 
             number, \emph{k} say, is given then a grid line is placed after every 
             \emph{k}th  row and column. If a vector of values is given then the 
             number of grid lines is the length of the vector and the spacing between 
             each is specified by the elements of the vector.}
  \item{factors.per.grid}{A \code{numeric} specifying the number of factors to include 
             within each grid of differences. The \code{gridspacing} will then be 
             computed based on the numbers of combinations observed within the 
             levels of the remaining factors in a single plot. The \code{gridspacing} 
             argument to this function will be ignored if \code{factors.per.grid} is 
             greater than zero. Grids are most useful when two or more factors index the 
             rows and columns of each plot.}
  \item{triangles}{A \code{character} indicating whether the plot should include the 
             \code{lower}, \code{upper} or \code{both} traingle(s).}
\item{title}{A \code{character} string giving the main title for the plot and to which 
              is appended the levels combination of the sectioning factors, if any, for each plot.}
 \item{axis.labels}{A \code{logical} indicating whether a label is to be added to the x- and y-axes. 
              If TRUE, the label is the comma-separated list of factors whose levels 
              combinations are involved in the prediction differences for which the LSD values 
              are calculated.}
 \item{axis.text.size}{A \code{\link{numeric}} giving the size of the \code{labels} on the axes   
                  of the \code{heatmap}.}
 \item{sep}{A \code{character} giving the characters separating the levels of different 
               factors in the row and column names of the \code{sed} component. }
 \item{colours}{A vector of of colours to be passed to the \code{ggplot} function 
                \code{scale\_colour\_gradientn}.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element. 
                It is passed to \code{ggplot} via \code{\link{plotLSDerrors.data.frame}} 
                to be applied in creating the \code{ggplot} object.}
 \item{printPlot}{A \code{logical} indicating whether or not the a plot is to be printed. This would 
                be used when just the returned \code{data.frame} is required.}
\item{sortFactor}{A \code{\link{character}} containing the name of the 
              \code{factor} that indexes the set of predicted values that determines 
               the sorting of the components. If there is only one variable in the 
              \code{classify} term then \code{sortFactor} can be \code{NULL} and 
              the order is defined by the complete set of predicted values. 
              If there is more than one variable in the \code{classify} term 
              then \code{sortFactor} must be set. In this case the \code{sortFactor} 
              is sorted in the same order within each combination of the values of 
              the \code{sortParallelToCombo} variables: the \code{classify} variables, excluding the 
              \code{sortFactor}. There should be only one predicted value for 
              each unique value of \code{sortFactor} within each set defined by a 
              combination of the values of the \code{classify} variables, excluding the 
              \code{sortFactor} \code{factor}. 
              The order to use is determined by either \code{sortParallelToCombo} or 
              \code{sortOrder}.}
 \item{sortParallelToCombo}{A \code{\link{list}} that specifies a combination of the values 
              of the \code{factor}s and \code{numeric}s, excluding \code{sortFactor}, that 
              are in \code{classify}. Each of the components of the supplied \code{\link{list}} 
              is named for a \code{classify} variable and specifies a single value for it. The 
              combination of this set of values will be used to define a subset of the predicted 
              values whose order will define the order of \code{sortFactor}. Each of the other 
              combinations of the values of the \code{factor}s and \code{numeric}s will be sorted 
              in parallel. If \code{sortParallelToCombo} is \code{NULL} then the first value of               
              each \code{classify} variable, except for the \code{sortFactor} \code{factor},  
              in the \code{predictions} component is used to define \code{sortParallelToCombo}. 
              If there is only one variable in the \code{classify} then 
              \code{sortParallelToCombo} is ignored.}
 \item{sortNestingFactor}{A \code{\link{character}} containing the name of the 
              \code{factor} that defines groups of the \code{sortFactor} within which the predicted 
              values are to be ordered. 
              If there is only one variable in the \code{classify} then 
              \code{sortNestingFactor} is ignored.}
 \item{sortOrder}{A \code{character vector} whose length is the same as the number 
               of levels for \code{sortFactor} in the \code{predictions} component of the 
              \code{\link{alldiffs.object}}. It specifies the desired order of the 
              levels in the reordered components of the \code{\link{alldiffs.object}}. 
              The argument \code{sortParallelToCombo} is ignored.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
 \item{decreasing}{A \code{logical} passed to \code{order} that detemines whether 
              the order for sorting the \code{\link{alldiffs.object}} components is for 
               increasing or decreasing magnitude of the predicted values.}
 \item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
 }
\value{A \code{\link{list}} with components named \code{LSDresults} and \code{plots}. 
       The \code{LSDresults} component contains the \code{data.frame} with the columns \code{Rows}, 
       \code{Columns}, \code{LSDresults}, \code{sections1} and \code{sections2}. This \code{data.frame} 
       is formed using the \code{LSD} and \code{sed} components of \code{object} and is used 
       by \code{\link{plotLSDerrors.data.frame}} in producng the plots. The \code{plots} 
       component contains a list of \code{ggplot} objects, one for each plot produced. 
       Multiple plots are stored in the \code{plots} component if the \code{sections} argument 
       is set and the plots are are named for the levels combinations of the sections.}
\author{Chris Brien}
\seealso{\code{\link{plotLSDerrors.alldiffs}}, \code{\link{plotLSDerrors.data.frame}}, 
        \code{\link{plotLSDs.data.frame}},  \cr
        \code{\link{exploreLSDs}}, \code{\link{sort.alldiffs}}, \code{\link{subset.alldiffs}}, 
         \code{ggplot}}
\examples{
##Subset WaterRunoff data to reduce time to execute
data(WaterRunoff.dat)
tmp <- subset(WaterRunoff.dat, Date == "05-18" & Benches != "3")

##Use asreml to get predictions and associated statistics

\dontrun{
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= tmp)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
TS.diffs <- predictPlus.asreml(classify = "Sources:Type", 
                               asreml.obj = current.asr, tables = "none", 
                               wald.tab = current.asrt$wald.tab, 
                               present = c("Type","Species","Sources"))
}

## Use lmeTest and emmmeans to get predictions and associated statistics

if (requireNamespace("lmerTest", quietly = TRUE) & 
    requireNamespace("emmeans", quietly = TRUE))
{
  m1.lmer <- lmerTest::lmer(pH ~ Benches + (Sources * (Type + Species)) + 
                              (1|Benches:MainPlots),
                            data=na.omit(WaterRunoff.dat))
  TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Type)
  TS.preds <- summary(TS.emm)
  den.df <- min(TS.preds$df, na.rm = TRUE)
  ## Modify TS.preds to be compatible with a predictions.frame
  TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                   se = "SE", interval.type = "CI", 
                                   interval.names = c("lower.CL", "upper.CL"))
   
  ## Form an all.diffs object and check its validity
  TS.vcov <- vcov(TS.emm)
  TS.diffs <- allDifferences(predictions = TS.preds, classify = "Sources:Type", 
                             vcov = TS.vcov, tdf = den.df)
  validAlldiffs(TS.diffs)
}  

## Plot LSD values for predictions obtained using asreml or lmerTest
if (exists("TS.diffs"))
{
  plotLSDerrors(TS.diffs, gridspacing = rep(c(3,4), c(4,2)))
\donttest{
  plotLSDerrors(TS.diffs, sections = "Sources", axis.labels = TRUE)
}  
}
}
\keyword{hplot}
\keyword{dplot}
\keyword{asreml}
