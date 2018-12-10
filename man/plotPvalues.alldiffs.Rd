\name{plotPvalues.alldiffs}
\alias{plotPvalues.alldiffs}
\alias{plotPvalues}
\title{Plots a heat map of p-values for pairwise differences between predictions.}
\description{Produces a heat-map plot of the p-values for pairwise differences between 
             predictions that is stored in the \code{p.differences} component of an 
             \code{all.diffs} object. This is generally a matrix whose rows and columns 
             are labelled by the levels of one or more factors, the set of labels being 
             the same for rows and columns. The \code{sections} argument allows multiple 
             plots to be produced, one for each combination of the levels of the factors 
             listed in \code{sections}. Otherwise, a single plot is produced for all 
             observed combinations of the levels of the factors in the \code{classify} 
             for the \code{alldiffs} object. The plots are produced using 
             \code{\link{plotPvalues.data.frame}}. The order of plotting the levels of 
             one of the factors indexing the predictions can be modified and is achieved 
             using \code{\link{sort.alldiffs}}.
}
\usage{plotPvalues(object, ...)
\method{plotPvalues}{alldiffs}(object, sections = NULL, 
            gridspacing = 0, factors.per.grid = 0, 
            show.sig = FALSE, triangles = "both", 
            title = NULL, axis.labels = TRUE, sep=",", 
            colours = RColorBrewer::brewer.pal(3, "Set2"), 
            ggplotFuncs = NULL, sortFactor = NULL, 
            sortWithinVals = NULL, sortOrder = NULL, 
            decreasing = FALSE, ...)}
\arguments{
 \item{object}{An \code{alldiffs} object with a \code{p.differences} component that is not 
             \code{NULL}. with either 3 or 4 columns. Only if there are 4 columns, the
             first should be a factor indexing sections for which separate variogram plots are 
             to be produced. In either case, the other 3 columns should be, in order, 
             (i) a factor indexing the x-direction, (ii) a factor indexing the y-direction, 
             and (iii) the residuals for the observed response.}
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
\item{show.sig}{A \code{logical} that specifies whether asterisks indicating the level 
             of significance are to be added to the plot. If they are then 
             `***' indicates that \eqn{p \leq 0.001}{p <= 0.001}, 
             `**' that \eqn{ 0.001 < p \leq 0.01}{0.001 < p <= 0.01}, 
             `*' that \eqn{ 0.01 < p \leq 0.05}{0.01 < p <= 0.05}
             `.' that \eqn{ 0.05 < p \leq 0.10}{0.05 < p <= 0.10}.}
  \item{triangles}{A \code{character} indicating whether the plot should include the 
             \code{lower}, \code{upper} or \code{both} traingle(s).}
\item{title}{A \code{character} string giving the main title for the plot and to which 
              is appended the levels combination of the sectioning factors, if any, for each plot.}
 \item{axis.labels}{A \code{logical} indicating whether a label is to be added to the x- and y-axes. 
              If TRUE, the label is the comma-separated list of factors whose levels 
              combinations are involved in the prediction differences for which the p-values 
              are calculated.}
 \item{colours}{A vector of of colours to be passed to the \code{ggplot} function 
                \code{scale\_colour\_gradientn}.}
 \item{sep}{A \code{character} giving the characters separating the levels of different 
               factors in the row and column names of the \code{p.differences} component. }
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element. 
                It is passed to \code{ggplot} via \code{\link{plotPvalues.data.frame}}.}
\item{sortFactor}{A \code{character} containing the name of the 
              \code{factor} that indexes the set of predicted values that determines 
              the sorting of the \code{alldiffs} components by 
              \code{\link{sort.alldiffs}}. If \code{NULL} then sorting is not carried 
              out. If there is more than one variable 
              in the \code{classify} term then \code{sortFactor} is sorted for the 
              predicted values within each combination of the values of the sortWithin 
              variables: the \code{classify} variables, excluding the 
              \code{sortFactor}. There should be only one predicted value for 
              each unique value of \code{sortFactor} within each set defined by a 
              combination of the values of the sortWithin variables.}
 \item{sortWithinVals}{A \code{list} with a component named for each \code{factor} and 
              \code{numeric} that is a \code{classify} variable for the predictions, 
              excluding \code{sortFactor}. Each component should contain a single 
              value that is a value of the variable. The combination of this set 
              of values will be used to define a subset of the predicted values whose 
              order will define the order of \code{sortFactor} to be used for all 
              combinations of the \code{sortWithinVals} variables. If 
              \code{sortWithinVals} is \code{NULL} then the first value of each 
              sortWithin variable in \code{predictions} component is used 
              to define \code{sortWithinVals}. If there is only one variable in the 
              classify then \code{sortWithinVals} is ignored.}
 \item{sortOrder}{A \code{character vector} whose length is the same as the number 
               of levels for \code{sortFactor} in the \code{predictions} component of the 
              \code{alldiffs} object. It specifies the desired order of the 
              levels in the reordered components of the \code{alldiffs} object. 
              The argument \code{sortWithinVals} is ignored.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
 \item{decreasing}{A \code{logical} passed to \code{order} that detemines whether 
              the order for sorting the \code{alldiffs} components is for 
               increasing or decreasing magnitude of the predicted values.}
 \item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
 }
\value{A \code{data.frame} with the columns \code{X1}, \code{X2}, \code{p}, \code{sections1} 
       and \code{sections2}. This {data.frame} is formed from the \code{p.differences} component 
       of \code{object} and is used in producng the plot.}
\author{Chris Brien}
\seealso{\code{\link{plotPvalues.data.frame}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{sort.alldiffs}}, \code{\link{subset.alldiffs}}, \code{ggplot}}
\examples{\dontrun{
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= WaterRunoff.dat))
current.asrt <- asrtests(current.asr, NULL, NULL)
diffs <- predictPlus.asreml(classify = "Sources:Type", 
                            asreml.obj = current.asr, tables = "none", 
                            wald.tab = current.asrt$wald.tab, 
                            present = c("Type","Species","Sources"))
  
plotPvalues(diffs, gridspacing = rep(c(3,4), c(4,2)), show.sig = TRUE)
plotPvalues(diffs, sections = "Sources", show.sig = TRUE, axis.labels = TRUE)
}}
\keyword{hplot}
\keyword{dplot}
\keyword{asreml}
