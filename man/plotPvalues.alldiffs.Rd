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
             for the \code{\link{alldiffs.object}}. The plots are produced using 
             \code{\link{plotPvalues.data.frame}}. The order of plotting the levels of 
             one of the factors indexing the predictions can be modified using 
             \code{\link{sort.alldiffs}}.
}
\usage{plotPvalues(object, ...)
\method{plotPvalues}{alldiffs}(object, sections = NULL, 
            gridspacing = 0, factors.per.grid = 0, 
            show.sig = FALSE, alpha = 0.10, 
            sig.size = 3, sig.colour = "black", 
            sig.face = "plain", sig.family = "",  
            triangles = "both", 
            title = NULL, axis.labels = TRUE, axis.text.size = 12, 
            sep=",", colours = RColorBrewer::brewer.pal(3, "Set2"), 
            ggplotFuncs = NULL, printPlot = TRUE, 
            sortFactor = NULL, sortParallelToCombo = NULL, 
            sortNestingFactor = NULL, sortOrder = NULL, 
            decreasing = FALSE, ...)}
\arguments{
 \item{object}{An \code{\link{alldiffs.object}} with a \code{p.differences} component that is not 
             \code{NULL}.}
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
             `.' that \eqn{ 0.05 < p \leq 0.10}{0.05 < p <= 0.10}. The last is only 
             included for \code{alpha} = 0.10.}
 \item{alpha}{A \code{\link{numeric}} giving the significance level for testing 
             pairwise differences; must be 0.05 or 0.10.}
 \item{sig.size}{A \code{\link{numeric}} specifying the size, in pts, of the 
             significance asterisks.}
 \item{sig.colour}{A \code{\link{character}} specifying the colour to use for the 
             significance asterisks.}
 \item{sig.face}{A \code{\link{character}} specifying the font face for the significance 
             asterisks (\code{"plain"}, \code{"italic"}, \code{"bold"}, \code{"bold.italic"}).}
\item{sig.family}{A \code{\link{character}} specifying the font family for the significance 
             asterisks. The font families that are available depends on the system. For font 
             families other than the basic Postscript fonts, see the \code{extrafont} package.}
  \item{triangles}{A \code{\link{character}} indicating whether the plot should include the 
             \code{lower}, \code{upper} or \code{both} traingle(s).}
\item{title}{A \code{\link{character}}  string giving the main title for the plot and to which 
              is appended the levels combination of the sectioning factors, if any, for each plot.}
 \item{axis.labels}{A \code{\link{logical}} indicating whether a label is to be added to the x- and y-axes. 
              If TRUE, the label is the comma-separated list of factors whose levels 
              combinations are involved in the prediction differences for which the p-values 
              are calculated.}
 \item{axis.text.size}{A \code{\link{numeric}} giving the size of the \code{labels} on the axes   
                  of the \code{heatmap}.}
 \item{sep}{A \code{\link{character}} giving the characters separating the levels of different 
               factors in the row and column names of the \code{p.differences} component. }
 \item{colours}{A vector of of colours to be passed to the \code{ggplot} function 
                \code{scale\_colour\_gradientn}.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{ggplot2} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{ggplot2} function call for each element. 
                It is passed to \code{ggplot} via \code{\link{plotPvalues.data.frame}} 
                to be applied in creating the \code{ggplot} object.}
 \item{printPlot}{A \code{\link{logical}} indicating whether or not the a plot is to be printed. This would 
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
\value{A \code{\link{list}} with components named \code{pvalues} and \code{plots}. 
       The \code{pvalues} component contains the \code{data.frame} with the columns \code{Rows}, 
       \code{Columns}, \code{p}, \code{sections1} and \code{sections2}. This \code{data.frame} is 
       formed using the \code{sed} component of \code{object} and is used by 
       \code{\link{plotPvalues.data.frame}} in producng the plot. The \code{plots} 
       component contains a list of \code{ggplot} objects, one for each plot produced. 
       Multiple plots are stored in the \code{plots} component if the \code{sections} argument 
       is set and the plots are are named for the levels combinations of the sections.}
\author{Chris Brien}
\seealso{\code{\link{plotPvalues.data.frame}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{sort.alldiffs}}, \code{\link{subset.alldiffs}}, \code{ggplot}}
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

## Plot p-values for predictions obtained using asreml or lmerTest
if (exists("TS.diffs"))
{
  plotPvalues(TS.diffs, gridspacing = rep(c(3,4), c(4,2)), show.sig = TRUE)
\donttest{
  plotPvalues(TS.diffs, sections = "Sources", show.sig = TRUE, axis.labels = TRUE)
}  
}
}
\keyword{hplot}
\keyword{dplot}
\keyword{asreml}
