\name{plotLSDerrors.data.frame}
\alias{plotLSDerrors.data.frame}
\title{Plots a map of the supplied errors that occur in using the computed LSD values for pairwise differences between predictions.}
\description{Produces a plot of the errors that have been supplied in a \code{data.frame}. 
             The \code{data.frame} includes two factors whose levels specify, 
             for each LSD result, which combinations of factor levels are being compared. The function
             \code{\link{plotLSDerrors.alldiffs}} produces such \code{data.frame}s.}
\usage{\method{plotLSDerrors}{data.frame}(object, LSDresults = "LSDresults", x, y, 
              alpha = 0.05, triangles = "both", 
              gridspacing = 0, title = NULL, 
              axis.labels = NULL, axis.text.size = 12, 
              colours = c("white","blue","red","grey"), 
              ggplotFuncs = NULL, printPlot = TRUE, ...)}
\arguments{
 \item{object}{A \code{\link{data.frame}} containing the three columns specified by 
             \code{LSDresults}, \code{x} and \code{y}.}
 \item{LSDresults}{A \code{\link{character}} giving the name of the column in \code{object} 
             that contains the LSDresults values to be plotted. The column should be a 
             \code{\link{character}} or \code{\link{factor}} with values or levels that are a 
             subset of \code{Ok},\code{FN}, \code{FP} and \code{na}.}
 \item{x}{A \code{\link{character}} giving the name of the column in \code{object} that contains 
             the factor whose levels index the LSD values that are to be plotted in the same column.}
 \item{y}{A \code{character} giving the name of the column in \code{object} that contains 
             the labels of the LSD values that are to be plotted as the rows.}
 \item{alpha}{A \code{\link{numeric}} giving the significance level for the LSD.}
 \item{triangles}{A \code{character} indicating whether the plot should include the 
             \code{lower}, \code{upper} or \code{both} triangle(s). Here it is only used 
              to adjust gridlines for the omission of the diagonal.}
 \item{gridspacing}{A \code{numeric} specifying the number(s) of rows and columns 
             that form groups in the grid of differences. This is most useful when  
             two or more factors index the rows and columns. If a single, nonzero 
             number, \emph{k} say, is given then a grid line is placed after every 
             \emph{k}th  row and column. If a vector of values is given then the 
             number of grid lines is the length of the vector and the spacing between 
             each is specified by the elements of the vector.}
 \item{title}{A \code{character} string giving the main title for the plot.}
 \item{axis.labels}{A \code{character} string giving the label to use for both the x- and y-axis.}
 \item{axis.text.size}{A \code{\link{numeric}} giving the size of the \code{labels} on the axes   
                  of the \code{heatmap}.}
 \item{colours}{A vector of colours to be passed to the \code{ggplot} function 
                \code{scale\_colour\_gradientn}.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{ggplot2} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{ggplot2} function call for each element. 
                These functions are applied in creating the \code{ggplot} object.}
 \item{printPlot}{A \code{logical} indicating whether or not the a plot is to be printed. This would 
                be used when just the returned \code{ggplot} object is required.}
 \item{\dots}{Provision for passing arguments to functions called internally - 
            not used at present.}
 }
\value{An object of class "\code{ggplot}", which can be plotted using \code{print} or otherwise manipulated.}
\author{Chris Brien}
\seealso{\code{\link{plotLSDs.data.frame}}, \code{\link{plotLSDs.alldiffs}}, 
         \code{\link{exploreLSDs}}, \code{ggplot}}
\examples{
##Subset WaterRunoff data to reduce time to execute
data(WaterRunoff.dat)
tmp <- subset(WaterRunoff.dat, Date == "05-18")

##Use asreml to get predictions and associated statistics

\dontrun{
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= tmp))
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
  LSDresults <- within(reshape2::melt(TS.diffs$p.differences), 
                      { 
                        Var1 <- factor(Var1, levels=dimnames(TS.diffs$p.differences)[[1]])
                        Var2 <- factor(Var2, levels=levels(Var1))
                      })
  names(LSDresults) <- c("Rows","Columns","LSDresults")
  plotLSDerrors(LSDresults, x = "Rows", y = "Columns", gridspacing = rep(c(3,4), c(4,2)))
}
}
\keyword{hplot}
\keyword{dplot}
\keyword{asreml}
