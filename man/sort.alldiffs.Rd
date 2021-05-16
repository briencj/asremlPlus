\name{sort.alldiffs}
\alias{sort.alldiffs}
\title{Sorts the components in an \code{\link{alldiffs.object}} according to the predicted values 
       associated with a factor.}
\description{Sorts the rows of the components in an \code{\link{alldiffs.object}} 
             (see \code{\link{as.alldiffs}}) that are \code{data.frames} and the rows and columns 
             of those that are \code{matrices} according to the predicted values in the 
             \code{predictions} component. These predicted values are generally obtained using 
             \code{predict.asreml} by specifying a \code{classify} term comprised of one or 
             more variables. Generally, the values associated with one variable are sorted in 
             parallel within each combination of values of the other variables. When there is more 
             than one variable in the \code{classify} term, the sorting is controlled using 
             one or more of \code{sortFactor}, \code{sortWithinVals} and \code{sortOrder}. 
             If there is only one variable in the \code{classify} then all components are sorted 
             according to the order of the complete set of predictions.

             Note that \code{\link{renewClassify.alldiffs}} is called after sorting to ensure that 
             the order of the rows and columns of the components is in standard order for the new 
             variable order.}
\usage{\method{sort}{alldiffs}(x, decreasing = FALSE, classify = NULL, 
     sortFactor = NULL, sortWithinVals = NULL, sortOrder = NULL, ...)}
\arguments{
 \item{x}{An \code{\link{alldiffs.object}}.}
 \item{decreasing}{A \code{\link{logical}} passed to \code{order} that detemines whether 
              the order is for increasing or decreasing magnitude of the predicted 
              values.}
\item{classify}{A \code{\link{character}} string giving the variables that 
                 define the margins of the multiway table that was predicted. 
                 Multiway tables are specified by forming an interaction type 
                 term from the classifying variables, that is, separating the 
                 variable names with the \code{:} operator. If \code{NULL}, 
                 it will be obtained from the \code{classify} attribute of the 
                 \code{\link{as.alldiffs}} object supplied through \code{x}.}
\item{sortFactor}{A \code{\link{character}} containing the name of the 
              \code{factor} that indexes the set of predicted values that determines 
               the sorting of the components. If there is only one variable in the 
              \code{classify} term then \code{sortFactor} can be \code{NULL} and 
              the order is defined by the complete set of predicted values. 
              If there is more than one variable in the \code{classify} term 
              then \code{sortFactor} must be set. In this case the \code{sortFactor} 
              is sorted in the same order within each combination of the values of 
              the sortWithin variables: the \code{classify} variables, excluding the 
              \code{sortFactor}. There should be only one predicted value for 
              each unique value of \code{sortFactor} within each set defined by a 
              combination of the values of the sortWithin variables. 
              The order to use is determined by either \code{sortWithinVals} or 
              \code{sortOrder}.}
 \item{sortWithinVals}{A \code{\link{list}} with a component named for each \code{factor} and 
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
 \item{sortOrder}{A \code{\link{character}} vector whose length is the same as the number 
               of levels for \code{sortFactor} in the \code{predictions} component of the 
              \code{\link{alldiffs.object}}. It specifies the desired order of the 
              levels in the reordered components of the \code{\link{alldiffs.object}}. 
              The argument \code{sortWithinVals} is ignored.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
 \item{\dots}{further arguments passed to or from other methods. Not used at present.}
}
\value{The \code{\link{alldiffs.object}} supplied with the following components, 
       if present, sorted: \code{predictions}, \code{vcov}, \code{backtransforms}, \code{differences}, 
       \code{p.differences} and \code{sed}. Also, the \code{sortFactor} and \code{sortOrder} 
       attributes are set.}
\details{The basic technique is to change the order of the levels of the \code{sortFactor} 
         within the \code{predictions} and, if present, \code{backtransforms} components so 
         that they are ordered for a subset of predicted values, one for each levels of the 
         \code{sortFactor}. When the \code{classify} term consists of more than one 
         variable then a subset of one combination of the values of variables other than
         the \code{sortFactor}, the sortWithin set, must be chosen for determining the 
         order of the \code{sortFactor} levels. Then the sorting of the rows (and columns) 
         will be in parallel within each combination of the values of sortWithin variables: 
         the \code{classify} term, excluding the \code{sortFactor}.}


\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{print.alldiffs}}, \cr 
         \code{\link{renewClassify.alldiffs}}, \code{\link{redoErrorIntervals.alldiffs}}, 
         \code{\link{recalcLSD.alldiffs}},  \cr
         \code{\link{predictPlus.asreml}}, \code{\link{predictPresent.asreml}}}
\examples{
##Halve WaterRunoff data to reduce time to execute
data(WaterRunoff.dat)
tmp <- subset(WaterRunoff.dat, Date == "05-18")

##Use asreml to get predictions and associated statistics

\dontrun{
#Analyse pH  
m1.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                 random = ~ Benches:MainPlots,
                 keep.order=TRUE, data= tmp)
current.asrt <- as.asrtests(m1.asr, NULL, NULL)
current.asrt <- as.asrtests(m1.asr)
current.asrt <- rmboundary(current.asrt)
m1.asr <- current.asrt$asreml.obj

#Get predictions and associated statistics  
TS.diffs <- predictPlus.asreml(classify = "Sources:Type", 
                               asreml.obj = m1.asr, tables = "none", 
                               wald.tab = current.asrt$wald.tab, 
                               present = c("Type","Species","Sources"))
  
#Use sort.alldiffs and save order for use with other response variables
TS.diffs.sort <- sort(diffs, sortFactor = "Sources", sortWithinVals = list(Type = "Control"))
sort.order <- attr(TS.diffs.sort, which = "sortOrder")
  
#Analyse Turbidity
m2.asr <- asreml(fixed = Turbidity ~ Benches + (Sources * (Type + Species)), 
                 random = ~ Benches:MainPlots,
                 keep.order=TRUE, data= tmp)
current.asrt <- as.asrtests(m2.asr)
#Use pH sort.order to sort Turbidity alldiffs object
diffs2.sort <- predictPlus(m2.asr, classify = "Sources:Type", 
                           pairwise = FALSE, error.intervals = "Stand", 
                           tables = "none", present = c("Type","Species","Sources"),
                           sortFactor = "Sources", 
                           sortOrder = sort.order)
}

## Use lmeTest and emmmeans to get predictions and associated statistics

if (requireNamespace("lmerTest", quietly = TRUE) & 
    requireNamespace("emmeans", quietly = TRUE))
{
  #Analyse pH
  m1.lmer <- lmerTest::lmer(pH ~ Benches + (Sources * (Type + Species)) + 
                              (1|Benches:MainPlots),
                            data=na.omit(tmp))
  TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Type)
  TS.preds <- summary(TS.emm)
  den.df <- min(TS.preds$df, na.rm = TRUE)
  ## Modify TS.preds to be compatible with a predictions.frame
  TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                   se = "SE", interval.type = "CI", 
                                   interval.names = c("lower.CL", "upper.CL"))
  
  ## Form an all.diffs object and check its validity
  TS.vcov <- vcov(TS.emm)
  TS.diffs <- allDifferences(predictions = TS.preds, 
                             classify = "Sources:Type", 
                             vcov = TS.vcov, tdf = den.df)
  validAlldiffs(TS.diffs)
    
  #Use sort.alldiffs and save order for use with other response variables
  TS.diffs.sort <- sort(TS.diffs, sortFactor = "Sources", 
                        sortWithinVals = list(Type = "Control"))
  sort.order <- attr(TS.diffs.sort, which = "sortOrder")
  
  #Analyse Turbidity
  m2.lmer <- lmerTest::lmer(Turbidity ~ Benches + (Sources * (Type + Species)) + 
                              (1|Benches:MainPlots),
                            data=na.omit(tmp))
  TS.emm <- emmeans::emmeans(m2.lmer, specs = ~ Sources:Type)
  TS.preds <- summary(TS.emm)
  den.df <- min(TS.preds$df, na.rm = TRUE)
  ## Modify TS.preds to be compatible with a predictions.frame
  TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                   se = "SE", interval.type = "CI", 
                                   interval.names = c("lower.CL", "upper.CL"))
    
  ## Form an all.diffs object, sorting it using the pH sort.order and check its validity
  TS.vcov <- vcov(TS.emm)
  TS.diffs2.sort <- allDifferences(predictions = TS.preds, 
                                   classify = "Sources:Type", 
                                   vcov = TS.vcov, tdf = den.df,
                                   sortFactor = "Sources", 
                                   sortOrder = sort.order)
  validAlldiffs(TS.diffs2.sort)
}  
}
\keyword{asreml}