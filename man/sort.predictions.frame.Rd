\name{sort.predictions.frame}
\alias{sort.predictions.frame}
\title{Sorts a \code{\link{predictions.frame}} according to the predicted values 
       associated with a factor.}
\description{Sorts the rows of a \code{\link{predictions.frame}} according to the predicted values 
             in the \code{predictions.frame}. These predicted values are generally obtained using 
             \code{predict.asreml} by specifying a \code{classify} term comprised of one or 
             more variables. Generally, the values associated with one variable are sorted in 
             parallel within each combination of values of the other variables. When there is more 
             than one variable in the \code{classify} term, the sorting is controlled using 
             one or more of \code{sortFactor}, \code{sortParallelToCombo} and \code{sortOrder}. 
             If there is only one variable in the \code{classify} then the \code{\link{predictions.frame}} 
             is sorted according to the order of the complete set of predictions.}
\usage{\method{sort}{predictions.frame}(x, decreasing = FALSE, classify, sortFactor = NULL, 
     sortParallelToCombo = NULL, sortNestingFactor = NULL, 
     sortOrder = NULL, ...)}
\arguments{
 \item{x}{A \code{\link{predictions.frame}}.}
 \item{decreasing}{A \code{\link{logical}} passed to \code{order} that detemines whether 
              the order is for increasing or decreasing magnitude of the predicted 
              values.}
\item{classify}{A \code{\link{character}} string giving the variables that 
                 define the margins of the multiway table that was predicted. 
                 Multiway tables are specified by forming an interaction type 
                 term from the classifying variables, that is, separating the 
                 variable names with the \code{:} operator.}
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
 \item{sortOrder}{A \code{\link{character}} vector whose length is the same as the number 
               of levels for \code{sortFactor} in the \code{predictions.frame}. 
              It specifies the desired order of the levels in the reordered 
              the \code{\link{predictions.frame}}. 
              The argument \code{sortParallelToCombo} is ignored.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
 \item{\dots}{further arguments passed to or from other methods. Not used at present.}
}
\value{The sorted \code{\link{predictions.frame}}. Also, the \code{sortFactor} and 
             \code{sortOrder} attributes are set.}
\details{The basic technique is to change the order of the levels of the \code{sortFactor} 
         within the \code{predictions.frame} so 
         that they are ordered for a subset of predicted values, one for each levels of the 
         \code{sortFactor}. When the \code{classify} term consists of more than one 
         variable then a subset of one combination of the values of variables other than
         the \code{sortFactor}, the \code{sortParallelToCombo} combination, must be chosen for determining the 
         order of the \code{sortFactor} levels. Then the sorting of the rows (and columns) 
         will be in parallel within each combination of the values of \code{sortParallelToCombo} variables: 
         the \code{classify} term, excluding the \code{sortFactor}.}


\author{Chris Brien}
\seealso{\code{\link{as.predictions.frame}}, \code{\link{print.predictions.frame}}, 
         \code{\link{sort.alldiffs}}, \cr 
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
  
#Use sort.predictions.frame and save order for use with other response variables
TS.preds <- TS.diffs$predictions
TS.preds.sort <- sort(TS.preds, sortFactor = "Sources", 
                      sortParallelToCombo = list(Type = "Control"))
sort.order <- attr(TS.preds.sort, which = "sortOrder")
  
#Analyse Turbidity
m2.asr <- asreml(fixed = Turbidity ~ Benches + (Sources * (Type + Species)), 
                 random = ~ Benches:MainPlots,
                 keep.order=TRUE, data= tmp)
current.asrt <- as.asrtests(m2.asr)
#Use pH sort.order to sort Turbidity alldiffs object
TS.diffs2 <- predictPlus(m2.asr, classify = "Sources:Type", 
                         pairwise = FALSE, error.intervals = "Stand", 
                         tables = "none", present = c("Type","Species","Sources"))
TS.preds2 <- TS.diffs2$predictions
TS.preds2.sort <- sort(TS.preds, sortFactor = "Sources", sortOder = sort.order)
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
  
  #Use sort.predictions.frame and save order for use with other response variables
  TS.preds.sort <- sort(TS.preds, classify = "Sources:Type", sortFactor = "Sources", 
                        sortParallelToCombo = list(Type = "Control"))
  sort.order <- attr(TS.preds.sort, which = "sortOrder")
  
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
}  
}
\keyword{asreml}