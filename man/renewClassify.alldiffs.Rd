\name{renewClassify.alldiffs}
\alias{renewClassify.alldiffs}
\alias{renewClassify}
\title{Renews the components in an \code{\link{alldiffs.object}} according to a new classify.}
\description{The \code{classify} is an attribute of an \code{\link{alldiffs.object}} and determines 
             the order within the components of an unsorted \code{\link{alldiffs.object}}. 
             This function resets the \code{classify} attribute and re-orders the components of 
             \code{\link{alldiffs.object}} to be in standard order for the variables in a 
             \code{newclassify}, using \code{\link{allDifferences.data.frame}}. The \code{newclassify} 
             may be just a re-ordering of the variable names in the previous \code{classify}, or be 
             based on a new set of variable names. The latter is particularly useful when 
             \code{\link{linTransform.alldiffs}} has been used with a \code{\link{matrix}} and it 
             is desired to replace the resulting \code{Combination} \code{classify} with a 
             \code{newclassify} comprised of a more meaningful set of variables that have 
             replaced \code{Combination} in the \code{predictions} component.}
\usage{\method{renewClassify}{alldiffs}(alldiffs.obj, newclassify, 
              sortFactor = NULL, sortWithinVals = NULL, 
              sortOrder = NULL, decreasing = FALSE, ...)}
\arguments{
\item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
\item{newclassify}{A \code{\link{character}} string giving the variables that 
                 define the margins of the multiway table that was predicted, 
                 but ordered so that the predictions are in the desired order when they are 
                 arranged in standard order for the \code{newclassify}. 
                 Multiway tables are specified by forming an interaction type 
                 term from the classifying variables, that is, separating the 
                 variable names with the \code{:} operator. The number of combined values of the 
                 set of variable name(s) must equal the number of rows in the \code{predictions} 
                 component.}
\item{sortFactor}{A \code{\link{character}} containing the name of the 
              \code{factor} that indexes the set of predicted values that determines 
               the sorting of the components. If there is only one variable in the 
              \code{classify} term then \code{sortFactor} can be \code{NULL} and 
              the order is defined by the 
              complete set of predicted values. If there is more than one variable 
              in the \code{classify} term then \code{sortFactor} must be set. 
              In this case the \code{sortFactor} is sorted for the predicted values 
              within each combination of the values of the sortWithin 
              variables: the \code{classify} variables, excluding the 
              \code{sortFactor}. There should be only one predicted value for 
              each unique value of \code{sortFactor} within each set defined by a 
              combination of the values of the sortWithin variables.}
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
 \item{decreasing}{A \code{\link{logical}} passed to \code{order} that detemines whether 
              the order is for increasing or decreasing magnitude of the predicted 
              values.}
 \item{\dots}{further arguments passed to \code{\link{allDifferences.data.frame}}.}
}
\value{The \code{\link{alldiffs.object}} supplied with the following components, 
       if present, sorted: \code{predictions}, \code{vcov}, \code{backtransforms}, \code{differences}, 
       \code{p.differences} and \code{sed}. Also, the \code{sortFactor} and \code{sortOrder} 
       attributes are set.}
\details{First, the components of the \code{\link{alldiffs.object}} is arranged in standard order for
         the \code{newclassify}. Then predictions are reordered according to the settings of 
         \code{sortFactor}, \code{sortWithinVals}, \code{sortOrder} and \code{decreasing} (see 
         \code{\link{sort.alldiffs}} for details).}


\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}}, \cr 
         \code{\link{redoErrorIntervals.alldiffs}}, \code{\link{recalcLSD.alldiffs}},  \cr
          \code{\link{predictPlus.asreml}}, \code{\link{predictPresent.asreml}}}
\examples{
data(WaterRunoff.dat)

##Use asreml to get predictions and associated statistics

\dontrun{
#Analyse pH  
m1.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                 random = ~ Benches:MainPlots,
                 keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- as.asrtests(m1.asr, NULL, NULL)
current.asrt <- as.asrtests(m1.asr)
current.asrt <- rmboundary(current.asrt)
m1.asr <- current.asrt$asreml.obj

#Get predictions and associated statistics  
TS.diffs <- predictPlus.asreml(classify = "Sources:Type", 
                               asreml.obj = m1.asr, tables = "none", 
                               wald.tab = current.asrt$wald.tab, 
                               present = c("Type","Species","Sources"))
}

## Use lmeTest and emmmeans to get predictions and associated statistics

if (requireNamespace("lmerTest", quietly = TRUE) & 
    requireNamespace("emmeans", quietly = TRUE))
{
  #Analyse pH
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
  TS.diffs <- allDifferences(predictions = TS.preds, 
                               classify = "Sources:Type", 
                               vcov = TS.vcov, tdf = den.df)
  validAlldiffs(TS.diffs)
}  

#Re-order predictions from asreml or lmerTest so all Sources for the same Type are together 
#for each combination of A and B
if (exists("TS.diffs"))
{
  TS.diffs.reord <- renewClassify(TS.diffs, newclassify = "Type:Sources")
  validAlldiffs(TS.diffs.reord)
}
}
\keyword{asreml}