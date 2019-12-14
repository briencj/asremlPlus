\name{subset.alldiffs}
\alias{subset.alldiffs}
\title{Subsets the components in an \code{\link{alldiffs.object}} according to the supplied condition.}
\description{Subsets each of the components of an \code{\link{alldiffs.object}}. The subset is 
             determined by applying the condition to the \code{prediction} component to 
             determine which of its rows are to be included in the subset. Then, if present, 
             this subset is applied to the rows of \code{backtransforms} and to the rows 
             and columns of \code{differences}, \code{p.differences} and \code{sed} 
             components. In addition, if \code{sed} is present, 
             \code{\link{recalcLSD.alldiffs}} is called to recalculate the values in the 
             \code{LSD} component, with any arguments supplied via the \code{\dots} 
             argument passed ot it.

             The \code{select} argument of \code{subset} is not implemented, but can be 
             achieved for variables in the \code{classify} using the \code{rmClassifyVars} 
             argument.}
\usage{\method{subset}{alldiffs}(x, subset = rep(TRUE, nrow(x$predictions)), 
       rmClassifyVars = NULL, ...)}
\arguments{
 \item{x}{An \code{\link{alldiffs.object}}.}
 \item{subset}{A \code{\link{logical}} that detemines rows of the \code{predictions}  
              component of \code{x} to be included in the subset. By default allws are 
              included.}
 \item{rmClassifyVars}{A \code{\link{character}} that contains the names of the 
              variables in the \code{classify} attribute of \code{x} that are to be 
              removed from the predictions \code{data.frame} and the names of the 
              dimensions of the other components of \code{x}. In doing this, 
              the combinations of the remaining \code{classify} variables must 
              uniquely index the predictions.}
 \item{\dots}{further arguments passed to \code{\link{recalcLSD.alldiffs}}.}
}
\value{An \code{\link{alldiffs.object}} with the following components of the supplied 
       \code{\link{alldiffs.object}} subsetted, if present in the original object: 
       \code{predictions}, \code{vcov}, \code{backtransforms}, \code{differences}, 
       \code{p.differences} and \code{sed}. In addition, if \code{sed} is present, the 
       \code{LSD} component will be recalculated.}

\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}}, \cr
         \code{\link{redoErrorIntervals.alldiffs}}, \code{\link{recalcLSD.alldiffs}}, \cr
          \code{\link{predictPlus.asreml}}, \code{\link{predictPresent.asreml}}}
\examples{
data(WaterRunoff.dat)

##Use asreml to get predictions and associated statistics

\dontrun{
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= WaterRunoff.dat)
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
  ##Use subset.alldiffs to select a subset of the alldiffs object
  TS.diffs.subs <- subset(TS.diffs, 
                          subset = grepl("R", Sources, fixed = TRUE) & 
                                   Type \%in\% c("Control","Medicinal"))
}
}
\keyword{asreml}