\name{recalcLSD.alldiffs}
\alias{recalcLSD.alldiffs}
\alias{recalcLSD}
\title{Adds or recalculates the LSD component of an \code{\link{alldiffs.object}}.}
\description{Given an \code{\link{alldiffs.object}}, adds or recalculate its LSD component.
}
\usage{
\method{recalcLSD}{alldiffs}(alldiffs.obj, meanLSD.type = "overall", LSDby = NULL, 
          alpha = 0.05, ...)}
\arguments{
\item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
 \item{meanLSD.type}{A \code{\link{character}} string determining whether the mean LSD stored is 
             (i) the \code{overall} mean, based on the square root of the mean of the variances of 
             all pairwise variances, (ii) the mean for each \code{factor.combination} of the 
             \code{\link{factor}s} specified by \code{LSDby}, which is based on the square root of 
             the mean of the variances for all pairwise differences for each factor combination, unless 
             there is only one predction for a \code{factor.combination}, when notional LSDs are 
             calculated that are based on the standard error of the prediction multiplied by the square 
              root of two, or 
             (iii) the \code{per.prediction} mean, based, for each prediction, 
             on the square root of the mean of the variances for all pairwise differences involving 
             that prediction. It also 
             determines, in conjunction with \code{avsed.tolerance}, which LSD will be used
             in calculating \code{error.intervals} and, hence, is used for plots.}
\item{LSDby}{A \code{\link{character}} (vector) of variables names, being the names of the 
             \code{\link{factor}s} or \code{\link{numeric}s} in the \code{classify} for each 
             combination of which a mean LSD, minLSD and max LSD is stored in the \code{LSD} 
             component of the \code{\link{alldiffs.object}} when \code{meanLSD.type} is 
             \code{factor.combinatons}.}
\item{alpha}{The significance level for an LSD to compare a pair of predictions.}
\item{\dots}{further arguments passed to \code{\link{allDifferences.data.frame}}.}
}
\value{An \code{\link{alldiffs.object}} with components 
       \code{predictions}, \code{vcov}, \code{differences}, \code{p.differences} 
       \code{sed}, \code{LSD} and, if present in \code{alldiffs.obj}, \code{backtransforms}.}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{as.alldiffs}}, \code{\link{sort.alldiffs}}, 
         \code{\link{subset.alldiffs}}, \code{\link{print.alldiffs}},  \cr 
         \code{\link{renewClassify.alldiffs}}, \code{\link{redoErrorIntervals.alldiffs}},  
         \code{\link{plotPredictions.data.frame}}, \cr 
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
TS.diffs <- predictPlus(classify = "Sources:Type", 
                        asreml.obj = current.asr, 
                        wald.tab = current.asrt$wald.tab, 
                        present = c("Sources", "Type", "Species"))
}

## Use lmeTest and emmmeans to get predictions and associated statistics

if (requireNamespace("lmerTest", quietly = TRUE) & 
    requireNamespace("emmeans", quietly = TRUE))
{
  m1.lmer <- lmerTest::lmer(pH ~ Benches + (Sources * (Type + Species)) + 
                              (1|Benches:MainPlots),
                            data=na.omit(WaterRunoff.dat))
  TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Species)
  TS.preds <- summary(TS.emm)
  den.df <- min(TS.preds$df, na.rm = TRUE)
  ## Modify TS.preds to be compatible with a predictions.frame
  TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                   se = "SE", interval.type = "CI", 
                                   interval.names = c("lower.CL", "upper.CL"))
  
  ## Form an all.diffs object and check its validity
  TS.vcov <- vcov(TS.emm)
  TS.diffs <- allDifferences(predictions = TS.preds, classify = "Sources:Species", 
                             vcov = TS.vcov, tdf = den.df)
  validAlldiffs(TS.diffs)
}  

## Plot p-values for predictions obtained using asreml or lmerTest
if (exists("TS.diffs"))
{
  ##Recalculate the LSD values for predictions obtained using asreml or lmerTest  
  TS.diffs <- recalcLSD.alldiffs(TS.diffs, meanLSD.type = "factor.combinations", 
                                 LSDby = "Sources")
}
}
\keyword{asreml}