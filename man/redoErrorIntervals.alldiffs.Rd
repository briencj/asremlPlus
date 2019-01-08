\name{redoErrorIntervals.alldiffs}
\alias{redoErrorIntervals.alldiffs}
\alias{redoErrorIntervals}
\title{Adds or replaces the error intervals stored in a prediction component 
       of an \code{\link{alldiffs.object}}.}
\description{Given an \code{\link{alldiffs.object}}, adds or replaces \code{error.intervals} 
             for its \code{prediction} component. If the \code{backtransforms} component is 
             present, \code{error.intervals} will be calculated for it too. 
}
\usage{
\method{redoErrorIntervals}{alldiffs}(alldiffs.obj, error.intervals = "Confidence", 
                   alpha = 0.05, avsed.tolerance = 0.25, 
                   meanLSD.type = "overall", LSDby = NULL, ...)
}
\arguments{
 \item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
 \item{error.intervals}{A \code{\link{character}} string indicating the type of error interval, if any, 
             to calculate in order to indicate uncertainty in the predicted values. 
             Possible values are \code{"none"}, \code{"StandardError"}, \code{"Confidence"} 
             and \code{"halfLeastSignificant"}. The default is for confidence limits to 
             be used. The \code{"halfLeastSignificant"} option results in half the 
             Least Significant Difference (LSD) being added and subtracted to the 
             predictions, the LSD being calculated using the square root of the mean of the 
             variances of all or a subset of pairwise differences between the predictions. 
             If \code{meanLSD.type} is set to \code{overall}, the \code{avsed.tolerance} is not 
             \code{NA} and the range of the SEDs divided by the average of the SEDs exceeds 
             \code{avsed.tolerance} then the \code{error.intervals} calculations and the plotting 
             will revert to confidence intervals. Also, half LSDs cannot be used for 
             backtansformed values and so confidence intervals will be used instead.}
\item{alpha}{The significance level for an LSD to compare a pair of predictions.}
\item{avsed.tolerance}{A \code{\link{numeric}} giving the value of the SED range, the range of the SEDs 
             divided by the square root of the mean of the variances of all or a subset of the 
             pairwise differences, that is considered reasonable in calculating 
             \code{error.intervals}. It should be a value between 0 and 1. The following rules apply:
             \enumerate{
             \item If \code{avsed.tolerance} is \code{NA} then mean LSDs of the type specified by 
             \code{meanLSD.type} are calculated and used in \code{error.intervals} and plots.
             \item Irrespective of the setting of \code{meanLSD.type}, if \code{avsed.tolerance} is not 
             exceeded then the mean LSDs are used in \code{error.intervals} and plots. 
             \item If \code{meanLSD.type} is set to \code{overall}, \code{avsed.tolerance} is not 
             \code{NA}, and \code{avsed.tolerance} is exceeded then  \code{error.intervals} and 
             plotting revert to confidence intervals. 
             \item If \code{meanLSD.type} is set to \code{factor.combinations} and \code{avsed.tolerance} 
             is not exceeded for any factor combination then the half LSDs are 
             used in \code{error.intervals} and plots; otherwise, \code{error.intervals} and 
             plotting revert to confidence intervals.
             \item If \code{meanLSD.type} is set to \code{per.prediction} and \code{avsed.tolerance} 
             is not exceeded for any prediction then the half LSDs are used in \code{error.intervals} 
             and plots; otherwise, \code{error.intervals} and plotting revert to confidence intervals.
             }}
\item{meanLSD.type}{A \code{\link{character}} string determining whether the mean LSD stored is 
             (i) the \code{overall} mean, based on the square root of the mean of the variances of 
             all pairwise variances, (ii) the mean for each \code{factor.combination} of the 
             \code{\link{factor}s} specified by \code{LSDby}, which is based on the square root of 
             the mean of the variances for all pairwise differences for each factor combination, or 
             (iii) the \code{per.prediction} mean, based, for each prediction, 
             on the square root of the mean of the variances for all pairwise differences 
             involving that prediction. It also 
             determines, in conjunction with \code{avsed.tolerance}, which LSD will be used
             in calculating \code{error.intervals} and, hence, is used for plots.}
\item{LSDby}{A \code{\link{character}} (vector) of factor names, being the names of the 
             \code{\link{factor}s} for each of whose levels combinations a mean LSD, 
             minLSD and max LSD is stored in the \code{LSD} component of the 
             \code{\link{alldiffs.object}} object when \code{meanLSD.type} is \code{factor.combinatons}.}
 \item{\dots}{provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{An \code{\link{alldiffs.object}} with components 
       \code{predictions}, \code{vcov}, \code{differences}, \code{p.differences} 
       \code{sed}, and \code{LSD}. 

           If \code{error.intervals} is not \code{"none"}, then the 
           \code{predictions} component and, if present, the 
           \code{backtransforms} component will contain columns for the lower 
           and upper values of the limits for the interval. The names of these 
           columns will consist of three parts separated by full stops: 
           1) the first part will be \code{lower} or \code{upper}; 
           2) the second part will be one of \code{Confidence}, 
              \code{StandardError} or \code{halfLeastSignificant}; 
           3) the third component will be \code{limits}.

           The name of the \code{response}, the \code{term}, the \code{classify} 
           and \code{tdf}, as well as the degrees of freedom 
           of the standard error, will be set as attributes to the object.}
\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}},  
         \code{\link{subset.alldiffs}}, \cr
         \code{\link{allDifferences.data.frame}}, 
         \code{\link{recalcLSD.alldiffs}},  \code{\link{predictPresent.asreml}},  \cr
         \code{\link{plotPredictions.data.frame}}, 
         \code{as.Date}, \code{predict.asreml}}
\examples{
data(WaterRunoff.dat)

##Use asreml to get predictions and associated statistics

\dontrun{
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- asrtests(current.asr, NULL, NULL)
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
  TS.diffs <- allDifferences(predictions = TS.preds, 
                             classify = "Sources:Species", 
                             vcov = TS.vcov, tdf = den.df)
  validAlldiffs(TS.diffs)
}  

## Plot p-values for predictions obtained using asreml or lmerTest
if (exists("TS.diffs"))
{
  ##Recalculate the LSD values for predictions obtained using asreml or lmerTest  
  TS.diffs <- redoErrorIntervals.alldiffs(TS.diffs, 
                                          error.intervals = "halfLeastSignificant")

}
}
\keyword{asreml}