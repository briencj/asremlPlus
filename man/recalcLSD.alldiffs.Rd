\name{recalcLSD.alldiffs}
\alias{recalcLSD.alldiffs}
\alias{recalcLSD}
\title{Adds or recalculates the \code{\link{LSD.frame}} that is a component of an 
             \code{\link{alldiffs.object}}.}
\description{Given an \code{\link{alldiffs.object}}, adds or recalculate its \code{\link{LSD.frame}}. 
\bold{N.B. No changes are made to the error.intervals --- use \code{\link{redoErrorIntervals.alldiffs}}
to modify both the error.intervals and the \code{\link{LSD.frame}}.}
}
\usage{
\method{recalcLSD}{alldiffs}(alldiffs.obj, LSDtype = "overall", LSDsupplied = NULL, 
          LSDby = NULL, LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
          alpha = 0.05, ...)}
\arguments{
\item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
\item{LSDtype}{A \code{\link{character}} string that can be \code{overall}, \code{factor.combinations}, 
              \code{per.prediction} or \code{supplied}. It determines whether the values stored in a row  
              of a \code{\link{LSD.frame}} are the values calculated 
             (i) \code{overall} from the LSD values for all pairwise comparison2, 
             (ii) the values calculated from the pairwise LSDs for the levels of each 
             \code{factor.combination}, unless there is only one prediction for a level  of the 
             \code{factor.combination}, when a notional LSD is calculated,  
             (iii) \code{per.prediction}, being based, for each prediction, on all pairwise differences 
             involving that prediction, or 
             (iv) as \code{supplied} values of the LSD, specified with the \code{LSDsupplied} argument; 
             these supplied values are to be placed in the \code{assignedLSD} column of the 
             \code{\link{LSD.frame}} stored in an \code{\link{alldiffs.object}} so that they can be used 
             in LSD calculations.
             
             See \code{\link{LSD.frame}} for further information on the values in a row of this 
             \code{data.frame} and how they are calculated.}
\item{LSDsupplied}{A \code{\link{data.frame}} or a named \code{\link{numeric}} containing a set of \code{LSD} 
             values that correspond to the observed combinations of the values of the \code{LSDby} variables 
             in the \code{\link{predictions.frame}} or a single LSD value that is an overall LSD. 
             If a \code{\link{data.frame}}, it may have (i) a column for the \code{LSDby} variable and a column 
             of \code{LSD} values or (ii) a single column of \code{LSD} values with rownames being the 
             combinations of the observed values of the \code{LSDby} variables. Any name can be used 
             for the column of \code{LSD} values; \code{assignedLSD} is sensible, but not obligatory. Otherwise, 
             a \code{\link{numeric}} containing the \code{LSD} values, each of which is named for the observed 
             combination of the values of the \code{LSDby} variables to which it corresponds. (Applying the 
             \code{function} \code{dae::fac.combine} to the \code{predictions} component is one way of 
             forming the required combinations for the (row) names.)  The values supplied 
             will be incorporated into \code{assignedLSD} column of the \code{\link{LSD.frame}} stored as the 
             \code{LSD} component of the \code{\link{alldiffs.object}}.}
\item{LSDby}{A \code{\link{character}} (vector) of variables names, being the names of the 
             \code{\link{factor}s} or \code{\link{numeric}s} in the \code{classify}; for each 
             combination of their levels and values, there will be or is a row in the \code{\link{LSD.frame}} 
             stored in the \code{LSD} component of the \code{\link{alldiffs.object}} when \code{LSDtype} is 
             \code{factor.combinatons}.}
\item{LSDstatistic}{A \code{\link{character}} nominating one or more of \code{minimum}, \code{q10}, \code{q25},  
             \code{mean}, \code{median}, \code{q75}, \code{q90} or \code{maximum} as the value(s) to be 
             stored in the \code{assignedLSD} column in an \code{\link{LSD.frame}}; the values in the 
             \code{assignedLSD} column are used in computing \code{halfLeastSignificant} \code{error.intervals}. 
             Here \code{q10},  \code{q25}, \code{q75} and \code{q90} indicate the sample quantiles corresponding 
             to probabilities of 0.1, 0.25, 0.75 and 0.9 for the group of LSDs from which a single LSD value 
             is calculated. The function \code{\link{quantile}} is used to obtain them. The \code{mean} LSD is 
             calculated as the square root of the mean of the squares of the LSDs for the group. The 
             \code{median} is calculated using the \code{\link{median}} function. Multiple values are only 
             produced for \code{LSDtype} set to \code{factor.combination}, in which case \code{LSDby} must 
             not be \code{NULL} and the number of values must equal the number of observed combinations of 
             the values of the variables specified by \code{LSDby}.  If \code{LSDstatistic} is \code{NULL}, 
             it is reset to \code{mean}.}
\item{LSDaccuracy}{A \code{\link{character}} nominating one of \code{maxAbsDeviation}, \code{maxDeviation}, 
             \code{q90Deviation} or \code{RootMeanSqDeviation} as the statistic to be calculated as a measure 
             of the accuracy of \code{assignedLSD}. The option \code{q90Deviation} produces the sample quantile 
             corresponding to a probability of 0.90. The deviations are the differences between the LSDs used in 
             calculating the LSD statistics and each assigned LSD and the accuracy is expressed as a 
             proportion of the assigned LSD value. The calculated values are stored in the column named 
             \code{accuracyLSD} in an \code{\link{LSD.frame}}.}
\item{alpha}{The significance level for an LSD to compare a pair of predictions.
              It is stored as an attribute to the \code{\link{alldiffs.object}}.}

\item{\dots}{further arguments passed to \code{\link{allDifferences.data.frame}}; 
              attributes \code{transform.power}, \code{offset} and \code{scale} cannot be passed.}
}
\value{An \code{\link{alldiffs.object}} with components 
       \code{predictions}, \code{vcov}, \code{differences}, \code{p.differences} 
       \code{sed}, \code{LSD} and, if present in \code{alldiffs.obj}, \code{backtransforms}.}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{as.alldiffs}}, \code{\link{sort.alldiffs}}, 
         \code{\link{subset.alldiffs}}, \code{\link{print.alldiffs}},  \cr 
         \code{\link{renewClassify.alldiffs}}, \code{\link{exploreLSDs.alldiffs}}, 
         \code{\link{pickLSDstatistics.alldiffs}}, \cr
          \code{\link{redoErrorIntervals.alldiffs}},   
         \code{\link{plotPredictions.data.frame}}, 
         \code{\link{predictPlus.asreml}}, \cr
         \code{\link{predictPresent.asreml}}}
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
  TS.diffs <- recalcLSD.alldiffs(TS.diffs, LSDtype = "factor.combinations", 
                                 LSDby = "Sources")
}
}
\keyword{asreml}