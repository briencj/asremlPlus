\name{exploreLSDs.alldiffs}
\alias{exploreLSDs.alldiffs}
\alias{exploreLSDs}
\title{Explores the computed LSD values for pairwise differences between predictions.}
\description{Given an \code{\link{alldiffs.object}} with an \code{sed} component, the LSDs are 
             calculated for all pairwise comparisons of predictions. It then calculates 
             (i) a table of frequencies of the LSD values, (ii) the distinct values of the LSDs 
             after rounding, (iii) various statistics from the LSD values, (iv) a measure of 
             the accuracy of each of the LSD statistics, (v) the accuracy of each statistic in 
             representing the LSD values for each prediction and (vi) a matrix containing 
             the LSD values for comparinbg each pair of predictions. Histograms of the 
             frequencies can also be produced.}
\usage{
\method{exploreLSDs}{alldiffs}(alldiffs.obj, LSDtype = "overall", LSDby = NULL, 
            LSDaccuracy = "maxAbsDeviation", alpha = 0.05, digits = 3, 
            retain.zeroLSDs = FALSE, 
            zero.tolerance = .Machine$double.eps ^ 0.5, 
            plotHistogram = FALSE, ...)}
\arguments{
\item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
\item{LSDtype}{A \code{\link{character}} string that can be \code{overall} or \code{factor.combinations}. 
              It determines whether the LSD values that are investigated and stored are  
              (i) the \code{overall} minimum, quantile10, quantile25, mean, median, quantile75, quantile90, 
              and maximum of all pairwise LSDs, or 
              (ii) the minimum, quantile10, quantile25, mean, median, quantile75, quantile90, and maximum 
               for the pairwise LSDs for each \code{factor.combination}, unless there is only one prediction 
              for a \code{factor.combination}, when notional LSDs are calculated. The \code{LSDtype} 
              specified here does not have to match that used in the creating the \code{\link{alldiffs.object}}.}
\item{LSDby}{A \code{\link{character}} (vector) of variables names, being the names of the 
             \code{\link{factor}s} or \code{\link{numeric}s} in the \code{classify} for each 
             combination of which LSD statistics and accuracy are to be computed, as well histograms 
             plotted, when \code{LSDtype} is \code{factor.combinatons}. The \code{LSDby} specified here does 
             not have to match that used in the creating the \code{\link{alldiffs.object}}.}
\item{LSDaccuracy}{A \code{\link{character}} nominating one of \code{maxAbsDeviation}, \code{maxDeviation}, 
             \code{q90Deviation} or \code{RootMeanSqDeviation} as the statistic to be calculated as a measure 
             of the accuracy of an LSD statistic. The option \code{q90Deviation} produces the sample quantile 
             corresponding to a probability of 0.90. The deviations are the differences between the LSDs used in 
             calculating an LSD statistic and the LSD statistic and the accuracy is expressed as a 
             proportion of the value of the LSD statistic.}
\item{alpha}{A \code{\link{numeric}} specifying the significance level for an LSD to compare a pair of 
            predictions.}
\item{digits}{A \code{\link{numeric}} specifying the number of significant digits to retain in rounding 
            the LSDs before determining the disinct rounded LSDs.}
\item{retain.zeroLSDs}{A \code{\link{logical}} indicating whether to retain or omit LSDs that are zero when 
             calculating the summaries of LSDs.}
\item{zero.tolerance}{A \code{\link{numeric}} specifying the value such that if an LSD is less than it, it will be 
            considered to be zero.}
\item{plotHistogram}{A \code{\link{logical}} indicating whether or not histograms of the LSD values are to 
            be plotted. The \code{LSDtype} argument controls whether one histogram of all LSD values is 
            plotted or histograms are plotted for each combination of the levels of the factors specified 
            by the \code{LSDby} argument.}
\item{\dots}{Provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{A \code{\link{list}} with components 
       \code{frequencies}, \code{distinct.vals}, \code{statistics}, \code{accuracy}, \code{per.pred.accuracy} 
       and \code{LSD}:  \code{frequencies} is a \code{data.frame} with the frequency distribution of the LSD 
       values; \code{distinct.vals} is a \code{list}, each component of which contains the distinct values 
       of the LSDs after rounding; \code{statistics} is a \code{data.frame} with the minimum, quantile10, 
       quantile25, mean, median, quantile75, quantile90, and maximum of the LSD values; \code{accuracy} is 
       a \code{data.frame} with the accuracies of the minimum, quantile10, quantile25, mean, median, 
       quantile75, quantile90, and maximum of the LSD values with respect to the values from which these 
       statistics are calculated; \code{per.pred.accuracy} is a \code{data.frame} with the accuracies of the 
       minimum, quantile10, quantile25, mean, median, quantile75, quantile90, and maximum of the LSD values 
       for a set of predictions when these statistics are used to represent the LSDs for the comparisons 
       amongst the set of predictions; \code{LSD} is a square \code{matrix} containing the LSD values for 
       all pairwise comparisons of the predictions contained in the supplied \code{alldiffs.obj}. 
       The accuracy measure is specified by the \code{LSDaccuracy} argument.}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{plotLSDs.data.frame}}, \code{\link{plotLSDs.alldiffs}},  
         \code{\link{recalcLSD.alldiffs}}, \cr 
         \code{\link{redoErrorIntervals.alldiffs}}}
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
  TS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Type)
  TS.preds <- summary(TS.emm)
  den.df <- min(TS.preds$df, na.rm = TRUE)
  ## Modify TS.preds to be compatible with a predictions.frame
  TS.preds <- as.predictions.frame(TS.preds, predictions = "emmean", 
                                   se = "SE", interval.type = "CI", 
                                   interval.names = c("lower.CL", "upper.CL"))
   
  ## Form an all.diffs object and check its validity
  els <- as.numeric(rownames(TS.preds))
  TS.vcov <- vcov(TS.emm)[els,els]
  TS.diffs <- allDifferences(predictions = TS.preds, classify = "Sources:Type", 
                             vcov = TS.vcov, tdf = den.df)
  validAlldiffs(TS.diffs)
}  

## Plot p-values for predictions obtained using asreml or lmerTest
if (exists("TS.diffs"))
{
  ##Explore the LSD values for predictions obtained using asreml or lmerTest  
  LSDstat <- exploreLSDs(TS.diffs, LSDtype = "factor.combinations", 
                         LSDby = "Sources")
}
}
\keyword{asreml}