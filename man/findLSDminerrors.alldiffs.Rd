\name{findLSDminerrors.alldiffs}
\alias{findLSDminerrors.alldiffs}
\alias{findLSDminerrors}
\title{Find LSD values that minimize the number of errors in pairwise comparisons of predictions.}
\description{Given an \code{\link{alldiffs.object}} with an \code{sed} component, 
             a search is made of a set of equally spaced values between the minimum and 
             maximum values of the LSDs, calculated from the \code{sed} component of the 
             \code{\link{alldiffs.object}}, to identify LSD values that minimize the number 
             of errors made in deciding on the significance of pairs of predicted values stored 
             in the \code{\link{alldiffs.object}}. If \code{LSDtype} is set to \code{overall}, 
             a search is made over the range of LSD values for all pairwise comparisons for a 
             single LSD value; if \code{LSDtype} is set to \code{factor.combinations}, 
             a separate search is made over the LSD values for the set of pairwise comparisons 
             for each \code{factor.combination} in order to identify a single value for each set. 
             The number of values used in the search is controlled by the argument \code{nvalues}. 
             For each value in the search, the numbers of false positives and false negatives 
             resulting from employing it as the LSD for each set of pairwise comparisons is calculated. 
             A criterion that combines the false positives and negative is calculated using the 
             \code{false.pos.wt}, the criterion being the number of false postives times the 
             \code{false.pos.wt} plus the number of false negatives. The value chosen for the LSD 
             is the smallest value from amongst those with the minimum value of the criterion and 
             the least number of false positives. A secondary search with 10 equally spaced values 
             is made of the interval below the chosen value and the search value immediately below 
             it to check whether the chosen grid value can be further reduced without changing the 
             value of either its criterion or the number of false positives.
             
             The primary options for changing the numbers of errors associated with the values 
             resulting from the searching is to manipulate the \code{LSDby} and/or 
             \code{false.pos.wt} arguments.}
\usage{
\method{findLSDminerrors}{alldiffs}(alldiffs.obj, 
                 LSDtype = "overall", LSDby = NULL, 
                 alpha = 0.05, 
                 false.pos.wt = 10, nvalues = 100,
                 retain.zeroLSDs = FALSE, 
                 zero.tolerance = .Machine$double.eps ^ 0.5, 
                 trace = FALSE, ...)}
\arguments{
\item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
\item{LSDtype}{A \code{\link{character}} string that can be \code{overall} or \code{factor.combinations}. 
              It determines whether the  minimum LSD values that are obtained are those for   
              (i) all pairwise comparisons of the predicted values, unless there is only one, 
              possibly repeated, prediction, when a notional LSD is calculated, or 
              (ii) for each combination of the values of the \code{\link{factor}s} and 
              \code{\link{numeric}s} named in \code{LSDby}, unless there is only one prediction for a 
              combination, when notional LSDs are calculated. The \code{LSDtype} specified here does not 
              have to match that used in the creating the \code{\link{alldiffs.object}}.}
\item{LSDby}{A \code{\link{character}} (vector) of variables names, being the names of some of the 
             \code{\link{factor}s} or \code{\link{numeric}s} in the \code{classify}; for each 
             combination of the values of the \code{\link{factor}s} and \code{\link{numeric}s}, the 
             LSD errors are to be computed when \code{LSDtype} is \code{factor.combinatons}. The \code{LSDby} 
             specified here does not have to match that used in the creating the 
             \code{\link{alldiffs.object}}.}
\item{alpha}{A \code{\link{numeric}} specifying the significance level for an LSD to compare a pair of 
            predictions.}
\item{false.pos.wt}{A \code{\link{numeric}} that specifies the 
            weight (e.g. 3 so that a false positive is considered to be equivalent to three false negatives) 
            to apply to the number of false positives in calculating the weighted sums of the numbers of 
            false positives and negatives that is used as the criterion to be minimized in selecting 
            the LSD value that results in the minimum number of errors. If, for a particular LSD  value, 
            the number of false positives is \emph{p}, the number of false negatives is \emph{n} and 
            \emph{w} the value of \code{false.pos.wt}, then then criterion for that LSD value is 
            \eqn{(w p) + n}{(w p) + n}.
            
            The default of 10 for \code{false.pos.wt} greatly favours false negatives; a value with 
            one false positive and no false negative will only be chosen over a value with 
            no false positive when the latter value has more than 10 false negatives.}
\item{nvalues}{A \code{\link{numeric}} specifying the number of equally spaced LSD values, between the 
            minimum and maximum LSD values for a set of paired comparisons, to be used in the search 
            for the LSD value with the minimum number of errors.}
\item{retain.zeroLSDs}{A \code{\link{logical}} indicating whether to retain or omit LSDs that are zero 
             when calculating the numbers of errors.}
\item{zero.tolerance}{A \code{\link{numeric}} specifying the value such that if an LSD is less than it, 
             the LSD will be considered to be zero.}
\item{trace}{A \code{\link{logical}} indicating whether details of the searching are to be output.}
\item{\dots}{Provision for passing arguments to functions called internally - 
            not used at present.}
}
\value{A \code{\link{data.frame}} containing the chosen LSD(s), its(their) numbers of false positives 
            and negatives and the value(s) of the false criterion.
}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{exploreLSDs.alldiffs}} \code{\link{plotLSDs.data.frame}}, 
         \code{\link{plotLSDs.alldiffs}}, \cr 
         \code{\link{plotLSDerrors.alldiffs}}, \code{\link{plotLSDerrors.data.frame}}, \code{\link{recalcLSD.alldiffs}}, \cr 
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
  TS.vcov <- vcov(TS.emm)
  TS.diffs <- allDifferences(predictions = TS.preds, classify = "Sources:Type", 
                             vcov = TS.vcov, tdf = den.df)
  validAlldiffs(TS.diffs)
}  

## Choose LSD values with the minimum mumber of error for pairwise comparisons of 
##   the predictions obtained using asreml or lmerTest
if (exists("TS.diffs"))
{
  ##Pick the LSD values for predictions obtained using asreml or lmerTest  
  minLSD <- findLSDminerrors(TS.diffs)
  TS.diffs <- redoErrorIntervals(TS.diffs, LSDtype = "supplied", LSDsupplied = minLSD["LSD"])
  TS.diffs$LSD
  minLSDs <- findLSDminerrors(TS.diffs, LSDtype = "factor.combinations", 
                              LSDby = "Sources")
  TS.diffs <- redoErrorIntervals(TS.diffs, LSDtype = "supplied", 
                                 LSDby = "Sources", LSDsupplied = minLSDs["LSD"])
  TS.diffs$LSD
}
}
\keyword{asreml}