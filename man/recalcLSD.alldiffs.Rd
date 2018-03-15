\name{recalcLSD.alldiffs}
\alias{recalcLSD.alldiffs}
\alias{recalcLSD}
\title{Adds or recalculates the LSD component of an alldiffs object.}
\description{Given an alldiffs object, adds or recalculate its LSD component.
}
\usage{
\method{recalcLSD}{alldiffs}(alldiffs.obj, meanLSD.type = "overall", LSDby = NULL, 
          alpha = 0.05, ...)}
\arguments{
 \item{alldiffs.obj}{An alldiffs object that is a list with components 
       \code{predictions} containing the predictions and 
       their standard errors, \code{vcov} containing the variance matrix of the 
       \code{predictions}, \code{differences} containing all pairwise 
       differences between the predictions, \code{p.differences} containing 
       p-values for all pairwise differences between the predictions, 
       \code{sed} containing the standard errors of all pairwise differences 
       between the predictions, and an \code{LSD} containing the mean, minimum 
       and maximum LSDs.}
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
             minLSD and max LSD is stored in the \code{LSD} component of the \code{alldiffs} 
             object when \code{meanLSD.type} is \code{factor.combinatons}.}
\item{alpha}{The significance level for an LSD to compare a pair of predictions.}
\item{\dots}{further arguments passed to \code{\link{allDifferences.data.frame}}.}
}
\value{An alldiffs object that is a list with components 
       \code{predictions} containing the predictions and 
       their standard errors, \code{vcov} containing the variance matrix of the 
       \code{predictions}, \code{differences} containing all pairwise 
       differences between the predictions, \code{p.differences} containing 
       p-values for all pairwise differences between the predictions, 
       \code{sed} containing the standard errors of all pairwise differences 
       between the predictions, and an \code{LSD} containing the mean, minimum 
       and maximum LSDs.}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{as.alldiffs}}, \code{\link{sort.alldiffs}}, 
         \code{\link{subset.alldiffs}}, \code{\link{print.alldiffs}},  \cr 
         \code{\link{redoErrorIntervals.alldiffs}},  \code{\link{plotPredictions.data.frame}}, 
         \code{\link{predictPlus.asreml}}, \cr
         \code{\link{predictPresent.asreml}}}
\examples{\dontrun{
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE) #required for asreml4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- asrtests(current.asr, NULL, NULL)
diffs <- predictPlus(classify = "Sources:Type", 
                     asreml.obj = current.asr, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("Sources", "Type", "Species"))
diffs <- recalcLSD.alldiffs(diffs, meanLSD.type = "factor.combinations", LSDby = "Sources")
}}
\keyword{asreml}