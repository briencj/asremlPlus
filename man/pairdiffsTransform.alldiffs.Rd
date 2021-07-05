\name{pairdiffsTransform.alldiffs}
\alias{pairdiffsTransform.alldiffs}
\alias{pairdiffsTransform}
\title{Calculates the differences between nominated pairs of predictions stored in 
       an \code{\link{alldiffs.object}}.}
\description{Predictions of differences and their error intervals are formed for two levels of 
             a factor, the \code{pairs.factor}. For each pair of a level of the 
             \code{pairs.factor} in \code{numerator.levels} with a level in 
             \code{denominator.levels}, an \code{\link{alldiffs.object}} is formed that 
             contains the differences between predictions with this pair of levels for all of 
             the combinations of the levels of the other factors in the \code{classify} of the 
             \code{\link{alldiffs.object}}. These prediction differences are obtained using 
             \code{\link{linTransform}} by forming a suitable contrast matrix to specfy 
             the \code{linear.transformation}. This function has the advantage that the 
             factors indexing the differences are included in the components of the 
             \code{\link{alldiffs.object}}s.

             If \code{pairwise = TRUE}, all pairwise differences between the 
             linear transforms of the \code{predictions}, their standard errors, 
             p-values and LSD statistics are computed as using 
             \code{\link{allDifferences.data.frame}}. 
             This adds them to the \code{\link{alldiffs.object}} as additional 
             \code{list} components named \code{differences}, \code{sed}, 
             \code{p.differences} and \code{LSD}.

             The printing of the components produced is controlled by the 
             \code{tables} argument. The order of plotting the levels of 
             one of the factors indexing the predictions can be modified 
             and is achieved using \code{\link{sort.alldiffs}}.
}
\usage{
\method{pairdiffsTransform}{alldiffs}(alldiffs.obj, pairs.factor, first.levels, second.levels, 
                   Vmatrix = FALSE, 
                   error.intervals = "Confidence", avsed.tolerance = 0.25, 
                   meanLSD.type = "overall", LSDby = NULL, 
                   response = NULL, response.title = NULL, tables = "all", 
                   pairwise = TRUE, alpha = 0.05, ...)
}
\arguments{
 \item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
 \item{pairs.factor}{A \code{\link{character}} string giving the name of the factor 
                 for whose levels the differences are to be calculated.}
 \item{first.levels}{A \code{\link{character}} string containing the levels of the \code{pairs.factor}
               whose predictions are those subtracted from.}
 \item{second.levels}{A \code{\link{character}} string containing the levels of the \code{pairs.factor}
               whose predictions are those that are subtracted.}
 \item{Vmatrix}{A \code{\link{logical}} indicating whether the variance matrix of the 
             \code{predictions} will be stored as a component of the \code{\link{alldiffs.object}} 
             that is returned.}
 \item{error.intervals}{A \code{\link{character}} string indicating the type of error interval, if any, 
             to calculate in order to indicate uncertainty in the results. 
             Possible values are \code{"none"}, \code{"StandardError"}, \code{"Confidence"} 
             and \code{"halfLeastSignificant"}. The default is for confidence limits to 
             be used. The \code{"halfLeastSignificant"} option results in half the 
             Least Significant Difference (LSD) being added and subtracted to the 
             predictions, the LSD being calculated using the square root of the mean of the 
             variances of all or a subset of pairwise differences between the predictions. 
             If the LSD is zero, as can happen when predictions are constrained to be equal,
             then the limits of the error intervals are set to \code{NA}. 
             If \code{meanLSD.type} is set to \code{overall}, the \code{avsed.tolerance} is not 
             \code{NA} and the range of the SEDs divided by the average of the SEDs exceeds 
             \code{avsed.tolerance} then the \code{error.intervals} calculations and the plotting 
             will revert to confidence intervals.}
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
             involving that prediction, unless there is only one prediction for a 
             \code{factor.combination}, when notional LSDs are calculated that are based on the 
             standard error of the prediction multiplied by the square root of two. It also 
             determines, in conjunction with \code{avsed.tolerance}, which LSD will be used
             in calculating \code{error.intervals} and, hence, is used for plots.}
\item{LSDby}{A \code{\link{character}} (vector) of variables names, being the names of the 
             \code{\link{factor}s} or \code{\link{numeric}s} in the \code{classify} for each 
             combination of which a mean LSD, minLSD and max LSD is stored in the \code{LSD} 
             component of the \code{\link{alldiffs.object}} when \code{meanLSD.type} is 
             \code{factor.combinatons}.}
 \item{response}{A \code{character} specifying the response variable for the 
                    predictions. It is stored as an attribute to the \code{\link{alldiffs.object}} .}
 \item{response.title}{A \code{character} specifying the title for the response variable 
                    for the predictions. It is stored as an attribute to the 
                    \code{\link{alldiffs.object}}.}
 \item{tables}{A \code{\link{character}} vector containing a combination of 
              \code{none}, \code{predictions}, \code{vcov}, \code{backtransforms}, \code{differences}, 
              \code{p.differences}, \code{sed}, \code{LSD} and \code{all}.  
              These nominate which components of the \code{\link{alldiffs.object}} to print.
              Only \code{none}, \code{predictions} and \code{all} apply to \code{method} 
              \code{Fieller}.}
\item{pairwise}{A \code{\link{logical}} indicating whether all pairwise differences of the
             \code{predictions} and their standard errors and p-values are to be 
             computed and stored. If \code{tables} is equal to 
             \code{"differences"} or \code{"all"} or \code{error.intervals} is 
             equal to \code{"halfLeastSignificant"}, they will be stored 
             irrespective of the value of \code{pairwise}.}
 \item{alpha}{A \code{\link{numeric}} giving the significance level for LSDs or one minus 
              the confidence level for confidence intervals.}
 \item{\dots}{further arguments passed to \code{\link{linTransform.alldiffs}}.}
}
\value{A \code{list} of \code{\link{alldiffs.object}}s with a component for each combination 
       of a \code{first.levels} with a \code{second.levels}. The name of a component will be 
       a level from \code{first.levels} combined with a level from \code{second.levels}, 
       separated by a comma. If the predictions in the supplied \code{\link{alldiffs.object}} 
       are based on a \code{response} that was transformed, each \code{\link{alldiffs.object}} 
       in the \code{list} will include a \code{backtransforms} component that contains   
       a column labelled \code{backtransformed.predictions}, along with the backtransforms of 
       the nominated \code{error.intervals}. If the transformation was the logarithmic 
       transformation, these \code{backtransformed.predictions} are predicted ratios of the 
       untransformed \code{response}. 
       
       If \code{sortFactor} attribute is set and is not the 
       \code{ratio.factor}, the predictions and, if present, their backtransforms will be sorted using 
       the \code{sortOrder} attribute of the  \code{\link{alldiffs.object}},
       and both \code{sortFactor} and \code{sortOrder} will be set as attributes to the object.
}
\author{Chris Brien}

\seealso{\code{\link{linTransform}}, \code{\link{ratioTransform}}, \code{\link{predictPlus.asreml}}, 
         \code{\link{as.alldiffs}}, \code{\link{print.alldiffs}}, \cr 
         \code{\link{sort.alldiffs}},  \code{\link{subset.alldiffs}}, 
         \code{\link{allDifferences.data.frame}}, \cr
         \code{\link{redoErrorIntervals.alldiffs}},  
         \code{\link{recalcLSD.alldiffs}},  \code{\link{predictPresent.asreml}}, \cr
         \code{\link{plotPredictions.data.frame}}, \cr
         \code{as.Date}, \code{predict.asreml}}
\examples{
#### Form the differences for log(RGR) for Salinity
load(system.file("extdata", "testDiffs.rda", package = "asremlPlus", mustWork = TRUE))
#### For the ratios for Cl per WU Temperature - use backtransforms of log-predictions
Preds.ratio.ClUp <- pairdiffsTransform(diffs.ClUp, method = "log",
                                       pairs.factor = "Temperature", 
                                       first.levels = "Hot",
                                       second.levels = "Cool",
                                       error.intervals = "halfLeast",
                                       tables = "backtransforms") #BAcktransforms are ratios

#### Form the differences for Nitrogen compared to no Nitrogen                                  
data("Oats.dat")
\dontrun{
m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                 random=~Blocks/Wplots,
                 data=Oats.dat)
current.asrt <- as.asrtests(m1.asr)
wald.tab <-  current.asrt$wald.tab
Var.diffs <- predictPlus(m1.asr, classify="Nitrogen:Variety", pairwise = TRUE,
                         Vmatrix = TRUE, error.intervals = "halfLeast",
                         meanLSD.type = "factor", LSDby = "Variety",
                         wald.tab = wald.tab)
}

 ## Use lme4 and emmmeans to get predictions and associated statistics
if (requireNamespace("lmerTest", quietly = TRUE) & 
    requireNamespace("emmeans", quietly = TRUE))
{
  m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
  ## Set up a wald.tab
  int <- as.data.frame(rbind(rep(NA,4)))
  rownames(int) <- "(Intercept)"
  wald.tab <- anova(m1.lmer, ddf = "Kenward", type = 1)[,3:6]
  names(wald.tab) <- names(int) <- c("Df", "denDF", "F.inc", "Pr")
  wald.tab <- rbind(int, wald.tab)
  #Get predictions
  Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
  Var.preds <- summary(Var.emm)
  ## Modify Var.preds to be compatible with a predictions.frame
  Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                    se = "SE", interval.type = "CI", 
                                    interval.names = c("lower.CL", "upper.CL"))
  Var.vcov <- vcov(Var.emm)
  Var.sed <- NULL
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  
  #Create alldiffs object
  Var.diffs <- as.alldiffs(predictions = Var.preds, 
                           sed = Var.sed, vcov = Var.vcov, 
                           classify = "Nitrogen:Variety", response = "Yield", tdf = den.df)
} 

if (exists("Var.diffs"))
  Preds.diffs.OatsN <- pairdiffsTransform(alldiffs.obj = Var.diffs,
                                          pairs.factor = "Nitrogen", 
                                          first.levels = c("0.2","0.4","0.6"),
                                          second.levels = "0", error.intervals = "halfLeast",
                                          tables = "none")
}  
 
\keyword{asreml}