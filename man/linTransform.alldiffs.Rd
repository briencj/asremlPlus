\name{linTransform.alldiffs}
\alias{linTransform.alldiffs}
\alias{linTransform}
\title{Calculates a linear transformation of the predictions stored in an 
       \code{\link{alldiffs.object}}.}
\description{Effects the linear transformation of the predictions in the 
             supplied \code{\link{alldiffs.object}}, the transformation being specified 
             by a \code{\link{matrix}} or a \code{\link{formula}}. The values of 
             the transformed values are stored in an \code{\link{alldiffs.object}}. 
             A \code{\link{matrix}} might be a contrast \code{\link{matrix}} or 
             a \code{\link{matrix}} of weights for the levels of a 
             \code{\link{factor}} used to obtain the weighted average over  
             the levels of that \code{\link{factor}}. A \code{\link{formula}} gives 
             rise to a projection \code{\link{matrix}} that linearly transforms 
             the predictions so that they conform to the model specified by the 
             \code{\link{formula}}, this model being a submodel of that inherent 
             in the \code{classify}. 

             If \code{pairwise = TRUE}, all pairwise differences between the 
             linear transforms of the \code{predictions}, their standard errors, 
             p-values and LSD statistics are computed as using 
             \code{\link{allDifferences.data.frame}}. 
             This adds them to the \code{\link{alldiffs.object}} as additional 
             \code{list} components named \code{differences}, \code{sed}, 
             \code{p.differences} and \code{LSD}.

             If a transformation has been applied (any one of 
             \code{transform.power} is not one, \code{scale} is  not one and 
             \code{offset} is nonzero), the backtransforms of the transformed 
             values and of the lower and upper limits of their \code{error.intervals} 
             are added to a \code{data.frame} that is consistent with a 
             \code{\link{predictions.frame}}.  If \code{transform.power} is other than 
             one, the \code{standard.error} column of the \code{data.frame} 
             is set to \code{NA}. This \code{data.frame} is added to the 
             \code{\link{alldiffs.object}} as a \code{list} component called 
             \code{backtransforms}.

             The printing of the components produced is controlled by the 
             \code{tables} argument. The order of plotting the levels of 
             one of the factors indexing the predictions can be modified 
             and is achieved using \code{\link{sort.alldiffs}}.
}
\usage{
\method{linTransform}{alldiffs}(alldiffs.obj, classify = NULL, term = NULL, 
             linear.transformation = NULL, EGLS.linTransform = TRUE, 
             Vmatrix = FALSE, error.intervals = "Confidence", 
             avsed.tolerance = 0.25, accuracy.threshold = NA, 
             LSDtype = "overall", LSDsupplied = NULL, 
             LSDby = NULL, LSDstatistic = "mean", 
             LSDaccuracy = "maxAbsDeviation", 
             zero.tolerance = .Machine$double.eps ^ 0.5, 
             response = NULL, response.title = NULL, 
             x.num = NULL, x.fac = NULL, 
             tables = "all", level.length = NA, 
             pairwise = TRUE, alpha = 0.05,
             inestimable.rm = TRUE, ...)
}
\arguments{
\item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
\item{classify}{A \code{\link{character}} string giving the variables that 
                 define the margins of the multiway table corresponding to the 
                 \code{predictions} in \code{alldiffs.obj}. 
                 Multiway tables are specified by forming an interaction type 
                 term from the classifying variables, that is, separating the 
                 variable names with the \code{:} operator.}
 \item{term}{A \code{\link{character}} string giving the variables that define the term 
               that was fitted using \code{asreml} and that corresponds to \code{classify}. 
               It only needs to be specified when it is different to \code{classify}; it 
               is stored as an attribute of the \code{\link{alldiffs.object}}. 
               It is likely to be needed when the fitted model includes terms that involve  
               both a \code{\link{numeric}} covariate and a \code{\link{factor}} that 
               parallel each other; the \code{classify} would include the covariate and 
               the \code{term} would include the \code{factor}.}
 \item{linear.transformation}{A \code{\link{formula}} or a \code{\link{matrix}}. 
               If a \code{\link{formula}} is given then it is taken to be a submodel of 
               a model term corresponding to the \code{classify}. The projection matrix 
               that transforms the \code{predictions} so that they conform to the submodel 
               is obtained; the submodel does not have to involve variables in the 
               \code{classify}, but the variables must be columns in the \code{predictions} 
               component of \code{alldiffs.obj} and the space for the submodel must be a 
               subspace of the space for the term specified by the \code{classify}. 
               For example, for \code{classify} set to \code{"A:B"}, the submodel 
               \code{~ A + B} will result in the \code{predictions} for the combinations of 
               \code{A} and \code{B} being made additive for the \code{\link{factor}s} 
               \code{A} and \code{B}. The submodel space corresponding to \code{A + B} is 
               a subspace of the space \code{A:B}. In this case both the submodel and the 
               classify involve only the factors A and B. To fit an intercept-only 
               submodel, specify \code{linear.transformation} to be the formula \code{~1}.
               
               If a \code{\link{matrix}} is provided then it will be 
               used to apply the linear transformation to the \code{predictions}. 
               The number of rows in the \code{\link{matrix}} should equal the 
               number of linear combinations of the \code{predictions} desired and 
               the number of columns should equal the number of \code{predictions}.
               
               In either case, as well as the values of the linear combinations, 
               their standard errors, pairwise differences and associated statistics 
               are returned.}
 \item{EGLS.linTransform}{A \code{\link{logical}} indicating whether or not the 
              \code{linear.transformation} of the predictions stored in an 
              \code{\link{alldiffs.object}} by fitting a submodel supplied in a 
              \code{\link{formula}} is to take into account the variance of the 
              predictions using a Estimated Generalized Least Squares (EGLS) approach. 
              This is likely to be appropriate when the variance matrix of the predictions 
              is not compound symmetric i.e. when not all the variances are equal or not 
              all the covariances are equal. If the variance matrix is compund symmetric, 
              then the setting of \code{EGLS.linTransform} will not affect the transformed 
              predictions.} 
\item{Vmatrix}{A \code{\link{logical}} indicating whether the variance matrix of the 
             \code{predictions} will be stored as a component of the \code{\link{alldiffs.object}} 
             that is returned. If \code{linear.transformation} is set, it will 
             be stored irrespective of the value of \code{Vmatrix}.}
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
             If \code{LSDtype} is set to \code{overall}, the \code{avsed.tolerance} is not 
             \code{NA} and the range of the SEDs divided by the average of the SEDs exceeds 
             \code{avsed.tolerance} then the \code{error.intervals} calculations and the plotting 
             will revert to confidence intervals.}
\item{avsed.tolerance}{A \code{\link{numeric}} giving the value of the SED range, the range of the SEDs 
             divided by the square root of the mean of the variances of all or a subset of the 
             pairwise differences, that is considered reasonable in calculating 
             \code{error.intervals}. To have it ignored, set it to \code{NA}. It should be a 
             value between 0 and 1. The following rules apply:
             \enumerate{
             \item If \code{avsed.tolerance} is \code{NA} then mean LSDs of the type specified by 
             \code{LSDtype} are calculated and used in \code{error.intervals} and plots.
             \item Irrespective of the setting of \code{LSDtype}, if \code{avsed.tolerance} is not 
             exceeded then the mean LSDs are used in \code{error.intervals} and plots. 
             \item If \code{LSDtype} is set to \code{overall}, \code{avsed.tolerance} is not 
             \code{NA}, and \code{avsed.tolerance} is exceeded then  \code{error.intervals} and 
             plotting revert to confidence intervals. 
             \item If \code{LSDtype} is set to \code{factor.combinations} and \code{avsed.tolerance} 
             is not exceeded for any factor combination then the half LSDs are 
             used in \code{error.intervals} and plots; otherwise, \code{error.intervals} and 
             plotting revert to confidence intervals. 
             \item If \code{LSDtype} is set to \code{per.prediction} and \code{avsed.tolerance} 
             is not exceeded for any prediction then the half LSDs are used in \code{error.intervals} 
             and plots; otherwise, \code{error.intervals} and plotting revert to confidence intervals.
             }}
\item{accuracy.threshold}{A \code{\link{numeric}} specifying the value of the LSD accuracy measure,  
             which measure is specified by \code{LSDaccuracy}, as a threshold value in determining whether the 
             \code{hallfLeastSignificant} \code{error.interval} for a predicted value is a reasonable  
             approximation; this will be the case if the LSDs across all pairwise comparisons for which 
             the interval's LSD was computed, as specified by \code{LSDtype} and \code{LSDby}, 
             are similar enough to the interval's LSD, as measured by \code{LSDaccuracy}.  
             If it is \code{NA}, it will be ignored. If it is 
             not \code{NA}, a column of \code{\link{logical}s} named \code{LSDwarning}  will be added 
             to the \code{predictions} component of the \code{\link{alldiffs.object}}. The value of 
             \code{LSDwarning} for a \code{predicted.value} will be \code{TRUE} if the value of the 
             \code{LSDaccuracy} measure computed from the LSDs for differences between this 
             \code{predicted.value} and the other \code{predicted.values} as compared to its 
             \code{assignedLSD} exceeds the value of \code{accuracy.threshold}. Otherwise, the 
             value of \code{LSDwarning} for a \code{predicted.value} will be \code{FALSE}.}
\item{LSDtype}{A \code{\link{character}} string that can be \code{overall}, \code{factor.combinations}, 
              \code{per.prediction} or \code{supplied}. It determines whether the values stored in a row  
              of a \code{\link{LSD.frame}} are the values calculated 
             (i) \code{overall} from the LSD values for all pairwise comparisons, unless there is only one, 
              possibly repeated, prediction, when a notional LSD is calculated, 
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
\item{zero.tolerance}{A \code{\link{numeric}} specifying the value such that if a \code{predicted.value}, 
             its variance-covariance, or an \code{LSD} is less than it, the LSD will be 
             considered to be zero.}
\item{response}{A \code{character} specifying the response variable for the 
                    predictions. It is stored as an attribute to the \code{\link{alldiffs.object}} .}
\item{response.title}{A \code{character} specifying the title for the response variable 
                    for the predictions. It is stored as an attribute to the 
                    \code{\link{alldiffs.object}}.}
\item{x.num}{A \code{\link{character}} string giving the name of the numeric covariate that 
              (i) is potentially included in terms in the fitted model and (ii) is the 
              x-axis variable for plots. Its values will not be converted to a \code{\link{factor}}.}
 \item{x.fac}{A \code{\link{character}} string giving the name of the factor that 
              (i) corresponds to \code{x.num} and (ii) is potentially included in 
               terms in the fitted model. It should have the same number of levels as the 
              number of unique values in \code{x.num}. The levels of 
              \code{x.fac} must be in the order in which they are to be plotted 
              - if they are dates, then they should be in the form 
              yyyymmdd, which can be achieved using \code{as.Date}. However, the levels
              can be non-numeric in nature, provided that \code{x.num} is also set.}
\item{tables}{A \code{\link{character}} vector containing a combination of 
              \code{none}, 
              \code{predictions}, \code{vcov}, \code{backtransforms}, \code{differences}, 
              \code{p.differences}, \code{sed}, \code{LSD} and \code{all}.
              These nominate which components of the \code{\link{alldiffs.object}} to print.}
\item{level.length}{The maximum number of characters from the levels of 
            factors to use in the row and column labels of the tables of 
            pairwise differences and their p-values and standard errors.}
\item{pairwise}{A \code{\link{logical}} indicating whether all pairwise differences of the
             \code{predictions} and their standard errors and p-values are to be 
             computed and stored. If \code{tables} is equal to 
             \code{"differences"} or \code{"all"} or \code{error.intervals} is 
             equal to \code{"halfLeastSignificant"}, they will be stored 
             irrespective of the value of \code{pairwise}.}
\item{alpha}{A \code{\link{numeric}} giving the significance level for LSDs or one minus 
              the confidence level for confidence intervals.
              It is stored as an attribute to the \code{\link{alldiffs.object}}.}
\item{inestimable.rm}{A \code{logical} indicating whether rows for  predictions 
               that are not estimable are to be removed from the components of 
               the \code{\link{alldiffs.object}}.}
\item{\dots}{further arguments passed to \code{\link{redoErrorIntervals.alldiffs}}.}
}
\value{A \code{\link{alldiffs.object}} with the linear transformation of the predictions 
       and their standard errors and all pairwise differences between the linear 
       transforms of their predictions, their standard errors and p-values 
       and LSD statistics.

       If the supplied \code{\link{alldiffs.object}} contained a \code{backtransforms} 
       component, then the returned \code{\link{alldiffs.object}} will contain 
       a \code{backtransforms} component with the backtransformed linear transformation 
       of the predictions. The backtransformation will, after backtransforming for any 
       power transformation, subtract the \code{offset} and then divide by the \code{scale}.
           
           If \code{error.intervals} is not \code{"none"}, then the 
           \code{predictions} component and, if present, the 
           \code{backtransforms} component will contain columns for the lower 
           and upper values of the limits for the interval. The names of these 
           columns will consist of three parts separated by full stops: 
           1) the first part will be \code{lower} or \code{upper}; 
           2) the second part will be one of \code{Confidence}, 
              \code{StandardError} or \code{halfLeastSignificant}; 
           3) the third component will be \code{limits}.

           The name of the \code{response}, the \code{response.title}, 
           the \code{term}, the \code{classify}, \code{tdf}, \code{alpha}, \code{sortFactor} 
           and the \code{sortOrder} will be set as attributes to the object.
           Also, if \code{error.intervals} is \code{"halfLeastSignificant"}, then those of 
           \code{LSDtype}, \code{LSDby} and \code{LSDstatistic} that are not \code{NULL} 
           will be added as attributes of the object and of the \code{predictions frame}; 
           additionally, \code{LSDvalues} will be added as attribute of the 
           \code{predictions frame}, \code{LSDvalues} being the LSD values used in 
           calculating the \code{error.intervals}. 


}
\details{For a matrix \eqn{\mathbf{L}}{L}, vector of 
             predictions \eqn{\mathbf{p}}{p} and variance matrix of the predictions 
             \eqn{\mathbf{V}_p}{V_p}, the linear transformed predictions are given 
             by \eqn{\mathbf{Lp}}{Lp} with variance matrix \eqn{\mathbf{LV}_p\mathbf{L}^\mathrm{T}}{LVL^T}. 
             The last matrix is used to compute the variance of pairwise 
             differences between the transformed values. 

             If \code{linear.transformation} is a \code{\link{matrix}}, \eqn{\mathbf{M}}{M} say, then the 
             linear-transformation \code{\link{matrix}}, \eqn{\mathbf{L}}{L}, is just the supplied 
             \code{\link{matrix}} \eqn{\mathbf{M}}{M}. 
             
             If \code{linear.transformation} is a \code{\link{formula}} and \code{EGLS.linTransform} 
             is \code{TRUE}, then a matrix \eqn{\mathbf{M}}{M} is obtained that is the design matrix 
             for all of the terms in the \code{\link{formula}}. Using \eqn{\mathbf{M}}{M}, the 
             linear-transformation \code{\link{matrix}}, \eqn{\mathbf{L}}{L}, is formed as 
             \eqn{\mathbf{M} (\mathbf{M}^\top \widehat{\mathbf{V}}^- \mathbf{M})^- (\mathbf{M}^\top \widehat{\mathbf{V}}^-)}{M [t(M) V^-  M]^- t(M) V^-}.

             On the other hand, for \code{linear.transformation} a \code{\link{formula}} 
             and \code{EGLS.linTransform} set to \code{FALSE}, \eqn{\mathbf{L}}{L} is 
             formed as the sum of the orthogonal projection matrices obtained using 
             \code{pstructure.formula} from the package \code{dae}; \code{grandMean} 
             is set to \code{TRUE} and \code{orthogonalize} to \code{"eigenmethods"}.
}
\author{Chris Brien}
\seealso{\code{\link{linTransform}}, \code{\link{predictPlus.asreml}}, \code{\link{as.alldiffs}}, 
         \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}},  \cr
         \code{\link{subset.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{redoErrorIntervals.alldiffs}},  \cr 
         \code{\link{recalcLSD.alldiffs}}, \code{\link{pickLSDstatistics.alldiffs}},  
         \code{\link{predictPresent.asreml}}, \cr
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
current.asrt <- as.asrtests(current.asr, NULL, NULL)
#Get additive predictions directly using predictPlus
diffs.sub <- predictPlus.asreml(classify = "Sources:Species", Vmatrix = TRUE, 
                                linear.transformation = ~ Sources + Species,
                                asreml.obj = current.asr, tables = "none", 
                                wald.tab = current.asrt$wald.tab, 
                                present = c("Type","Species","Sources"))
}

## Use lmeTest and emmmeans to get predictions and associated statistics

if (requireNamespace("lmerTest", quietly = TRUE) & 
    requireNamespace("emmeans", quietly = TRUE))
{
  m1.lmer <- lmerTest::lmer(pH ~ Benches + (Sources * Species) + 
                              (1|Benches:MainPlots),
                            data=na.omit(WaterRunoff.dat))
  SS.emm <- emmeans::emmeans(m1.lmer, specs = ~ Sources:Species)
  SS.preds <- summary(SS.emm)
  den.df <- min(SS.preds$df, na.rm = TRUE)
  ## Modify SS.preds to be compatible with a predictions.frame
  SS.preds <- as.predictions.frame(SS.preds, predictions = "emmean", 
                                   se = "SE", interval.type = "CI", 
                                   interval.names = c("lower.CL", "upper.CL"))
  
  ## Form an all.diffs object and check its validity
  SS.vcov <- vcov(SS.emm)
  SS.diffs <- allDifferences(predictions = SS.preds, classify = "Sources:Species", 
                             vcov = SS.vcov, tdf = den.df)
  validAlldiffs(SS.diffs)

  #Get additive predictions
  diffs.sub <- linTransform(SS.diffs, classify = "Sources:Species", 
                            linear.transformation = ~ Sources + Species,
                            Vmatrix = TRUE, tables = "none")
}  
 
##Calculate contrasts from prediction obtained using asreml or lmerTest 
if (exists("diffs.sub"))
{ 
  #Contrast matrix for differences between each species and non-planted for the last source
  L <- cbind(matrix(rep(0,7*32), nrow = 7, ncol = 32),
             diag(1, nrow = 7), 
             matrix(rep(-1, 7), ncol = 1))
  rownames(L) <- as.character(diffs.sub$predictions$Species[33:39])
  diffs.L <- linTransform(diffs.sub, 
                          classify = "Sources:Species",
                          linear.transformation = L,
                          tables = "predictions")
}
}
\keyword{asreml}