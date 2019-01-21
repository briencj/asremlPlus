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
             values and their lower and upper confidence intervals are added 
             to a \code{data.frame} that is consistent with a 
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
             linear.transformation = NULL, Vmatrix = FALSE, 
             error.intervals = "Confidence", avsed.tolerance = 0.25, 
             meanLSD.type = "overall", LSDby = NULL, 
             response = NULL, response.title = NULL, 
             x.num = NULL, x.fac = NULL, 
             tables = "all", level.length = NA, 
             pairwise = TRUE, alpha = 0.05,
             transform.power = 1, offset = 0, scale = 1, 
             inestimable.rm = TRUE, ...)
}
\arguments{
\item{alldiffs.obj}{An \code{\link{alldiffs.object}}.}
\item{classify}{A \code{\link{character}} string giving the variables that 
                 define the margins of the multiway table to be predicted. 
                 Multiway tables are specified by forming an interaction type 
                 term from the classifying variables, that is, separating the 
                 variable names with the \code{:} operator.}
 \item{term}{A \code{\link{character}} string giving the variables that define the term 
               that was fitted using \code{asreml} and that corresponds 
               to \code{classify}. It only needs to be specified when 
               it is different to \code{classify}.}
 \item{linear.transformation}{A \code{\link{formula}} or a \code{\link{matrix}}. 
               If a \code{\link{formula}} is given then it is taken to be a submodel of 
               the model term corresponding to the \code{classify}. The projection matrix 
               that transforms the \code{predictions} so that they conform to the submodel 
               is obtained; the submodel should involving the variables in the 
               \code{classify}. For example, 
               for \code{classify} set to \code{"A:B"}, the submodel \code{~ A + B} will 
               result in the \code{predictions} for the combinations of 
               \code{A} and \code{B} being made additive for the \code{\link{factor}s} 
               \code{A} and \code{B}.
               
               If a \code{\link{matrix}} is provided then it will be 
               used to apply the linear transformation to the \code{predictions}. 
               The number of rows in the \code{\link{matrix}} should equal the 
               number of linear combinations of the \code{predictions} desired and 
               the number of columns should equal the number of \code{predictions}.
               
               In either case, as well as the values of the linear combinations, 
               their standard errors, pairwise differences and associated statistics 
               are returned.}
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
             If \code{meanLSD.type} is set to \code{overall}, the \code{avsed.tolerance} is not 
             \code{NA} and the range of the SEDs divided by the average of the SEDs exceeds 
             \code{avsed.tolerance} then the \code{error.intervals} calculations and the plotting 
             will revert to confidence intervals. Also, half LSDs cannot be used for 
             backtansformed values and so confidence intervals will be used instead.}
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
             \code{\link{alldiffs.object}} when \code{meanLSD.type} is \code{factor.combinatons}.}
 \item{response}{A \code{character} specifying the response variable for the 
                    predictions. It is stored as an attribute to the \code{\link{alldiffs.object}} .}
 \item{response.title}{A \code{character} specifying the title for the response variable 
                    for the predictions. It is stored as an attribute to the 
                    \code{\link{alldiffs.object}}.}
 \item{x.num}{A \code{\link{character}} string giving the name of the numeric covariate that 
              (i) corresponds to \code{x.fac}, (ii) is potentially included in 
              terms in the fitted model, and (iii) which corresponds to the 
              x-axis variable. It should have the same number of unique values 
              as the number of levels in \code{x.fac}.}
 \item{x.fac}{A \code{\link{character}} string giving the name of the factor that 
              (i) corresponds to \code{x.num}, (ii) is potentially included in 
               terms in the fitted model, and (iii) which corresponds to the 
              x-axis variable. It should have the same number of levels as the 
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
 \item{level.length}{The maximum number of characters from the the levels of 
            factors to use in the row and column labels of the tables of 
            pairwise differences and their p-values and standard errors.}
\item{pairwise}{A \code{\link{logical}} indicating whether all pairwise differences of the
             \code{predictions} and their standard errors and p-values are to be 
             computed and stored. If \code{tables} is equal to 
             \code{"differences"} or \code{"all"} or \code{error.intervals} is 
             equal to \code{"halfLeastSignificant"}, they will be stored 
             irrespective of the value of \code{pairwise}.}
 \item{alpha}{A \code{\link{numeric}} giving the significance level for LSDs or one minus 
              the confidence level for confidence intervals.}
 \item{transform.power}{A \code{\link{numeric}} specifying the power of a transformation, if 
              one has been applied to the response variable. Unless it is equal 
              to 1, the default, back-transforms of the predictions will be 
              obtained and presented in tables or graphs as appropriate. 
              The back-transformation raises the predictions to the power equal 
              to the reciprocal of \code{transform.power}, unless it equals 0 in 
              which case the exponential of the predictions is taken.} 
 \item{offset}{A \code{\link{numeric}} that has been added to each value of the 
               response after any scaling and before applying any power transformation.}
 \item{scale}{A \code{\link{numeric}} by which each value of the response has been multiplied 
               before adding any offset and applying any power transformation.}
 \item{inestimable.rm}{A \code{logical} indicating whether rows for  predictions 
               that are not estimable are to be removed from the components of 
               the \code{\link{alldiffs.object}}.}
 \item{\dots}{provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{A \code{\link{alldiffs.object}} with the linear transformation of the predictions 
       and their standard errors and all pairwise differences between the linear 
       transforms of their predictions, their standard errors and p-values 
       and LSD statistics.

           If \code{transform.power} or \code{scale} is not one or \code{offset} 
           is not zero, it will contain a \code{data.frame} with the backtransformed 
           linear transformation of the predictions. The backtransformation will, after 
           backtransforming for any power transformation, subtract the \code{offset} 
           and then divide by the \code{scale}.
           
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
           the \code{term}, the \code{classify}, \code{tdf}, \code{sortFactor} 
           and the \code{sortOrder} will be set as attributes to the object.
}
\details{For a matrix \eqn{\mathbf{L}}{L}, vector of 
             predictions \eqn{\mathbf{p}}{p} and variance matrix of the predictions 
             \eqn{\mathbf{V}_p}{V_p}, the linear transformed predictions are given 
             by \eqn{\mathbf{Lp}}{Lp} with variance matrix \eqn{\mathbf{LV}_p\mathbf{L}^\mathrm{T}}{LVL^T}. 
             The last matrix is used to compute the variance of pairwise 
             differences between the transformed values. 

             The \code{\link{matrix}} \eqn{\mathbf{L}}{L} is directly specified by setting 
             \code{linear.transformation} to it. If \code{linear.transformation} is a 
             \code{\link{formula}} then \eqn{\mathbf{L}}{L} is formed as the sum of the 
             orthogonal projection matrices obtained using \code{pstructure.formula} 
             from the package \code{dae}; \code{grandMean} is set to \code{TRUE} and 
             \code{orthogonalize} to \code{"eigenmethods"}.
}
\author{Chris Brien}
\seealso{\code{\link{predictPlus.asreml}}, \code{\link{as.alldiffs}}, \code{\link{print.alldiffs}}, 
         \code{\link{sort.alldiffs}},  
         \code{\link{subset.alldiffs}}, \cr 
         \code{\link{allDifferences.data.frame}}, 
         \code{\link{redoErrorIntervals.alldiffs}},  \code{\link{recalcLSD.alldiffs}},  \cr
         \code{\link{predictPresent.asreml}}, 
         \code{\link{plotPredictions.data.frame}}, \code{as.Date}, \code{predict.asreml}}
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