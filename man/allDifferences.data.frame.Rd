\name{allDifferences.data.frame}
\alias{allDifferences.data.frame}
\alias{allDifferences}
\title{Using supplied predictions and standard errors of pairwise differences or the 
       variance matrix of predictions, forms all pairwise differences between the 
       set of predictions, and p-values for the differences.}
\description{Uses supplied predictions and standard errors of pairwise differences, 
       or the variance matrix of predictions to form, in an 
       \code{\link{alldiffs.object}}, for those components not already present, 
       (i) a table of all pairwise differences of the predictions, 
       (ii) the p-value of each pairwise difference, and 
       (iii) the minimum, mean, maximum and accuracy of LSD values. 
       Predictions that are aliased (or nonestimable) are removed from the 
       \code{predictions} component of the \code{\link{alldiffs.object}} and 
       standard errors of differences involving them are removed from the \code{sed} 
       component. 

       If necessary, the order of the columns of the variables in the \code{predictions} 
       component are changed to be the initial columns of the \code{\link{predictions.frame}} 
       and to match their order in the \code{classify}.  Also, the rows of \code{predictions} 
       component are ordered so that they are in standard order for the variables in the 
       \code{classify}. That is, the values of the last variable change with every row, 
       those of the second-last variable only change after all the values of the last 
       variable have been traversed; in general, the values of a variable are the same for 
       all the combinations of the values to the variables to its right in the 
       \code{classify}. The \code{sortFactor} or \code{sortOrder} arguments can be used 
       to order of the values for the \code{classify} variables, which is achieved using 
       \code{\link{sort.alldiffs}}. 

       Each p-value is computed as the probability of a t-statistic as large as or larger 
       than the absolute value of the observed difference divided by its standard error. The 
       p-values are stored in the \code{p.differences} component. The degrees of freedom of 
       the t-distribution is the degrees of freedom stored in the \code{tdf} attribute of 
       the \code{\link{alldiffs.object}}. This t-distibrution is also used in calculating 
       the LSD statistics stored in the \code{LSD} component of the \code{\link{alldiffs.object}}.
}
\usage{
\method{allDifferences}{data.frame}(predictions, classify, vcov = NULL, 
               differences = NULL, p.differences = NULL, sed = NULL, 
               LSD = NULL, LSDtype = "overall", LSDsupplied = NULL, 
               LSDby = NULL, LSDstatistic = "mean", 
               LSDaccuracy = "maxAbsDeviation", 
               retain.zeroLSDs = FALSE, 
               zero.tolerance = .Machine$double.eps ^ 0.5, 
               backtransforms = NULL, 
               response = NULL, response.title = NULL, 
               term = NULL, tdf = NULL,  
               x.num = NULL, x.fac = NULL,  
               level.length = NA, 
               pairwise = TRUE, alpha = 0.05,
               transform.power = 1, offset = 0, scale = 1, 
               transform.function = "identity", 
               inestimable.rm = TRUE,
               sortFactor = NULL, sortParallelToCombo = NULL, 
               sortNestingFactor = NULL, sortOrder = NULL, 
               decreasing = FALSE, ...)}
\arguments{
 \item{predictions}{A \code{\link{predictions.frame}}, or a \code{data.frame}, beginning 
                    with the variables classifying the predictions and also containing columns 
                    named \code{predicted.value}, \code{standard.error} and \code{est.status}; 
                    each row contains a single predicted value. It may also contain columns 
                    for the lower and upper limits of error intervals for the predictions.
                    Note that the names \code{standard.error} and 
                    \code{est.status} have been changed to \code{std.error} and \code{status} 
                    in the \code{pvals} component produced by \code{asreml-R4}; if the new names 
                    are in the \code{data.frame} supplied to \code{predictions}, they will be 
                    returned to the previous names. }
 \item{classify}{A \code{character} string giving the variables that define the margins 
                    of the multiway table that has been predicted. Multiway tables are 
                    specified by forming an interaction type term from the 
                    classifying variables, that is, separating the variable names 
                    with the \code{:} operator.}
\item{vcov}{A \code{matrix} containing the variance matrix of the predictions; it is used in 
                    computing the variance of linear transformations of the predictions.}
\item{differences}{A \code{matrix} containing all pairwise differences between
                    the predictions; it should have the same number of rows and columns as there are
                    rows in \code{predictions}.}
 \item{p.differences}{A \code{matrix} containing p-values for all pairwise differences 
                    between the predictions; each p-value is computed as the probability of a t-statistic 
                    as large as or larger than the observed difference divided by its standard error. 
                    The degrees of freedom of the t distribution for computing it are computed as 
                    the denominator degrees of freedom of the F value for the fixed term, if available; 
                    otherwise, the degrees of freedom stored in the attribute \code{tdf} are used; 
                    the matrix should be of the same size as that for \code{differences}.}
 \item{sed}{A \code{matrix} containing the standard errors of all pairwise differences 
                    between the predictions; they are used in computing the p-values.}
\item{LSD}{An \code{\link{LSD.frame}} containing the mean, minimum and maximum LSD for determining             
                    the significance of pairwise differences, as well as an assigned LSD and a measure 
                    of the accuracy of the LSD. If \code{LSD} is \code{NULL} then the \code{\link{LSD.frame}} 
                     stored in the \code{LSD} component will be calculated and 
                    the values of \code{LSDtype}, \code{LSDby} and \code{LSDstatistic} added as attributes 
                    of the \code{\link{alldiffs.object}}. The LSD for a single prediction 
                    assumes that any predictions to be compared are independent; this is not the case if 
                    residual errors are correlated.}
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
\item{LSDstatistic}{A \code{\link{character}} nominating one or more of \code{minmum}, \code{q10}, \code{q25},  
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
\item{retain.zeroLSDs}{A \code{\link{logical}} indicating whether to retain or omit LSDs that are zero when 
             calculating the summaries of LSDs.}
\item{zero.tolerance}{A \code{\link{numeric}} specifying the value such that if an LSD is less than it, it will be 
            considered to be zero.}
\item{backtransforms}{A \code{data.frame} containing the backtransformed values of the predicted 
                    values  that is consistent with the \code{predictions} component, except 
                    that the column named \code{predicted.value} is replaced by one called 
                    \code{backtransformed.predictions}. Any \code{error.interval} values will also 
                    be the backtransformed values. Each row contains a single predicted value.}
\item{response}{A \code{character} specifying the response variable for the 
                    predictions. It is stored as an attribute to the \code{\link{alldiffs.object}}.}
\item{response.title}{A \code{character} specifying the title for the response variable 
                    for the predictions. It is stored as an attribute to the \code{\link{alldiffs.object}}.}
 \item{term}{A \code{\link{character}} string giving the variables that define the term 
               that was fitted using \code{asreml} and that corresponds to \code{classify}. 
               It only needs to be specified when it is different to \code{classify}; it 
               is stored as an attribute of the \code{\link{alldiffs.object}}. 
               It is likely to be needed when the fitted model includes terms that involve  
               both a \code{\link{numeric}} covariate and a \code{\link{factor}} that 
               parallel each other; the \code{classify} would include the covariate and 
               the \code{term} would include the \code{factor}.}
\item{tdf}{an \code{integer} specifying the degrees of freedom of the standard error. It is used as 
                the degrees of freedom for the t-distribution on which p-values and confidence 
                intervals are based. 
                It is stored as an attribute to the \code{\link{alldiffs.object}}.}
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
\item{level.length}{The maximum number of characters from the levels of 
            factors to use in the row and column labels of the tables of 
            pairwise differences and their p-values and standard errors.}
\item{pairwise}{A logical indicating whether all pairwise differences of the
                \code{predictions} and their standard errors and p-values are to be 
                computed and stored. If \code{FALSE}, the components \code{differences} 
                and \code{p.differences} will be \code{NULL} in the returned 
                \code{\link{alldiffs.object}}.}
\item{alpha}{A \code{\link{numeric}} giving the significance level for LSDs or one minus 
              the confidence level for confidence intervals.
              It is stored as an attribute to the \code{\link{alldiffs.object}}.}
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
 \item{transform.function}{A \code{\link{character}} giving the name of a function that 
              specifies the scale on which the predicted values are defined. This may be the 
              result of a transformation of the data using the function or the use of the 
              function as a link function in the fitting of a generalized linear (mixed) 
              model (GL(M)M). The possible \code{transform.function}s are 
              \code{identity}, \code{log}, \code{inverse}, \code{sqrt}, \code{logit}, 
              \code{probit}, and \code{cloglog}.  The \code{predicted.values} and 
              \code{error.intervals}, if not \code{StandardError} intervals, will be 
              back-transformed using the inverse fucntion of the \code{transform.function}. 
              The \code{standard.error} column will be set to \code{NA}, unless (i) 
              \code{asreml} returns columns named \code{transformed.value} and 
              \code{approx.se}, as well as those called \code{predicted.values} and 
              \code{standard.error} (such as when a GLM is fitted) and 
              (ii) the values in \code{transformed.value} are equal to those obtained by 
              backtransforming the \code{predicted.value}s using the inverse function 
              of the \code{transform.function}. Then, the \code{approx.se} values will be 
              saved in the \code{standard.error} column of the \code{backtransforms} 
              component of the returned \code{alldiffs.obj}. Also, the 
              \code{transformed.value} and \code{approx.se} columns are removed from both 
              the \code{predictions} and \code{backtransforms} components of the 
              \code{alldiffs.obj}. Note that the values that end up in the \code{standard errors} 
              column are approximate for the backtransformed values and are not used in 
              calculating \code{error.intervals}.}
\item{inestimable.rm}{A \code{logical} indicating whether rows for  predictions that 
                       are not estimable are to be removed from the components of the
                       \code{\link{alldiffs.object}}.}
\item{sortFactor}{A \code{\link{character}} containing the name of the 
              \code{factor} that indexes the set of predicted values that determines 
               the sorting of the components. If there is only one variable in the 
              \code{classify} term then \code{sortFactor} can be \code{NULL} and 
              the order is defined by the complete set of predicted values. 
              If there is more than one variable in the \code{classify} term 
              then \code{sortFactor} must be set. In this case the \code{sortFactor} 
              is sorted in the same order within each combination of the values of 
              the \code{sortParallelToCombo} variables: the \code{classify} variables, excluding the 
              \code{sortFactor}. There should be only one predicted value for 
              each unique value of \code{sortFactor} within each set defined by a 
              combination of the values of the \code{classify} variables, excluding the 
              \code{sortFactor} \code{factor}. 
              The order to use is determined by either \code{sortParallelToCombo} or 
              \code{sortOrder}.}
\item{sortParallelToCombo}{A \code{\link{list}} that specifies a combination of the values 
              of the \code{factor}s and \code{numeric}s, excluding \code{sortFactor}, that 
              are in \code{classify}. Each of the components of the supplied \code{\link{list}} 
              is named for a \code{classify} variable and specifies a single value for it. The 
              combination of this set of values will be used to define a subset of the predicted 
              values whose order will define the order of \code{sortFactor}. Each of the other 
              combinations of the values of the \code{factor}s and \code{numeric}s will be sorted 
              in parallel. If \code{sortParallelToCombo} is \code{NULL} then the first value of               
              each \code{classify} variable, except for the \code{sortFactor} \code{factor},  
              in the \code{predictions} component is used to define \code{sortParallelToCombo}. 
              If there is only one variable in the \code{classify} then 
              \code{sortParallelToCombo} is ignored.}
\item{sortNestingFactor}{A \code{\link{character}} containing the name of the 
              \code{factor} that defines groups of the \code{sortFactor} within which the predicted 
              values are to be ordered. 
              If there is only one variable in the \code{classify} then 
              \code{sortNestingFactor} is ignored.}
\item{sortOrder}{A \code{character vector} whose length is the same as the number 
               of levels for \code{sortFactor} in the \code{predictions} component of the 
              \code{\link{alldiffs.object}}. It specifies the desired order of the 
              levels in the reordered components of the \code{\link{alldiffs.object}}. 
              The argument \code{sortParallelToCombo} is ignored.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}:  \cr
              \code{levs <- levels(f)[order(x)]}.}
\item{decreasing}{A \code{logical} passed to \code{order} that detemines whether 
              the order for sorting the components of the \code{\link{alldiffs.object}} 
              is for increasing or decreasing magnitude of the predicted values.}
\item{\dots}{provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{An \code{\link{alldiffs.object}} with components 
       \code{predictions}, \code{vcov}, \code{differences}, \code{p.differences} 
       \code{sed}, and \code{LSD}. 

       The name of the \code{response}, the \code{response.title}, 
       the \code{term}, the \code{classify}, \code{tdf}, \code{alpha}, \code{sortFactor} 
       and the \code{sortOrder} will be set as attributes to the object. 
       Note that the \code{classify} in an \code{\link{alldiffs.object}} is based on the 
       variables indexing the predictions, which may differ from the 
       \code{classify} used to obtain the original predictions (for example, 
       when the \code{\link{alldiffs.object}}s stores a linear transformation of predictions.

       Also, see \code{\link{predictPlus.asreml}} for more information.
}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{as.alldiffs}},  \code{\link{as.predictions.frame}}, 
         \code{\link{sort.alldiffs}}, \code{\link{subset.alldiffs}},  \cr 
         \code{\link{print.alldiffs}},  \code{\link{renewClassify.alldiffs}},  
         \code{\link{redoErrorIntervals.alldiffs}},  \cr
         \code{\link{recalcLSD.alldiffs}},  \code{\link{pickLSDstatistics.alldiffs}}, 
         \code{\link{plotPredictions.data.frame}}, \cr 
         \code{\link{predictPlus.asreml}}, \code{\link{predictPresent.asreml}}}
\examples{
  data(Oats.dat)
  
  ## Use asreml to get predictions and associated statistics

  \dontrun{
  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  current.asrt <- as.asrtests(m1.asr)
  Var.pred <- asreml::predict.asreml(m1.asr, classify="Nitrogen:Variety", 
                                      sed=TRUE)
  if (getASRemlVersionLoaded(nchar = 1) == "3")
    Var.pred <- Var.pred$predictions
  Var.preds <- Var.pred$pvals
  Var.sed <- Var.pred$sed
  Var.vcov <- NULL
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  }

  ## Use lmerTest and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
    Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
    Var.preds <- summary(Var.emm)
    den.df <- min(Var.preds$df)
    ## Modify Var.preds to be compatible with a predictions.frame
    Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
    Var.vcov <- vcov(Var.emm)
    Var.sed <- NULL
  }

  ## Use the predictions obtained with either asreml or lmerTest
  if (exists("Var.preds"))
  {
    ## Order the Varieties in decreasing order for the predictions values in the 
    ## first N level 
    Var.diffs <- allDifferences(predictions = Var.preds, 
                                classify = "Nitrogen:Variety", 
                                sed = Var.sed, vcov = Var.vcov, tdf = den.df,
                                sortFactor = "Variety", decreasing = TRUE)
    print.alldiffs(Var.diffs, which="differences")
  
    ## Change the order of the factors in the alldiffs object and reorder components
    Var.reord.diffs <- allDifferences(predictions = Var.preds,
                                classify = "Variety:Nitrogen", 
                                sed = Var.sed, vcov = Var.vcov, tdf = den.df)
    print.alldiffs(Var.reord.diffs, which="predictions")
  }
}
\keyword{asreml}