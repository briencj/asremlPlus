\name{predictPlus.asreml}
\alias{predictPlus.asreml}
\alias{predictPlus}
\title{Forms the predictions for a term, their pairwise differences and 
       associated statistics. A factor having parallel values may occur in 
       the model and a linear transformation of the predictions can be specified.
       It results in an object of class alldifffs. }
\description{This function forms the predictions for \code{term} using 
             \code{classify} and the supplied \code{asreml} object and stores 
             them in an \code{\link{alldiffs.object}}. If \code{x.num} is 
             supplied, the predictions will be obtained for the values supplied 
             in \code{x.pred.values} and, if supplied, \code{x.plot.values} will 
             replace them in the \code{\link{alldiffs.object}} that is returned. 
             If \code{x.fac}, but not \code{x.num}, is specified, predictions 
             will involve it and, if supplied, \code{x.plot.values} will replace 
             the levels of \code{x.fac} in the \code{\link{alldiffs.object}} 
             that is returned. In order to get the correct predictions you may 
             need to supply additional arguments to \code{predict.asreml} 
             through \code{\dots} e.g. \code{present}, \code{parallel}, \code{levels}. 
             Any aliased predictions will be removed, as 
             will any standard error of pairwise differences involving them.

             Also calculated are the approximate degrees of freedom of the 
             standard errors of the predictions. If the denominator degrees of 
             freedom for \code{term} are available in \code{wald.tab}, they are 
             used. Otherwise the residual degrees of freedom or the maximum of 
             the denominator degrees in \code{wald.tab}, excluding the  
             Intercept, are used. Which is used depends on the setting of 
             \code{dDF.na}. These degrees of freedom are used for the 
             t-distribution on which p-values and confidence intervals are 
             based. It is stored as an  attribute to the \code{\link{alldiffs.object}}. 
             The degrees of freedom are also used in calculating the minimum, 
             mean and maximum LSD for comparing pairs of predictions, which are 
             also stored in the \code{\link{alldiffs.object}}.

             If \code{pairwise = TRUE}, all pairwise differences between the 
             \code{predictions}, their standard errors, p-values and LSD 
             statistics are computed using \code{\link{allDifferences.data.frame}}. 
             This adds them to the \code{\link{alldiffs.object}} as additional 
             \code{list} components named \code{differences}, \code{sed}, 
             \code{p.differences} and \code{LSD}.

            If a linear transformation of the predictions is specified then the 
            values of this linear transformation are returned, 
            instead of the original predictions, along with their standard errors 
            and the pairwise differences and associated statistics. 

             If a transformation has been applied in the analysis (any one of 
             \code{transform.power} is not one, \code{scale} is  not one and 
             \code{offset} is nonzero), the backtransforms of the transformed 
             values and their lower and upper error intervals are added 
             to a \code{data.frame} that is consistent with the \code{predictions} 
             \code{data.frame}. 
              If \code{transform.power} is other than 
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
\method{predictPlus}{asreml}(asreml.obj, classify, term = NULL, 
            inestimable.rm = TRUE, 
            linear.transformation = NULL, EGLS.linTransform = TRUE, 
            error.intervals = "Confidence", alpha = 0.05, 
            wald.tab = NULL, dDF.na = "residual",  dDF.values = NULL, 
            pairwise = TRUE, Vmatrix = FALSE, 
            avsed.tolerance = 0.25, accuracy.threshold = NA, 
            LSDtype = "overall", LSDsupplied = NULL, LSDby = NULL, 
            LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
            x.num = NULL, x.fac = NULL,  
            x.pred.values = NULL, x.plot.values = NULL, 
            titles = NULL,  tables = "all" , level.length = NA, 
            transform.power = 1, offset = 0, scale = 1, 
            transform.function = "identity", 
            sortFactor = NULL, sortParallelToCombo = NULL, 
            sortNestingFactor = NULL, sortOrder = NULL, 
            decreasing = FALSE, trace = FALSE, ...)
}
\arguments{
 \item{asreml.obj}{\code{asreml} object for a fitted model.}
 \item{classify}{A \code{\link{character}} string giving the variables that 
                 define the margins of the multiway table to be predicted. 
                 Multiway tables are specified by forming an interaction type 
                 term from the classifying variables, that is, separating the 
                 variable names with the \code{:} operator. To predict the overall 
                 mean, set the \code{classify} to "\code{(Intercept)}".}
 \item{term}{A \code{\link{character}} string giving the variables that define the term 
               that was fitted using \code{asreml} and that corresponds to \code{classify}. 
               It only needs to be specified when it is different to \code{classify}; it 
               is stored as an attribute of the \code{\link{alldiffs.object}}. 
               It is likely to be needed when the fitted model includes terms that involve  
               both a \code{\link{numeric}} covariate and a \code{\link{factor}} that 
               parallel each other; the \code{classify} would include the covariate and 
               the \code{term} would include the \code{factor}.}
\item{inestimable.rm}{A \code{logical} indicating whether rows for  predictions 
               that are not estimable are to be removed from the components of 
               the \code{\link{alldiffs.object}}.}
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
               It might be a contrast \code{\link{matrix}} or a \code{\link{matrix}} of 
               weights for a factor used to obtain the weighted average over that factor. 
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
\item{alpha}{A \code{\link{numeric}} giving the significance level for LSDs or one minus 
              the confidence level for confidence intervals.
              It is stored as an attribute to the \code{\link{alldiffs.object}}.}
\item{wald.tab}{A \code{\link{data.frame}} containing the pseudo-anova table for the 
               fixed terms produced by a call to \code{wald.asreml}. The main 
               use of it here is in determining the degrees of freedom for 
               calculating confidence or half-LSD \code{error.intervals} and p-values, 
               the latter to be stored in the \code{p.differences} component of the 
               \code{alldiffs.object} that is created.}
\item{dDF.na}{A \code{\link{character}} specifying the method to use to 
               obtain approximate denominator degrees of 
               freedom. when the numeric or algebraic methods produce an 
               \code{NA}. Consistent with when no denDF are available, the 
               default is \code{"residual"} and so the residual degrees of 
               freedom from \code{asreml.obj$nedf} are used. 
               If \code{dDF.na = "none"}, no substitute denominator degrees of 
               freedom are employed; if \code{dDF.na = "maximum"}, the maximum 
               of those denDF that are available, excluding that for the 
               Intercept, is used; if all denDF are \code{NA}, \code{asreml.obj$nedf} is used. If 
              \code{dDF.na = "supplied"}, a \code{vector} of values for the 
              denominator degrees of freedom is to be supplied in \code{dDF.values}. 
              Any other setting is ignored and a warning message produced. 
              Generally, substituting these degrees of freedom is 
              anticonservative in that it is likely that the degrees of freedom 
              used will be too large.}
\item{dDF.values}{A \code{\link{vector}} of values to be used when 
              \code{dDF.na = "supplied"}. Its values will be used when 
              \code{denDF} in a test for a fixed effect is \code{NA}. 
              This vector must be the same length as the number of 
              fixed terms, including (Intercept) whose value could be 
              \code{NA}.}
\item{pairwise}{A \code{\link{logical}} indicating whether all pairwise differences of the
             \code{predictions} and their standard errors and p-values are to be 
             computed and stored. If \code{tables} is equal to 
             \code{"differences"} or \code{"all"} or \code{error.intervals} is 
             equal to \code{"halfLeastSignificant"}, they will be stored 
             irrespective of the value of \code{pairwise}.}
\item{Vmatrix}{A \code{\link{logical}} indicating whether the variance matrix of the 
             \code{predictions} will be stored as a component of the 
             \code{\link{alldiffs.object}} that is returned. 
              If \code{linear.transformation} is set, it will 
             be stored irrespective of the value of \code{Vmatrix}.}
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
 \item{titles}{A \code{\link{list}}, each component of which is named for a column in 
               the \code{data.frame} for \code{asreml.obj} and contains a 
               \code{character string} giving a title to use 
               in output (e.g. tables and graphs). Here they will 
               be used for table headings.}
\item{tables}{A \code{\link{character}} vector containing a combination of 
              \code{none}, 
              \code{predictions}, \code{vcov}, \code{backtransforms}, \code{differences}, 
              \code{p.differences}, \code{sed}, \code{LSD} and \code{all}.
              These nominate which components of the \code{\link{alldiffs.object}} to print.}
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
 \item{x.pred.values}{The values of \code{x.num} for which predicted values are 
              required. If \code{levels} is set for passing to \code{predict.asreml}, 
              \code{x.pred.values} is ignored. Note that while \code{levels} is an
              alternative to \code{x.pred.values}, \code{x.pred.values} allows more 
              general setting of the levels to be predicted.}
 \item{x.plot.values}{The actual values to be plotted on the x axis. They are 
               needed when values different to those in \code{x.num} are to be 
               plotted or \code{x.fac} is to be plotted because there is no 
               \code{x.num} term corresponding to the same term with \code{x.fac}.} 
\item{level.length}{The maximum number of characters from the levels of 
            factors to use in the row and column labels of the tables of 
            pairwise differences and their p-values and standard errors.}
\item{transform.power}{A \code{\link{numeric}} specifying the power of a transformation, if 
              one has been applied to the response variable. Unless it is equal 
              to 1, the default, back-transforms of the predictions will be 
              obtained and stored in the \code{backtransforms} component of the 
              \code{\link{alldiffs.object}}. 
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
              back-transformed using the inverse function of the \code{transform.function}. 
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
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
\item{decreasing}{A \code{logical} passed to \code{order} that detemines whether 
              the order for sorting the components of the \code{\link{alldiffs.object}} is for 
               increasing or decreasing magnitude of the predicted values.}
\item{trace}{A \code{\link{logical}} that control output from ASReml-R. 
              If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
\item{\dots}{further arguments passed to \code{predict.asreml}.}
}
\value{For \code{linear.transformations} set to \code{NULL}, an S3-class 
       \code{\link{alldiffs.object}} with predictions and their standard 
       errors and, depending on the settings of the arguments, all pairwise 
       differences between predictions, their standard errors and p-values 
       and LSD statistics.  Also, unless the \code{sortFactor} or \code{sortOrder} 
       arguments are invoked, the rows of \code{predictions} component are ordered 
       so that they are in standard order for the variables in the \code{classify}. 
       That is, the values of the last variable change with every row, those of the 
       second-last variable only change after all the values of the last variable have 
       been traversed; in general, the values of a variable are the same for all the 
       combinations of the values to the variables to its right in the \code{classify}. 
       In addition, if necessary, the order of the columns of the variables in the 
       \code{predictions} component are changed to match their order in the \code{classify}.

           If \code{transform.power} or \code{scale} is not one or \code{offset} 
           is not zero, it will contain a \code{data.frame} with the backtransformed 
           linear transformation of the predictions. The backtransformation will, after 
           backtransforming for any power transformation, subtract the \code{offset} 
           and then divide by the \code{scale}.
           
           If \code{error.intervals} is not \code{"none"}, then the 
           \code{predictions} component and, if present, the 
           \code{backtransforms} component will contain columns for the lower 
           and upper values of the limits for the interval.

           The name of the \code{response}, the \code{response.title}, 
           the \code{term}, the \code{classify}, \code{tdf}, \code{sortFactor} 
           and the \code{sortOrder} will be set as attributes to the object. 
           Also, if \code{error.intervals} is \code{"halfLeastSignificant"}, then those of 
           \code{LSDtype}, \code{LSDby} and \code{LSDstatistic} that are not \code{NULL} 
           will be added as attributes of the object and of the \code{predictions frame}; 
           additionally, \code{LSDvalues} will be added as attribute of the 
           \code{predictions frame}, \code{LSDvalues} being the LSD values used in 
           calculating the \code{error.intervals}. 
           Note that the \code{classify} in an \code{\link{alldiffs.object}} is based on the 
           variables indexing the predictions, which may differ from the 
           \code{classify} used to obtain the original predictions (for example, 
           when the \code{\link{alldiffs.object}}s stores a linear transformation of predictions.
           
           For \code{linear.transformations} set to other than \code{NULL}, 
           an \code{\link{alldiffs.object}} with the \code{linear.transformation} 
           applied to the \code{predictions} and their standard errors and, 
           depending on the settings of the arguments, all pairwise 
           differences between the linearly transformed predictions, their 
           standard errors and p-values and LSD statistics. 
           (See also \code{\link{linTransform.alldiffs}}.)}
\author{Chris Brien}
\seealso{\code{\link{alldiffs.object}}, \code{\link{as.alldiffs}}, \code{\link{print.alldiffs}}, 
         \code{\link{linTransform.alldiffs}}, \code{\link{sort.alldiffs}},  \cr 
         \code{\link{subset.alldiffs}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{redoErrorIntervals.alldiffs}},  \cr 
         \code{\link{recalcLSD.alldiffs}},  \code{\link{exploreLSDs.alldiffs}}, 
         \code{\link{pickLSDstatistics.alldiffs}}, \cr
         \code{\link{predictPresent.asreml}}, 
         \code{\link{plotPredictions.data.frame}}, \code{as.Date}, \code{predict.asreml}}
\examples{\dontrun{
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
diffs <- predictPlus(classify = "Sources:Type", 
                     asreml.obj = current.asr, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("Sources", "Type", "Species"))
}}
\keyword{asreml}