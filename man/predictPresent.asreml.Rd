\name{predictPresent.asreml}
\alias{predictPresent.asreml}
\alias{predictPresent}
\title{Forms the predictions for each of one or more terms and presents 
       them in tables and/or graphs.}
\description{This function forms the predictions for each term in \code{terms} 
             using a supplied \code{asreml} object and 
             \code{\link{predictPlus.asreml}}. 
             Tables are produced using \code{\link{predictPlus.asreml}}, 
             in conjunction with \cr 
             \code{\link{allDifferences.data.frame}}, 
             with the argument \code{tables} specifying which tables are printed. 
             The argument \code{plots}, along with \code{transform.power}, 
             controls which plots are produced. The plots are 
             produced using \code{\link{plotPredictions.data.frame}}, with 
             line plots produced when variables involving \code{x.num} or \code{x.fac} 
             are involved in \code{classify} for the predictions and bar charts 
             otherwise.  
             In order to get the correct predictions you may 
             need to supply additional arguments to \code{predict.asreml} 
             through \code{\dots} e.g. \code{present}, \code{parallel}, \code{levels}.

             The order of plotting the levels of 
             one of the factors indexing the predictions can be modified and is achieved 
             using \code{\link{sort.alldiffs}}.
}
\usage{
\method{predictPresent}{asreml}(asreml.obj, terms, inestimable.rm = TRUE, 
               linear.transformation = NULL, EGLS.linTransform, 
               error.intervals = "Confidence", alpha = 0.05, 
               wald.tab = NULL, dDF.fault = "residual", dDF.values = NULL, 
               pairwise = TRUE, Vmatrix = FALSE, 
               avsed.tolerance = 0.25, accuracy.threshold = NA, 
               LSDtype = "overall", LSDsupplied = NULL, LSDby = NULL, 
               LSDstatistic = "mean", LSDaccuracy = "maxAbsDeviation", 
               x.num = NULL, x.fac = NULL, nonx.fac.order = NULL,  
               x.pred.values = NULL, x.plot.values = NULL, 
               plots = "predictions", panels = "multiple", 
               graphics.device = NULL, interval.annotate = TRUE,
               titles = NULL, colour.scheme = "colour", save.plots = FALSE, 
               transform.power = 1, offset = 0, scale = 1, 
               transform.function = "identity", 
               tables = "all", level.length = NA, 
               sortFactor = NULL, sortParallelToCombo = NULL, 
               sortNestingFactor = NULL, sortOrder = NULL, 
               decreasing = FALSE, 
               trace = FALSE, ggplotFuncs = NULL, ...)
}
\arguments{
 \item{asreml.obj}{\code{asreml} object for a fitted model.}
 \item{terms}{A \code{character vector} giving the terms for which predictions 
                 are required.}
\item{inestimable.rm}{A \code{logical} indicating whether rows for  predictions that 
                       are not estimable are to be removed from the components of the
                       \code{\link{alldiffs.object}}.}
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
               are returned in the \code{\link{alldiffs.object}}.}
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
 \item{error.intervals}{A \code{character} string indicating the type of error interval, if any, 
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
 \item{dDF.fault}{A \code{\link{character}} specifying the method to use to obtain substitute 
              denominator degrees of freedom. when the numeric or algebraic methods produce 
              faulty values, viz. \code{NA}, \code{Inf} or less than 0.01. Consistent with 
              when no denDF are available, the default is \code{"residual"} and so the 
              residual degrees of freedom from \code{asreml.obj$nedf} are used. 
              If \code{dDF.fault = "none"}, no substitute denominator degrees of freedom are 
              employed; if \code{dDF.fault = "residual"}; if \code{dDF.fault = "maximum"}, 
              the maximum of those denDF that are available, excluding that for the 
              Intercept, is used; if all denDF are faulty, \code{asreml.obj$nedf} is used. 
              If \code{dDF.fault = "supplied"}, a \code{vector} of values for the denominator 
              degrees of freedom is to be supplied in \code{dDF.values}. Any other setting 
              is ignored and a warning message produced. Generally, substituting these 
              degrees of freedom is anticonservative in that it is likely that the degrees 
              of freedom used will be too large.}
 \item{dDF.values}{A \code{vector} of values to be used when \code{dDF.fault = "supplied"}. 
              Its values will be used when \code{denDF} in a test for a fixed effect
              is \code{NA}. This vector must be the same length as the number of 
              fixed terms, including (Intercept) whose value could be \code{NA}.}
 \item{pairwise}{A logical indicating whether all pairwise differences of the
                \code{predictions} and their standard errors and p-values are to be 
                computed and stored. If \code{tables} is equal to \code{"differences"} 
                or \code{"all"} or \code{error.intervals} is equal to 
                \code{"halfLeastSignificant"}, they will be stored irrespective of the 
                value of \code{pairwise}.}
\item{Vmatrix}{A \code{\link{logical}} indicating whether the variance matrix of the 
             \code{predictions} will be stored as a component of the \code{\link{alldiffs.object}} 
             that is returned. If \code{linear.transformation} is set, it will 
             be stored irrespective of the value of \code{Vmatrix}.}
\item{avsed.tolerance}{A \code{\link{numeric}} giving the value of the SED range, the range of the SEDs 
             divided by the square root of the mean of the variances of all or a subset of the 
             pairwise differences, that is considered reasonable in calculating 
             \code{error.intervals}. It should be a value between 0 and 1. The following rules apply:
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
 \item{nonx.fac.order}{A \code{character vector} giving the order in which factors other 
              than \code{x.fac} are to be plotted in plots with multiple panels 
              (i.e. where the number of non-x factors is greater than 1). 
              The first factor in the vector 
              will be plotted on the X axis (if there is no \code{x.num} or 
              \code{x.fac}. Otherwise, the order of plotting the factors is in 
              columns (X facets) and then rows (Y facets). By default the order is 
              in decreasing order for the numbers of levels of the non x factors.}
 \item{x.pred.values}{The values of \code{x.num} for which predicted values are required.}
 \item{x.plot.values}{The actual values to be plotted on the x axis or in the labels of 
               tables. They are 
               needed when values different to those in \code{x.num} are to be 
               plotted or \code{x.fac} is to be plotted because there is no 
               \code{x.num} term corresponding to the same term with \code{x.fac}.} 
 \item{plots}{Possible values are \code{"none"}, \code{"predictions"}, 
              \code{"backtransforms"} and \code{"both"}. Plots are not produced 
              if the value is \code{"none"}. If data are not transformed for 
              analysis (\code{transform.power} = 1), a plot of the predictions 
              is produced provided \code{plots} is not \code{"none"}. If the 
              data are transformed, the value of \code{plots} determines what 
              is produced.}
 \item{panels}{Possible values are \code{"single"} and \code{"multiple"}. 
              When line plots are to be produced, because variables involving 
              \code{x.num} or \code{x.fac} are involved in \code{classify} for 
              the predictions, \code{panels} determines whether or not a single 
              panel or multiple panels in a single window are produced. The 
              \code{panels} argument is ignored for bar charts.}
 \item{graphics.device}{A \code{character} specifying a graphics device for plotting. 
              The default is \cr\code{graphics.device = NULL}, which will result 
              in plots being produced on the current graphics device. Setting it to 
              \code{"windows"}, for example,  will result in a windows graphics 
              device being opened.}
 \item{interval.annotate}{A \code{logical} indicating whether the plot annotation indicating the 
             type of \code{error.interval} is to be included in the plot.}
\item{titles}{A \code{list}, each component of which is named for a column in 
               the \code{data.frame} for \code{asreml.obj} and contains a 
               \code{character string} giving a title to use in output (e.g. tables and 
               graphs). Here they will be used for axis labels.}
\item{colour.scheme}{A character string specifying the colour scheme for the plots. 
              The default is \code{"colour"} which produces coloured lines and bars, 
              a grey background and white gridlines. A value of \code{"black"} 
              results in black lines, grey bars and gridlines and a white background.}
\item{save.plots}{A \code{logical} that determines whether any plots will be saved. 
              If they are to be saved, a file name will be generated that consists of the 
              following elements separated by full stops: the response variable name with 
              \code{.back} if backtransformed values are being plotted, 
              the classify term, \code{Bar} or \code{Line} and, if \code{error.intervals} 
              is not \code{"none"}, one of \code{SE}, \code{CI} or \code{LSI}. The 
              file will be saved as a `png' file in the current work directory.}
\item{transform.power}{A \code{\link{numeric}} specifying the power of a transformation, if 
              one has been applied to the response variable. Unless it is equal 
              to 1, the default, back-transforms of the predictions will be 
              obtained and stored in the \code{backtransforms} component of the 
              \code{\link{alldiffs.object}}. The \code{plots} and \code{tables} arguments 
              control the plotting and output of the \code{predictions} and 
              \code{backtransforms}. 
              The back-transformation raises the predictions to the power equal 
              to the reciprocal of \code{transform.power}, unless it equals 0 in 
              which case the exponential of the predictions is taken.} 
\item{offset}{A number that has been added to each value of the response after any scaling 
               and before applying any power transformation. Unless it is equal to 0, the 
               default, back-transforms of the predictions will be obtained and stored in 
               the \code{backtransforms} component of the \code{\link{alldiffs.object}}. 
               The \code{plots} and \code{tables} arguments 
              control the plotting and output of the \code{predictions} and 
              \code{backtransforms}. The backtransformation will, after 
               backtransforming for any power transformation, subtract the \code{offset}.}
\item{scale}{A number by which each value of the response has been multiply before adding 
               any offset and applying any power transformation. Unless it is equal to 1, the
               default, back-transforms of the predictions will be obtained and stored in 
               the \code{backtransforms} component of the \code{\link{alldiffs.object}}. 
               The \code{plots} and \code{tables} arguments 
              control the plotting and output of the \code{predictions} and 
              \code{backtransforms}. The backtransformation will, after backtransforming 
               for any power transformation and then subtracting the offset, divide by the \code{scale}.}
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
\item{tables}{A \code{character vector} containing a combination of 
              \code{predictions}, \code{vcov}, \code{backtransforms}, 
              \code{differences}, \code{p.differences}, \code{sed}, 
              \code{LSD} and \code{all}.
              These nominate which components of the \code{\link{alldiffs.object}} 
              to print.}
\item{level.length}{The maximum number of characters from the levels of 
            factors to use in the row and column labels of the tables produced by 
            \code{\link{allDifferences.data.frame}}.}
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
\item{trace}{If TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
\item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{ggplot2} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{ggplot2} function call for each element. 
                It is passed to \code{\link{plotPredictions.data.frame}}.}
\item{\dots}{further arguments passed to \code{predict.asreml} via 
            \code{\link{predictPlus.asreml}} and to \code{ggplot} 
            via \code{\link{plotPredictions.data.frame}}.}
}
\value{A \code{\link{list}} containing an \code{\link{alldiffs.object}} for each \code{term} for 
             which tables are produced. The names of the components of this list are 
             the \code{terms} with full-stops (\code{.}) replacing colons (\code{:}). 
             Plots are also preduced depending on the setting of the \code{plot} 
             argument.}

\author{Chris Brien}
\seealso{\code{\link{predictPlus.asreml}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{sort.alldiffs}}, \code{\link{subset.alldiffs}}, \cr 
         \code{\link{redoErrorIntervals.alldiffs}},  \code{\link{recalcLSD.alldiffs}},  
         \code{\link{pickLSDstatistics.alldiffs}}, \cr 
         \code{\link{plotPredictions.data.frame}}, 
         \code{\link{print.alldiffs}}, \code{as.Date}, \code{Devices}}
\examples{\dontrun{
data(WaterRunoff.dat)
titles <- list("Days since first observation", "Days since first observation", 
               "pH", "Turbidity (NTU)")
names(titles) <- names(WaterRunoff.dat)[c(5,7,11:12)]
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = log.Turbidity ~ Benches + Sources + Type + Species + 
                                 Sources:Type + Sources:Species + Sources:Species:xDay + 
                                 Sources:Species:Date, 
                      data = WaterRunoff.dat, keep.order = TRUE)
current.asrt <- as.asrtests(current.asr, NULL, NULL)

#### Get the observed combinations of the factors and variables in classify
class.facs <- c("Sources","Species","Date","xDay")
levs <- as.data.frame(table(WaterRunoff.dat[class.facs]))
levs <- levs[do.call(order, levs), ]
levs <- as.list(levs[levs$Freq != 0, class.facs])
levs$xDay <- as.numfac(levs$xDay)
  
#### parallel and levels are arguments from predict.asreml
diff.list <- predictPresent.asreml(asreml.obj = current.asrt$asreml.obj, 
                                   terms = "Date:Sources:Species:xDay",
                                   x.num = "xDay", x.fac = "Date", 
                                   parallel = TRUE, levels = levs, 
                                   wald.tab = current.asrt$wald.tab, 
                                   plots = "predictions", 
                                   error.intervals = "StandardError", 
                                   titles = titles, 
                                   transform.power = 0, 
                                   present = c("Type","Species","Sources"), 
                                   tables = "none", 
                                   level.length = 6)
}}
\keyword{asreml}