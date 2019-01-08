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
             standard errors of the predictions. If the deominator degrees of 
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
             \code{data.frame}. If the \code{error.intervals} are 
             \code{halfLeastSignificant}, then backtransforms of the Confidence 
             limits are stored, instead of the \code{halfLeastSignificant} 
             intervals. If \code{transform.power} is other than 
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
            linear.transformation = NULL, titles = NULL,  
            x.num = NULL, x.fac = NULL,  
            x.pred.values = NULL, x.plot.values = NULL, 
            error.intervals = "Confidence", avsed.tolerance = 0.25, 
            meanLSD.type = "overall", LSDby = NULL, 
            pairwise = TRUE, Vmatrix = FALSE, 
            tables = "all" , level.length = NA, 
            transform.power = 1, offset = 0, scale = 1, 
            inestimable.rm = TRUE, 
            sortFactor = NULL, sortWithinVals = NULL, 
            sortOrder = NULL, decreasing = FALSE, 
            wald.tab = NULL, alpha = 0.05, 
            dDF.na = "residual",  dDF.values = NULL, 
            trace = FALSE, ...)
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
               It might be a contrast \code{\link{matrix}} or a \code{\link{matrix}} of 
               weights for a factor used to obtain the weighted average over that factor. 
               The number of rows in the \code{\link{matrix}} should equal the 
               number of linear combinations of the \code{predictions} desired and 
               the number of columns should equal the number of \code{predictions}.
               
               In either case, as well as the values of the linear combinations, 
               their standard errors, pairwise differences and associated statistics 
               are returned.}
 \item{titles}{A \code{\link{list}}, each component of which is named for a column in 
               the \code{data.frame} for \code{asreml.obj} and contains a 
               \code{character string} giving a title to use 
               in output (e.g. tables and graphs). Here they will 
               be used for table headings.}
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
 \item{x.pred.values}{The values of \code{x.num} for which predicted values are 
              required. If \code{levels} is set for passing to \code{predict.asreml}, 
              \code{x.pred.values} is ignored. Note that while \code{levels} is and
              alternative to \code{x.pred.values}, it allows more general setting 
              of the levels to be predicted.}
 \item{x.plot.values}{The actual values to be plotted on the x axis. They are 
               needed when values different to those in \code{x.num} are to be 
               plotted or \code{x.fac} is to be plotted because there is no 
               \code{x.num} term corresponding to the same term with \code{x.fac}.} 
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
\item{pairwise}{A \code{\link{logical}} indicating whether all pairwise differences of the
             \code{predictions} and their standard errors and p-values are to be 
             computed and stored. If \code{tables} is equal to 
             \code{"differences"} or \code{"all"} or \code{error.intervals} is 
             equal to \code{"halfLeastSignificant"}, they will be stored 
             irrespective of the value of \code{pairwise}.}
\item{Vmatrix}{A \code{\link{logical}} indicating whether the variance matrix of the 
             \code{predictions} will be stored as a component of the \code{\link{alldiffs.object}} 
             that is returned. If \code{linear.transformation} is set, it will 
             be stored irrespective of the value of \code{Vmatrix}.}
 \item{tables}{A \code{\link{character}} vector containing a combination of 
              \code{none}, 
              \code{predictions}, \code{vcov}, \code{backtransforms}, \code{differences}, 
              \code{p.differences}, \code{sed}, \code{LSD} and \code{all}.
              These nominate which components of the \code{\link{alldiffs.object}} to print.}
 \item{level.length}{The maximum number of characters from the the levels of 
            factors to use in the row and column labels of the tables of 
            pairwise differences and their p-values and standard errors.}
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
\item{sortFactor}{A \code{character} containing the name of the 
              \code{factor} that indexes the set of predicted values that determines 
              the sorting of the components of the \code{\link{alldiffs.object}} by 
              \code{\link{sort.alldiffs}}. If \code{NULL} then sorting is not carried 
              out. If there is more than one variable 
              in the \code{classify} term then \code{sortFactor} is sorted for the 
              predicted values within each combination of the values of the sortWithin 
              variables: the \code{classify} variables, excluding the 
              \code{sortFactor}. There should be only one predicted value for 
              each unique value of \code{sortFactor} within each set defined by a 
              combination of the values of the sortWithin variables.}
 \item{sortWithinVals}{A \code{list} with a component named for each \code{factor} and 
              \code{numeric} that is a \code{classify} variable for the predictions, 
              excluding \code{sortFactor}. Each component should contain a single 
              value that is a value of the variable. The combination of this set 
              of values will be used to define a subset of the predicted values whose 
              order will define the order of \code{sortFactor} to be used for all 
              combinations of the \code{sortWithinVals} variables. If 
              \code{sortWithinVals} is \code{NULL} then the first value of each 
              sortWithin variable in \code{predictions} component is used 
              to define \code{sortWithinVals}. If there is only one variable in the 
              classify then \code{sortWithinVals} is ignored.}
 \item{sortOrder}{A \code{character vector} whose length is the same as the number 
               of levels for \code{sortFactor} in the \code{predictions} component of the 
              \code{\link{alldiffs.object}}. It specifies the desired order of the 
              levels in the reordered components of the \code{\link{alldiffs.object}}. 
              The argument \code{sortWithinVals} is ignored.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
 \item{decreasing}{A \code{logical} passed to \code{order} that detemines whether 
              the order for sorting the components of the \code{\link{alldiffs.object}} is for 
               increasing or decreasing magnitude of the predicted values.}
 \item{wald.tab}{A \code{\link{data.frame}} containing the pseudo-anova table for the 
                 fixed terms produced by a call to \code{wald.asreml}. The main 
                 use of it here is in determinining the degrees of freedom of the standard errors
                 of the predictions. denominator degrees of freedom when
                 p-values or confidence intervals are to be calculated.}
 \item{alpha}{A \code{\link{numeric}} giving the significance level for LSDs or one minus 
              the confidence level for confidence intervals.}
 \item{dDF.na}{A \code{\link{character}} specifying the method to use to 
               obtain approximate denominator degrees of 
               freedom. when the numeric or algebraic methods produce an 
               \code{NA}. Consistent with when no denDF are available, the 
               default is \code{"residual"} and so the residual degrees of 
               freedom from \code{asreml.obj$nedf} are used. 
               If \code{dDF.na = "none"}, no subtitute denominator degrees of 
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
         \code{\link{recalcLSD.alldiffs}},  \code{\link{predictPresent.asreml}}, 
         \code{\link{plotPredictions.data.frame}}, \code{as.Date}, \cr
         \code{predict.asreml}}
\examples{\dontrun{
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- asrtests(current.asr, NULL, NULL)
diffs <- predictPlus(classify = "Sources:Type", 
                     asreml.obj = current.asr, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("Sources", "Type", "Species"))
}}
\keyword{asreml}