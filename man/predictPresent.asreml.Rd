\name{predictPresent.asreml}
\alias{predictPresent.asreml}
\alias{predictPresent}
\title{Forms the predictions for each of one or more terms and presents 
       them in tables and/or graphs.}
\description{This function forms the predictions for each term in \code{terms} 
             using a supplied \code{asreml} object and 
             \code{\link{predictPlus.asreml}}. 
             Tables are produced using \code{\link{predictPlus.asreml}}, 
             in conjunction with \code{\link{allDifferences.data.frame}}, 
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
\method{predictPresent}{asreml}(asreml.obj, terms, 
               linear.transformation = NULL, 
               wald.tab = NULL, dDF.na = "residual", dDF.values = NULL, 
               x.num = NULL, x.fac = NULL, nonx.fac.order = NULL,  
               x.pred.values = NULL, x.plot.values = NULL, 
               plots = "predictions", panels = "multiple", 
               graphics.device = NULL, 
               error.intervals = "Confidence", meanLSD.type = "overall", 
               LSDby = NULL, avsed.tolerance = 0.25, titles = NULL, 
               colour.scheme = "colour", save.plots = FALSE, 
               transform.power = 1, offset = 0, scale = 1, 
               pairwise = TRUE, Vmatrix = FALSE, 
               tables = "all", level.length = NA, 
               alpha = 0.05, inestimable.rm = TRUE, 
               sortFactor = NULL, sortWithinVals = NULL, 
               sortOrder = NULL, decreasing = FALSE, 
               trace = FALSE, ggplotFuncs = NULL, ...)
}
\arguments{
 \item{asreml.obj}{\code{asreml} object for a fitted model.}
 \item{terms}{a \code{character vector} giving the terms for which predictions 
                 are required.}
 \item{linear.transformation}{a \code{\link{formula}} or a \code{\link{matrix}} specifying 
               a linear transformation to be applied to the predictions. 
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
               are returned in the \code{alldiffs} object.}
 \item{wald.tab}{a \code{data frame} containing the pseudo-anova table for the 
                 fixed terms produced by a call to \code{wald.asreml}. The main 
                 use of it here is in getting denominator degrees of freedom when
                 confidence intervals are to be plotted.}
 \item{dDF.na}{the method to use to obtain approximate denominator degrees of freedom. 
              when the numeric or algebraic methods produce an \code{NA}. Consistent with 
              when no denDF are available, the default is \code{"residual"} and so the residual 
              degrees of freedom from \code{asreml.obj$nedf} are used. If 
              \code{dDF.na = "none"}, no subtitute denominator degrees of freedom 
              are employed; if \code{dDF.na = "maximum"}, the maximum of those denDF 
              that are available, excluding that for the Intercept, is used; if all 
              denDF are \code{NA}, \code{asreml.obj$nedf} is used. If 
              \code{dDF.na = "supplied"}, a \code{vector} of values for the 
              denominator degrees of freedom is to be supplied in \code{dDF.values}. 
              Any other setting is ignored and a warning message produced. Generally, 
              substituting these degrees of freedom is anticonservative in that it 
              is likely that the degrees of freedom used will be too large.}
 \item{dDF.values}{A \code{vector} of values to be used when \code{dDF.na = "supplied"}. 
              Its values will be used when \code{denDF} in a test for a fixed effect
              is \code{NA}. This vector must be the same length as the number of 
              fixed terms, including (Intercept) whose value could be \code{NA}.}
 \item{x.num}{A character string giving the name of the numeric covariate that 
              corresponds to \code{x.fac}, is potentially included in terms in the 
              fitted model and which corresponds to the x-axis variable. It should 
              have the same number of unique values as the number of levels in 
              \code{x.fac}.}
 \item{x.fac}{A \code{character} string giving the name of the factor that corresponds to 
              \code{x.num}, is potentially included in terms in the fitted model and 
              which corresponds to the x-axis variable. It should have the same 
              number of levels as the number of unique values in \code{x.num}. 
              The levels of \code{x.fac} must be in the order in which they are to 
              be plotted - if they are dates, then they should be in the form 
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
              \code{panels} argument is ignored for for bar charts.}
 \item{graphics.device}{A \code{character} specifying a graphics device for plotting. 
              The default is \cr\code{graphics.device = NULL}, which will result 
              in plots being produced on the current graphics device. Setting it to 
              \code{"windows"}, for example,  will result in a windows graphics 
              device being opened.}
 \item{error.intervals}{A \code{character} string indicating the type of error interval, if any, 
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
             minLSD and max LSD is stored in the \code{LSD} component of the \code{alldiffs} 
             object when \code{meanLSD.type} is \code{factor.combinatons}.}
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
 \item{transform.power}{A number specifying the power of a transformation, if one 
              has been applied to the response variable. Unless it is equal to 1, the default, 
              back-transforms of the predictions will be obtained and presented in tables 
              or graphs as appropriate. The back-transform will raise the predictions to the 
              power equal to the reciprocal of \code{transform.power}, unless it equals 0 in 
              which case the exponential will be taken. Any scaling and offsetting will also be 
              taken into account in the backtransformation.} 
 \item{offset}{A number that has been added to each value of the response after any scaling 
               and before applying any power transformation. Unless it is equal to 0, the 
               default, back-transforms of the predictions will be obtained and presented in 
               tables or graphs as appropriate. The backtransformation will, after 
               backtransforming for any power transformation, subtract the \code{offset}.}
 \item{scale}{A number by which each value of the response has been multiply before adding 
               any offset and applying any power transformation. Unless it is equal to 1, the
               default, back-transforms of the predictions will be obtained and presented in 
               tables or graphs as appropriate. The backtransformation will, after backtransforming 
               for any power transformation and then subtracting the offset, divide by the \code{scale}.}
 \item{pairwise}{A logical indicating whether all pairise differences of the
                \code{predictions} and their standard errors and p-values are to be 
                computed and stored. If \code{tables} is equal to \code{"differences"} 
                or \code{"all"} or \code{error.intervals} is equal to 
                \code{"halfLeastSignificant"}, they will be stored irrespective of the 
                value of \code{pairwise}.}
\item{Vmatrix}{A \code{\link{logical}} indicating whether the variance matrix of the 
             \code{predictions} will be stored as a component of the \code{alldiffs} 
             object that is returned. If \code{linear.transformation} is set, it will 
             be stored irrespective of the value of \code{Vmatrix}.}
 \item{tables}{A \code{character vector} containing a combination of 
              \code{predictions}, \code{vcov}, \code{backtransforms}, 
              \code{differences}, \code{p.differences}, \code{sed}, 
              \code{LSD} and \code{all}.
              These nominate which components of the \code{alldiffs} 
              object to print.}
 \item{level.length}{The maximum number of characters from the the levels of 
            factors to use in the row and column labels of the tables produced by 
            \code{\link{allDifferences.data.frame}}.}
 \item{alpha}{the significance level for LSDs or 1 - alpha is the confidence level 
              for confidence intervals.}
 \item{inestimable.rm}{A \code{logical} indicating whether rows for  predictions that 
                       are not estimable are to be removed from the components of the
                       \code{alldiffs} object.}
\item{sortFactor}{A \code{character} containing the name of the 
              \code{factor} that indexes the set of predicted values that determines 
              the sorting of the \code{alldiffs} components by 
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
              \code{alldiffs} object. It specifies the desired order of the 
              levels in the reordered components of the \code{alldiffs} object. 
              The argument \code{sortWithinVals} is ignored.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
 \item{decreasing}{A \code{logical} passed to \code{order} that detemines whether 
              the order for sorting the \code{alldiffs} components is for 
               increasing or decreasing magnitude of the predicted values.}
 \item{trace}{if TRUE then partial iteration details are displayed when ASReml-R 
              functions are invoked; if FALSE then no output is displayed.}
 \item{ggplotFuncs}{A \code{\link{list}}, each element of which contains the 
                results of evaluating a \code{\link{ggplot}} function. 
                It is created by calling the \code{\link{list}} function with 
                a \code{\link{ggplot}} function call for each element. 
                It is passed to \code{\link{plotPredictions.data.frame}}.}
 \item{\dots}{further arguments passed to \code{predict.asreml} via 
            \code{\link{predictPlus.asreml}} and to \code{ggplot} 
            via \code{\link{plotPredictions.data.frame}}.}
}
\value{a \code{\link{list}} containing an \code{alldiffs} object for each \code{term} for 
             which tables are produced. The names of the components of this list are 
             the \code{terms} with full-stops (\code{.}) replacing colons (\code{:}). 
             Plots are also preduced depending on the setting of the \code{plot} 
             argument.}

\author{Chris Brien}
\seealso{\code{\link{predictPlus.asreml}}, \code{\link{allDifferences.data.frame}}, 
         \code{\link{sort.alldiffs}}, \code{\link{subset.alldiffs}}, \cr 
         \code{\link{redoErrorIntervals.alldiffs}},  \code{\link{recalcLSD.alldiffs}},  
         \code{\link{plotPredictions.data.frame}}, \cr 
         \code{\link{print.alldiffs}}, \code{as.Date}, \code{Devices}}
\examples{\dontrun{
data(WaterRunoff.dat)
titles <- list("Days since first observation", "Days since first observation", 
               "pH", "Turbidity (NTU)")
names(titles) <- names(WaterRunoff.dat)[c(5,7,11:12)]
asreml.options(keep.order = TRUE) #required for asreml4 only
current.asr <- asreml(fixed = log.Turbidity ~ Benches + Sources + Type + Species + 
                                 Sources:Type + Sources:Species + Sources:Species:xDay + 
                                 Sources:Species:Date, 
                      data = WaterRunoff.dat, keep.order = TRUE)
current.asrt <- asrtests(current.asr, NULL, NULL)

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