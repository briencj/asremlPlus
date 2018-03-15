\name{as.alldiffs}
\alias{as.alldiffs}
\title{Forms an object of S3-class \code{alldiffs} that stores the predictions for a model 
       fitted using asreml, along with supplied statistics for all pairwise differences.}
\description{Creates an object of S3-class \code{alldiffs} that consists of a \code{list} 
             containing the following components: \code{predictions}, \code{differences}, 
             \code{p.differences}, \code{sed}, \code{LSD} and \code{backtransforms}. 
             Predictions must be supplied to the functions while the others will be set 
             only if they are supplied; those not supplied are set to \code{NULL}. 
             It also has attributes \code{response}, \code{response.title}, \code{term}, 
             \code{classify}, \code{tdf}, \code{sortFactor} and \code{sortOrder}.
             which will be set to the values supplied or \code{NULL} if none are supplied.}

\usage{as.alldiffs(predictions, vcov = NULL, differences = NULL, 
            p.differences = NULL, sed = NULL, LSD = NULL, 
            backtransforms = NULL, 
            response = NULL, response.title = NULL, 
            term = NULL, classify = NULL, tdf = NULL, 
            sortFactor = NULL, sortOrder = NULL)}
\arguments{
 \item{predictions}{a \code{data.frame} containing the predicted values that is consistent with an 
                    object of class \code{asremlPredict} such as 
                    is stored in the \code{pvals} component of the \code{prediction} component 
                    of the value produced by \code{predict.asreml}. That is, in addition to 
                    variables classifying the predictions, it will include columns named 
                    \code{predicted.value}, \code{standard.error} and \code{est.status}; 
                    each row contains a single predicted value. Note that the names 
                    \code{standard.error} and \code{est.status} have been changed to 
                    \code{std.error} and \code{status} in asreml4; if the new names are 
                    present in \code{predictions}, they will be returned to the previous 
                    names. It may also contain
                    columns for the lower and upper confidence limits for the predictions. 
                    If \code{LSD} is not \code{NULL} then the mean LSD will be added as an 
                    attribute named \code{meanLSD}, it being calculated using the
                    square root of the mean of the variances of pairwise differences.}
 \item{differences}{a \code{matrix} containing all pairwise differences between
                    the predictions; it should have the same number of rows and columns as there are
                    rows in \code{predictions}.}
 \item{p.differences}{a \code{matrix} containing p-values for all pairwise differences 
                    between the predictions; each p-value is computed as the probability of a t-statistic 
                    as large as or larger than the observed difference divided by its standard error. 
                    The degrees of freedom of the t distribution for computing it are computed as 
                    the denominator degrees of freedom of the F value for the fixed term, if available; 
                    otherwise, the degrees of freedom stored in the attribute \code{tdf} are used; 
                    the matrix should be of the same size as that for \code{differences}.}
 \item{sed}{a \code{matrix} containing the standard errors of all pairwise differences 
                    between the predictions; they are used in computing the p-values.}
 \item{vcov}{a \code{matrix} containing the variance matrix of the predictions; it is used in 
                    computing the variance of linear transformations of the predictions.}
\item{LSD}{a \code{data.frame} containing the mean, minimum and maximum LSD for determining             
                    the significance of pairwise differences, the mean LSD being calculated using 
                    the square root of the mean of the variances of pairwise differences.}
 \item{backtransforms}{a \code{data.frame} containing the backtransformed values of the predicted values 
                    that is consistent with an object of class \code{asremlPredict} such as 
                    is stored in the \code{pvals} component of the \code{prediction} component 
                    of the value produced by \code{predict.asreml}. That is, in addition to 
                    variables classifying the predictions, it will include columns named 
                    \code{backtransformed.predictions} and \code{est.status}; it may also contain 
                    columns for the lower and upper confidence limits; 
                    each row contains a single predicted value.}
 \item{response}{a \code{character} specifying the response variable for the 
                    predictions. It is stored as an attribute to the \code{alldiffs} 
                    object.}
 \item{response.title}{a \code{character} specifying the title for the response variable 
                    for the predictions. It is stored as an attribute to the \code{alldiffs} 
                    object.}
 \item{term}{a character string giving the variables that define the term 
                 that was fitted using \code{asreml} and that corresponds 
                 to \code{classify}. It is often the same as \code{classify}. 
                 It is stored as an attribute to the \code{alldiffs} object.}
 \item{classify}{a character string giving the variables that define the margins
                 of the multiway table used in the prediction. Multiway tables are 
                 specified by forming an interaction type term from the 
                 classifying variables, that is, separating the variable names 
                 with the \code{:} operator. It is stored as an attribute to the 
                 \code{alldiffs} object.}
 \item{tdf}{an \code{integer} specifying the degrees of freedom of the standard error. It is used as 
                the degrees of freedom for the t-distribution on which p-values and confidence 
                intervals are based. 
                It is stored as an attribute to the \code{alldiffs} object.}
 \item{sortFactor}{A \code{character} containing the name of the 
              \code{factor} that indexes the set of predicted values that 
              determined the sorting of the components.}
 \item{sortOrder}{A \code{character} vector that is the same length as the number of levels for 
              \code{sortFactor} in the \code{predictions} component of the 
              \code{alldiffs} object. It specifies the order of the 
              levels in the reordered components of the \code{alldiffs} object.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
}
\value{An object of S3-class \code{alldiffs} that consists of a \code{list} 
                containing the following components as desribed in the arguments above: 
                \code{predictions}, \code{vcov}, \code{differences}, \code{p.differences}, \code{sed}, 
                \code{LSD} and \code{backtransforms} and with attributes \code{names}, 
                \code{class}, \code{response}, \code{response.title}, \code{term}, 
                \code{classify}, \code{tdf}, \code{meanLSD}, \code{sortFactor} and \code{sortOrder}. 
                The last two attributes are also set by \code{\link{sort.alldiffs}}.

                The name of the \code{response}, the \code{response.title}, 
                the \code{term}, the \code{classify}, \code{tdf}, \code{sortFactor} 
                and the \code{sortOrder} will be set as attributes to the object. 
                Note that the \code{classify} in an \code{alldiffs} object is based on the 
                variables indexing the predictions, which may differ from the 
                \code{classify} used to obtain the original predictions (for example, 
                when the \code{alldiffs} objects stores a linear transformation of predictions.

                See \code{\link{predictPlus.asreml}} for mnore information.
}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{print.alldiffs}}, \code{\link{sort.alldiffs}}, 
         \code{\link{subset.alldiffs}},  \cr
         \code{\link{allDifferences.data.frame}}, \code{\link{redoErrorIntervals.alldiffs}},  
         \code{\link{recalcLSD.alldiffs}}, \cr
         \code{\link{predictPlus.asreml}}, \code{\link{plotPredictions.data.frame}}, 
        \code{\link{predictPresent.asreml}}}
\examples{\dontrun{
  Var.pred <- predict(current.asr, classify="Variety", sed=TRUE)$predictions
  wald.tab <-  current.asrt$wald.tab
  den.df <- wald.tab[match("Variety", rownames(wald.tab)), "denDF"]
  Var.diffs <- as.alldiffs(predictions = Var.pred$pvals, 
                           sed = Var.pred$sed, 
                           tdf = den.df)
}}
\keyword{asreml}
\keyword{htest}