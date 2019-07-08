\name{alldiffs.object}
\alias{alldiffs.object}
\alias{alldiffs-class}
\title{Description of an alldiffs object}
\description{An object of S3-class \code{alldiffs} that that stores the predictions for a model, 
             along with supplied statistics for all pairwise differences. While 
             \code{\link{alldiffs.object}} can be constructed by defining a \code{list} with 
             the appropriate components, it can be formed by passing the components to 
             \code{\link{as.alldiffs}}, or from a \code{predictions} \code{data.frame} using 
             \code{\link{allDifferences.data.frame}}.

             \code{\link{as.alldiffs}} is function that assembles an object of this class from 
              supplied components.

             \code{\link{is.alldiffs}} is the membership function for this class; it tests 
             that an object is of class \code{alldiffs}.

             \code{\link{validAlldiffs}(object)} can be used to test the validity of an object 
             with this class. 

             \code{\link{allDifferences.data.frame}} is the function that constructs an 
             object of this class by calculating components from statistics supplied via 
             its arguments and then using \code{\link{as.alldiffs}} to make the object.
}
\value{A \code{list} of class \code{alldiffs} containing the following components: 
             \code{predictions}, \code{vcov}, \code{differences}, 
             \code{p.differences}, \code{sed}, \code{LSD} and \code{backtransforms}. 
             Except for \code{predictions}, the components are optional and can be set 
             to \code{NULL}. 
             
             An \code{\link{alldiffs.object}} also has attributes \code{response}, 
             \code{response.title}, \code{term}, \code{classify}, \code{tdf}, 
             \code{sortFactor} and \code{sortOrder}, which may be set to \code{NULL}.

            The details of the components are as follows:

       \enumerate{
          \item \code{predictions}: A \code{\link{predictions.frame}}, being a \code{data.frame} 
                    beginning with the variables classifying the predictions, in the same order 
                    as in the \code{classify}, and also containing columns named 
                    \code{predicted.value}, \code{standard.error} and \code{est.status}; 
                    each row contains a single predicted value. The number of rows should equal the 
                    number of unique combinations of the \code{classify} variables and will be in 
                    standard order for the \code{classify} variables. That is, the values of the 
                    last variable change with every row, those of the second-last variable only 
                    change after all the values of the last variable have been traversed; 
                    in general, the values of a variable are the same for all the 
                    combinations of the values to the variables to its right in the \code{classify}.

                    The \code{data.frame} may also include columns for the lower and upper 
                    values of error intervals, either standard error, confidence or half-LSD 
                    intervals.  The names of these columns will consist of three parts 
                    separated by full stops: 
                     1) the first part will be \code{lower} or \code{upper}; 
                     2) the second part will be one of \code{Confidence}, 
                        \code{StandardError} or \code{halfLeastSignificant}; 
                     3) the third component will be \code{limits}.

                    Note that the names \code{standard.error} and 
                    \code{est.status} have been changed to \code{std.error} and \code{status} 
                    in the \code{pvals} component produced by \code{asreml-R4}; if the new names 
                    are in the \code{data.frame} supplied to \code{predictions}, they will be 
                    returned to the previous names. 
                    
          \item{differences}{A \code{matrix} containing all pairwise differences between
                    the predictions; it should have the same number of rows and columns as there are
                    rows in \code{predictions}.}
          \item \code{p.differences}: A \code{matrix} containing p-values for all pairwise differences 
                    between the predictions; each p-value is computed as the probability of a t-statistic 
                    as large as or larger than the observed difference divided by its standard error. 
                    The degrees of freedom of the t distribution for computing it are computed as 
                    the denominator degrees of freedom of the F value for the fixed term, if available; 
                    otherwise, the degrees of freedom stored in the attribute \code{tdf} are used; 
                    the matrix should be of the same size as that for \code{differences}.
          \item \code{sed}: A \code{matrix} containing the standard errors of all pairwise differences 
                    between the predictions; they are used in computing the p-values in \code{p.differences}.
          \item \code{vcov}: A \code{matrix} containing the variance matrix of the predictions; it is 
                    used in computing the variance of linear transformations of the predictions.
          \item \code{LSD}: A \code{data.frame} containing the mean, minimum and maximum LSD for determining             
                    the significance of pairwise differences, the mean LSD being calculated using 
                    the square root of the mean of the variances of pairwise differences. 
                    If \code{factor.combination} was specified for \code{meanLSD.type} when the 
                    LSDs were being calculated, then \code{LSD} contains an LSD for each 
                    \code{factor.combination} of the \code{\link{factor}s} specified by \code{LSDby}.  
                    Each LSD is calculated from the square root of the mean of the variances for all 
                    pairwise differences for each factor combination, unless there is only one prediction 
                    for a \code{factor.combination}, when notional LSDs are calculated that are based 
                    on the standard error of the prediction multiplied by the square root of two. 
                    If \code{LSD} is not \code{NULL} then the overall mean LSD will be added as 
                    an attribute named \code{meanLSD} of the \code{\link{alldiffs.object}}, as will 
                    the values of \code{meanLSD.type} and \code{LSDby}. The LSD for a single prediction 
                    assumes that any predictions to be compared are independent; this is not the case if 
                    residual errors are correlated.
\item \code{backtransforms}: When the response values have been transformed for analysis, 
                    a \code{data.frame} containing the backtransformed values of the 
                    predicted values is added to the \code{alldiffs.object}. This \code{data.frame} 
                    is consistent with the \code{predictions} component, except 
                    that the column named \code{predicted.value} is replaced by one called 
                    \code{backtransformed.predictions}. Any \code{error.interval} values will also 
                    be the backtransformed values. Each row contains a single predicted value.
}

          The details of the \code{attributes} of an \code{\link{alldiffs.object}} are:
    \enumerate{
          \item \code{response}: A \code{character} specifying the response variable for the 
                    predictions.
          \item \code{response.title}: A \code{character} specifying the title for the response variable 
                    for the predictions.
          \item \code{term}: A character string giving the variables that define the term 
                 that was fitted using \code{asreml} and that corresponds 
                 to \code{classify}. It is often the same as \code{classify}.
          \item \code{classify}: A character string giving the variables that define the margins
                 of the multiway table used in the prediction. Multiway tables are 
                 specified by forming an interaction type term from the 
                 classifying variables, that is, separating the variable names 
                 with the \code{:} operator.
          \item \code{tdf}: An \code{integer} specifying the degrees of freedom of the standard error. It is used as 
                the degrees of freedom for the t-distribution on which p-values and confidence 
                intervals are based.
          \item \code{meanLSD}: If the \code{LSD} component is not \code{NULL} then the mean LSD is added as an 
                attribute, calculated using the square root of the mean of the variances of pairwise differences.
          \item \code{meanLSD.type}: If the \code{LSD} component is not \code{NULL} then \code{meanLSD.type} is 
                added as an attribute.
          \item \code{LSDby}: If the \code{LSD} component is not \code{NULL} then \code{LSDby} is added as an attribute.
          \item \code{sortFactor}: A \code{character} containing the name of the 
              \code{factor} that indexes the set of predicted values that 
              determined the sorting of the components.
          \item \code{sortOrder}: A \code{character} vector that is the same length as the number of levels for 
              \code{sortFactor} in the \code{predictions} component of the 
              \code{\link{alldiffs.object}}. It specifies the order of the 
              levels in the reordered components of the \code{\link{alldiffs.object}}.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: \cr
              \code{levs <- levels(f)[order(x)]}.
}
     See \code{\link{predictPlus.asreml}} for more information.
}
\author{Chris Brien}
\seealso{\code{\link{is.alldiffs}}, \code{\link{as.alldiffs}}, \code{\link{validAlldiffs}}, \code{\link{allDifferences.data.frame}}}
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
    ## Form an all.diffs object
     Var.diffs <- as.alldiffs(predictions = Var.preds, classify = "Nitrogen:Variety", 
                              sed = Var.sed, vcov = Var.vcov, tdf = den.df)

    ## Check the class and validity of the alldiffs object
    is.alldiffs(Var.diffs)
    validAlldiffs(Var.diffs)
  }
}
\keyword{asreml}
\keyword{htest}