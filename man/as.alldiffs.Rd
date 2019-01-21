\name{as.alldiffs}
\alias{as.alldiffs}
\title{Forms an \code{\link{alldiffs.object}} from the supplied predictions, along with those statistics, 
       associated with the predictions and their pairwise differences, that have been supplied.}
\description{Creates an \code{\link{alldiffs.object}} that consists of a \code{list} 
             containing the following components: \code{predictions}, \code{vcov}, \code{differences}, 
             \code{p.differences}, \code{sed}, \code{LSD} and \code{backtransforms}. 
             Predictions must be supplied to the function while the others will be set 
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
 \item{predictions}{A \code{\link{predictions.frame}}, being a \code{data.frame} beginning 
                    with the variables classifying the predictions and also containing columns 
                    named \code{predicted.value}, \code{standard.error} and \code{est.status}; 
                    each row contains a single predicted value. It may also contain columns 
                    for the lower and upper limits of error intervals for the predictions.
                    Note that the names \code{standard.error} and 
                    \code{est.status} have been changed to \code{std.error} and \code{status} 
                    in the \code{pvals} component produced by \code{asreml-R4}; if the new names 
                    are in the \code{data.frame} supplied to \code{predictions}, they will be 
                    returned to the previous names.}
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
 \item{vcov}{A \code{matrix} containing the variance matrix of the predictions; it is used in 
                    computing the variance of linear transformations of the predictions.}
\item{LSD}{A \code{data.frame} containing the mean, minimum and maximum LSD for determining             
                    the significance of pairwise differences, the mean LSD being calculated using 
                    the square root of the mean of the variances of pairwise differences. 
                    If \code{factor.combination} was specified for \code{meanLSD.type} when the 
                    LSDs were being calculated, then \code{LSD} contains an LSD for each 
                    \code{factor.combination} of the \code{\link{factor}s} specified by \code{LSDby}; 
                    each LSD is the square root of the mean of the variances for all 
                    pairwise differences for each factor combination. 
                    If \code{LSD} is not \code{NULL} then the mean LSD will be added as an attribute 
                    named \code{meanLSD} of the \code{\link{alldiffs.object}}.}
 \item{backtransforms}{A \code{data.frame} containing the backtransformed values of the predicted 
                    values  that is consistent with the \code{predictions} component, except 
                    that the column named \code{predicted.value} is replaced by one called 
                    \code{backtransformed.predictions}. Any \code{error.interval} values will also 
                    be the backtransformed values. Each row contains a single predicted value.}
 \item{response}{A \code{character} specifying the response variable for the 
                    predictions. It is stored as an attribute to the \code{\link{alldiffs.object}}.}
 \item{response.title}{A \code{character} specifying the title for the response variable 
                    for the predictions. It is stored as an attribute to the \code{\link{alldiffs.object}}.}
 \item{term}{A character string giving the variables that define the term 
                 that was fitted using \code{asreml} and that corresponds 
                 to \code{classify}. It is often the same as \code{classify}. 
                 It is stored as an attribute to the \code{\link{alldiffs.object}}.}
 \item{classify}{A character string giving the variables that define the margins
                 of the multiway table used in the prediction. Multiway tables are 
                 specified by forming an interaction type term from the 
                 classifying variables, that is, separating the variable names 
                 with the \code{:} operator. It is stored as an attribute to the 
                 \code{\link{alldiffs.object}}.}
 \item{tdf}{an \code{integer} specifying the degrees of freedom of the standard error. It is used as 
                the degrees of freedom for the t-distribution on which p-values and confidence 
                intervals are based. 
                It is stored as an attribute to the \code{\link{alldiffs.object}}.}
 \item{sortFactor}{A \code{character} containing the name of the 
              \code{factor} that indexes the set of predicted values that 
              determined the sorting of the components.}
 \item{sortOrder}{A \code{character} vector that is the same length as the number of levels for 
              \code{sortFactor} in the \code{predictions} component of the 
              \code{\link{alldiffs.object}}. It specifies the order of the 
              levels in the reordered components of the \code{\link{alldiffs.object}}.

              The following creates a \code{sortOrder} vector \code{levs} for factor 
              \code{f} based on the values in \code{x}: 
              \code{levs <- levels(f)[order(x)]}.}
}
\value{An S3-class \code{\link{alldiffs.object}}. Also, see \code{\link{predictPlus.asreml}} for 
              more information.}


\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{alldiffs.object}},  \code{\link{is.alldiffs}}, \code{\link{as.alldiffs}}, 
         \code{\link{print.alldiffs}}, \cr
         \code{\link{sort.alldiffs}}, \code{\link{subset.alldiffs}}, \code{\link{allDifferences.data.frame}}, \cr
         \code{\link{reorderClassify.alldiffs}}, \code{\link{redoErrorIntervals.alldiffs}},  \code{\link{recalcLSD.alldiffs}}, \cr
         \code{\link{predictPlus.asreml}}, \code{\link{plotPredictions.data.frame}},  \code{\link{predictPresent.asreml}}}
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