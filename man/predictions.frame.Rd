\name{predictions.frame}
\alias{predictions.frame}
\alias{predictions.frame-class}
\title{Description of a predictions object}
\description{A \code{data.frame} of S3-class \code{predictions.frame} that stores the 
             predictions for a fitted model.

             \code{\link{as.predictions.frame}} is function that converts a 
             \code{\link{data.frame}} to an object of this class.

             \code{\link{is.predictions.frame}} is the membership function for this class; it tests 
             that an object has class \code{predictions.frame}.

             \code{\link{validPredictionsFrame}} can be used to test the validity of a 
             \code{predictions.frame}.
}
\value{A \code{data.frame} being a \code{data.frame} beginning with the variables 
                    classifying the predictions, in the same order as in the 
                    \code{classify}, and also containing columns named \code{standard.error} 
                    and \code{est.status}; each row contains a single predicted value. 
                    The usual name of the column containing the predictions is  
                    \code{predicted.value} or \code{backtransformed.predictions}. 
                    The number of rows should equal the number of unique combinations 
                    of the classifying variables. While such a \code{data.frame} 
                    can be constructed from the beginning, the \code{pvals} component 
                    of the value produced by \code{predict.asreml} is a suitable value to 
                    supply for this argument. Note that the names \code{standard.error} and 
                    \code{est.status} have been changed to \code{std.error} and \code{status} 
                    in the \code{pvals} component produced by \code{asreml-R4}; if the new names 
                    are in the \code{data.frame} supplied to \code{predictions}, they will be 
                    returned to the previous names. 

                    The \code{data.frame} may also include columns for the lower and upper 
                    values of error intervals, either standard error, confidence or half-LSD 
                    intervals.  The names of these columns will consist of three parts 
                    separated by full stops: 
                     1) the first part will be \code{lower} or \code{upper}; 
                     2) the second part will be one of \code{Confidence}, 
                        \code{StandardError} or \code{halfLeastSignificant}; 
                     3) the third component will be \code{limits}.
                    
     See \code{\link{predictPlus.asreml}} for more information.
}
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
  Var.preds <- as.predictions.frame(Var.pred$pvals, se = "std.error", 
                                    est.status = "status")
  }
  
  ## Use lmerTest and emmmeans to get predictions and associated statistics
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
    Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
    Var.preds <- summary(Var.emm)
    Var.preds <- as.predictions.frame(Var.preds, predictions = "emmean", 
                                      se = "SE", interval.type = "CI", 
                                      interval.names = c("lower.CL", "upper.CL"))
  }

  if (exists("Var.preds"))
  {
    ## Check the class and validity of the alldiffs object
    is.predictions.frame(Var.preds)
    validPredictionsFrame(Var.preds)
  }
}

\author{Chris Brien}
\seealso{\code{\link{predictPlus.asreml}}, \code{\link{is.predictions.frame}}, 
          \code{\link{as.predictions.frame}}, \code{\link{validPredictionsFrame}}}
\keyword{asreml}
\keyword{htest}