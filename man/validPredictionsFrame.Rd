\name{validPredictionsFrame}
\alias{validPredictionsFrame}
\title{Checks that an object is a valid predictions.frame.}
\description{Checks that an \code{object} is a \code{\link{predictions.frame}} 
             of S3-class \code{data.frame} that contains the columns 
             predicted.value, standard.error and est.status.}

\usage{
validPredictionsFrame(object)}
\arguments{
 \item{object}{an \code{\link{predictions.frame}}.}
}
\value{\code{TRUE} or a \code{character} describing why the \code{object} 
       is not a valid \code{\link{predictions.frame}}.}


\author{Chris Brien}
\seealso{\code{\link{predictions.frame}}, \code{\link{validAsrtests}}, \code{\link{validAlldiffs}}}
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
    ## Check the class and validity of the predictions.frame
    is.predictions.frame(Var.preds)
    validPredictionsFrame(Var.preds)
  }
}
\keyword{asreml}
\keyword{htest}