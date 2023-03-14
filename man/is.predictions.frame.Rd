\name{is.predictions.frame}
\alias{is.predictions.frame}
\title{Tests whether an object is of class predictions.frame}
\description{A single-line \code{function} that tests whether an object is of \code{class} \code{predictions.frame}.}
\usage{is.predictions.frame(object)}
\arguments{
 \item{object}{An \code{object} to be tested.}
}
\value{A \code{logical}.}
\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{predictions.frame}},  \code{\link{validPredictionsFrame}}, \code{\link{as.predictions.frame}}}
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
  }
}
\keyword{manip}