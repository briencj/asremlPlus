\name{is.alldiffs}
\alias{is.alldiffs}
\title{Tests whether an object is of class alldiffs}
\description{A single-line \code{function} that tests whether an object is of \code{class} \code{alldiffs}.}
\usage{is.alldiffs(object)}
\arguments{
 \item{object}{An \code{object} to be tested.}
}
\value{A \code{logical}.}
\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{alldiffs.object}},  \code{\link{is.alldiffs}}, \code{\link{as.alldiffs}}}
\examples{
  data(Oats.dat)
  
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

  ## Form an all.diffs object
   Var.diffs <- as.alldiffs(predictions = Var.preds, classify = "Nitrogen:Variety", 
                            sed = Var.sed, vcov = Var.vcov, tdf = den.df)

  ## check the class of Var.diffs
  is.alldiffs(Var.diffs)
  }
}
\keyword{manip}