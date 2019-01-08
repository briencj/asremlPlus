\name{validAlldiffs}
\alias{validAlldiffs}
\title{Checks that an object is a valid alldiffs object.}
\description{Checks that an \code{object} is an \code{\link{alldiffs.object}} of S3-class 
             \code{alldiffs} containing the components \code{asreml.obj}, 
             \code{wald.tab} and \code{test.summary}.}

\usage{
validAlldiffs(object)}
\arguments{
 \item{object}{an \code{\link{alldiffs.object}}.}
}
\value{\code{TRUE} or a \code{character} describing why the \code{object} 
       is not a valid \code{\link{alldiffs.object}}.}


\author{Chris Brien}
\seealso{\code{\link{alldiffs.object}}, \code{\link{is.alldiffs}}, \code{\link{as.alldiffs}}, \cr 
         \code{\link{validPredictionsFrame}}, \code{\link{validAsrtests}}}
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

    ## check the validity of Var.diffs
    validAlldiffs(Var.diffs)
  }
}
\keyword{asreml}
\keyword{htest}