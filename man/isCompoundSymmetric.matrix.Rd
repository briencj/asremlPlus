\name{isCompoundSymmetric.matrix}
\alias{isCompoundSymmetric.matrix}
\alias{isCompoundSymmetric}
\title{Tests whether an object of class matrix is compound symmetric}
\description{Tests whether an object of class \code{\link{matrix}} is compound symmetric by checking whether all diagonal elements are equal and all off-diagonal elements are equal.
}
\usage{
\method{isCompoundSymmetric}{matrix}(object, tol = 100 * .Machine$double.eps, ...) 
}

\arguments{
 \item{object}{An \code{object} to be tested.}
 \item{tol}{a \code{\link{numeric}} scalar > 0 specifying that values smaller than it 
            are considered to be zero.}
 \item{\dots}{Provision for passing arguments to functions called internally - 
            not used at present.}
}
\value{A \code{logical}.}
\author{Chris Brien}
\seealso{\code{\link{isSymmetric}}}
\examples{
  data(Oats.dat)
  
  \dontrun{
  ## Use asreml to get the variance matrix of a set of predictions

  m1.asr <- asreml(Yield ~ Nitrogen*Variety, 
                   random=~Blocks/Wplots,
                   data=Oats.dat)
  current.asrt <- as.asrtests(m1.asr)
  Var.pred <- asreml::predict.asreml(m1.asr, classify="Nitrogen:Variety", 
                                      vcov = TRUE)
                                    est.status = "status")
  isCompoundSymmetrix(Var.pred$vcov)
}
  
  ## Use lmerTest and emmmeans to get the variance matrix of a set of predictions
  if (requireNamespace("lmerTest", quietly = TRUE) & 
      requireNamespace("emmeans", quietly = TRUE))
  {
    m1.lmer <- lmerTest::lmer(Yield ~ Nitrogen*Variety + (1|Blocks/Wplots),
                              data=Oats.dat)
    Var.emm <- emmeans::emmeans(m1.lmer, specs = ~ Nitrogen:Variety)
    Var.vcov <- vcov(Var.emm)
  }

  if (exists("Var.vcov"))
  {
    ## Check the variance matrix of the predictions is compound symmetric
    isCompoundSymmetric(Var.vcov)
  }
}
\keyword{manip}