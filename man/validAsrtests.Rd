\name{validAsrtests}
\alias{validAsrtests}
\title{Checks that an object is a valid asrtests object.}
\description{Checks that an \code{object} is an \code{\link{asrtests.object}} of S3-class 
             \code{asrtests} containing the components \code{asreml.obj}, 
             \code{wald.tab} and \code{test.summary}.}

\usage{
validAsrtests(object)}
\arguments{
 \item{object}{an \code{\link{asrtests.object}}.}
}
\value{\code{TRUE} or a \code{character} describing why the \code{object} 
       is not a valid \code{\link{asrtests.object}}.}


\author{Chris Brien}
\seealso{\code{\link{asrtests.object}}, \code{\link{is.asrtests}}, \code{\link{as.asrtests}}, \cr 
         \code{\link{validPredictionsFrame}}, \code{\link{validAlldiffs}}}
\examples{\dontrun{
library(dae)
library(asreml)
library(asremlPlus)
## use ?Wheat.dat for data set details
data(Wheat.dat)

# Fit initial model
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)

# Load current fit into an asrtests object
current.asrt <- asrtests(current.asr, NULL, NULL)

# check validity of current.asrt
validAsrtests(current.asrt)
}}
\keyword{asreml}
\keyword{htest}