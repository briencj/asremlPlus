\name{is.asrtests}
\alias{is.asrtests}
\title{Tests whether an object is of class asrtests}
\description{A single-line \code{function} that tests whether an object is of \code{class} \code{asrtests}.}
\usage{is.asrtests(object)}
\arguments{
 \item{object}{An \code{object} to be tested.}
}
\value{A \code{logical}.}
\author{Chris Brien}
\seealso{\code{\link{asremlPlus-package}}, \code{\link{asrtests.object}},  \code{\link{is.asrtests}}, \code{\link{as.asrtests}}}
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

# check the class of current.asrt
is.asrtests(current.asrt)
}}
\keyword{manip}