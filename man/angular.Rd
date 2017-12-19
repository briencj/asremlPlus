\name{angular}
\alias{angular}
\title{Applies the angular transformation to proportions.}
\description{Applies the angular transformation to numeric values.
  It is given by \eqn{\sin^{-1}(\sqrt{proportions})}{arcsin( sqrt( proportions ) )}}
\usage{angular(proportions, n)}
\arguments{
 \item{proportions}{The proportions.}
 \item{n}{The divisor(s) for each proportion}
}
\value{A numeric.}
\author{Chris Brien}
\seealso{\code{\link{angular.mod}}, \code{\link{powerTransform}}.}
\examples{
n <-25
y <- rbinom(10, n, 0.5)
y <- c(y,0,n)
p <- y/n
p.ang <- angular(p, n)
}
\keyword{manip}