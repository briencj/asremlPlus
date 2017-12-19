\name{angular.mod}
\alias{angular.mod}
\title{Applies the modified angular transformation to a vector of counts.}
\description{Applies the angular transformation to a vector of counts. A 
  modified transformation is used that is appropriate when N < 50 and the 
  proportion is not between 0.3 and 0.7. The transformation is given by 
  \eqn{\sin^{-1}{\frac{count + 0.375}{n + 0.75}}}%
  {arcsin( sqrt( (count + 0.375) / (n + 0.75) ) )}.}
\usage{angular.mod(count, n)}
\arguments{
 \item{count}{The numeric vector of counts.}
 \item{n}{The number(s) of observations from which the count(s) were obtained.}
}
\value{A numeric vector.}
\author{Chris Brien}
\seealso{\code{\link{angular}}, \code{\link{powerTransform}}.}
\examples{
n <-25
y <- rbinom(10, n, 0.5)
y <- c(y,0,n)
p.ang.mod <- angular.mod(y, n)
}
\keyword{manip}