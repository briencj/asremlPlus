\name{permute.to.zero.lowertri}
\alias{permute.to.zero.lowertri}
\title{Permutes a square matrix until all the lower triangular elements are zero.}
\description{Permutes a square matrix until all the lower triangular elements are zero.}
\usage{permute.to.zero.lowertri(x)}
\arguments{
 \item{x}{A square matrix of order n with at least n*(n-1)/2 zero elements.}
}
\value{A square \code{matrix}.}
\author{Chris Brien}
\seealso{\code{\link{permute.square}}}
\examples{
  terms.marginality <-  matrix(c(1,0,0,0,0,  0,1,0,0,0, 0,1,1,0,0, 
                                 1,1,1,1,0, 1,1,1,1,1), nrow=5)
  terms.marginality <- permute.to.zero.lowertri(terms.marginality)
}
\keyword{array}
