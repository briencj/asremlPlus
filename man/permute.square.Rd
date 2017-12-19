\name{permute.square}
\alias{permute.square}
\title{Permutes the rows and columns of a square matrix.}
\description{Permutes the rows and columns of a square matrix.}
\usage{permute.square(x, permutation)}
\arguments{
 \item{x}{A square matrix.}
 \item{permutation}{A vector specifying the new order of rows and columns.}
}
\value{A square \code{matrix}.}
\author{Chris Brien}
\seealso{\code{\link{permute.to.zero.lowertri}}}
\examples{
  terms.marginality <-  matrix(c(1,0,0,0,0,  0,1,0,0,0, 0,1,1,0,0, 
                                 1,1,1,1,0, 1,1,1,1,1), nrow=5)
  permtn <- c(1,3,2,4,5)
  terms.marginality <- permute.square(terms.marginality, permtn)
}
\keyword{array}
