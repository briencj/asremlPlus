\name{num.recode}
\alias{num.recode}
\title{Recodes the unique values of a vector using the values in a new vector.}
\description{Recodes the unique values of a variate using the value in position
        i of the \code{new.values} vector to replace the ith 
        sorted \code{unique} values of x. The new levels do
             not have to be unique.}
\usage{num.recode(x, new.values)}
\arguments{
 \item{x}{The \code{vector} to be recoded.}
 \item{new.values}{A \code{vector} of length \code{unique(x)} containing values to use in 
                  the recoding.}
}
\value{A \code{vector}.}
\author{Chris Brien}
\seealso{\code{dae::fac.recast}.}
\examples{
## set up a factor with labels
x <- rep(c(-42, -14, 14, 42), 4)
 
## recode x
b <- num.recode(x, c(0, 28, 56, 84))
}
\keyword{manip}