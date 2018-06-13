\name{getASRemlVersionLoaded}
\alias{getASRemlVersionLoaded}
\title{Finds the version of asreml that is loaded and returns the initial characters in version.}
\description{Checks that \code{asreml} is loaded and, if it is, returns the first \code{nchar} characters of the version that is loaded.}
\usage{getASRemlVersionLoaded(nchar = NULL, notloaded.fault = FALSE)}
\arguments{
 \item{nchar}{The number of characters in the asreml version to get.}
 \item{notloaded.fault}{A \code{\link{logical}} indicating whether a fault is to occur 
       if \code{asreml} is not loaded.}
}
\value{A character, being the first \code{nchar} characters of the version of \code{asreml} 
       that is loaded.}
\author{Chris Brien}
\seealso{\code{\link{loadASRemlVersion}}.}
\examples{\dontrun{
getASRemlVersionLoaded()}}
\keyword{asreml}