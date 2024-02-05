\name{loadASRemlVersion}
\alias{loadASRemlVersion}
\title{Ensures that a specific version of \code{asreml} is loaded.}
\description{Loads the specified version of \code{asreml}, provided that it is not already loaded. If the version of \code{asreml} is not the required version, then the loaded version is unloaded first.}
\usage{loadASRemlVersion(version = 4, ...)}
\arguments{
 \item{version}{The \code{version} that is to be loaded, the version consisting of just the initial 
                characters that are significant in the version that should be loaded. 
                For example, the default value of 4 implies that any version that 
                begins with \code{"4"} is acceptable. It is used to check that the 
                required version is loaded.}
 \item{\dots}{Other library/require arguments that are needed to load the specified \code{version} of \code{asreml}.}
}
\value{A character, being all characters in the version of asreml that is loaded on exit from the function.}
\author{Chris Brien}
\seealso{\code{\link{getASRemlVersionLoaded}}.}
\examples{\dontrun{
loadASRemlVersion(3, lib.loc = "D:\\Analyses\\R asreml3")}}
\keyword{asreml}