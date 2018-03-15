\name{print.alldiffs}
\alias{print.alldiffs}
\title{Prints the values in an alldiffs object in a nice format.}
\description{Prints the predictions and standard errors from using asreml 
  to fit models in the same way as \code{asreml} prints them. Also prints out all 
  pairwise differences between the predictions to 2 signigicant figures,  
  along with their p-values and standard errors to 4 decimal places. If LSDs are 
  requested the mean, minimum and maximum LSDs will be printed.}
\usage{\method{print}{alldiffs}(x, which = "all", ...)}
\arguments{
 \item{x}{An object of S3-class \code{alldiffs}.}
 \item{which}{A \code{character vector} containing a combination of 
              \code{predictions}, \code{vcov}, \code{backtransforms}, 
              \code{differences}, \code{p.differences}, \code{sed}, 
              \code{LSD} and \code{all}.
              These nominate which components of the \code{alldiffs} 
              object to print.}
 \item{\dots}{further arguments passed to or from other methods.}
}
\value{No value is returned, but the components of \code{x} are printed.}


\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}}
\examples{\dontrun{
print.alldiffs(diffs, which = "predictions")
}}
\keyword{asreml}