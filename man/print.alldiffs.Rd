\name{print.alldiffs}
\alias{print.alldiffs}
\title{Prints the values in an \code{\link{alldiffs.object}} in a nice format.}
\description{Prints the predictions and standard errors from a fitted model, 
  including the attributes of the \code{\link{predictions.frame}}. Also prints out all 
  pairwise differences between the predictions to 2 signigicant figures,  
  along with their p-values and standard errors to 4 decimal places. If LSDs are 
  requested the mean, minimum and maximum LSDs will be printed.}
\usage{\method{print}{alldiffs}(x, which = "all", colourise = FALSE, ...)}
\arguments{
 \item{x}{An \code{\link{alldiffs.object}}.}
 \item{which}{A \code{character vector} containing a combination of 
              \code{predictions}, \code{vcov}, \code{backtransforms}, 
              \code{differences}, \code{p.differences}, \code{sed}, 
              \code{LSD} and \code{all}.
              These nominate which components of the \code{\link{alldiffs.object}} 
              to print.}
 \item{colourise}{A \code{\link{logical}} which, if \code{TRUE}, results in the header 
                  text produced by \code{predict.asreml} being displayed in a different colour, if 
                  supported by the output terminal device. It overrides the \code{TRUE} 
                  setting of the \code{colourise} argument of \code{asreml::asreml.options}.}
 \item{\dots}{further arguments passed to or from other methods.}
}
\value{No value is returned, but the components of \code{x} are printed.}


\author{Chris Brien}
\seealso{\code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}}
\examples{\dontrun{
print.alldiffs(diffs, which = "predictions")
}}
\keyword{asreml}