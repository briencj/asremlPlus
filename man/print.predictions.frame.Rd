\name{print.predictions.frame}
\alias{print.predictions.frame}
\title{Prints the values in a \code{\link{predictions.frame}}, with or without title and heading.}
\description{Prints the predictions from a fitted model, along with their standard errors and, if present, their error intervals, with or without title and headings.}
\usage{\method{print}{predictions.frame}(x, title = NULL,  
      which.predictions = c("title", "heading", "table"), 
      colourise = FALSE, ...)}
\arguments{
 \item{x}{An object that, ideally, is of class \code{\link{predictions.frame}}.}
 \item{title}{A \code{\link{character}} giving a title to be printed out before the 
              \code{heading} and \code{table} for the \code{\link{predictions.frame}},}
 \item{which.predictions}{what Which elements of the \code{\link{predictions.frame}} 
             to print. Possible values are some combination of \code{title}, 
             \code{heading}, \code{table} and \code{all}. The \code{heading} is an 
             attribute of \code{x}.}
 \item{colourise}{A \code{\link{logical}} which, if \code{TRUE}, results in the header 
                  text produced by \code{predict.asreml} being displayed in a different colour, if 
                  supported by the output terminal device. It overrides the \code{TRUE} 
                  setting of the \code{colourise} argument of \code{asreml::asreml.options}, 
                  but is only operational when the \code{table} is also printed.}
 \item{\dots}{further arguments passed to \code{print.predictions.frame}.}
}
\value{No value is returned, but the components of \code{x} are printed.}


\author{Chris Brien}
\seealso{\code{\link{print.alldiffs}}, \code{\link{as.alldiffs}}, \code{\link{allDifferences.data.frame}}}
\examples{\dontrun{
print.predictions.frame(diffs$predictions, which = "all")
}}
\keyword{asreml}