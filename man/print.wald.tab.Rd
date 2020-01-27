\name{print.wald.tab}
\alias{print.wald.tab}
\title{Prints a data.frame containing a Wald or pseudoanova table.}
\description{Prints a \code{wald.tab} with or without title and/or heading. 
             The printing of the p-values is limited to 4 digits.}
\usage{\method{print}{wald.tab}(x, which.wald = c("title", "heading", "table"), 
      colourise = FALSE, ...)}
\arguments{
 \item{x}{An object that, ideally, is of class \code{wald.tab}.}
 \item{which.wald}{Which elements of the \code{wald.tab} to print. Possible
             values are some combination of \code{title}, \code{heading}, 
             \code{table} and \code{all}. The \code{heading} is an 
             attribute of \code{x}.}
 \item{colourise}{A \code{\link{logical}} which, if \code{TRUE}, results in the header 
                  text produced by \code{wald.asreml} being displayed in a different colour, if 
                  supported by the output terminal device. It overrides the \code{TRUE} 
                  setting of the \code{colourise} argument of \code{asreml::asreml.options}.}
 \item{\dots}{further arguments passed to \code{print} and \code{print.wald.tab}, 
                  but is only operational when the \code{table} is also printed.}
}
\value{No value is returned, but \code{x} is printed as specified in \code{which.wald}.}


\author{Chris Brien}
\seealso{\code{\link{print.test.summary}}, \code{\link{print.asrtests}}, \code{\link{as.asrtests}}, \code{\link{asremlPlus-package}}}
\examples{\dontrun{
data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
current.asrt <- rmboundary.asrtests(current.asrt)
# Test Row autocorrelation
current.asrt <- testresidual(current.asrt, "~ Row:ar1(Column)", 
                             label="Row autocorrelation", simpler=TRUE)
print(current.asrt$wald.tab)
}}
\keyword{asreml}
\keyword{htest}