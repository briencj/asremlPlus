\name{print.asrtests}
\alias{print.asrtests}
\title{Prints the values in an \code{\link{asrtests.object}}}
\description{Prints a \code{summary} of the asreml object, the \code{pseusdanova} and 
             the test.summary \code{data.frame} that are stored in the 
             \code{\link{asrtests.object}}.}
\usage{\method{print}{asrtests}(x, which = "all", colourise = FALSE, ...)}
\arguments{
 \item{x}{An \code{\link{asrtests.object}}.}
 \item{which}{Which elements of the \code{\link{asrtests.object}} to print. Possible
             values are some combination of \code{asremlsummary}, \code{vparametersummary},
             \code{pseudoanova}, \code{wald.tab}, \code{testsummary} and \code{key} or 
             \code{all}. The option \code{wald.tab} is a synonym for \code{pseudoanova}. 
             The options \code{key} and \code{all} are mutually exclusive; \code{key} 
             includes \code{vparametersummary}, but not the rest of \code{asremlsummary}, while 
             \code{all} includes the full \code{asremlsummary} that includes the 
             \code{vparametersummary}.}
 \item{colourise}{A \code{\link{logical}} which, if \code{TRUE}, results in the header 
                  text produced by \code{wald.asreml} being displayed in a different colour, if 
                  supported by the output terminal device. It overrides the \code{TRUE} 
                  setting of the \code{colourise} argument of \code{asreml::asreml.options}.}
 \item{\dots}{further arguments passed to \code{print} and \code{print.wald.tab}.}
}
\value{No value is returned, but the elements of the list in \code{x} are 
       printed.}


\author{Chris Brien}
\seealso{\code{\link{print.wald.tab}}, \code{\link{as.asrtests}}, \code{\link{asremlPlus-package}}}
\examples{\dontrun{
data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
current.asrt <- rmboundary.asrtests(current.asrt)
# Test Row autocorrelation
current.asrt <- testresidual(current.asrt, , "~ Row:ar1(Column)", 
                             label="Row autocorrelation", simpler=TRUE)
print(current.asrt)
}}
\keyword{asreml}
\keyword{htest}