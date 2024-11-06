\name{print.test.summary}
\alias{print.test.summary}
\title{Prints a data.frame containing a test.summary.}
\description{Prints a \code{test.summary} (also a \code{choose.summary}) with or without a title and with p-values limited to 4-digits.}
\usage{\method{print}{test.summary}(x, which.print = c("title", "table"), 
      omit.columns = NULL, response = NULL, ...)}
\arguments{
 \item{x}{A object that, ideally, is of class \code{test.summary}.}
 \item{which.print}{A \code{character} specifying the aspects of the \code{test.summary} to print. 
              Possible values are some combination of \code{title}, 
             \code{table} and \code{all}.}
 \item{omit.columns}{A \code{character} specifying the columns of the \code{test.summary} table 
              to be omitted from the print. If \code{NULL}, none are omitted.}
 \item{response}{A \code{character} specifying the name of the response that the 
                 \code{test.summary} table is based on and is to be printed in the 
                 heading for the table. If \code{NULL}, no response name is printed.}
 \item{\dots}{further arguments passed to \code{print}, 
                  but is only operational when the \code{table} is also printed.}
}
\value{No value is returned, but \code{x} is printed, possibly with a title.}


\author{Chris Brien}
\seealso{\code{\link{print.wald.tab}}, \code{\link{print.asrtests}}, \code{\link{as.asrtests}}, \code{\link{asremlPlus-package}}}
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
print(current.asrt$test.summary)
}}
\keyword{asreml}
\keyword{htest}