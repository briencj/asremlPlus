\name{getTestPvalue.asrtests}
\alias{getTestPvalue.asrtests}
\alias{getTestPvalue}
\title{Gets the p-value for a test recorded in the test.summary data.frame of an \code{\link{asrtests.object}}}
\description{Matches the \code{label} in the \code{term} column of the supplied 
  \code{\link{asrtests.object}} and extracts the its p-value. It only matches 
  the first occurence of \code{label}.}
\usage{
\method{getTestPvalue}{asrtests}(asrtests.obj, label, ...)}
\arguments{
 \item{asrtests.obj}{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj},
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{label}{A \code{character} specifying the \code{label} of the test for which  
              the p-value is required. If \code{\link{testranfix.asrtests}} was used 
              for the test of interest, then the \code{label} will be the value of 
              the \code{term} argument supplied to \code{\link{testranfix.asrtests}}.}
 \item{\dots}{provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{An \code{numeric} containing the p-value. It can be \code{NA}, for example when a p-value could not be calculated.
}

\author{Chris Brien}
\seealso{\code{\link{as.asrtests}}, \code{\link{testranfix.asrtests}}, 
         \code{\link{testswapran.asrtests}}, \code{\link{testresidual.asrtests}}, \cr    
         \code{\link{changeTerms.asrtests}}, \code{\link{chooseModel.asrtests}}}
\examples{\dontrun{
data(Wheat.dat)
current.asr <- asreml(yield ~ Rep + WithinColPairs + Variety, 
                      random = ~ Row + Column + units,
                      residual = ~ ar1(Row):ar1(Column), 
                      data=Wheat.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
current.asrt <- rmboundary(current.asrt)
# Test nugget term
current.asrt <- testranfix(current.asrt, "units", positive=TRUE)
getTestPvalue(current.asrt, label = "units")
}}
\keyword{asreml}
\keyword{htest}