\name{getTestEntry.asrtests}
\alias{getTestEntry.asrtests}
\alias{getTestEntry}
\title{Gets the entry for a test recorded in the test.summary data.frame of an \code{\link{asrtests.object}}}
\description{Matches the \code{label} in the \code{term} column of the \code{test.summary} \code{data.frame} 
  in the supplied \code{\link{asrtests.object}} and extracts the line for it. It only matches 
  the last occurence of \code{label}.}
\usage{
\method{getTestEntry}{asrtests}(asrtests.obj, label, ...)}
\arguments{
 \item{asrtests.obj}{An \code{\link{asrtests.object}} containing the components (i) \code{asreml.obj},
       (ii) \code{wald.tab}, and (iii) \code{test.summary}.}
 \item{label}{A \code{character} specifying the \code{label} of the test for which  
              the entry is required. If \code{\link{testranfix.asrtests}} was used 
              for the test of interest, then the \code{label} will be the value of 
              the \code{term} argument supplied to \code{\link{testranfix.asrtests}}. 
              For \code{\link{changeModelOnIC.asrtests}}, the \code{label} will be the value of 
              the \code{label} argument. Other arguments will be relevant for other 
              test and change functions.}
 \item{\dots}{provision for passsing arguments to functions called internally - 
            not used at present.}
}
\value{A one-line \code{data.frame} containing the entry.
}

\author{Chris Brien}
\seealso{\code{\link{getTestPvalue.asrtests}}, \code{\link{as.asrtests}}, \cr
         \code{\link{testranfix.asrtests}}, 
         \code{\link{testswapran.asrtests}}, \code{\link{testresidual.asrtests}}, \cr    
         \code{\link{changeModelOnIC.asrtests}}, \code{\link{changeTerms.asrtests}}, 
         \code{\link{chooseModel.asrtests}}}
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
getTestEntry(current.asrt, label = "units")
}}
\keyword{asreml}
\keyword{htest}