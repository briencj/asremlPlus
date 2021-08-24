\name{print.LSDdata}
\alias{print.LSDdata}
\title{Prints the components of a list containing data on the LSDs for all pairwise differences of predictions.}
\description{Prints the components of an \code{LSDdata} \code{\link{list}} created by \code{\link{exploreLSDs}}, 
             that contains data on the LSDs for all pairwise differences of predictions stored in an 
             \code{\link{alldiffs.object}}.}
\usage{\method{print}{LSDdata}(x, which.print = c("statistics", "accuracy"), ...)}
\arguments{
 \item{x}{An object that, ideally, is of class \code{LSDdata}.}
 \item{which.print}{Which components of the \code{LSDdata} \code{\link{list}} to print. 
                    Possible values are some combination of \code{frequencies}, 
                    \code{distinct.vals}, \code{statistics},   \code{accuracy}, 
                    \code{per.pred.accuracy}, \code{LSD}, \code{summary} and \code{all}. 
                    Only one of \code{summary} and \code{all} is allowed. 
                    The option \code{summary} results in the printing of 
                    \code{distinct.vals}, \code{statistics}, \code{accuracy}. 
                    For a description of these components, see \code{\link{alldiffs.object}}.}
 \item{\dots}{further arguments passed to \code{print}.}
}
\value{No value is returned, but \code{x} is printed as specified in \code{which.print}.}


\author{Chris Brien}
\seealso{\code{\link{exploreLSDs.alldiffs}}, \code{\link{alldiffs.object}}}
\examples{\dontrun{
data(WaterRunoff.dat)
asreml.options(keep.order = TRUE) #required for asreml-R4 only
current.asr <- asreml(fixed = pH ~ Benches + (Sources * (Type + Species)), 
                      random = ~ Benches:MainPlots,
                      keep.order=TRUE, data= WaterRunoff.dat)
current.asrt <- as.asrtests(current.asr, NULL, NULL)
diffs <- predictPlus(classify = "Sources:Type", 
                     asreml.obj = current.asr, 
                     wald.tab = current.asrt$wald.tab, 
                     present = c("Sources", "Type", "Species"))
LSDdata <- exploreLSDs(diffs, LSDtype = "factor.combinations", LSDby = "Sources")
print(LSDdata)
}}
\keyword{asreml}
\keyword{htest}